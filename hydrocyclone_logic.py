import numpy as np
from scipy.optimize import minimize
from typing import List, Tuple, Optional
from hydrocyclone_models import (
    HydrocycloneAnalysisRequest, HydrocycloneAnalysisResponse, 
    PartitionCurvePoint, GranulometryPoint, BalanceRow, WaterBalance,
    TrompParameters, GlobalBalance, FlowData
)

def _interpolate_size(sizes, efficiencies, target_eff):
    """Interpola un tamaño para una eficiencia dada en escala logarítmica."""
    try:
        log_sizes = np.log10(sizes)
        # Asegurar que esté ordenado para interp
        s_idx = np.argsort(efficiencies)
        eff_s = efficiencies[s_idx]
        log_s_s = log_sizes[s_idx]
        
        if np.max(eff_s) >= target_eff >= np.min(eff_s):
            return 10**np.interp(target_eff, eff_s, log_s_s)
        return 0.0
    except:
        return 0.0

def _calculate_d50c(sizes, efficiencies):
    return _interpolate_size(sizes, efficiencies, 0.5)

def analyze_hydrocyclone(request: HydrocycloneAnalysisRequest) -> HydrocycloneAnalysisResponse:
    # 1. Preparación de Datos Experimentales
    sizes = np.array([s.mesh_size for s in request.sieves])
    # Ordenar por tamaño descendente para coherencia
    sort_idx = np.argsort(sizes)[::-1]
    sizes = sizes[sort_idx]
    
    w_f = np.array([request.sieves[i].weight_feed for i in sort_idx])
    w_o = np.array([request.sieves[i].weight_overflow for i in sort_idx])
    w_u = np.array([request.sieves[i].weight_underflow for i in sort_idx])
    
    # Incluir fondo (Pan) al final
    f_raw = np.append(w_f, request.pan_feed)
    o_raw = np.append(w_o, request.pan_overflow)
    u_raw = np.append(w_u, request.pan_underflow)
    
    tot_f, tot_o, tot_u = np.sum(f_raw), np.sum(o_raw), np.sum(u_raw)
    if tot_f <= 0: return _empty_response()

    # Porcentajes parciales experimentales
    p_f_exp = (f_raw / tot_f) * 100
    p_o_exp = (o_raw / tot_o) * 100
    p_u_exp = (u_raw / tot_u) * 100

    # 2. Optimización de la Recuperación de Sólidos (S)
    def obj_s(S):
        error = 0
        for f, o, u in zip(p_f_exp, p_o_exp, p_u_exp):
            balance_calc = (1 - S) * o + S * u
            error += (f - balance_calc)**2
        return error

    res_s = minimize(obj_s, 0.5, bounds=[(0.01, 0.99)])
    S_opt = res_s.x[0]

    # 3. Reconciliación de Porcentajes por Malla
    p_f_adj, p_o_adj, p_u_adj = [], [], []
    for f_e, o_e, u_e in zip(p_f_exp, p_o_exp, p_u_exp):
        def obj_mesh(x):
            f, o, u = x
            return (f - f_e)**2 + (o - o_e)**2 + (u - u_e)**2
        
        cons = {'type': 'eq', 'fun': lambda x: x[0] - ((1 - S_opt) * x[1] + S_opt * x[2])}
        res_m = minimize(obj_mesh, [f_e, o_e, u_e], constraints=cons, bounds=[(0, 100), (0, 100), (0, 100)])
        p_f_adj.append(res_m.x[0]), p_o_adj.append(res_m.x[1]), p_u_adj.append(res_m.x[2])

    p_f_adj, p_o_adj, p_u_adj = np.array(p_f_adj), np.array(p_o_adj), np.array(p_u_adj)
    p_f_adj = (p_f_adj / np.sum(p_f_adj)) * 100
    p_o_adj = (p_o_adj / np.sum(p_o_adj)) * 100
    p_u_adj = (p_u_adj / np.sum(p_u_adj)) * 100

    # 4. Cálculo de Eficiencias y Curvas
    E_a_adj_all = np.clip(S_opt * p_u_adj / p_f_adj, 0, 1)
    Rf = float(E_a_adj_all[-1])
    E_c_adj_all = (E_a_adj_all - Rf) / (1 - Rf) if Rf < 1 else np.zeros_like(E_a_adj_all)
    E_c_adj_all = np.clip(E_c_adj_all, 0, 1)

    d50c_exp = _calculate_d50c(sizes, (S_opt * p_u_exp[:-1] / p_f_exp[:-1] - Rf)/(1-Rf))
    d50c_adj = _calculate_d50c(sizes, E_c_adj_all[:-1])
    
    # 4.1 Parámetros de Tromp
    d25c = _interpolate_size(sizes, E_c_adj_all[:-1], 0.25)
    d75c = _interpolate_size(sizes, E_c_adj_all[:-1], 0.75)
    tromp = TrompParameters(
        d25c=float(d25c), d50c=float(d50c_adj), d75c=float(d75c),
        imperfection=float((d75c - d25c) / (2 * d50c_adj)) if d50c_adj > 0 else 0
    )

    # 5. Balance de Agua y Diagnóstico
    # ... (código existente hasta el diagnóstico)
    Rw = None
    E_pan = float(E_a_adj_all[-1])
    consistency_err = None
    diag_msg = "Datos procesados correctamente."
    diag_level = "success"

    # --- NUEVAS TABLAS SOLICITADAS ---
    
    # 1. TABLA BALANCE POR PESOS (O + U)
    # F = O + U en cada malla
    w_f_calc_weights = o_raw + u_raw
    tot_f_calc_weights = np.sum(w_f_calc_weights)
    S_weights = np.sum(u_raw) / tot_f_calc_weights if tot_f_calc_weights > 0 else 0
    
    p_f_weights = (w_f_calc_weights / tot_f_calc_weights) * 100
    p_o_weights = (o_raw / np.sum(o_raw)) * 100
    p_u_weights = (u_raw / np.sum(u_raw)) * 100
    
    balance_weights_table = []
    for i in range(len(sizes)):
        rec_malla = (S_weights * p_u_weights[i] / p_f_weights[i]) if p_f_weights[i] > 0 else 0
        balance_weights_table.append(BalanceRow(
            size=f"{sizes[i]} µm",
            feed_pct=float(p_f_weights[i]), overflow_pct=float(p_o_weights[i]), underflow_pct=float(p_u_weights[i]),
            recovery_underflow=float(rec_malla)
        ))
    balance_weights_table.append(BalanceRow(
        size="Fondo (Pan)",
        feed_pct=float(p_f_weights[-1]), overflow_pct=float(p_o_weights[-1]), underflow_pct=float(p_u_weights[-1]),
        recovery_underflow=float(S_weights * p_u_weights[-1] / p_f_weights[-1] if p_f_weights[-1] > 0 else 0)
    ))

    # 2. TABLA BALANCE POR % SOLIDOS (F recalcuado)
    balance_solids_table = None
    if request.feed_p_solids and request.overflow_p_solids and request.underflow_p_solids:
        cp, co, cu = request.feed_p_solids/100, request.overflow_p_solids/100, request.underflow_p_solids/100
        # Fórmula de Recuperación de Sólidos S (Material Balance de 2 productos)
        if (cu - co) > 0 and cp > 0:
            S_solids = (cu * (cp - co)) / (cp * (cu - co)) # Simplifica a (cp-co)/(cu-co) * (cu/cp)
            # Asegurar S_solids entre 0 y 1 por errores de muestreo
            S_solids = max(0.01, min(0.99, S_solids))
            
            # Recalcular alimento granulométrico: F_i = (1-S)*O_i + S*U_i
            p_o_norm = (o_raw / np.sum(o_raw)) * 100
            p_u_norm = (u_raw / np.sum(u_raw)) * 100
            p_f_solids_calc = (1 - S_solids) * p_o_norm + S_solids * p_u_norm
            
            balance_solids_table = []
            for i in range(len(sizes)):
                rec_malla = (S_solids * p_u_norm[i] / p_f_solids_calc[i]) if p_f_solids_calc[i] > 0 else 0
                balance_solids_table.append(BalanceRow(
                    size=f"{sizes[i]} µm",
                    feed_pct=float(p_f_solids_calc[i]), overflow_pct=float(p_o_norm[i]), underflow_pct=float(p_u_norm[i]),
                    recovery_underflow=float(rec_malla)
                ))
            balance_solids_table.append(BalanceRow(
                size="Fondo (Pan)",
                feed_pct=float(p_f_solids_calc[-1]), overflow_pct=float(p_o_norm[-1]), underflow_pct=float(p_u_norm[-1]),
                recovery_underflow=float(S_solids * p_u_norm[-1] / p_f_solids_calc[-1] if p_f_solids_calc[-1] > 0 else 0)
            ))

    # --- FIN NUEVAS TABLAS ---

    if request.feed_p_solids and request.overflow_p_solids and request.underflow_p_solids:
        cp, co, cu = request.feed_p_solids/100, request.overflow_p_solids/100, request.underflow_p_solids/100
        Rw = S_opt * ((1 - cu) / cu) * (cp / (1 - cp))
        consistency_err = abs(E_pan - Rw)
        
        if consistency_err < 0.02:
            diag_msg = "Excelente consistencia entre el fondo (pan) y el balance de agua."
            diag_level = "success"
        elif consistency_err < 0.05:
            diag_msg = "Consistencia aceptable."
            diag_level = "warning"
        else:
            diag_msg = "Atención: Discrepancia significativa (>5%) entre el fondo y el agua."
            diag_level = "error"

    # 5.1 Balance Global Absoluto
    global_bal = None
    if request.feed_flow_rate and request.feed_p_solids:
        cp = request.feed_p_solids / 100
        co = (request.overflow_p_solids or 0) / 100
        cu = (request.underflow_p_solids or 0) / 100
        
        # Calcular masa de sólidos en alimento (tph)
        if request.feed_flow_unit == "tph":
            f_mass_s = request.feed_flow_rate
        else: # m3h de pulpa
            factor = (1 / request.solid_density) + ((1 - cp) / cp)
            f_mass_s = request.feed_flow_rate / factor
            
        def get_flow_data(mass_s, p_sol):
            p_sol = max(p_sol, 0.001)
            mass_w = mass_s * (1 - p_sol) / p_sol
            vol_s = mass_s / request.solid_density
            vol_w = mass_w / 1.0
            return FlowData(
                mass_solids=float(mass_s), mass_water=float(mass_w),
                vol_solids=float(vol_s), vol_water=float(vol_w),
                vol_pulp=float(vol_s + vol_w), p_solids=float(p_sol * 100)
            )

        global_bal = GlobalBalance(
            feed=get_flow_data(f_mass_s, cp),
            overflow=get_flow_data(f_mass_s * (1 - S_opt), co or cp),
            underflow=get_flow_data(f_mass_s * S_opt, cu or cp)
        )

    water_bal = WaterBalance(
        solids_recovery_S=float(S_opt * 100),
        water_recovery_Rw=float(Rw * 100) if Rw is not None else None,
        bypass_Rf=float(Rf * 100),
        pan_recovery_Epan=float(E_pan * 100),
        consistency_error=float(consistency_err * 100) if consistency_err is not None else None,
        feed_flow=100.0,
        overflow_flow=float((1 - S_opt) * 100),
        underflow_flow=float(S_opt * 100),
        global_balance=global_bal
    )

    # 6. Construcción de Respuesta
    f_pass_exp, o_pass_exp, u_pass_exp = 100 - np.cumsum(p_f_exp[:-1]), 100 - np.cumsum(p_o_exp[:-1]), 100 - np.cumsum(p_u_exp[:-1])
    f_pass_adj, o_pass_adj, u_pass_adj = 100 - np.cumsum(p_f_adj[:-1]), 100 - np.cumsum(p_o_adj[:-1]), 100 - np.cumsum(p_u_adj[:-1])

    granulometry_pts = [
        GranulometryPoint(
            size=float(sizes[i]), 
            feed_passing=float(f_pass_exp[i]), overflow_passing=float(o_pass_exp[i]), underflow_passing=float(u_pass_exp[i]),
            feed_passing_adj=float(f_pass_adj[i]), overflow_passing_adj=float(o_pass_adj[i]), underflow_passing_adj=float(u_pass_adj[i])
        ) for i in range(len(sizes))
    ]

    partition_pts = [
        PartitionCurvePoint(
            size=float(sizes[i]), 
            actual_recovery=float(S_opt * p_u_exp[i] / p_f_exp[i]), 
            corrected_recovery=float((S_opt * p_u_exp[i] / p_f_exp[i] - Rf)/(1-Rf)),
            adjusted_recovery=float(E_a_adj_all[i])
        ) for i in range(len(sizes))
    ]

    balance_table = []
    for i in range(len(sizes)):
        balance_table.append(BalanceRow(
            size=f"{sizes[i]} µm",
            feed_pct=float(p_f_exp[i]), overflow_pct=float(p_o_exp[i]), underflow_pct=float(p_u_exp[i]),
            feed_pct_adj=float(p_f_adj[i]), overflow_pct_adj=float(p_o_adj[i]), underflow_pct_adj=float(p_u_adj[i]),
            recovery_underflow=float(E_a_adj_all[i])
        ))
    balance_table.append(BalanceRow(
        size="Fondo (Pan)",
        feed_pct=float(p_f_exp[-1]), overflow_pct=float(p_o_exp[-1]), underflow_pct=float(p_u_exp[-1]),
        feed_pct_adj=float(p_f_adj[-1]), overflow_pct_adj=float(p_o_adj[-1]), underflow_pct_adj=float(p_u_adj[-1]),
        recovery_underflow=float(E_a_adj_all[-1])
    ))

    return HydrocycloneAnalysisResponse(
        d50c_experimental=float(d50c_exp),
        d50c_adjusted=float(d50c_adj),
        tromp=tromp,
        diagnosis_message=diag_msg,
        diagnosis_level=diag_level,
        water_balance=water_bal,
        partition_curve=partition_pts,
        granulometry_curve=granulometry_pts,
        balance_table=balance_table,
        balance_solids_table=balance_solids_table,
        balance_weights_table=balance_weights_table,
        summary={"iterations": 1, "status": "Success"}
    )

def _empty_response():
    # Implementación básica para evitar errores
    return HydrocycloneAnalysisResponse(
        d50c_experimental=0, d50c_adjusted=0, 
        water_balance=WaterBalance(solids_recovery_S=0, bypass_Rf=0, feed_flow=100, overflow_flow=0, underflow_flow=0),
        partition_curve=[], granulometry_curve=[], balance_table=[], summary={}
    )
