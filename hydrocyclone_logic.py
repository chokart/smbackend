import numpy as np
from scipy.optimize import minimize
from typing import List, Tuple, Optional
from hydrocyclone_models import (
    HydrocycloneAnalysisRequest, HydrocycloneAnalysisResponse, 
    PartitionCurvePoint, GranulometryPoint, BalanceRow, WaterBalance,
    TrompParameters, GlobalBalance, FlowData, HydrocycloneMetrics
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

def analyze_hydrocyclone(request: HydrocycloneAnalysisRequest) -> HydrocycloneAnalysisResponse:
    # 1. Preparación de Datos Experimentales
    sizes = np.array([s.mesh_size for s in request.sieves])
    sort_idx = np.argsort(sizes)[::-1]
    sizes = sizes[sort_idx]
    
    w_f = np.array([request.sieves[i].weight_feed for i in sort_idx])
    w_o = np.array([request.sieves[i].weight_overflow for i in sort_idx])
    w_u = np.array([request.sieves[i].weight_underflow for i in sort_idx])
    
    f_raw = np.append(w_f, request.pan_feed)
    o_raw = np.append(w_o, request.pan_overflow)
    u_raw = np.append(w_u, request.pan_underflow)
    
    tot_f, tot_o, tot_u = np.sum(f_raw), np.sum(o_raw), np.sum(u_raw)
    if tot_f <= 0: return _empty_response()

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

    # 4. Cálculo de Eficiencias y Métricas (Reconciliado)
    E_a_adj_all = np.clip(S_opt * p_u_adj / p_f_adj, 0, 1)
    Rf_adj = float(E_a_adj_all[-1])
    E_c_adj_all = (E_a_adj_all - Rf_adj) / (1 - Rf_adj) if Rf_adj < 1 else np.zeros_like(E_a_adj_all)
    E_c_adj_all = np.clip(E_c_adj_all, 0, 1)

    d50_adj = _interpolate_size(sizes, E_a_adj_all[:-1], 0.5)
    d50c_adj = _interpolate_size(sizes, E_c_adj_all[:-1], 0.5)
    d25c_adj = _interpolate_size(sizes, E_c_adj_all[:-1], 0.25)
    d75c_adj = _interpolate_size(sizes, E_c_adj_all[:-1], 0.75)
    
    reconciled_metrics = HydrocycloneMetrics(
        d50=float(d50_adj), d50c=float(d50c_adj), bypass_rf=float(Rf_adj * 100),
        imperfection=float((d75c_adj - d25c_adj) / (2 * d50c_adj)) if d50c_adj > 0 else 0,
        solids_recovery_s=float(S_opt * 100)
    )

    # 5. Balance por % Sólidos
    balance_solids_table = None
    solids_metrics = None
    S_solids = 0
    if request.feed_p_solids and request.overflow_p_solids and request.underflow_p_solids:
        cp, co, cu = request.feed_p_solids/100, request.overflow_p_solids/100, request.underflow_p_solids/100
        if (cu - co) > 0 and cp > 0:
            S_solids = (cu * (cp - co)) / (cp * (cu - co))
            S_solids = max(0.01, min(0.99, S_solids))
            
            p_o_norm = (o_raw / np.sum(o_raw)) * 100
            p_u_norm = (u_raw / np.sum(u_raw)) * 100
            p_f_sol_calc = (1 - S_solids) * p_o_norm + S_solids * p_u_norm
            
            E_a_sol_all = np.clip(S_solids * p_u_norm / p_f_sol_calc, 0, 1)
            Rf_sol = float(E_a_sol_all[-1])
            E_c_sol_all = (E_a_sol_all - Rf_sol) / (1 - Rf_sol) if Rf_sol < 1 else np.zeros_like(E_a_sol_all)
            
            d50_sol = _interpolate_size(sizes, E_a_sol_all[:-1], 0.5)
            d50c_sol = _interpolate_size(sizes, E_c_sol_all[:-1], 0.5)
            d25c_sol = _interpolate_size(sizes, E_c_sol_all[:-1], 0.25)
            d75c_sol = _interpolate_size(sizes, E_c_sol_all[:-1], 0.75)

            solids_metrics = HydrocycloneMetrics(
                d50=float(d50_sol), d50c=float(d50c_sol), bypass_rf=float(Rf_sol * 100),
                imperfection=float((d75c_sol - d25c_sol) / (2 * d50c_sol)) if d50c_sol > 0 else 0,
                solids_recovery_s=float(S_solids * 100)
            )

            # Construir tabla por sólidos con pasantes
            f_pass_s = 100 - np.cumsum(p_f_sol_calc)
            o_pass_s = 100 - np.cumsum(p_o_norm)
            u_pass_s = 100 - np.cumsum(p_u_norm)
            f_pass_s[-1], o_pass_s[-1], u_pass_s[-1] = 0, 0, 0

            balance_solids_table = []
            for i in range(len(sizes)):
                balance_solids_table.append(BalanceRow(
                    size=f"{sizes[i]} µm",
                    feed_pct=float(p_f_sol_calc[i]), overflow_pct=float(p_o_norm[i]), underflow_pct=float(p_u_norm[i]),
                    feed_pass=float(f_pass_s[i]), overflow_pass=float(o_pass_s[i]), underflow_pass=float(u_pass_s[i]),
                    recovery_underflow=float(E_a_sol_all[i]), recovery_corrected=float(E_c_sol_all[i])
                ))
            balance_solids_table.append(BalanceRow(
                size="Fondo (Pan)",
                feed_pct=float(p_f_sol_calc[-1]), overflow_pct=float(p_o_norm[-1]), underflow_pct=float(p_u_norm[-1]),
                feed_pass=0.0, overflow_pass=0.0, underflow_pass=0.0,
                recovery_underflow=float(E_a_sol_all[-1]), recovery_corrected=0.0
            ))
            balance_solids_table.append(BalanceRow(
                size="TOTAL", feed_pct=100.0, overflow_pct=100.0, underflow_pct=100.0,
                recovery_underflow=float(S_solids), recovery_corrected=1.0
            ))

    # 6. Diagnóstico y Balance de Agua
    diag_msg, diag_level = "Datos procesados correctamente.", "success"
    Rw = None
    if request.feed_p_solids and request.overflow_p_solids and request.underflow_p_solids:
        cp, co, cu = request.feed_p_solids/100, request.overflow_p_solids/100, request.underflow_p_solids/100
        Rw = S_opt * ((1 - cu) / cu) * (cp / (1 - cp))
        err = abs(Rf_adj - Rw)
        if err < 0.02: diag_msg = "Excelente consistencia."
        elif err < 0.05: diag_msg, diag_level = "Consistencia aceptable.", "warning"
        else: diag_msg, diag_level = "Atención: Discrepancia significativa (>5%).", "error"

    # 7. Balance Global Absoluto
    global_bal, global_bal_sol = None, None
    if request.feed_flow_rate and request.feed_p_solids:
        cp = request.feed_p_solids / 100
        f_mass_s = request.feed_flow_rate if request.feed_flow_unit == "tph" else request.feed_flow_rate / ((1 / request.solid_density) + ((1 - cp) / cp))
        def get_flow(mass_s, p_sol):
            p_sol = max(p_sol, 0.001)
            mw = mass_s * (1 - p_sol) / p_sol
            vs, vw = mass_s / request.solid_density, mw / 1.0
            return FlowData(mass_solids=float(mass_s), mass_water=float(mw), vol_solids=float(vs), vol_water=float(vw), vol_pulp=float(vs+vw), p_solids=float(p_sol*100))
        
        global_bal = GlobalBalance(feed=get_flow(f_mass_s, cp), overflow=get_flow(f_mass_s * (1 - S_opt), co), underflow=get_flow(f_mass_s * S_opt, cu))
        if S_solids > 0:
            global_bal_sol = GlobalBalance(feed=get_flow(f_mass_s, cp), overflow=get_flow(f_mass_s * (1 - S_solids), co), underflow=get_flow(f_mass_s * S_solids, cu))

    # 8. Gráficos y Tabla Reconciliada Final
    f_pass_exp, o_pass_exp, u_pass_exp = 100 - np.cumsum(p_f_exp), 100 - np.cumsum(p_o_exp), 100 - np.cumsum(p_u_exp)
    f_pass_exp[-1], o_pass_exp[-1], u_pass_exp[-1] = 0, 0, 0
    f_pass_adj, o_pass_adj, u_pass_adj = 100 - np.cumsum(p_f_adj), 100 - np.cumsum(p_o_adj), 100 - np.cumsum(p_u_adj)
    f_pass_adj[-1], o_pass_adj[-1], u_pass_adj[-1] = 0, 0, 0

    balance_table = []
    for i in range(len(sizes)):
        balance_table.append(BalanceRow(
            size=f"{sizes[i]} µm", feed_w=float(f_raw[i]), overflow_w=float(o_raw[i]), underflow_w=float(u_raw[i]),
            feed_pct=float(p_f_exp[i]), overflow_pct=float(p_o_exp[i]), underflow_pct=float(p_u_exp[i]),
            feed_pass=float(f_pass_exp[i]), overflow_pass=float(o_pass_exp[i]), underflow_pass=float(u_pass_exp[i]),
            recovery_underflow=float(E_a_adj_all[i]), recovery_corrected=float(E_c_adj_all[i])
        ))
    balance_table.append(BalanceRow(size="Fondo (Pan)", feed_w=float(f_raw[-1]), overflow_w=float(o_raw[-1]), underflow_w=float(u_raw[-1]), feed_pct=float(p_f_exp[-1]), overflow_pct=float(p_o_exp[-1]), underflow_pct=float(p_u_exp[-1]), recovery_underflow=float(E_a_adj_all[-1]), recovery_corrected=0.0))
    balance_table.append(BalanceRow(size="TOTAL", feed_w=float(tot_f), overflow_w=float(tot_o), underflow_w=float(tot_u), feed_pct=100.0, overflow_pct=100.0, underflow_pct=100.0, recovery_underflow=float(S_opt), recovery_corrected=1.0))

    granulometry_pts = []
    for i in range(len(sizes)):
        g_pt = GranulometryPoint(size=float(sizes[i]), feed_passing=float(f_pass_exp[i]), overflow_passing=float(o_pass_exp[i]), underflow_passing=float(u_pass_exp[i]), feed_passing_adj=float(f_pass_adj[i]), overflow_passing_adj=float(o_pass_adj[i]), underflow_passing_adj=float(u_pass_adj[i]))
        if balance_solids_table: g_pt.feed_passing_sol = float(balance_solids_table[i].feed_pass)
        granulometry_pts.append(g_pt)

    partition_pts = []
    for i in range(len(sizes)):
        p_pt = PartitionCurvePoint(size=float(sizes[i]), actual_recovery=float(S_opt * p_u_exp[i] / p_f_exp[i]), corrected_recovery=float(E_c_adj_all[i]), adjusted_recovery=float(E_a_adj_all[i]))
        if solids_metrics: p_pt.solids_recovery = float(E_a_sol_all[i])
        partition_pts.append(p_pt)

    return HydrocycloneAnalysisResponse(
        reconciled_metrics=reconciled_metrics, solids_metrics=solids_metrics,
        d50c_experimental=float(d50_adj), d50c_adjusted=float(d50c_adj),
        tromp=TrompParameters(d25c=float(d25c_adj), d50c=float(d50c_adj), d75c=float(d75c_adj), imperfection=reconciled_metrics.imperfection),
        diagnosis_message=diag_msg, diagnosis_level=diag_level,
        water_balance=WaterBalance(solids_recovery_S=float(S_opt*100), water_recovery_Rw=float(Rw*100) if Rw else None, bypass_Rf=float(Rf_adj*100), pan_recovery_Epan=float(Rf_adj*100), feed_flow=100.0, overflow_flow=float((1-S_opt)*100), underflow_flow=float(S_opt*100), global_balance=global_bal, global_balance_solids=global_bal_sol),
        partition_curve=partition_pts, granulometry_curve=granulometry_pts, balance_table=balance_table, balance_solids_table=balance_solids_table, summary={"status": "Success"}
    )

def _empty_response():
    return HydrocycloneAnalysisResponse(d50c_experimental=0, d50c_adjusted=0, water_balance=WaterBalance(solids_recovery_S=0, bypass_Rf=0, feed_flow=100, overflow_flow=0, underflow_flow=0), partition_curve=[], granulometry_curve=[], balance_table=[], summary={})
