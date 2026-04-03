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
    if len(sizes) < 2: return 0.0
    try:
        log_sizes = np.log10(sizes)
        # Ordenar por eficiencia ascendente para np.interp
        s_idx = np.argsort(efficiencies)
        eff_s = efficiencies[s_idx]
        log_s_s = log_sizes[s_idx]
        
        if target_eff < np.min(eff_s) or target_eff > np.max(eff_s):
            # Extrapolación simple si está cerca de los bordes
            if target_eff < np.min(eff_s): return float(np.min(sizes) * 0.8)
            return float(np.max(sizes) * 1.2)
            
        return float(10**np.interp(target_eff, eff_s, log_s_s))
    except:
        return 0.0

def analyze_hydrocyclone(request: HydrocycloneAnalysisRequest) -> HydrocycloneAnalysisResponse:
    # 1. Preparación de Datos Experimentales
    sizes = np.array([s.mesh_size for s in request.sieves])
    if len(sizes) == 0: return _empty_response()
    
    sort_idx = np.argsort(sizes)[::-1]
    sizes = sizes[sort_idx]
    
    w_f = np.array([request.sieves[i].weight_feed for i in sort_idx])
    w_o = np.array([request.sieves[i].weight_overflow for i in sort_idx])
    w_u = np.array([request.sieves[i].weight_underflow for i in sort_idx])
    
    f_raw = np.append(w_f, request.pan_feed)
    o_raw = np.append(w_o, request.pan_overflow)
    u_raw = np.append(w_u, request.pan_underflow)
    
    tot_f = max(np.sum(f_raw), 0.001)
    tot_o = max(np.sum(o_raw), 0.001)
    tot_u = max(np.sum(u_raw), 0.001)

    p_f_exp = (f_raw / tot_f) * 100
    p_o_exp = (o_raw / tot_o) * 100
    p_u_exp = (u_raw / tot_u) * 100

    # 2. Optimización de S
    def obj_s(S):
        return np.sum((p_f_exp - ((1 - S) * p_o_exp + S * p_u_exp))**2)

    res_s = minimize(obj_s, 0.5, bounds=[(0.01, 0.99)])
    S_opt = float(res_s.x[0])

    # 3. Reconciliación
    p_f_adj, p_o_adj, p_u_adj = [], [], []
    for f_e, o_e, u_e in zip(p_f_exp, p_o_exp, p_u_exp):
        def obj_mesh(x):
            return (x[0] - f_e)**2 + (x[1] - o_e)**2 + (x[2] - u_e)**2
        cons = {'type': 'eq', 'fun': lambda x: x[0] - ((1 - S_opt) * x[1] + S_opt * x[2])}
        res_m = minimize(obj_mesh, [f_e, o_e, u_e], constraints=cons, bounds=[(0, 100), (0, 100), (0, 100)])
        p_f_adj.append(res_m.x[0]); p_o_adj.append(res_m.x[1]); p_u_adj.append(res_m.x[2])

    p_f_adj = (np.array(p_f_adj) / np.sum(p_f_adj)) * 100
    p_o_adj = (np.array(p_o_adj) / np.sum(p_o_adj)) * 100
    p_u_adj = (np.array(p_u_adj) / np.sum(p_u_adj)) * 100

    # 4. Métricas Reconciliado
    E_a_adj = np.clip(S_opt * p_u_adj / np.where(p_f_adj > 0, p_f_adj, 0.001), 0, 1)
    Rf_adj = float(E_a_adj[-1])
    E_c_adj = np.clip((E_a_adj - Rf_adj) / (1 - Rf_adj) if Rf_adj < 1 else 0, 0, 1)

    d50_adj = _interpolate_size(sizes, E_a_adj[:-1], 0.5)
    d50c_adj = _interpolate_size(sizes, E_c_adj[:-1], 0.5)
    d25c_adj = _interpolate_size(sizes, E_c_adj[:-1], 0.25)
    d75c_adj = _interpolate_size(sizes, E_c_adj[:-1], 0.75)
    
    reconciled_metrics = HydrocycloneMetrics(
        d50=d50_adj, d50c=d50c_adj, bypass_rf=Rf_adj * 100,
        imperfection=(d75c_adj - d25c_adj) / (2 * d50c_adj) if d50c_adj > 0 else 0,
        solids_recovery_s=S_opt * 100
    )

    # 5. Balance % Sólidos
    balance_solids_table = None
    solids_metrics = None
    S_sol = 0
    if request.feed_p_solids and request.overflow_p_solids and request.underflow_p_solids:
        cp, co, cu = request.feed_p_solids/100, request.overflow_p_solids/100, request.underflow_p_solids/100
        if (cu - co) > 0 and cp > 0:
            S_sol = max(0.01, min(0.99, (cu * (cp - co)) / (cp * (cu - co))))
            p_f_sol = (1 - S_sol) * p_o_exp + S_sol * p_u_exp
            E_a_sol = np.clip(S_sol * p_u_exp / np.where(p_f_sol > 0, p_f_sol, 0.001), 0, 1)
            Rf_sol = float(E_a_sol[-1])
            E_c_sol = np.clip((E_a_sol - Rf_sol) / (1 - Rf_sol) if Rf_sol < 1 else 0, 0, 1)

            solids_metrics = HydrocycloneMetrics(
                d50=_interpolate_size(sizes, E_a_sol[:-1], 0.5),
                d50c=_interpolate_size(sizes, E_c_sol[:-1], 0.5),
                bypass_rf=Rf_sol * 100,
                imperfection=0, # Simplificado
                solids_recovery_s=S_sol * 100
            )
            
            f_pass_s = 100 - np.cumsum(p_f_sol); f_pass_s[-1] = 0
            o_pass_s = 100 - np.cumsum(p_o_exp); o_pass_s[-1] = 0
            u_pass_s = 100 - np.cumsum(p_u_exp); u_pass_s[-1] = 0

            balance_solids_table = [BalanceRow(
                size=f"{sizes[i]} µm", feed_pct=float(p_f_sol[i]), overflow_pct=float(p_o_exp[i]), underflow_pct=float(p_u_exp[i]),
                feed_pass=float(f_pass_s[i]), overflow_pass=float(o_pass_s[i]), underflow_pass=float(u_pass_s[i]),
                recovery_underflow=float(E_a_sol[i]), recovery_corrected=float(E_c_sol[i])
            ) for i in range(len(sizes))]
            balance_solids_table.append(BalanceRow(size="Fondo (Pan)", feed_pct=float(p_f_sol[-1]), overflow_pct=float(p_o_exp[-1]), underflow_pct=float(p_u_exp[-1]), recovery_underflow=Rf_sol, recovery_corrected=0))
            balance_solids_table.append(BalanceRow(size="TOTAL", feed_pct=100, overflow_pct=100, underflow_pct=100, recovery_underflow=S_sol, recovery_corrected=1))

    # 6. Global y Respuesta
    global_bal, global_bal_sol = None, None
    if request.feed_flow_rate and request.feed_p_solids:
        cp = request.feed_p_solids / 100
        fm = request.feed_flow_rate if request.feed_flow_unit == "tph" else request.feed_flow_rate / ((1/request.solid_density) + ((1-cp)/cp))
        def get_f(m, p):
            mw = m*(1-p/100)/(p/100) if p>0 else 0
            vs, vw = m/request.solid_density, mw/1.0
            return FlowData(mass_solids=m, mass_water=mw, vol_solids=vs, vol_water=vw, vol_pulp=vs+vw, p_solids=p)
        global_bal = GlobalBalance(feed=get_f(fm, cp*100), overflow=get_f(fm*(1-S_opt), (request.overflow_p_solids or cp*100)), underflow=get_f(fm*S_opt, (request.underflow_p_solids or cp*100)))
        if S_sol > 0: global_bal_sol = GlobalBalance(feed=get_f(fm, cp*100), overflow=get_f(fm*(1-S_sol), (request.overflow_p_solids or cp*100)), underflow=get_f(fm*S_sol, (request.underflow_p_solids or cp*100)))

    f_pass_e = 100 - np.cumsum(p_f_exp); f_pass_e[-1] = 0
    o_pass_e = 100 - np.cumsum(p_o_exp); o_pass_e[-1] = 0
    u_pass_e = 100 - np.cumsum(p_u_exp); u_pass_e[-1] = 0

    balance_table = [BalanceRow(
        size=f"{sizes[i]} µm", feed_w=float(f_raw[i]), overflow_w=float(o_raw[i]), underflow_w=float(u_raw[i]),
        feed_pct=float(p_f_exp[i]), overflow_pct=float(p_o_exp[i]), underflow_pct=float(p_u_exp[i]),
        feed_pass=float(f_pass_e[i]), overflow_pass=float(o_pass_e[i]), underflow_pass=float(u_pass_e[i]),
        recovery_underflow=float(E_a_adj[i]), recovery_corrected=float(E_c_adj[i])
    ) for i in range(len(sizes))]
    balance_table.append(BalanceRow(size="Fondo (Pan)", feed_w=float(f_raw[-1]), overflow_w=float(o_raw[-1]), underflow_w=float(u_raw[-1]), feed_pct=float(p_f_exp[-1]), overflow_pct=float(p_o_exp[-1]), underflow_pct=float(p_u_exp[-1]), recovery_underflow=Rf_adj, recovery_corrected=0))
    balance_table.append(BalanceRow(size="TOTAL", feed_w=tot_f, overflow_w=tot_o, underflow_w=tot_u, feed_pct=100, overflow_pct=100, underflow_pct=100, recovery_underflow=S_opt, recovery_corrected=1))

    return HydrocycloneAnalysisResponse(
        reconciled_metrics=reconciled_metrics, solids_metrics=solids_metrics,
        d50c_experimental=d50_adj, d50c_adjusted=d50c_adj,
        tromp=TrompParameters(d25c=d25c_adj, d50c=d50c_adj, d75c=d75c_adj, imperfection=reconciled_metrics.imperfection),
        diagnosis_message="Cálculo exitoso", diagnosis_level="success",
        water_balance=WaterBalance(solids_recovery_S=S_opt*100, bypass_Rf=Rf_adj*100, feed_flow=100, overflow_flow=(1-S_opt)*100, underflow_flow=S_opt*100, global_balance=global_bal, global_balance_solids=global_bal_sol),
        partition_curve=[PartitionCurvePoint(size=float(sizes[i]), actual_recovery=float(S_opt*p_u_exp[i]/max(p_f_exp[i],0.1)), corrected_recovery=float(E_c_adj[i]), adjusted_recovery=float(E_a_adj[i])) for i in range(len(sizes))],
        granulometry_curve=[GranulometryPoint(size=float(sizes[i]), feed_passing=float(f_pass_e[i]), overflow_passing=float(o_pass_e[i]), underflow_passing=float(u_pass_e[i])) for i in range(len(sizes))],
        balance_table=balance_table, balance_solids_table=balance_solids_table, summary={"status": "Success"}
    )

def _empty_response():
    return HydrocycloneAnalysisResponse(d50c_experimental=0, d50c_adjusted=0, water_balance=WaterBalance(solids_recovery_S=0, bypass_Rf=0, feed_flow=100, overflow_flow=0, underflow_flow=0), partition_curve=[], granulometry_curve=[], balance_table=[], summary={})
