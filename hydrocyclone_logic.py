import numpy as np
from scipy.optimize import minimize
from typing import List, Tuple
from hydrocyclone_models import (
    HydrocycloneAnalysisRequest, HydrocycloneAnalysisResponse, 
    PartitionCurvePoint, GranulometryPoint, BalanceRow, WaterBalance
)

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
    # Buscamos S que minimice el error en: f_i = (1-S)*o_i + S*u_i
    def obj_s(S):
        error = 0
        for f, o, u in zip(p_f_exp, p_o_exp, p_u_exp):
            balance_calc = (1 - S) * o + S * u
            error += (f - balance_calc)**2
        return error

    res_s = minimize(obj_s, 0.5, bounds=[(0.01, 0.99)])
    S_opt = res_s.x[0]

    # 3. Reconciliación de Porcentajes por Malla
    # Para cada malla, ajustamos f_i, o_i, u_i para que el balance sea exacto
    # min sum((f_i-f_exp)^2 + (o_i-o_exp)^2 + (u_i-u_exp)^2) s.t. f_i = (1-S)*o_i + S*u_i
    p_f_adj, p_o_adj, p_u_adj = [], [], []
    for f_e, o_e, u_e in zip(p_f_exp, p_o_exp, p_u_exp):
        def obj_mesh(x):
            f, o, u = x
            return (f - f_e)**2 + (o - o_e)**2 + (u - u_e)**2
        
        cons = {'type': 'eq', 'fun': lambda x: x[0] - ((1 - S_opt) * x[1] + S_opt * x[2])}
        res_m = minimize(obj_mesh, [f_e, o_e, u_e], constraints=cons, bounds=[(0, 100), (0, 100), (0, 100)])
        p_f_adj.append(res_m.x[0]), p_o_adj.append(res_m.x[1]), p_u_adj.append(res_m.x[2])

    p_f_adj, p_o_adj, p_u_adj = np.array(p_f_adj), np.array(p_o_adj), np.array(p_u_adj)
    
    # Normalizar para que sumen 100% (pequeños ajustes residuales)
    p_f_adj = (p_f_adj / np.sum(p_f_adj)) * 100
    p_o_adj = (p_o_adj / np.sum(p_o_adj)) * 100
    p_u_adj = (p_u_adj / np.sum(p_u_adj)) * 100

    # 4. Cálculo de Eficiencias y Curvas
    # Ea = S * (u_i / f_i)
    E_a_adj_all = np.clip(S_opt * p_u_adj / p_f_adj, 0, 1)
    
    # Bypass (Rf) -> Usamos la eficiencia de la fracción más fina (fondo)
    Rf = float(E_a_adj_all[-1])
    
    # Eficiencia Corregida (Ec)
    E_c_adj_all = (E_a_adj_all - Rf) / (1 - Rf) if Rf < 1 else np.zeros_like(E_a_adj_all)
    E_c_adj_all = np.clip(E_c_adj_all, 0, 1)

    # d50c experimental vs ajustado
    d50c_exp = _calculate_d50c(sizes, (S_opt * p_u_exp[:-1] / p_f_exp[:-1] - Rf)/(1-Rf))
    d50c_adj = _calculate_d50c(sizes, E_c_adj_all[:-1])

    # 5. Balance de Agua
    Rw = None
    if request.feed_p_solids and request.overflow_p_solids and request.underflow_p_solids:
        # F = O + U (Sólidos)
        # Fw = Ow + Uw (Agua) -> F(1-Cp)/Cp = O(1-Co)/Co + U(1-Cu)/Cu
        cp, co, cu = request.feed_p_solids/100, request.overflow_p_solids/100, request.underflow_p_solids/100
        # Rw = Uw / Fw
        # Usando S = U/F -> U = S*F, O = (1-S)*F
        # Rw = [S*F*(1-cu)/cu] / [F*(1-cp)/cp] = S * (1-cu)/cu * cp/(1-cp)
        Rw = S_opt * ((1 - cu) / cu) * (cp / (1 - cp))
        # Si Rw es inconsistente, podríamos usar Rf como estimador o viceversa

    water_bal = WaterBalance(
        solids_recovery_S=float(S_opt * 100),
        water_recovery_Rw=float(Rw * 100) if Rw is not None else None,
        bypass_Rf=float(Rf * 100),
        feed_flow=100.0,
        overflow_flow=float((1 - S_opt) * 100),
        underflow_flow=float(S_opt * 100)
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
        water_balance=water_bal,
        partition_curve=partition_pts,
        granulometry_curve=granulometry_pts,
        balance_table=balance_table,
        summary={"iterations": 1, "status": "Success"}
    )

def _calculate_d50c(sizes, efficiencies):
    """Interpolación logarítmica para encontrar el tamaño donde la eficiencia es 0.5."""
    try:
        log_sizes = np.log10(sizes)
        # Asegurar que esté ordenado para interp
        s_idx = np.argsort(efficiencies)
        eff_s = efficiencies[s_idx]
        log_s_s = log_sizes[s_idx]
        
        if np.max(eff_s) >= 0.5 >= np.min(eff_s):
            log_d50c = np.interp(0.5, eff_s, log_s_s)
            return 10**log_d50c
        return 0.0
    except:
        return 0.0

def _empty_response():
    # Implementar según sea necesario para manejar errores
    pass
