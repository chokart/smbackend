import numpy as np
from typing import List, Tuple
from hydrocyclone_models import (
    RaoLynchRequest, RaoLynchResponse, PartitionCurvePoint, 
    GranulometryPoint, BalanceRow
)

def calculate_rao_lynch(request: RaoLynchRequest) -> RaoLynchResponse:
    # 1. Procesar granulometría experimental
    sizes = [s.mesh_size for s in request.sieves]
    w_f = np.array([s.weight_feed for s in request.sieves])
    w_o = np.array([s.weight_overflow for s in request.sieves])
    w_u = np.array([s.weight_underflow for s in request.sieves])
    
    # Agregar pesos del FONDO (PAN)
    w_f_all = np.append(w_f, request.pan_feed)
    w_o_all = np.append(w_o, request.pan_overflow)
    w_u_all = np.append(w_u, request.pan_underflow)
    
    # Totales y Porcentajes Retenidos Parciales
    tot_f = np.sum(w_f_all)
    tot_o = np.sum(w_o_all)
    tot_u = np.sum(w_u_all)
    
    pr_f_all = (w_f_all / tot_f) * 100 if tot_f > 0 else np.zeros_like(w_f_all)
    pr_o_all = (w_o_all / tot_o) * 100 if tot_o > 0 else np.zeros_like(w_o_all)
    pr_u_all = (w_u_all / tot_u) * 100 if tot_u > 0 else np.zeros_like(w_u_all)
    
    # 2. Balance de masa y Recuperación de sólidos (S)
    S_vals = []
    for f, o, u in zip(pr_f_all, pr_o_all, pr_u_all):
        if abs(u - o) > 0.1:
            S_est = (f - o) / (u - o)
            if 0 < S_est < 1: S_vals.append(S_est)
    S_solids = np.mean(S_vals) if S_vals else (tot_u / tot_f if tot_f > 0 else 0.5)
    
    # 3. Curva de Partición (Ea, Ec)
    E_a_all = []
    for f, u in zip(pr_f_all, pr_u_all):
        if f > 0:
            E_a_all.append(np.clip(S_solids * (u / f), 0, 1))
        else:
            E_a_all.append(0)
            
    Rf = float(E_a_all[-1]) # Bypass basado en el fondo
    E_a_sieves = np.array(E_a_all[:-1])
    E_c_sieves = (E_a_sieves - Rf) / (1 - Rf) if Rf < 1 else np.zeros_like(E_a_sieves)
    E_c_sieves = np.clip(E_c_sieves, 0, 1)
    
    # d50c Experimental (Interpolación Logarítmica)
    log_sizes = np.log10(sizes)
    try:
        if np.max(E_c_sieves) >= 0.5 >= np.min(E_c_sieves):
            log_d50c = np.interp(0.5, E_c_sieves[::-1], log_sizes[::-1])
            d50c_exp = 10**log_d50c
        else:
            d50c_exp = 0.0
    except:
        d50c_exp = 0.0
        
    # 4. Distribución Granulométrica (Pasante Acumulado)
    # Calculamos acumulados pasantes
    f_cum_pass = 100 - np.cumsum(pr_f_all)
    o_cum_pass = 100 - np.cumsum(pr_o_all)
    u_cum_pass = 100 - np.cumsum(pr_u_all)
    
    granulometry_points = []
    for i in range(len(sizes)):
        granulometry_points.append(GranulometryPoint(
            size=sizes[i],
            feed_passing=float(f_cum_pass[i]),
            overflow_passing=float(o_cum_pass[i]),
            underflow_passing=float(u_cum_pass[i])
        ))
        
    # 5. Tabla de Balance
    balance_rows = []
    for i in range(len(sizes)):
        balance_rows.append(BalanceRow(
            size=f"{sizes[i]} µm",
            feed_pct=float(pr_f_all[i]),
            overflow_pct=float(pr_o_all[i]),
            underflow_pct=float(pr_u_all[i]),
            recovery_underflow=float(E_a_all[i])
        ))
    # Agregar fila de fondo
    balance_rows.append(BalanceRow(
        size="Fondo (Pan)",
        feed_pct=float(pr_f_all[-1]),
        overflow_pct=float(pr_o_all[-1]),
        underflow_pct=float(pr_u_all[-1]),
        recovery_underflow=float(E_a_all[-1])
    ))
    
    # 6. CALIBRACIÓN Y PREDICCIÓN (RAO & LYNCH)
    # Factores geométricos fijos
    geom_factor_Q = (request.geometry.Do**0.73) * (request.geometry.Di**0.53) * (request.pressure**0.42)
    rho_s = request.solid_density
    rho_l = 1.0
    Cv = (request.feed_p_solids / rho_s) / (request.feed_p_solids / rho_s + (100 - request.feed_p_solids) / rho_l) * 100
    
    geom_factor_d50c = (request.geometry.Do**0.52) * (request.geometry.Di**-0.5) * (request.geometry.Du**-0.2) * (request.pressure**-0.3) * np.exp(Cv * 0.063)
    
    # Calibración de K1 (si hay caudal experimental)
    if request.exp_capacity_Q:
        k1_calc = request.exp_capacity_Q / geom_factor_Q
    else:
        k1_calc = 13.5 # Estándar
        
    # Calibración de K3 (si hay d50c experimental)
    if d50c_exp > 0:
        k3_calc = d50c_exp / geom_factor_d50c
    else:
        k3_calc = 35.0 # Estándar
        
    # Predicción (usando las constantes calibradas o las enviadas manualmente)
    k1_to_use = request.k1 or k1_calc
    k3_to_use = request.k3 or k3_calc
    
    Q_pred = k1_to_use * geom_factor_Q
    d50c_pred = k3_to_use * geom_factor_d50c
    
    # 7. Preparar Respuesta
    partition_points = []
    for i in range(len(sizes)):
        partition_points.append(PartitionCurvePoint(
            size=sizes[i],
            actual_recovery=float(E_a_all[i]),
            corrected_recovery=float(E_c_sieves[i])
        ))
        
    return RaoLynchResponse(
        d50c_experimental=float(d50c_exp),
        d50c_predicted=float(d50c_pred),
        water_bypass_Rf=float(Rf * 100),
        capacity_Q=float(Q_pred),
        k1_calculated=float(k1_calc),
        k3_calculated=float(k3_calc),
        partition_curve=partition_points,
        granulometry_curve=granulometry_points,
        balance_table=balance_rows,
        summary={
            "S_solids_recovery": float(S_solids * 100),
            "feed_total": float(tot_f),
            "overflow_total": float(tot_o),
            "underflow_total": float(tot_u),
            "vol_concentration_Cv": float(Cv)
        }
    )
