import numpy as np
from scipy.optimize import minimize
from typing import List, Tuple, Optional
from hydrocyclone_models import (
    HydrocycloneAnalysisRequest, HydrocycloneAnalysisResponse, 
    PartitionCurvePoint, GranulometryPoint, BalanceRow, WaterBalance,
    TrompParameters, GlobalBalance, FlowData, HydrocycloneMetrics,
    PlittGeometry, PlittParameters
)

def _interpolate_size(sizes, efficiencies, target_eff):
    if len(sizes) < 2: return 0.0
    try:
        log_sizes = np.log10(sizes)
        s_idx = np.argsort(efficiencies)
        eff_s = efficiencies[s_idx]; log_s_s = log_sizes[s_idx]
        if target_eff < np.min(eff_s) or target_eff > np.max(eff_s):
            if target_eff < np.min(eff_s): return float(np.min(sizes) * 0.8)
            return float(np.max(sizes) * 1.2)
        return float(10**np.interp(target_eff, eff_s, log_s_s))
    except: return 0.0

def _get_plitt_results(request: HydrocycloneAnalysisRequest, p_f_exp: np.ndarray, sizes: np.ndarray):
    """Calcula los parámetros y repartición según el modelo CIMM."""
    if not request.geometry: return None
    g = request.geometry
    p = request.plitt_params or PlittParameters()
    DC_in, h_in, DI_in, DO_in, DU_in = g.Dc/2.54, g.h/2.54, g.Di/2.54, g.Do/2.54, g.Du/2.54
    rho_s, rho_l = request.solid_density, request.liquid_density
    cp = (request.feed_p_solids or 30) / 100
    f = (cp/rho_s)/((cp/rho_s)+((1-cp)/rho_l))
    rho_p = 1/((cp/rho_s)+((1-cp)/rho_l))
    Q_m3h = request.feed_flow_rate if request.feed_flow_unit == "m3h" else (request.feed_flow_rate/rho_p/cp if cp>0 else 100)

    H_calc = p.a1 * (Q_m3h**1.46 * np.exp(-7.63*f + 10.79*f**2)) / (DC_in**0.2 * h_in**0.15 * DI_in**0.51 * DO_in**1.65 * DU_in**0.53)
    H_oper = (request.pressure/(rho_p*9.81))*3.28084 if request.pressure else H_calc
    d50c_calc = p.a2 * (DC_in**0.44 * DI_in**0.58 * DO_in**1.91 * np.exp(11.12*f)) / (DU_in**0.8 * h_in**0.37 * Q_m3h**0.44 * (rho_s-1)**0.5)
    S_v_calc = p.a3 * (h_in**0.19 * (DU_in/DO_in)**2.64 * np.exp(-4.33*f + 8.77*f**2)) / (H_oper**0.54 * DC_in**0.38)
    m_calc = np.exp(p.a4 - 1.58 * (S_v_calc/(S_v_calc+1))) * ((DC_in**2 * h_in)/Q_m3h)**0.15
    
    temp_ec = np.array([1.0 - np.exp(-0.693 * (d/d50c_calc)**m_calc) for d in sizes])
    Rs_c = np.sum((p_f_exp[:-1]/100) * temp_ec)
    Bpw = (S_v_calc/(S_v_calc+1) - f*Rs_c) / (1 - f*(1 - p.l_const*(1 - Rs_c)))
    Bpf = p.l_const * Bpw
    ea_sim = temp_ec * (1 - Bpf) + Bpf
    ea_sim = np.append(ea_sim, Bpf)
    
    return {"d50c": d50c_calc, "m": m_calc, "S": S_v_calc, "P_kPa": H_calc/3.28084*rho_p*9.81, "Bpf": Bpf, "ea": ea_sim}

def analyze_hydrocyclone(request: HydrocycloneAnalysisRequest) -> HydrocycloneAnalysisResponse:
    sizes = np.array([s.mesh_size for s in request.sieves])
    if len(sizes) == 0: return _empty_response()
    sort_idx = np.argsort(sizes)[::-1]; sizes = sizes[sort_idx]
    
    w_f = np.array([request.sieves[i].weight_feed for i in sort_idx])
    f_raw = np.append(w_f, request.pan_feed); tot_f = max(np.sum(f_raw), 0.001); p_f_exp = (f_raw / tot_f) * 100

    w_o = np.array([request.sieves[i].weight_overflow for i in sort_idx])
    w_u = np.array([request.sieves[i].weight_underflow for i in sort_idx])
    o_raw = np.append(w_o, request.pan_overflow); u_raw = np.append(w_u, request.pan_underflow)
    tot_o, tot_u = max(np.sum(o_raw), 0.001), max(np.sum(u_raw), 0.001)
    p_o_exp, p_u_exp = (o_raw / tot_o) * 100, (u_raw / tot_u) * 100

    plitt = _get_plitt_results(request, p_f_exp, sizes)

    # 1. BALANCE POR MALLAS (Optimización)
    def obj_s(S): return np.sum((p_f_exp - ((1 - S) * p_o_exp + S * p_u_exp))**2)
    S_mesh = float(minimize(obj_s, 0.5, bounds=[(0.01, 0.99)]).x[0])
    p_f_mesh, p_o_mesh, p_u_mesh = [], [], []
    for f_e, o_e, u_e in zip(p_f_exp, p_o_exp, p_u_exp):
        cons = {'type': 'eq', 'fun': lambda x: x[0] - ((1-S_mesh)*x[1] + S_mesh*x[2])}
        res = minimize(lambda x: (x[0]-f_e)**2 + (x[1]-o_e)**2 + (x[2]-u_e)**2, [f_e, o_e, u_e], constraints=cons, bounds=[(0,100),(0,100),(0,100)]).x
        p_f_mesh.append(res[0]); p_o_mesh.append(res[1]); p_u_mesh.append(res[2])
    p_f_mesh, p_o_mesh, p_u_mesh = np.array(p_f_mesh), np.array(p_o_mesh), np.array(p_u_mesh)
    ea_mesh = np.clip(S_mesh*p_u_mesh/np.where(p_f_mesh>0, p_f_mesh, 0.001), 0, 1)
    Rf_mesh = float(ea_mesh[-1])
    ec_mesh = np.clip((ea_mesh - Rf_mesh)/(1-Rf_mesh) if Rf_mesh<1 else 0, 0, 1)

    # 2. BALANCE POR SOLIDOS (Recalculando Alimento)
    cp_f = request.feed_p_solids / 100 if request.feed_p_solids else 0.3
    cp_o = request.overflow_p_solids / 100 if request.overflow_p_solids else cp_f * 0.8
    cp_u = request.underflow_p_solids / 100 if request.underflow_p_solids else cp_f * 1.5
    
    try:
        dil_f = (1 - cp_f) / cp_f if cp_f > 0 else 0
        dil_o = (1 - cp_o) / cp_o if cp_o > 0 else 0
        dil_u = (1 - cp_u) / cp_u if cp_u > 0 else 0
        S_solids = (dil_f - dil_o) / (dil_u - dil_o) if (dil_u - dil_o) != 0 else S_mesh
        S_solids = np.clip(S_solids, 0.01, 0.99)
    except:
        S_solids = S_mesh

    p_o_solids = p_o_exp
    p_u_solids = p_u_exp
    p_f_solids = (1 - S_solids) * p_o_solids + S_solids * p_u_solids
    ea_solids = np.clip(S_solids * p_u_solids / np.where(p_f_solids > 0, p_f_solids, 0.001), 0, 1)
    Rf_solids = float(ea_solids[-1])
    ec_solids = np.clip((ea_solids - Rf_solids)/(1-Rf_solids) if Rf_solids<1 else 0, 0, 1)

    # Simulación CIMM
    p_u_sim = (p_f_exp * plitt["ea"]) / (np.sum(p_f_exp * plitt["ea"])/100) if plitt else np.zeros_like(p_f_exp)
    p_o_sim = (p_f_exp * (1 - plitt["ea"])) / (1 - np.sum(p_f_exp * plitt["ea"])/100) if plitt else np.zeros_like(p_f_exp)

    comparison_table = []
    for i in range(len(sizes) + 1):
        idx = i if i < len(sizes) else -1
        sz_str = f"{sizes[i]} µm" if i < len(sizes) else "Fondo (Pan)"
        comparison_table.append(BalanceRow(
            size=sz_str,
            feed_pct_real=float(p_f_exp[idx]),
            overflow_pct_real=float(p_o_exp[idx]),
            underflow_pct_real=float(p_u_exp[idx]),
            
            feed_pct_mesh=float(p_f_mesh[idx]),
            overflow_pct_mesh=float(p_o_mesh[idx]),
            underflow_pct_mesh=float(p_u_mesh[idx]),
            recovery_mesh=float(ea_mesh[idx]),
            
            feed_pct_solids=float(p_f_solids[idx]),
            overflow_pct_solids=float(p_o_solids[idx]),
            underflow_pct_solids=float(p_u_solids[idx]),
            recovery_solids=float(ea_solids[idx]),
            
            feed_pct_sim=float(p_f_exp[idx]),
            overflow_pct_sim=float(p_o_sim[idx]),
            underflow_pct_sim=float(p_u_sim[idx]),
            recovery_sim=float(plitt["ea"][idx]) if plitt else 0
        ))

    metrics_mesh = HydrocycloneMetrics(d50=_interpolate_size(sizes, ea_mesh[:-1], 0.5), d50c=_interpolate_size(sizes, ec_mesh[:-1], 0.5), bypass_rf=Rf_mesh*100, solids_recovery_s=S_mesh*100)
    metrics_solids = HydrocycloneMetrics(d50=_interpolate_size(sizes, ea_solids[:-1], 0.5), d50c=_interpolate_size(sizes, ec_solids[:-1], 0.5), bypass_rf=Rf_solids*100, solids_recovery_s=S_solids*100)
    
    if plitt:
        metrics_mesh.m_plitt = metrics_solids.m_plitt = plitt["m"]
        metrics_mesh.s_plitt = metrics_solids.s_plitt = plitt["S"]
        metrics_mesh.p_plitt = metrics_solids.p_plitt = plitt["P_kPa"]

    part_pts = []
    for i in range(len(sizes)):
        part_pts.append(PartitionCurvePoint(
            size=float(sizes[i]),
            adjusted_recovery_mesh=float(ea_mesh[i]),
            adjusted_recovery_solids=float(ea_solids[i]),
            plitt_recovery=float(plitt["ea"][i]) if plitt else None
        ))

    fe, oe, ue = 100-np.cumsum(p_f_exp), 100-np.cumsum(p_o_exp), 100-np.cumsum(p_u_exp); fe[-1]=oe[-1]=ue[-1]=0
    gran_pts = [GranulometryPoint(size=float(sizes[i]), feed_passing=float(fe[i]), overflow_passing=float(oe[i]), underflow_passing=float(ue[i])) for i in range(len(sizes))]

    return HydrocycloneAnalysisResponse(
        metrics_mesh=metrics_mesh,
        metrics_solids=metrics_solids,
        d50c_adjusted=plitt["d50c"] if plitt else 0,
        tromp_mesh=TrompParameters(d25c=_interpolate_size(sizes, ec_mesh[:-1], 0.25), d50c=metrics_mesh.d50c, d75c=_interpolate_size(sizes, ec_mesh[:-1], 0.75)),
        water_balance=WaterBalance(solids_recovery_S=S_mesh*100, bypass_Rf=Rf_mesh*100, feed_flow=100, overflow_flow=(1-S_mesh)*100, underflow_flow=S_mesh*100),
        partition_curve=part_pts,
        granulometry_curve=gran_pts,
        comparison_table=comparison_table,
        summary={"status": "Success"}
    )

def calibrate_hydrocyclone(request: HydrocycloneAnalysisRequest) -> dict:
    res = analyze_hydrocyclone(request); m_exp = res.metrics_mesh # Calibramos usando el balance por mallas
    if not m_exp or not request.geometry: return {"error": "No hay datos suficientes"}
    g = request.geometry; rho_s, rho_l = request.solid_density, request.liquid_density
    cp = (request.feed_p_solids or 30)/100; f = (cp/rho_s)/((cp/rho_s)+((1-cp)/rho_l)); rho_p = 1/((cp/rho_s)+((1-cp)/rho_l))
    Q_m3h = request.feed_flow_rate if request.feed_flow_unit == "m3h" else (request.feed_flow_rate/rho_p/cp if cp>0 else 100)
    DC_in, h_in, DI_in, DO_in, DU_in = g.Dc/2.54, g.h/2.54, g.Di/2.54, g.Do/2.54, g.Du/2.54
    
    H_target = (request.pressure/(rho_p*9.81))*3.28084 if request.pressure else 10.0
    a1_new = H_target / ((Q_m3h**1.46 * np.exp(-7.63*f + 10.79*f**2)) / (DC_in**0.2 * h_in**0.15 * DI_in**0.51 * DO_in**1.65 * DU_in**0.53))
    a2_new = m_exp.d50c / ((DC_in**0.44 * DI_in**0.58 * DO_in**1.91 * np.exp(11.12*f)) / (DU_in**0.8 * h_in**0.37 * Q_m3h**0.44 * (rho_s-1)**0.5))
    S_target = m_exp.solids_recovery_s / (100.0 - m_exp.solids_recovery_s + 0.001)
    a3_new = S_target / ((h_in**0.19 * (DU_in / DO_in)**2.64 * np.exp(-4.33*f + 8.77*f**2)) / (H_target**0.54 * DC_in**0.38))
    
    def fit_m(m_test): return np.sum([(1.0-np.exp(-0.693*(pt.size/m_exp.d50c)**m_test) - np.clip((pt.adjusted_recovery_mesh - m_exp.bypass_rf/100)/(1-m_exp.bypass_rf/100),0,1))**2 for pt in res.partition_curve])
    m_real = float(minimize(fit_m, 1.5, bounds=[(0.5, 5.0)]).x[0])
    a4_new = np.log(max(m_real / (((DC_in**2 * h_in) / Q_m3h)**0.15), 0.001)) + 1.58 * (S_target / (S_target + 1))
    
    def fit_l(l_test):
        Rs_c = np.sum((np.array([r.feed_pct_mesh for r in res.comparison_table[:-1]])/100) * np.array([1.0-np.exp(-0.693*(float(r.size.split()[0])/m_exp.d50c)**m_real) for r in res.comparison_table[:-1]]))
        bpw = (S_target/(S_target+1)-f*Rs_c)/(1-f*(1-l_test*(1-Rs_c))); return (l_test*bpw - m_exp.bypass_rf/100)**2
    l_new = float(minimize(fit_l, 1.0, bounds=[(0.1, 2.0)]).x[0])
    return {"a1": round(a1_new, 3), "a2": round(a2_new, 3), "a3": round(a3_new, 3), "a4": round(a4_new, 3), "l_const": round(l_new, 3)}

def _empty_response():
    return HydrocycloneAnalysisResponse(water_balance=WaterBalance(solids_recovery_S=0, bypass_Rf=0, feed_flow=100, overflow_flow=0, underflow_flow=0), partition_curve=[], granulometry_curve=[], comparison_table=[], summary={})