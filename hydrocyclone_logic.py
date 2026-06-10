import numpy as np
from scipy.optimize import minimize
from typing import List, Tuple, Optional
from hydrocyclone_models import (
    HydrocycloneAnalysisRequest, HydrocycloneAnalysisResponse, 
    PartitionCurvePoint, GranulometryPoint, BalanceRow, WaterBalance,
    TrompParameters, GlobalBalance, FlowData, HydrocycloneMetrics,
    FlowComparison, GlobalFlowBalance
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
            recovery_solids=float(ea_solids[idx])
        ))

    # Flow Calculations
    rho_s = request.solid_density
    rho_l = request.liquid_density
    
    if request.feed_flow_unit == "tph":
        Ms_f = request.feed_flow_rate if request.feed_flow_rate else 100.0
    else: # m3h
        Vp_f = request.feed_flow_rate if request.feed_flow_rate else 100.0
        factor = (1/rho_s + (1-cp_f)/(cp_f*rho_l)) if cp_f > 0 else 1.0
        Ms_f = Vp_f / factor if factor > 0 else 0

    def calc_flows(Ms, cp):
        if cp <= 0: return 0.0, 0.0, 0.0
        Mw = Ms * (1 - cp) / cp
        Mp = Ms + Mw
        return Mp, Ms, Mw

    Mp_f, _, Mw_f = calc_flows(Ms_f, cp_f)

    # Mesh Flow Balance
    Ms_u_m = Ms_f * S_mesh
    Ms_o_m = Ms_f * (1 - S_mesh)
    Mp_u_m, _, Mw_u_m = calc_flows(Ms_u_m, cp_u)
    Mp_o_m, _, Mw_o_m = calc_flows(Ms_o_m, cp_o)
    
    # Cálculo detallado de flujos globales usando S_solids para consistencia física
    def get_flow_data(Ms, cp, rho_s, rho_l):
        if cp <= 0: return FlowData(mass_solids=0, mass_water=0, vol_solids=0, vol_water=0, vol_pulp=0, p_solids=0)
        Mw = Ms * (1 - cp) / cp
        Vs = Ms / rho_s
        Vw = Mw / rho_l
        Vp = Vs + Vw
        return FlowData(
            mass_solids=float(Ms),
            mass_water=float(Mw),
            vol_solids=float(Vs),
            vol_water=float(Vw),
            vol_pulp=float(Vp),
            p_solids=float(cp * 100)
        )

    # Balance global (Consistente con % Sólidos)
    Ms_u_s = Ms_f * S_solids
    Ms_o_s = Ms_f * (1 - S_solids)
    
    global_balance = GlobalBalance(
        feed=get_flow_data(Ms_f, cp_f, rho_s, rho_l),
        overflow=get_flow_data(Ms_o_s, cp_o, rho_s, rho_l),
        underflow=get_flow_data(Ms_u_s, cp_u, rho_s, rho_l)
    )

    global_flow_balance = GlobalFlowBalance(
        feed=FlowComparison(
            pulp_mesh=Mp_f, solids_mesh=Ms_f, water_mesh=Mw_f,
            pulp_solids=global_balance.feed.vol_pulp, solids_solids=Ms_f, water_solids=global_balance.feed.mass_water
        ),
        overflow=FlowComparison(
            pulp_mesh=Mp_o_m, solids_mesh=Ms_o_m, water_mesh=Mw_o_m,
            pulp_solids=global_balance.overflow.vol_pulp, solids_solids=Ms_o_s, water_solids=global_balance.overflow.mass_water
        ),
        underflow=FlowComparison(
            pulp_mesh=Mp_u_m, solids_mesh=Ms_u_m, water_mesh=Mw_u_m,
            pulp_solids=global_balance.underflow.vol_pulp, solids_solids=Ms_u_s, water_solids=global_balance.underflow.mass_water
        )
    )

    metrics_mesh = HydrocycloneMetrics(d50=_interpolate_size(sizes, ea_mesh[:-1], 0.5), d50c=_interpolate_size(sizes, ec_mesh[:-1], 0.5), bypass_rf=Rf_mesh*100, solids_recovery_s=S_mesh*100)
    metrics_solids = HydrocycloneMetrics(d50=_interpolate_size(sizes, ea_solids[:-1], 0.5), d50c=_interpolate_size(sizes, ec_solids[:-1], 0.5), bypass_rf=Rf_solids*100, solids_recovery_s=S_solids*100)

    part_pts = []
    for i in range(len(sizes)):
        part_pts.append(PartitionCurvePoint(
            size=float(sizes[i]),
            adjusted_recovery_mesh=float(ea_mesh[i]),
            adjusted_recovery_solids=float(ea_solids[i])
        ))

    fe, oe, ue = 100-np.cumsum(p_f_exp), 100-np.cumsum(p_o_exp), 100-np.cumsum(p_u_exp); fe[-1]=oe[-1]=ue[-1]=0
    gran_pts = [GranulometryPoint(size=float(sizes[i]), feed_passing=float(fe[i]), overflow_passing=float(oe[i]), underflow_passing=float(ue[i])) for i in range(len(sizes))]

    return HydrocycloneAnalysisResponse(
        metrics_mesh=metrics_mesh,
        metrics_solids=metrics_solids,
        tromp_mesh=TrompParameters(d25c=_interpolate_size(sizes, ec_mesh[:-1], 0.25), d50c=metrics_mesh.d50c, d75c=_interpolate_size(sizes, ec_mesh[:-1], 0.75)),
        water_balance=WaterBalance(
            solids_recovery_S=S_solids*100, 
            bypass_Rf=Rf_solids*100, 
            feed_flow=Ms_f, 
            overflow_flow=Ms_o_s, 
            underflow_flow=Ms_u_s,
            global_balance=global_balance
        ),
        partition_curve=part_pts,
        granulometry_curve=gran_pts,
        comparison_table=comparison_table,
        global_flow_balance=global_flow_balance,
        summary={"status": "Success"}
    )

def _empty_response():
    return HydrocycloneAnalysisResponse(water_balance=WaterBalance(solids_recovery_S=0, bypass_Rf=0, feed_flow=100, overflow_flow=0, underflow_flow=0), partition_curve=[], granulometry_curve=[], comparison_table=[], summary={})