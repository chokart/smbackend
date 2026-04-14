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

def calculate_plitt_recovery(size: float, d50c: float, m: float) -> float:
    """Calcula la recuperación corregida usando el modelo de Plitt: E_c = 1 - exp(-0.693 * (d/d50c)^m)"""
    if d50c <= 0: return 0.0
    return 1.0 - np.exp(-0.693 * (size / d50c)**m)

def analyze_plitt(request: HydrocycloneAnalysisRequest, results: HydrocycloneAnalysisResponse, S_opt: float):
    """Integra el modelamiento de Plitt en los resultados existentes."""
    if not request.geometry:
        return
    
    g = request.geometry
    p = request.plitt_params or PlittParameters()
    
    # phi es el % de sólidos en volumen en la alimentación
    rho_s = request.solid_density
    rho_l = request.liquid_density
    cp = (request.feed_p_solids or 0) / 100
    phi = (cp / rho_s) / ((cp / rho_s) + ((1 - cp) / rho_l)) * 100

    # Determinar Q (Flujo de alimentación en l/min)
    # Si tenemos tph, convertir a l/min usando densidad de pulpa
    Q_lmin = 0
    if request.feed_flow_rate:
        if request.feed_flow_unit == "tph":
            # fm (tph) -> l/min
            rho_p = 1 / ((cp / rho_s) + ((1 - cp) / rho_l))
            Q_lmin = (request.feed_flow_rate * 1000 / 60) / rho_p
        else: # m3/h
            Q_lmin = request.feed_flow_rate * 1000 / 60
    
    if Q_lmin <= 0: return

    # 1. Ecuación de Presión (P en kPa)
    # P = [F2 * Q^1.78 * exp(0.0055 * phi)] / [Dc^0.37 * Di^0.94 * h^0.28 * (Du^2 + Do^2)^0.87]
    P_calc = (p.F2 * (Q_lmin**1.78) * np.exp(0.0055 * phi)) / \
             (g.Dc**0.37 * g.Di**0.94 * g.h**0.28 * (g.Du**2 + g.Do**2)**0.87)
    
    # 2. Ecuación de d50c (en micrones)
    # d50c = [F1 * Dc^0.46 * Di^0.6 * Do^1.21 * exp(0.063 * phi)] / [Du^0.71 * h^0.38 * Q^0.45 * (rho_s - rho_l)^0.5]
    d50c_calc = (p.F1 * g.Dc**0.46 * g.Di**0.6 * g.Do**1.21 * np.exp(0.063 * phi)) / \
                (g.Du**0.71 * g.h**0.38 * Q_lmin**0.45 * (rho_s - rho_l)**0.5)
    
    # 3. Ecuación de S (Relación de división de volumen pulpa S = Vu/Vo)
    # S = [F3 * (Du/Do)^3.31 * h^0.54 * (Du^2 + Do^2)^0.36 * exp(0.0054 * phi)] / [Dc^1.11 * P^0.24]
    S_calc = (p.F3 * (g.Du / g.Do)**3.31 * g.h**0.54 * (g.Du**2 + g.Do**2)**0.36 * np.exp(0.0054 * phi)) / \
             (g.Dc**1.11 * P_calc**0.24)
    
    # 4. Ecuación de m (Agudeza)
    # m = F4 * exp(-1.58 * Rv) * (Dc^2 * h / Q)^0.15 donde Rv = S/(1+S)
    Rv = S_calc / (1 + S_calc)
    m_calc = p.F4 * np.exp(-1.58 * Rv) * ((g.Dc**2 * g.h) / Q_lmin)**0.15

    # Actualizar métricas
    if results.reconciled_metrics:
        results.reconciled_metrics.m_plitt = float(m_calc)
        results.reconciled_metrics.s_plitt = float(S_calc)
        results.reconciled_metrics.p_plitt = float(P_calc)
    
    results.d50c_adjusted = float(d50c_calc)

    # Calcular curva de Plitt
    Rf = results.reconciled_metrics.bypass_rf / 100 if results.reconciled_metrics else 0
    for pt in results.partition_curve:
        # Recuperación Corregida de Plitt
        ec_plitt = calculate_plitt_recovery(pt.size, d50c_calc, m_calc)
        # Recuperación Real de Plitt (Ea = Ec*(1-Rf) + Rf)
        pt.plitt_recovery = float(ec_plitt * (1 - Rf) + Rf)

def analyze_hydrocyclone(request: HydrocycloneAnalysisRequest) -> HydrocycloneAnalysisResponse:
    sizes = np.array([s.mesh_size for s in request.sieves])
    if len(sizes) == 0: return _empty_response()
    sort_idx = np.argsort(sizes)[::-1]; sizes = sizes[sort_idx]
    
    # Pesos de alimento (siempre requeridos)
    w_f = np.array([request.sieves[i].weight_feed for i in sort_idx])
    f_raw = np.append(w_f, request.pan_feed)
    tot_f = max(np.sum(f_raw), 0.001)
    p_f_exp = (f_raw / tot_f) * 100

    # Inicializar variables de resultados
    S_opt = 0.5
    p_f_adj, p_o_adj, p_u_adj = p_f_exp, np.zeros_like(p_f_exp), np.zeros_like(p_f_exp)
    
    if request.mode == "simulation":
        # --- MODO SIMULACIÓN ---
        # 1. Calcular parámetros de Plitt primero
        g = request.geometry
        p = request.plitt_params or PlittParameters()
        if not g: raise ValueError("La geometría es obligatoria para el modo simulación")

        # Densidades y concentraciones
        rho_s = request.solid_density
        rho_l = request.liquid_density
        cp = (request.feed_p_solids or 30) / 100
        phi = (cp / rho_s) / ((cp / rho_s) + ((1 - cp) / rho_l)) * 100

        # Flujo Q
        Q_lmin = 100.0
        if request.feed_flow_rate:
            if request.feed_flow_unit == "tph":
                rho_p = 1 / ((cp / rho_s) + ((1 - cp) / rho_l))
                Q_lmin = (request.feed_flow_rate * 1000 / 60) / rho_p
            else: Q_lmin = request.feed_flow_rate * 1000 / 60
        
        # Ecuaciones de Plitt
        P_calc = (p.F2 * (Q_lmin**1.78) * np.exp(0.0055 * phi)) / (g.Dc**0.37 * g.Di**0.94 * g.h**0.28 * (g.Du**2 + g.Do**2)**0.87)
        d50c_calc = (p.F1 * g.Dc**0.46 * g.Di**0.6 * g.Do**1.21 * np.exp(0.063 * phi)) / (g.Du**0.71 * g.h**0.38 * Q_lmin**0.45 * (rho_s - rho_l)**0.5)
        S_v_calc = (p.F3 * (g.Du / g.Do)**3.31 * g.h**0.54 * (g.Du**2 + g.Do**2)**0.36 * np.exp(0.0054 * phi)) / (g.Dc**1.11 * P_calc**0.24)
        Rv = S_v_calc / (1 + S_v_calc) # Fracción de volumen al UF
        m_calc = p.F4 * np.exp(-1.58 * Rv) * ((g.Dc**2 * g.h) / Q_lmin)**0.15

        # Asumimos Rf (Bypass) = Rv (Aproximación de Plitt para simulación)
        Rf_sim = Rv
        
        # Calcular recuperaciones por malla
        for i, d in enumerate(sizes):
            ec = 1.0 - np.exp(-0.693 * (d / d50c_calc)**m_calc)
            ea = ec * (1 - Rf_sim) + Rf_sim
            p_u_adj[i] = p_f_exp[i] * ea
            p_o_adj[i] = p_f_exp[i] * (1 - ea)
        
        # Pan (Fondo)
        p_u_adj[-1] = p_f_exp[-1] * Rf_sim
        p_o_adj[-1] = p_f_exp[-1] * (1 - Rf_sim)

        S_opt = np.sum(p_u_adj) / 100.0
        p_u_adj = (p_u_adj / max(np.sum(p_u_adj), 0.001)) * 100
        p_o_adj = (p_o_adj / max(np.sum(p_o_adj), 0.001)) * 100
        E_a_adj = np.clip(S_opt * p_u_adj / np.where(p_f_adj > 0, p_f_adj, 0.001), 0, 1)
        Rf_adj = Rf_sim
        E_c_adj = np.clip((E_a_adj - Rf_adj) / (1 - Rf_adj) if Rf_adj < 1 else 0, 0, 1)

    else:
        # --- MODO ANÁLISIS (RECONCILIACIÓN) ---
        w_o = np.array([request.sieves[i].weight_overflow for i in sort_idx])
        w_u = np.array([request.sieves[i].weight_underflow for i in sort_idx])
        o_raw = np.append(w_o, request.pan_overflow); u_raw = np.append(w_u, request.pan_underflow)
        tot_o = max(np.sum(o_raw), 0.001); tot_u = max(np.sum(u_raw), 0.001)
        p_o_exp = (o_raw / tot_o) * 100; p_u_exp = (u_raw / tot_u) * 100

        def obj_s(S): return np.sum((p_f_exp - ((1 - S) * p_o_exp + S * p_u_exp))**2)
        S_opt = float(minimize(obj_s, 0.5, bounds=[(0.01, 0.99)]).x[0])

        p_f_adj, p_o_adj, p_u_adj = [], [], []
        for f_e, o_e, u_e in zip(p_f_exp, p_o_exp, p_u_exp):
            def obj_m(x): return (x[0]-f_e)**2 + (x[1]-o_e)**2 + (x[2]-u_e)**2
            cons = {'type': 'eq', 'fun': lambda x: x[0] - ((1-S_opt)*x[1] + S_opt*x[2])}
            res = minimize(obj_m, [f_e, o_e, u_e], constraints=cons, bounds=[(0,100),(0,100),(0,100)]).x
            p_f_adj.append(res[0]); p_o_adj.append(res[1]); p_u_adj.append(res[2])

        p_f_adj = (np.array(p_f_adj)/np.sum(p_f_adj))*100; p_o_adj = (np.array(p_o_adj)/np.sum(p_o_adj))*100; p_u_adj = (np.array(p_u_adj)/np.sum(p_u_adj))*100
        E_a_adj = np.clip(S_opt*p_u_adj/np.where(p_f_adj>0, p_f_adj, 0.001), 0, 1)
        Rf_adj = float(E_a_adj[-1]); E_c_adj = np.clip((E_a_adj - Rf_adj)/(1-Rf_adj) if Rf_adj<1 else 0, 0, 1)

    reconciled_metrics = HydrocycloneMetrics(
        d50=_interpolate_size(sizes, E_a_adj[:-1], 0.5), d50c=_interpolate_size(sizes, E_c_adj[:-1], 0.5),
        bypass_rf=Rf_adj*100, solids_recovery_s=S_opt*100
    )

    mass_basis = request.feed_flow_rate if (request.feed_flow_rate and request.feed_flow_unit == "tph") else tot_f

    # Balance % Sólidos
    balance_solids_table = None; solids_metrics = None; S_sol = 0; p_f_sol = None
    if request.feed_p_solids and request.overflow_p_solids and request.underflow_p_solids:
        cp, co, cu = request.feed_p_solids/100, request.overflow_p_solids/100, request.underflow_p_solids/100
        if (cu - co) > 0 and cp > 0:
            S_sol = max(0.01, min(0.99, (cu*(cp-co))/(cp*(cu-co))))
            p_f_sol = (1-S_sol)*p_o_exp + S_sol*p_u_exp
            E_a_sol = np.clip(S_sol*p_u_exp/np.where(p_f_sol>0, p_f_sol, 0.001), 0, 1); Rf_sol = float(E_a_sol[-1]); E_c_sol = np.clip((E_a_sol-Rf_sol)/(1-Rf_sol) if Rf_sol<1 else 0, 0, 1)
            solids_metrics = HydrocycloneMetrics(d50=_interpolate_size(sizes, E_a_sol[:-1], 0.5), d50c=_interpolate_size(sizes, E_c_sol[:-1], 0.5), bypass_rf=Rf_sol*100, solids_recovery_s=S_sol*100)
            fs, os, us = 100-np.cumsum(p_f_sol), 100-np.cumsum(p_o_exp), 100-np.cumsum(p_u_exp); fs[-1]=os[-1]=us[-1]=0
            w_f_s, w_o_s, w_u_s = mass_basis*(p_f_sol/100), mass_basis*(1-S_sol)*(p_o_exp/100), mass_basis*S_sol*(p_u_exp/100)
            balance_solids_table = [BalanceRow(size=f"{sizes[i]} µm", feed_pct=float(p_f_sol[i]), overflow_pct=float(p_o_exp[i]), underflow_pct=float(p_u_exp[i]), feed_w_sol=float(w_f_s[i]), overflow_w_sol=float(w_o_s[i]), underflow_w_sol=float(w_u_s[i]), feed_pass=float(fs[i]), overflow_pass=float(os[i]), underflow_pass=float(us[i]), recovery_underflow=float(E_a_sol[i]), recovery_corrected=float(E_c_sol[i])) for i in range(len(sizes))]
            balance_solids_table.append(BalanceRow(size="Fondo (Pan)", feed_pct=float(p_f_sol[-1]), overflow_pct=float(p_o_exp[-1]), underflow_pct=float(p_u_exp[-1]), feed_w_sol=float(w_f_s[-1]), overflow_w_sol=float(w_o_s[-1]), underflow_w_sol=float(w_u_s[-1]), recovery_underflow=Rf_sol, recovery_corrected=0))
            balance_solids_table.append(BalanceRow(size="TOTAL", feed_pct=100, overflow_pct=100, underflow_pct=100, feed_w_sol=mass_basis, overflow_w_sol=mass_basis*(1-S_sol), underflow_w_sol=mass_basis*S_sol, recovery_underflow=S_sol, recovery_corrected=1))

    # Global
    global_bal, global_bal_sol = None, None
    if request.feed_flow_rate and request.feed_p_solids:
        cp = request.feed_p_solids / 100
        fm = request.feed_flow_rate if request.feed_flow_unit == "tph" else request.feed_flow_rate / ((1/request.solid_density) + ((1-cp)/cp))
        def get_f(m, p): mw = m*(1-p/100)/(p/100) if p>0 else 0; vs, vw = m/request.solid_density, mw/1.0; return FlowData(mass_solids=m, mass_water=mw, vol_solids=vs, vol_water=vw, vol_pulp=vs+vw, p_solids=p)
        global_bal = GlobalBalance(feed=get_f(fm, cp*100), overflow=get_f(fm*(1-S_opt), (request.overflow_p_solids or cp*100)), underflow=get_f(fm*S_opt, (request.underflow_p_solids or cp*100)))
        if S_sol > 0: global_bal_sol = GlobalBalance(feed=get_f(fm, cp*100), overflow=get_f(fm*(1-S_sol), (request.overflow_p_solids or cp*100)), underflow=get_f(fm*S_sol, (request.underflow_p_solids or cp*100)))

    # Curvas
    fe, oe, ue = 100-np.cumsum(p_f_exp), 100-np.cumsum(p_o_exp), 100-np.cumsum(p_u_exp); fe[-1]=oe[-1]=ue[-1]=0
    fa, oa, ua = 100-np.cumsum(p_f_adj), 100-np.cumsum(p_o_adj), 100-np.cumsum(p_u_adj); fa[-1]=oa[-1]=ua[-1]=0
    w_f_a, w_o_a, w_u_a = mass_basis*(p_f_adj/100), mass_basis*(1-S_opt)*(p_o_adj/100), mass_basis*S_opt*(p_u_adj/100)
    
    balance_table = [BalanceRow(size=f"{sizes[i]} µm", feed_w=float(f_raw[i]), overflow_w=float(o_raw[i]), underflow_w=float(u_raw[i]), feed_w_adj=float(w_f_a[i]), overflow_w_adj=float(w_o_a[i]), underflow_w_adj=float(w_u_a[i]), feed_pct=float(p_f_exp[i]), overflow_pct=float(p_o_exp[i]), underflow_pct=float(p_u_exp[i]), feed_pass=float(fe[i]), overflow_pass=float(oe[i]), underflow_pass=float(ue[i]), feed_pct_adj=float(p_f_adj[i]), overflow_pct_adj=float(p_o_adj[i]), underflow_pct_adj=float(p_u_adj[i]), feed_pass_adj=float(fa[i]), overflow_pass_adj=float(oa[i]), underflow_pass_adj=float(ua[i]), recovery_underflow=float(E_a_adj[i]), recovery_corrected=float(E_c_adj[i])) for i in range(len(sizes))]
    balance_table.append(BalanceRow(size="Fondo (Pan)", feed_w=float(f_raw[-1]), overflow_w=float(o_raw[-1]), underflow_w=float(u_raw[-1]), feed_w_adj=float(w_f_a[-1]), overflow_w_adj=float(w_o_a[-1]), underflow_w_adj=float(w_u_a[-1]), feed_pct=float(p_f_exp[-1]), overflow_pct=float(p_o_exp[-1]), underflow_pct=float(p_u_exp[-1]), feed_pct_adj=float(p_f_adj[-1]), overflow_pct_adj=float(p_o_adj[-1]), underflow_pct_adj=float(p_u_adj[-1]), recovery_underflow=Rf_adj, recovery_corrected=0))
    balance_table.append(BalanceRow(size="TOTAL", feed_w=tot_f, overflow_w=tot_o, underflow_w=tot_u, feed_w_adj=mass_basis, overflow_w_adj=mass_basis*(1-S_opt), underflow_w_adj=mass_basis*S_opt, feed_pct=100, overflow_pct=100, underflow_pct=100, feed_pct_adj=100, overflow_pct_adj=100, underflow_pct_adj=100, recovery_underflow=S_opt, recovery_corrected=1))

    gran_pts = []
    for i in range(len(sizes)):
        pt = GranulometryPoint(size=float(sizes[i]), feed_passing=float(fe[i]), overflow_passing=float(oe[i]), underflow_passing=float(ue[i]), feed_passing_adj=float(fa[i]), overflow_passing_adj=float(oa[i]), underflow_passing_adj=float(ua[i]))
        if p_f_sol is not None: pt.feed_passing_sol = float(fs[i]); pt.overflow_passing_sol = float(os[i]); pt.underflow_passing_sol = float(us[i])
        gran_pts.append(pt)

    part_pts = [PartitionCurvePoint(size=float(sizes[i]), actual_recovery=float(S_opt*p_u_exp[i]/max(p_f_exp[i],0.1)), corrected_recovery=float(E_c_adj[i]), adjusted_recovery=float(E_a_adj[i])) for i in range(len(sizes))]
    if solids_metrics:
        for i, pt in enumerate(part_pts): pt.solids_recovery = float(E_a_sol[i])

    res = HydrocycloneAnalysisResponse(
        reconciled_metrics=reconciled_metrics, solids_metrics=solids_metrics, d50c_experimental=reconciled_metrics.d50, d50c_adjusted=reconciled_metrics.d50c,
        tromp=TrompParameters(d25c=_interpolate_size(sizes, E_c_adj[:-1], 0.25), d50c=reconciled_metrics.d50c, d75c=_interpolate_size(sizes, E_c_adj[:-1], 0.75)),
        water_balance=WaterBalance(solids_recovery_S=S_opt*100, bypass_Rf=Rf_adj*100, feed_flow=100, overflow_flow=(1-S_opt)*100, underflow_flow=S_opt*100, global_balance=global_bal, global_balance_solids=global_bal_sol),
        partition_curve=part_pts, granulometry_curve=gran_pts, balance_table=balance_table, balance_solids_table=balance_solids_table, summary={"status": "Success"}
    )
    
    # Integrar modelamiento de Plitt si hay geometría disponible
    analyze_plitt(request, res, S_opt)
    
    return res

def _empty_response():
    return HydrocycloneAnalysisResponse(d50c_experimental=0, d50c_adjusted=0, water_balance=WaterBalance(solids_recovery_S=0, bypass_Rf=0, feed_flow=100, overflow_flow=0, underflow_flow=0), partition_curve=[], granulometry_curve=[], balance_table=[], summary={})
