import math
import numpy as np
from typing import List, Tuple, Dict, Any, Optional
from scipy.optimize import curve_fit

from cyclone_models import (
    CycloneBalanceRequest, CycloneBalanceResponse, ReconciledSieveRow, StreamProperties,
    CycloneGeometry, CycloneCalibrationRequest, CycloneCalibrationResponse,
    CycloneSimulationRequest, CycloneSimulationResponse
)

# ==============================================================================
# UTILIDADES MATEMÁTICAS Y DE DISTRIBUCIÓN GRANULOMÉTRICA
# ==============================================================================

def calculate_percentile(sizes: List[float], passings: List[float], target_pct: float = 80.0) -> float:
    """
    Calcula el tamaño característico (ej. P80, P50) mediante interpolación log-lineal estándar.
    sizes: lista de tamaños en micrones
    passings: lista de porcentajes pasantes (0 a 100)
    """
    if not sizes or not passings:
        return 0.0
    
    # Ordenar por tamaño ascendente
    sorted_pairs = sorted(zip(sizes, passings), key=lambda x: x[0])
    s_arr = [p[0] for p in sorted_pairs if p[0] > 0]
    p_arr = [p[1] for p in sorted_pairs if p[0] > 0]
    
    if target_pct <= p_arr[0]:
        return float(s_arr[0])
    if target_pct >= p_arr[-1]:
        return float(s_arr[-1])
        
    for i in range(len(p_arr) - 1):
        p1, p2 = p_arr[i], p_arr[i+1]
        s1, s2 = s_arr[i], s_arr[i+1]
        if (p1 <= target_pct <= p2) or (p2 <= target_pct <= p1):
            if abs(p2 - p1) < 1e-6:
                return float(s1)
            # Interpolación log-lineal: log(size) vs pasante
            frac = (target_pct - p1) / (p2 - p1)
            log_s = math.log10(s1) + frac * (math.log10(s2) - math.log10(s1))
            return round(10.0 ** log_s, 2)
            
    return float(s_arr[-1])


# ==============================================================================
# SUBMÓDULO 1: BALANCE DE MATERIA Y RECONCILIACIÓN GRANULOMÉTRICA
# ==============================================================================

def solve_cyclone_balance(req: CycloneBalanceRequest) -> CycloneBalanceResponse:
    """
    Resuelve la reconciliación granulométrica de dos productos (Feed = Over + Under)
    usando Mínimos Cuadrados Ponderados (WLS) y ajusta la curva de partición Tromp.
    """
    # 1. Ordenar mallas de mayor a menor tamaño (coarsest to finest)
    raw_sieves = sorted(req.sieves, key=lambda s: s.size_microns, reverse=True)
    n = len(raw_sieves)
    if n < 3:
        raise ValueError("Se requieren al menos 3 mallas para el balance granulométrico.")

    # 2. Convertir % pasantes acumulados a fracciones retenidas individuales (f_i, o_i, u_i)
    # Incluye las mallas y el fondo (pan)
    sizes = [s.size_microns for s in raw_sieves]
    f_pass = [s.feed_passing for s in raw_sieves]
    o_pass = [s.overflow_passing for s in raw_sieves]
    u_pass = [s.underflow_passing for s in raw_sieves]

    f_frac, o_frac, u_frac = [], [], []
    for i in range(n):
        prev_f = 100.0 if i == 0 else f_pass[i - 1]
        prev_o = 100.0 if i == 0 else o_pass[i - 1]
        prev_u = 100.0 if i == 0 else u_pass[i - 1]

        f_frac.append(max(0.0, (prev_f - f_pass[i]) / 100.0))
        o_frac.append(max(0.0, (prev_o - o_pass[i]) / 100.0))
        u_frac.append(max(0.0, (prev_u - u_pass[i]) / 100.0))

    # Fracción en el fondo (pan / finos finales)
    f_pan = max(0.0, f_pass[-1] / 100.0)
    o_pan = max(0.0, o_pass[-1] / 100.0)
    u_pan = max(0.0, u_pass[-1] / 100.0)
    f_frac.append(f_pan)
    o_frac.append(o_pan)
    u_frac.append(u_pan)

    # 3. Solución de Mínimos Cuadrados Ponderados para Beta (recuperación de sólidos al underflow)
    # f_i = (1 - beta) * o_i + beta * u_i = o_i - beta * (o_i - u_i)
    # beta * (o_i - u_i) = o_i - f_i
    sum_num = 0.0
    sum_den = 0.0

    for i in range(len(f_frac)):
        diff_ou = o_frac[i] - u_frac[i]
        diff_of = o_frac[i] - f_frac[i]
        
        # Ponderación basada en varianza experimental estimada
        weight = 1.0 / (max(0.015, f_frac[i] * 0.10 + 0.01) ** 2)
        
        sum_num += weight * diff_ou * diff_of
        sum_den += weight * (diff_ou ** 2)

    if sum_den <= 1e-9:
        beta = 0.50
    else:
        beta = sum_num / sum_den

    # Acotar beta a un rango físicamente posible
    beta = float(np.clip(beta, 0.05, 0.95))
    alpha = 1.0 - beta

    # 4. Reconciliación analítica simultánea por mínimos cuadrados para las 3 corrientes
    reconciled_f_frac = []
    reconciled_o_frac = []
    reconciled_u_frac = []

    for i in range(len(f_frac)):
        fi = f_frac[i]
        oi = o_frac[i]
        ui = u_frac[i]

        # Varianza de medición estimada (1% base + 5% del valor relativo)
        v_f = (0.01 + 0.05 * fi) ** 2
        v_o = (0.01 + 0.05 * oi) ** 2
        v_u = (0.01 + 0.05 * ui) ** 2

        # Residuo de balance de masa: r = f - (alpha * o + beta * u)
        r = fi - (alpha * oi + beta * ui)
        # Varianza del residuo: Var(r) = V_f + alpha^2 * V_o + beta^2 * V_u
        v_r = v_f + (alpha ** 2) * v_o + (beta ** 2) * v_u

        if v_r > 1e-12:
            adj_f = fi - (v_f / v_r) * r
            adj_o = oi + (alpha * v_o / v_r) * r
            adj_u = ui + (beta * v_u / v_r) * r
        else:
            adj_f = alpha * oi + beta * ui
            adj_o = oi
            adj_u = ui

        reconciled_f_frac.append(max(0.0, adj_f))
        reconciled_o_frac.append(max(0.0, adj_o))
        reconciled_u_frac.append(max(0.0, adj_u))

    # Normalizar cada corriente para que la suma de fracciones sea exactamente 1.0
    sum_f = sum(reconciled_f_frac)
    sum_o = sum(reconciled_o_frac)
    sum_u = sum(reconciled_u_frac)
    if sum_f > 0: reconciled_f_frac = [v / sum_f for v in reconciled_f_frac]
    if sum_o > 0: reconciled_o_frac = [v / sum_o for v in reconciled_o_frac]
    if sum_u > 0: reconciled_u_frac = [v / sum_u for v in reconciled_u_frac]

    # Reconstruir pasantes acumulados reconciliados para Alimentación, Rebose y Descarga
    rec_f_passing, rec_o_passing, rec_u_passing = [], [], []
    cum_f, cum_o, cum_u = 100.0, 100.0, 100.0
    for i in range(n):
        cum_f -= reconciled_f_frac[i] * 100.0
        cum_o -= reconciled_o_frac[i] * 100.0
        cum_u -= reconciled_u_frac[i] * 100.0
        rec_f_passing.append(round(float(np.clip(cum_f, 0.0, 100.0)), 2))
        rec_o_passing.append(round(float(np.clip(cum_o, 0.0, 100.0)), 2))
        rec_u_passing.append(round(float(np.clip(cum_u, 0.0, 100.0)), 2))

    # 5. Balances de Masa de Sólidos, Agua y Pulpa
    feed_solids = req.feed_solids_tph or 500.0
    underflow_solids = feed_solids * beta
    overflow_solids = feed_solids * alpha

    rho_s = req.solids_density
    rho_l = req.liquid_density

    cp_f = max(5.0, min(80.0, req.feed_cp or 55.0))
    cp_o = max(5.0, min(70.0, req.overflow_cp or 35.0))
    cp_u = max(40.0, min(88.0, req.underflow_cp or 75.0))

    # Flujos de agua
    feed_water = feed_solids * (100.0 - cp_f) / cp_f
    underflow_water = underflow_solids * (100.0 - cp_u) / cp_u
    overflow_water = max(0.0, feed_water - underflow_water)
    
    # Recalcular cp_o si es necesario para mantener balance estricto de agua
    overflow_pulp = overflow_solids + overflow_water
    if overflow_pulp > 0:
        cp_o = round((overflow_solids / overflow_pulp) * 100.0, 2)

    # Recuperación de agua al underflow (Rw)
    rw_pct = float(np.clip((underflow_water / feed_water) * 100.0 if feed_water > 0 else 25.0, 5.0, 60.0))
    rf = rw_pct / 100.0

    # 6. Cálculo de la Curva de Partición Real y Eficiencia Corregida Tromp
    actual_partition = []

    for i in range(n):
        # Usar fracciones reconciliadas balanceadas para consistencia estricta
        f_i = reconciled_f_frac[i]
        u_i = reconciled_u_frac[i]

        if f_i > 1e-6:
            y_act = (beta * u_i) / f_i
        else:
            y_act = 1.0 if sizes[i] > 200 else (rw_pct / 100.0)

        actual_partition.append(float(np.clip(y_act, 0.0, 1.0)))

    # Bypass efectivo de finos Rf: debe ser consistente con la asíntota inferior de la curva
    min_fine_recovery = min(actual_partition[-1], (beta * reconciled_u_frac[-1]) / max(1e-6, reconciled_f_frac[-1]))
    rf = float(np.clip(min(rw_pct / 100.0, min_fine_recovery), 0.01, 0.60))

    reduced_partition = []
    for y_act in actual_partition:
        if y_act <= rf:
            y_red = 0.0
        else:
            y_red = (y_act - rf) / max(0.001, 1.0 - rf)
        reduced_partition.append(float(np.clip(y_red, 0.0, 1.0)))

    # 7. Ajuste del Modelo de Lynch-Rao para d50c y parámetro de agudeza m
    # y'_ci = 1 - exp(-0.693 * (x_i / d50c)^m)
    d50c_est = 100.0
    m_est = 2.4
    r_squared = 0.95

    valid_pts = [(sizes[i], reduced_partition[i]) for i in range(n) if 0.04 < reduced_partition[i] < 0.96]

    if len(valid_pts) >= 2:
        try:
            x_vals = np.array([p[0] for p in valid_pts], dtype=float)
            y_vals = np.array([p[1] for p in valid_pts], dtype=float)

            # Regresión no lineal
            def lynch_model(x, d50c, m):
                return 1.0 - np.exp(-0.693 * ((x / np.maximum(d50c, 1.0)) ** np.clip(m, 0.5, 6.0)))

            popt, _ = curve_fit(lynch_model, x_vals, y_vals, p0=[100.0, 2.4], bounds=([5.0, 0.8], [600.0, 5.0]), maxfev=2000)
            d50c_est = float(popt[0])
            m_est = float(popt[1])

            # R^2
            y_pred = lynch_model(x_vals, d50c_est, m_est)
            ss_res = np.sum((y_vals - y_pred) ** 2)
            ss_tot = np.sum((y_vals - np.mean(y_vals)) ** 2)
            r_squared = float(max(0.0, 1.0 - (ss_res / (ss_tot + 1e-9))))
        except Exception:
            # Estimación por interpolación simple
            d50c_est = calculate_percentile(sizes, [y * 100.0 for y in reduced_partition], 50.0)
            m_est = 2.4
            r_squared = 0.90
    else:
        d50c_est = calculate_percentile(sizes, [y * 100.0 for y in reduced_partition], 50.0)

    # Valores ajustados por modelo para cada malla
    fitted_partition = []
    for s in sizes:
        y_fit = 1.0 - math.exp(-0.693 * ((s / max(1.0, d50c_est)) ** m_est))
        y_fit_actual = rf + (1.0 - rf) * y_fit
        fitted_partition.append(round(float(y_fit_actual * 100.0), 2))

    # Construir filas reconciliadas
    reconciled_rows = []
    residuals = []
    for i in range(n):
        s_obj = raw_sieves[i]
        rec_f = rec_f_passing[i]
        residuals.append((s_obj.feed_passing - rec_f) ** 2)
        reconciled_rows.append(ReconciledSieveRow(
            size_microns=s_obj.size_microns,
            feed_passing=s_obj.feed_passing,
            overflow_passing=s_obj.overflow_passing,
            underflow_passing=s_obj.underflow_passing,
            reconciled_feed_passing=rec_f,
            reconciled_overflow_passing=rec_o_passing[i],
            reconciled_underflow_passing=rec_u_passing[i],
            partition_actual=round(actual_partition[i] * 100.0, 2),
            partition_reduced=round(reduced_partition[i] * 100.0, 2),
            partition_fitted=fitted_partition[i]
        ))

    mse = float(np.mean(residuals)) if residuals else 0.0

    # Propiedades de las corrientes
    feed_pulp = feed_solids + feed_water
    underflow_pulp = underflow_solids + underflow_water
    
    feed_p80 = calculate_percentile(sizes, rec_f_passing, 80.0)
    over_p80 = calculate_percentile(sizes, rec_o_passing, 80.0)
    under_p80 = calculate_percentile(sizes, rec_u_passing, 80.0)

    feed_p50 = calculate_percentile(sizes, rec_f_passing, 50.0)
    over_p50 = calculate_percentile(sizes, rec_o_passing, 50.0)
    under_p50 = calculate_percentile(sizes, rec_u_passing, 50.0)

    # Densidades de pulpa
    feed_rho_p = 100.0 / ((cp_f / rho_s) + ((100.0 - cp_f) / rho_l))
    over_rho_p = 100.0 / ((cp_o / rho_s) + ((100.0 - cp_o) / rho_l))
    under_rho_p = 100.0 / ((cp_u / rho_s) + ((100.0 - cp_u) / rho_l))

    # Concentraciones volumétricas (Cv)
    feed_cv = (cp_f / rho_s) / ((cp_f / rho_s) + ((100.0 - cp_f) / rho_l)) * 100.0
    over_cv = (cp_o / rho_s) / ((cp_o / rho_s) + ((100.0 - cp_o) / rho_l)) * 100.0
    under_cv = (cp_u / rho_s) / ((cp_u / rho_s) + ((100.0 - cp_u) / rho_l)) * 100.0

    feed_prop = StreamProperties(
        solids_tph=round(feed_solids, 2),
        water_tph=round(feed_water, 2),
        pulp_tph=round(feed_pulp, 2),
        pulp_flow_m3h=round(feed_pulp / feed_rho_p, 2),
        solids_cp=round(cp_f, 2),
        solids_cv=round(feed_cv, 2),
        p80_microns=feed_p80,
        p50_microns=feed_p50
    )

    over_prop = StreamProperties(
        solids_tph=round(overflow_solids, 2),
        water_tph=round(overflow_water, 2),
        pulp_tph=round(overflow_pulp, 2),
        pulp_flow_m3h=round(overflow_pulp / over_rho_p, 2),
        solids_cp=round(cp_o, 2),
        solids_cv=round(over_cv, 2),
        p80_microns=over_p80,
        p50_microns=over_p50
    )

    under_prop = StreamProperties(
        solids_tph=round(underflow_solids, 2),
        water_tph=round(underflow_water, 2),
        pulp_tph=round(underflow_pulp, 2),
        pulp_flow_m3h=round(underflow_pulp / under_rho_p, 2),
        solids_cp=round(cp_u, 2),
        solids_cv=round(under_cv, 2),
        p80_microns=under_p80,
        p50_microns=under_p50
    )

    return CycloneBalanceResponse(
        solids_split_underflow_pct=round(beta * 100.0, 2),
        solids_split_overflow_pct=round(alpha * 100.0, 2),
        water_recovery_underflow_pct=round(rw_pct, 2),
        bypass_rf_pct=round(rf * 100.0, 2),
        d50c_microns=round(d50c_est, 1),
        sharpness_m=round(m_est, 2),
        r_squared=round(r_squared, 4),
        feed=feed_prop,
        overflow=over_prop,
        underflow=under_prop,
        reconciled_sieves=reconciled_rows,
        mean_squared_error=round(mse, 4)
    )


# ==============================================================================
# SUBMÓDULO 2: CALIBRACIÓN DEL MODELO DE NAGESWARARAO (PARAMETER ESTIMATION)
# ==============================================================================

def calibrate_nageswararao(req: CycloneCalibrationRequest) -> CycloneCalibrationResponse:
    """
    Despeja y calibra las 4 constantes adimensionales maestras de Nageswararao (KD0, KQ0, KW0, KS0)
    a partir de las variables geométricas y operacionales medidas en planta.
    """
    g = 9.81
    geo = req.geometry

    # Convertir dimensiones geométricas a metros
    Dc = geo.dc_mm / 1000.0
    Di = geo.di_mm / 1000.0
    Do = geo.do_mm / 1000.0
    Du = geo.du_mm / 1000.0
    Lc = geo.lc_mm / 1000.0
    theta = max(5.0, geo.theta_deg)

    # Densidades y concentraciones
    rho_s = req.solids_density
    rho_l = req.liquid_density
    cp_f = req.feed_cp
    cp_u = req.underflow_cp

    # Densidad de pulpa de alimentación (t/m³)
    rho_p = 100.0 / ((cp_f / rho_s) + ((100.0 - cp_f) / rho_l))

    # Concentración volumétrica de sólidos (Cv como fracción 0 a 1)
    cv_f = (cp_f / rho_s) / ((cp_f / rho_s) + ((100.0 - cp_f) / rho_l))
    cv_u = (cp_u / rho_s) / ((cp_u / rho_s) + ((100.0 - cp_u) / rho_l))

    # Caudal de pulpa por ciclón (Q en m³/h)
    num_cyc = max(1, req.num_cyclones_operating)
    if req.feed_slurry_flow_m3h is not None and req.feed_slurry_flow_m3h > 0:
        total_q = req.feed_slurry_flow_m3h
    elif req.feed_solids_tph is not None and req.feed_solids_tph > 0:
        pulp_tph = req.feed_solids_tph / (cp_f / 100.0)
        total_q = pulp_tph / rho_p
    else:
        # Caudal nominal por defecto
        total_q = 600.0 * num_cyc

    Q_per_cyc = total_q / num_cyc  # m³/h

    # Parámetro reológico de finos (lambda)
    # lambda = (Cv / (1 - Cv)) * (d_vf / 100)^(-0.5)
    fines_pct = max(5.0, req.feed_fines_pct_minus_38um or 25.0)
    lambda_param = (cv_f / max(0.01, 1.0 - cv_f)) * ((fines_pct / 100.0) ** -0.5)

    # Presión adimensional (P / (rho_p * g * Dc))
    # P en kPa -> 1000 Pa, rho_p en t/m³ -> 1000 kg/m³ => factor de 1000 se cancela
    P_kpa = max(10.0, req.pressure_kpa)
    P_adim = P_kpa / (rho_p * g * Dc)

    # -------------------------------------------------------------------------
    # 1. Calibración de KD0 (Ecuación de d50c)
    # d50c / Dc = KD0 * (Di/Dc)^0.45 * (Do/Dc)^0.68 * (Du/Dc)^-1.60 * (Lc/Dc)^0.20 * theta^-0.15 * P_adim^-0.22 * lambda^0.93
    # -------------------------------------------------------------------------
    d50c_m = (req.d50c_microns) * 1e-6
    d50c_over_dc = d50c_m / Dc

    term_d = (
        ((Di / Dc) ** 0.45) *
        ((Do / Dc) ** 0.68) *
        ((Du / Dc) ** -1.60) *
        ((Lc / Dc) ** 0.20) *
        (theta ** -0.15) *
        (P_adim ** -0.22) *
        (lambda_param ** 0.93)
    )
    kd0 = d50c_over_dc / max(1e-12, term_d)

    # -------------------------------------------------------------------------
    # 2. Calibración de KQ0 (Ecuación de Capacidad)
    # Q / (Dc^2 * sqrt(P / rho_p)) = KQ0 * (Di/Dc)^0.45 * (Do/Dc)^0.85 * (Lc/Dc)^0.20 * theta^-0.10 * lambda^0.10
    # -------------------------------------------------------------------------
    capacity_lhs = Q_per_cyc / ((Dc ** 2) * math.sqrt(P_kpa / rho_p))
    term_q = (
        ((Di / Dc) ** 0.45) *
        ((Do / Dc) ** 0.85) *
        ((Lc / Dc) ** 0.20) *
        (theta ** -0.10) *
        (lambda_param ** 0.10)
    )
    kq0 = capacity_lhs / max(1e-12, term_q)

    # -------------------------------------------------------------------------
    # 3. Calibración de KW0 (Ecuación de Recuperación de Agua)
    # Rw = KW0 * (Do/Dc)^-1.19 * (Du/Dc)^2.40 * (Lc/Dc)^0.22 * theta^-0.24 * P_adim^-0.53 * lambda^0.27
    # -------------------------------------------------------------------------
    Rw_frac = req.water_recovery_underflow_pct / 100.0
    term_w = (
        ((Do / Dc) ** -1.19) *
        ((Du / Dc) ** 2.40) *
        ((Lc / Dc) ** 0.22) *
        (theta ** -0.24) *
        (P_adim ** -0.53) *
        (lambda_param ** 0.27)
    )
    kw0 = Rw_frac / max(1e-12, term_w)

    # -------------------------------------------------------------------------
    # 4. Calibración de KS0 (Ecuación de Concentración Volumétrica en Underflow)
    # Cvu = KS0 * (Di/Dc)^-0.50 * (Do/Dc)^-0.50 * (Du/Dc)^0.75 * theta^-0.15 * Cvf^0.65
    # -------------------------------------------------------------------------
    term_s = (
        ((Di / Dc) ** -0.50) *
        ((Do / Dc) ** -0.50) *
        ((Du / Dc) ** 0.75) *
        (theta ** -0.15) *
        (cv_f ** 0.65)
    )
    ks0 = cv_u / max(1e-12, term_s)

    # Evaluación de riesgo de Roping (Acodalamiento)
    is_roping = cv_u >= 0.56
    if cv_u >= 0.58:
        roping_msg = f"PELIGRO CRÍTICO: Concentración volumétrica en underflow de {round(cv_u*100, 1)}% supera el límite de roping (56-58%). El hidrociclón opera en chorro cilíndrico colapsado."
    elif cv_u >= 0.54:
        roping_msg = f"PRECAUCIÓN: Concentración volumétrica de {round(cv_u*100, 1)}% se encuentra en zona de transición hacia roping."
    else:
        roping_msg = f"ÓPTIMO: Operación en régimen abanico (Spray) estable con {round(cv_u*100, 1)}% de sólidos en volumen."

    return CycloneCalibrationResponse(
        kd0=float(kd0),
        kq0=float(kq0),
        kw0=float(kw0),
        ks0=float(ks0),
        sharpness_m=round(req.sharpness_m or 2.4, 2),
        lambda_param=round(lambda_param, 4),
        feed_cv_pct=round(cv_f * 100.0, 2),
        underflow_cv_pct=round(cv_u * 100.0, 2),
        q_per_cyclone_m3h=round(Q_per_cyc, 2),
        roping_risk=is_roping,
        roping_message=roping_msg
    )


# ==============================================================================
# SUBMÓDULO 3: SIMULACIÓN / MODELAMIENTO PREDICTIVO (FORWARD SIMULATION)
# ==============================================================================

def simulate_nageswararao(req: CycloneSimulationRequest) -> CycloneSimulationResponse:
    """
    Ejecuta la simulación predictiva de clasificación del hidrociclón utilizando
    las constantes de Nageswararao calibradas y las nuevas condiciones operacionales.
    """
    g = 9.81
    geo = req.geometry

    Dc = geo.dc_mm / 1000.0
    Di = geo.di_mm / 1000.0
    Do = geo.do_mm / 1000.0
    Du = geo.du_mm / 1000.0
    Lc = geo.lc_mm / 1000.0
    theta = max(5.0, geo.theta_deg)

    rho_s = req.solids_density
    rho_l = req.liquid_density
    cp_f = req.feed_cp
    feed_solids = req.feed_solids_tph
    num_cyc = max(1, req.num_cyclones_operating)

    # Cálculos reológicos y volumétricos
    rho_p = 100.0 / ((cp_f / rho_s) + ((100.0 - cp_f) / rho_l))
    cv_f = (cp_f / rho_s) / ((cp_f / rho_s) + ((100.0 - cp_f) / rho_l))

    # Flujo másico total y volumétrico
    feed_water = feed_solids * (100.0 - cp_f) / cp_f
    feed_pulp_tph = feed_solids + feed_water
    total_q_m3h = feed_pulp_tph / rho_p
    Q_per_cyc = total_q_m3h / num_cyc

    # % Finos en la nueva alimentación
    sorted_feed_psd = sorted(req.feed_psd, key=lambda s: s.size_microns, reverse=True)
    fines_pct = 25.0
    for item in sorted_feed_psd:
        if item.size_microns <= 45.0:
            fines_pct = item.passing_pct
            break

    lambda_param = (cv_f / max(0.01, 1.0 - cv_f)) * ((max(5.0, fines_pct) / 100.0) ** -0.5)

    # 1. Predicción de la Presión (o Caudal si se fijó la presión)
    term_q = (
        ((Di / Dc) ** 0.45) *
        ((Do / Dc) ** 0.85) *
        ((Lc / Dc) ** 0.20) *
        (theta ** -0.10) *
        (lambda_param ** 0.10)
    )

    if req.fixed_pressure_kpa and req.fixed_pressure_kpa > 0:
        pred_P_kpa = req.fixed_pressure_kpa
    else:
        # Invertir ecuación de capacidad: Q / (Dc^2 * sqrt(P/rho_p)) = KQ0 * term_q
        denom = req.kq0 * (Dc ** 2) * term_q
        if denom > 0:
            sqrt_term = Q_per_cyc / denom
            pred_P_kpa = (sqrt_term ** 2) * rho_p
        else:
            pred_P_kpa = 100.0

    pred_P_kpa = float(np.clip(pred_P_kpa, 15.0, 500.0))
    P_adim = pred_P_kpa / (rho_p * g * Dc)

    # 2. Predicción de d50c
    term_d = (
        ((Di / Dc) ** 0.45) *
        ((Do / Dc) ** 0.68) *
        ((Du / Dc) ** -1.60) *
        ((Lc / Dc) ** 0.20) *
        (theta ** -0.15) *
        (P_adim ** -0.22) *
        (lambda_param ** 0.93)
    )
    d50c_m = req.kd0 * Dc * term_d
    pred_d50c_microns = float(np.clip(d50c_m * 1e6, 15.0, 800.0))

    # 3. Predicción de Recuperación de Agua (Rw)
    term_w = (
        ((Do / Dc) ** -1.19) *
        ((Du / Dc) ** 2.40) *
        ((Lc / Dc) ** 0.22) *
        (theta ** -0.24) *
        (P_adim ** -0.53) *
        (lambda_param ** 0.27)
    )
    pred_rw_frac = req.kw0 * term_w
    pred_rw_pct = float(np.clip(pred_rw_frac * 100.0, 8.0, 55.0))
    rf = pred_rw_pct / 100.0

    # 4. Predicción de Cvu y Chequeo de Roping
    term_s = (
        ((Di / Dc) ** -0.50) *
        ((Do / Dc) ** -0.50) *
        ((Du / Dc) ** 0.75) *
        (theta ** -0.15) *
        (cv_f ** 0.65)
    )
    pred_cv_u = float(np.clip(req.ks0 * term_s, 0.20, 0.65))
    is_roping = pred_cv_u >= 0.56
    if pred_cv_u >= 0.58:
        roping_severity = "Acodalamiento Severo (Roping)"
    elif pred_cv_u >= 0.54:
        roping_severity = "Transición Crítica"
    else:
        roping_severity = "Normal (Spray)"

    # Convertir Cvu a Cpu (% sólidos en peso)
    # Cv = (Cp/rho_s) / [Cp/rho_s + (1-Cp)/rho_l] => Cp = (Cv * rho_s) / [Cv * rho_s + (1 - Cv) * rho_l]
    pred_cp_u = (pred_cv_u * rho_s) / ((pred_cv_u * rho_s) + ((1.0 - pred_cv_u) * rho_l)) * 100.0

    # 5. Partición Granulométrica (Curva Tromp)
    m = req.sharpness_m or 2.4
    sizes = [item.size_microns for item in sorted_feed_psd]
    f_pass = [item.passing_pct for item in sorted_feed_psd]
    n = len(sizes)

    # Fracciones retenidas de alimentación
    f_frac = []
    for i in range(n):
        prev_f = 100.0 if i == 0 else f_pass[i - 1]
        f_frac.append(max(0.0, (prev_f - f_pass[i]) / 100.0))
    f_pan = max(0.0, f_pass[-1] / 100.0)
    f_frac.append(f_pan)

    # Calcular recuperación por tamaño
    u_mass_fractions = []
    o_mass_fractions = []
    predicted_sieves = []

    for i in range(n):
        s = sizes[i]
        # Eficiencia reducida
        y_red = 1.0 - math.exp(-0.693 * ((s / max(1.0, pred_d50c_microns)) ** m))
        # Eficiencia real
        y_act = rf + (1.0 - rf) * y_red
        y_act = float(np.clip(y_act, rf, 1.0))

        u_mass = f_frac[i] * y_act
        o_mass = f_frac[i] * (1.0 - y_act)

        u_mass_fractions.append(u_mass)
        o_mass_fractions.append(o_mass)

    # Fondo (pan)
    y_act_pan = rf  # Los ultrafinos siguen exclusivamente al agua
    u_mass_pan = f_pan * y_act_pan
    o_mass_pan = f_pan * (1.0 - y_act_pan)
    u_mass_fractions.append(u_mass_pan)
    o_mass_fractions.append(o_mass_pan)

    total_u_solids_frac = sum(u_mass_fractions)
    total_o_solids_frac = sum(o_mass_fractions)

    beta = total_u_solids_frac
    alpha = total_o_solids_frac

    under_solids_tph = feed_solids * beta
    over_solids_tph = feed_solids * alpha

    under_water_tph = feed_water * rf
    over_water_tph = feed_water - under_water_tph

    under_pulp_tph = under_solids_tph + under_water_tph
    over_pulp_tph = over_solids_tph + over_water_tph

    calc_cp_u = (under_solids_tph / under_pulp_tph * 100.0) if under_pulp_tph > 0 else pred_cp_u
    calc_cp_o = (over_solids_tph / over_pulp_tph * 100.0) if over_pulp_tph > 0 else 30.0

    # Reconstruir pasantes acumulados para underflow y overflow
    u_pass_pred = []
    o_pass_pred = []
    cum_u = 100.0
    cum_o = 100.0

    for i in range(n):
        u_norm = (u_mass_fractions[i] / total_u_solids_frac) if total_u_solids_frac > 0 else 0.0
        o_norm = (o_mass_fractions[i] / total_o_solids_frac) if total_o_solids_frac > 0 else 0.0

        cum_u -= u_norm * 100.0
        cum_o -= o_norm * 100.0

        u_pass_pred.append(round(float(np.clip(cum_u, 0.0, 100.0)), 2))
        o_pass_pred.append(round(float(np.clip(cum_o, 0.0, 100.0)), 2))

        # Partición
        s = sizes[i]
        y_red = 1.0 - math.exp(-0.693 * ((s / max(1.0, pred_d50c_microns)) ** m))
        y_act = rf + (1.0 - rf) * y_red

        predicted_sieves.append({
            "size_microns": s,
            "feed_passing": f_pass[i],
            "overflow_passing": o_pass_pred[-1],
            "underflow_passing": u_pass_pred[-1],
            "partition_actual": round(y_act * 100.0, 2),
            "partition_reduced": round(y_red * 100.0, 2)
        })

    # P80 y P50
    feed_p80 = calculate_percentile(sizes, f_pass, 80.0)
    over_p80 = calculate_percentile(sizes, o_pass_pred, 80.0)
    under_p80 = calculate_percentile(sizes, u_pass_pred, 80.0)

    feed_p50 = calculate_percentile(sizes, f_pass, 50.0)
    over_p50 = calculate_percentile(sizes, o_pass_pred, 50.0)
    under_p50 = calculate_percentile(sizes, u_pass_pred, 50.0)

    over_rho_p = 100.0 / ((calc_cp_o / rho_s) + ((100.0 - calc_cp_o) / rho_l))
    under_rho_p = 100.0 / ((calc_cp_u / rho_s) + ((100.0 - calc_cp_u) / rho_l))

    feed_prop = StreamProperties(
        solids_tph=round(feed_solids, 2),
        water_tph=round(feed_water, 2),
        pulp_tph=round(feed_pulp_tph, 2),
        pulp_flow_m3h=round(total_q_m3h, 2),
        solids_cp=round(cp_f, 2),
        solids_cv=round(cv_f * 100.0, 2),
        p80_microns=feed_p80,
        p50_microns=feed_p50
    )

    over_prop = StreamProperties(
        solids_tph=round(over_solids_tph, 2),
        water_tph=round(over_water_tph, 2),
        pulp_tph=round(over_pulp_tph, 2),
        pulp_flow_m3h=round(over_pulp_tph / over_rho_p, 2),
        solids_cp=round(calc_cp_o, 2),
        solids_cv=round(((calc_cp_o / rho_s) / ((calc_cp_o / rho_s) + ((100.0 - calc_cp_o) / rho_l))) * 100.0, 2),
        p80_microns=over_p80,
        p50_microns=over_p50
    )

    under_prop = StreamProperties(
        solids_tph=round(under_solids_tph, 2),
        water_tph=round(under_water_tph, 2),
        pulp_tph=round(under_pulp_tph, 2),
        pulp_flow_m3h=round(under_pulp_tph / under_rho_p, 2),
        solids_cp=round(calc_cp_u, 2),
        solids_cv=round(pred_cv_u * 100.0, 2),
        p80_microns=under_p80,
        p50_microns=under_p50
    )

    return CycloneSimulationResponse(
        predicted_pressure_kpa=round(pred_P_kpa, 1),
        total_slurry_flow_m3h=round(total_q_m3h, 1),
        q_per_cyclone_m3h=round(Q_per_cyc, 1),
        d50c_microns=round(pred_d50c_microns, 1),
        water_recovery_underflow_pct=round(pred_rw_pct, 2),
        solids_split_underflow_pct=round(beta * 100.0, 2),
        solids_split_overflow_pct=round(alpha * 100.0, 2),
        underflow_cv_pct=round(pred_cv_u * 100.0, 2),
        underflow_cp_pct=round(calc_cp_u, 2),
        overflow_cp_pct=round(calc_cp_o, 2),
        is_roping=is_roping,
        roping_severity=roping_severity,
        feed=feed_prop,
        overflow=over_prop,
        underflow=under_prop,
        predicted_sieves=predicted_sieves
    )
