import numpy as np
from typing import List, Tuple
from hydrocyclone_models import RaoLynchRequest, RaoLynchResponse, PartitionCurvePoint

def calculate_rao_lynch(request: RaoLynchRequest) -> RaoLynchResponse:
    # 1. Procesar granulometría experimental
    sizes = [s.mesh_size for s in request.sieves]
    w_feed = np.array([s.weight_feed for s in request.sieves])
    w_overflow = np.array([s.weight_overflow for s in request.sieves])
    w_underflow = np.array([s.weight_underflow for s in request.sieves])
    
    # Porcentajes de peso retenido
    total_f = np.sum(w_feed)
    total_o = np.sum(w_overflow)
    total_u = np.sum(w_underflow)
    
    pr_f = (w_feed / total_f) * 100
    pr_o = (w_overflow / total_o) * 100
    pr_u = (w_underflow / total_u) * 100
    
    # 2. Balance de masa simple para encontrar la recuperación de sólidos (S)
    # F * f_i = O * o_i + U * u_i  => (O + U) * f_i = O * o_i + U * u_i
    # (O/U + 1) * f_i = (O/U) * o_i + u_i
    # Relación U/F (S) estimada por mallas: S = (f_i - o_i) / (u_i - o_i)
    # Usaremos un promedio ponderado o mínimos cuadrados para S
    S_values = []
    for f, o, u in zip(pr_f, pr_o, pr_u):
        if abs(u - o) > 0.01:
            S_i = (f - o) / (u - o)
            if 0 < S_i < 1:
                S_values.append(S_i)
    
    S_solids = np.mean(S_values) if S_values else 0.5
    
    # 3. Calcular la curva de partición experimental (E_a)
    # E_a_i = (U * u_i) / (F * f_i) = S * (u_i / f_i)
    E_a = []
    for f, u in zip(pr_f, pr_u):
        if f > 0:
            E_a.append(S_solids * (u / f))
        else:
            E_a.append(0)
    E_a = np.array(E_a)
    
    # 4. Estimar Rf (Bypass de agua/finos) y d50c experimental
    # Rf es usualmente el valor de la curva E_a para tamaños muy finos
    Rf = E_a[-1] # Simplificación: último punto (finos)
    
    # Eficiencia corregida: Ec = (Ea - Rf) / (1 - Rf)
    E_c = (E_a - Rf) / (1 - Rf)
    E_c = np.clip(E_c, 0, 1) # Asegurar rango 0-1
    
    # d50c experimental: Tamaño donde Ec = 0.5 (Interpolación)
    d50c_experimental = np.interp(0.5, E_c[::-1], sizes[::-1])
    
    # 5. Modelo Rao & Lynch (Predicción)
    # Constantes estándar si no se proveen
    K1 = request.k1 or 13.5
    K2 = request.k2 or 0.05
    K3 = request.k3 or 35.0
    
    # Ecuaciones Rao & Lynch
    # Capacidad Q (m3/h)
    Q = K1 * (request.geometry.Do**0.73) * (request.geometry.Di**0.53) * (request.pressure**0.42)
    
    # Bypass de agua Rf (%) - Simplificado
    # Rf_pred = K2 * (request.geometry.Du / request.geometry.Do)**0.8 * (request.pressure**-0.1)
    
    # d50c corregido predicho (micrones)
    # d50c = K3 * Do^0.52 * Di^-0.5 * Du^-0.2 * P^-0.3 * exp(Cv * a)
    # Cv: % sólidos en volumen
    rho_s = request.solid_density
    rho_l = 1.0 # Densidad agua
    Cv = (request.feed_p_solids / rho_s) / (request.feed_p_solids / rho_s + (100 - request.feed_p_solids) / rho_l) * 100
    
    d50c_predicted = K3 * (request.geometry.Do**0.52) * (request.geometry.Di**-0.5) * (request.geometry.Du**-0.2) * (request.pressure**-0.3) * np.exp(Cv * 0.063)
    
    # 6. Preparar respuesta
    partition_points = []
    for i in range(len(sizes)):
        partition_points.append(PartitionCurvePoint(
            size=sizes[i],
            actual_recovery=float(E_a[i]),
            corrected_recovery=float(E_c[i])
        ))
        
    return RaoLynchResponse(
        d50c_predicted=float(d50c_predicted),
        d50c_experimental=float(d50c_experimental),
        water_bypass_Rf=float(Rf),
        capacity_Q=float(Q),
        partition_curve=partition_points,
        granulometry_report={
            "feed_total": float(total_f),
            "overflow_total": float(total_o),
            "underflow_total": float(total_u),
            "S_solids": float(S_solids)
        }
    )
