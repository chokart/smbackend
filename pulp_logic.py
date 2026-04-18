from pulp_models import PulpCalculationRequest, PulpCalculationResponse

def calculate_pulp_properties(req: PulpCalculationRequest) -> PulpCalculationResponse:
    # 1. Parámetros base
    gs = req.rho_solid or 2.7
    rl = req.rho_liquid or 1.0
    
    # Tolerancia para comparaciones matemáticas
    tol = 1e-6

    # Extraer entradas
    i_ws = req.ton_solid
    i_ww = req.ton_liquid
    i_wp = req.ton_pulp
    i_vs = req.vol_solid
    i_vw = req.vol_liquid
    i_vp = req.vol_pulp
    i_rp = req.rho_pulp
    i_cw = req.percent_solid_w
    i_cv = req.percent_solid_v

    # Nuestro objetivo es encontrar Ws (ton_solid) y Ww (ton_liquid).
    ws, ww = None, None

    # Intentaremos reducir cualquier par de variables a Ws y Ww
    
    # -----------------------------
    # GRUPO 1: Conversiones directas de Volumen a Masa
    # -----------------------------
    if i_vs is not None:
        if i_ws is None: i_ws = i_vs * gs
    if i_vw is not None:
        if i_ww is None: i_ww = i_vw * rl

    # -----------------------------
    # GRUPO 2: Resolución de ecuaciones
    # -----------------------------

    if i_ws is not None and i_ww is not None:
        ws, ww = i_ws, i_ww

    elif i_wp is not None and i_cw is not None:
        ws = i_wp * (i_cw / 100.0)
        ww = i_wp - ws

    elif i_wp is not None and i_ws is not None:
        ws = i_ws
        ww = i_wp - ws

    elif i_wp is not None and i_ww is not None:
        ww = i_ww
        ws = i_wp - ww

    elif i_wp is not None and i_rp is not None:
        # wp = ws + ww
        # vp = wp / rp = ws/gs + ww/rl
        vp = i_wp / i_rp if i_rp > 0 else 0
        factor = (1.0 - rl/gs)
        if abs(factor) > tol:
            ws = (i_wp - vp * rl) / factor
            ww = i_wp - ws

    elif i_vp is not None and i_rp is not None:
        wp = i_vp * i_rp
        factor = (1.0 - rl/gs)
        if abs(factor) > tol:
            ws = (wp - i_vp * rl) / factor
            ww = wp - ws

    elif i_vp is not None and i_cw is not None:
        # vp = ws/gs + ww/rl
        # ws = wp * cw/100 -> ww = wp * (1 - cw/100)
        # vp = wp * [ (cw/100)/gs + (1 - cw/100)/rl ]
        cw_frac = i_cw / 100.0
        factor_vp = (cw_frac / gs) + ((1.0 - cw_frac) / rl)
        if factor_vp > tol:
            wp = i_vp / factor_vp
            ws = wp * cw_frac
            ww = wp - ws

    elif i_vp is not None and i_cv is not None:
        cv_frac = i_cv / 100.0
        vs = i_vp * cv_frac
        vw = i_vp - vs
        ws = vs * gs
        ww = vw * rl

    elif i_vp is not None and i_ws is not None:
        ws = i_ws
        vs = ws / gs
        vw = i_vp - vs
        ww = vw * rl

    elif i_vp is not None and i_ww is not None:
        ww = i_ww
        vw = ww / rl
        vs = i_vp - vw
        ws = vs * gs

    elif i_cw is not None and i_rp is not None:
        vp_base = 1.0
        wp_base = vp_base * i_rp
        ws_base = wp_base * (i_cw / 100.0)
        ww_base = wp_base - ws_base
        
        scale = 1.0
        if i_wp is not None: scale = i_wp / wp_base
        elif i_vp is not None: scale = i_vp / vp_base
        elif i_ws is not None: scale = i_ws / ws_base if ws_base > 0 else 1.0
        elif i_ww is not None: scale = i_ww / ww_base if ww_base > 0 else 1.0
            
        ws = ws_base * scale
        ww = ww_base * scale

    elif i_cv is not None and i_rp is not None:
        vp_base = 1.0
        wp_base = vp_base * i_rp
        vs_base = vp_base * (i_cv / 100.0)
        ws_base = vs_base * gs
        ww_base = wp_base - ws_base
        
        scale = 1.0
        if i_wp is not None: scale = i_wp / wp_base
        elif i_vp is not None: scale = i_vp / vp_base
        elif i_ws is not None: scale = i_ws / ws_base if ws_base > 0 else 1.0
        
        ws = ws_base * scale
        ww = ww_base * scale

    elif i_cw is not None and i_ws is not None:
        ws = i_ws
        if i_cw > 0:
            wp = ws / (i_cw / 100.0)
            ww = wp - ws

    elif i_cw is not None and i_ww is not None:
        ww = i_ww
        if i_cw < 100:
            wp = ww / (1.0 - i_cw / 100.0)
            ws = wp - ww

    elif i_cv is not None and i_ws is not None:
        ws = i_ws
        vs = ws / gs
        if i_cv > 0:
            vp = vs / (i_cv / 100.0)
            vw = vp - vs
            ww = vw * rl

    elif i_rp is not None and i_ws is not None:
        ws = i_ws
        vs = ws / gs
        factor = (i_rp / rl - 1.0)
        if abs(factor) > tol:
            ww = (ws - i_rp * vs) / factor

    elif i_rp is not None and i_ww is not None:
        ww = i_ww
        vw = ww / rl
        factor = (1.0 - i_rp / gs)
        if abs(factor) > tol:
            ws = (i_rp * vw - ww) / factor

    # Si después de todas las reglas no logramos definir Ws y Ww, lanzamos error
    if ws is None or ww is None:
        raise ValueError("Datos insuficientes o inconsistentes. Ingrese al menos dos variables principales válidas (ej. Flujo Másico y Densidad, o %Sólidos y Densidad).")

    # Si Ws o Ww dan resultados negativos por incongruencias matemáticas en los datos del usuario:
    if ws < -tol or ww < -tol:
        raise ValueError("Los datos ingresados generan resultados físicamente imposibles (valores negativos). Revise las variables ingresadas.")

    # Asegurarse de que no sean "ligeramente negativos" por errores de flotante
    ws = max(0.0, ws)
    ww = max(0.0, ww)

    # 3. Calcular el resto de variables
    vs = ws / gs
    vw = ww / rl
    wp = ws + ww
    vp = vs + vw
    
    rp = wp / vp if vp > tol else 0
    cw = (ws / wp * 100.0) if wp > tol else 0
    cv = (vs / vp * 100.0) if vp > tol else 0
    ls_ratio = ww / ws if ws > tol else 0
    
    return PulpCalculationResponse(
        ton_solid=round(ws, 4),
        ton_liquid=round(ww, 4),
        ton_pulp=round(wp, 4),
        vol_solid=round(vs, 4),
        vol_liquid=round(vw, 4),
        vol_pulp=round(vp, 4),
        rho_solid=gs,
        rho_liquid=rl,
        rho_pulp=round(rp, 4),
        percent_solid_w=round(cw, 4),
        percent_solid_v=round(cv, 4),
        liquid_solid_ratio=round(ls_ratio, 4)
    )
