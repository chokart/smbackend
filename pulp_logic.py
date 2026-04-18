from pulp_models import PulpCalculationRequest, PulpCalculationResponse

def calculate_pulp_properties(req: PulpCalculationRequest) -> PulpCalculationResponse:
    # 1. Parámetros base
    gs = req.rho_solid or 2.7
    rl = req.rho_liquid or 1.0
    
    # Variables a encontrar
    ws, ww = None, None
    
    # Extraer entradas
    input_ws = req.ton_solid
    input_ww = req.ton_liquid
    input_wp = req.ton_pulp
    input_rp = req.rho_pulp
    input_vp = req.vol_pulp
    input_cw = req.percent_solid_w
    
    # Lógica de resolución (Intentar encontrar ws y ww)
    
    # Caso A: Tenemos Ton Solid y Ton Liquid directamente
    if input_ws is not None and input_ww is not None:
        ws, ww = input_ws, input_ww
        
    # Caso B: Tenemos Ton Pulpa y % Sólidos en peso
    elif input_wp is not None and input_cw is not None:
        ws = input_wp * (input_cw / 100.0)
        ww = input_wp - ws
        
    # Caso C: Tenemos Densidad de Pulpa y Flujo Volumétrico
    elif input_rp is not None and input_vp is not None:
        wp = input_vp * input_rp
        # Resolver para ws: wp = ws + (vp - ws/gs)*rl
        # ws(1 - rl/gs) = wp - vp*rl
        factor = (1.0 - rl/gs)
        if abs(factor) > 1e-6:
            ws = (wp - input_vp * rl) / factor
            ww = wp - ws
            
    # Caso D: Tenemos % Sólidos y Densidad de Pulpa
    elif input_cw is not None and input_rp is not None:
        # Usamos una base unitaria de 1 m3 de pulpa para hallar proporciones
        vp_base = 1.0
        wp_base = vp_base * input_rp
        ws_base = wp_base * (input_cw / 100.0)
        # Si el usuario ingresó algún tonelaje o volumen, escalamos. Si no, devolvemos proporciones
        scale = input_wp / wp_base if input_wp else (input_vp / vp_base if input_vp else 1.0)
        ws, ww = ws_base * scale, (wp_base - ws_base) * scale

    # Caso E: % Sólidos y Ton Solid
    elif input_cw is not None and input_ws is not None:
        ws = input_ws
        if input_cw > 0:
            wp = ws / (input_cw / 100.0)
            ww = wp - ws
    
    # Caso F: Densidad de pulpa y Ton Solid
    elif input_rp is not None and input_ws is not None:
        ws = input_ws
        vs = ws / gs
        # rp = (ws + ww) / (vs + ww/rl) -> rp*vs + rp*ww/rl = ws + ww
        # ww(rp/rl - 1) = ws - rp*vs
        factor = (input_rp / rl - 1.0)
        if abs(factor) > 1e-6:
            ww = (ws - input_rp * vs) / factor

    # Si no se pudo resolver con las combinaciones anteriores, lanzar error
    if ws is None or ww is None:
        # Intento de emergencia: si solo hay tonelaje de sólido, asumir pulpa al 100% o algo similar?
        # Mejor pedir más datos.
        raise ValueError("Se requieren al menos dos variables (ej. %Sólidos y Densidad, o Ton Sólido y %Sólidos) para calcular la pulpa.")

    # 3. Calcular el resto de variables a partir de ws y ww
    vs = ws / gs
    vw = ww / rl
    wp = ws + ww
    vp = vs + vw
    rp = wp / vp if vp > 0 else 0
    cw = (ws / wp * 100.0) if wp > 0 else 0
    cv = (vs / vp * 100.0) if vp > 0 else 0
    ls_ratio = ww / ws if ws > 0 else 0
    
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
