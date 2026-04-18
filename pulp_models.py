from pydantic import BaseModel
from typing import Optional

class PulpCalculationRequest(BaseModel):
    # Variables de entrada (todas opcionales)
    ton_solid: Optional[float] = None      # Ws (tph)
    ton_liquid: Optional[float] = None     # Ww (tph)
    ton_pulp: Optional[float] = None       # Wp (tph)
    vol_solid: Optional[float] = None      # Vs (m3h)
    vol_liquid: Optional[float] = None     # Vw (m3h)
    vol_pulp: Optional[float] = None       # Vp (m3h)
    rho_solid: Optional[float] = 2.7       # Gs (Gravedad específica)
    rho_liquid: Optional[float] = 1.0      # rl (Densidad del líquido)
    rho_pulp: Optional[float] = None       # rp (Densidad de pulpa)
    percent_solid_w: Optional[float] = None # %Cw (Concentración en peso)
    percent_solid_v: Optional[float] = None # %Cv (Concentración en volumen)

class PulpCalculationResponse(BaseModel):
    # Resultado completo
    ton_solid: float
    ton_liquid: float
    ton_pulp: float
    vol_solid: float
    vol_liquid: float
    vol_pulp: float
    rho_solid: float
    rho_liquid: float
    rho_pulp: float
    percent_solid_w: float
    percent_solid_v: float
    liquid_solid_ratio: float # Relación L/S
