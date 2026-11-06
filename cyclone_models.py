from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any

# ==============================================================================
# SUBMÓDULO 1: BALANCE DE MATERIA Y GRANULOMÉTRICO (TWO-PRODUCT RECONCILIATION)
# ==============================================================================

class SieveDataPoint(BaseModel):
    size_microns: float = Field(..., description="Abertura de malla en micrómetros")
    feed_passing: float = Field(..., description="% pasante acumulado en la alimentación")
    overflow_passing: float = Field(..., description="% pasante acumulado en el rebose (over)")
    underflow_passing: float = Field(..., description="% pasante acumulado en la descarga (under)")

class StreamProperties(BaseModel):
    solids_tph: Optional[float] = None
    water_tph: Optional[float] = None
    pulp_tph: Optional[float] = None
    pulp_flow_m3h: Optional[float] = None
    solids_cp: Optional[float] = None  # % sólidos en peso
    solids_cv: Optional[float] = None  # % sólidos en volumen
    p80_microns: Optional[float] = None
    p50_microns: Optional[float] = None

class ReconciledSieveRow(BaseModel):
    size_microns: float
    feed_passing: float
    overflow_passing: float
    underflow_passing: float
    reconciled_feed_passing: float
    reconciled_overflow_passing: float
    reconciled_underflow_passing: float
    partition_actual: float      # Eficiencia real hacia el underflow (%)
    partition_reduced: float     # Eficiencia corregida hacia el underflow (%)
    partition_fitted: float      # Curva ajustada por modelo (Lynch-Rao) (%)

class CycloneBalanceRequest(BaseModel):
    sieves: List[SieveDataPoint]
    feed_solids_tph: Optional[float] = Field(default=500.0, description="Flujo de sólidos secos de alimentación (t/h)")
    feed_cp: Optional[float] = Field(default=55.0, description="% de sólidos en peso de alimentación")
    overflow_cp: Optional[float] = Field(default=35.0, description="% de sólidos en peso del rebose")
    underflow_cp: Optional[float] = Field(default=75.0, description="% de sólidos en peso de la descarga")
    solids_density: float = Field(default=2.7, description="Gravedad específica / densidad de sólidos (t/m³)")
    liquid_density: float = Field(default=1.0, description="Densidad del líquido portador (t/m³)")

class CycloneBalanceResponse(BaseModel):
    solids_split_underflow_pct: float   # % de sólidos a la descarga (beta * 100)
    solids_split_overflow_pct: float    # % de sólidos al rebose (alpha * 100)
    water_recovery_underflow_pct: float # Recuperación de agua al underflow (Rw * 100)
    bypass_rf_pct: float                # Bypass de ultrafinos (Rf * 100)
    d50c_microns: float                 # Tamaño de corte corregido experimental
    sharpness_m: float                  # Parámetro de agudeza de separación m
    r_squared: float                    # Calidad de ajuste de la curva de partición
    
    feed: StreamProperties
    overflow: StreamProperties
    underflow: StreamProperties
    
    reconciled_sieves: List[ReconciledSieveRow]
    mean_squared_error: float


# ==============================================================================
# SUBMÓDULO 2: CALIBRACIÓN DEL MODELO DE NAGESWARARAO (PARAMETER ESTIMATION)
# ==============================================================================

class CycloneGeometry(BaseModel):
    dc_mm: float = Field(default=500.0, description="Diámetro interno del hidrociclón Dc (mm)")
    di_mm: float = Field(default=125.0, description="Diámetro equivalente de entrada Di (mm)")
    do_mm: float = Field(default=175.0, description="Diámetro del vortex finder Do (mm)")
    du_mm: float = Field(default=100.0, description="Diámetro del apex / tobera Du (mm)")
    lc_mm: float = Field(default=1200.0, description="Longitud de sección cilíndrica Lc (mm)")
    theta_deg: float = Field(default=20.0, description="Ángulo del cono en grados (°)")

class CycloneCalibrationRequest(BaseModel):
    geometry: CycloneGeometry
    num_cyclones_operating: int = Field(default=4, description="Número de ciclones operando en batería")
    
    # Datos de operación conocidos / medidos
    pressure_kpa: float = Field(..., description="Presión de alimentación en manifold (kPa)")
    feed_slurry_flow_m3h: Optional[float] = Field(default=None, description="Caudal volumétrico total de pulpa (m³/h)")
    feed_solids_tph: Optional[float] = Field(default=None, description="Flujo másico total de sólidos (t/h)")
    feed_cp: float = Field(..., description="% de sólidos en peso de la alimentación (Cp)")
    solids_density: float = Field(default=2.7, description="Gravedad específica de sólidos (t/m³)")
    liquid_density: float = Field(default=1.0, description="Densidad del líquido (t/m³)")
    
    # Resultados del balance (Submódulo 1 o medidos)
    d50c_microns: float = Field(..., description="Tamaño de corte corregido d50c (µm)")
    water_recovery_underflow_pct: float = Field(..., description="Recuperación de agua al underflow Rw (%)")
    underflow_cp: float = Field(..., description="% sólidos en peso del underflow")
    sharpness_m: float = Field(default=2.4, description="Parámetro de agudeza m")
    feed_fines_pct_minus_38um: Optional[float] = Field(default=25.0, description="% pasante de 38 µm en la alimentación")

class CycloneCalibrationResponse(BaseModel):
    kd0: float = Field(..., description="Constante de corte d50c KD0")
    kq0: float = Field(..., description="Constante de capacidad KQ0")
    kw0: float = Field(..., description="Constante de agua al underflow KW0")
    ks0: float = Field(..., description="Constante de concentración underflow KS0")
    sharpness_m: float
    
    lambda_param: float = Field(..., description="Parámetro reológico / corrección de finos lambda")
    feed_cv_pct: float = Field(..., description="Concentración volumétrica de sólidos alimentación (%)")
    underflow_cv_pct: float = Field(..., description="Concentración volumétrica de sólidos underflow (%)")
    q_per_cyclone_m3h: float = Field(..., description="Caudal de pulpa por ciclón individual (m³/h)")
    roping_risk: bool = Field(..., description="Indica si el underflow está en riesgo de acodalamiento")
    roping_message: str


# ==============================================================================
# SUBMÓDULO 3: SIMULACIÓN / MODELAMIENTO PREDICTIVO (FORWARD SIMULATION)
# ==============================================================================

class FeedPsdItem(BaseModel):
    size_microns: float
    passing_pct: float

class CycloneSimulationRequest(BaseModel):
    geometry: CycloneGeometry
    num_cyclones_operating: int = Field(default=4, description="Número de ciclones operando en batería")
    
    # Parámetros calibrados de Nageswararao (del Submódulo 2)
    kd0: float = Field(default=0.00045, description="Constante de corte KD0")
    kq0: float = Field(default=500.0, description="Constante de capacidad KQ0")
    kw0: float = Field(default=10.0, description="Constante de recuperación de agua KW0")
    ks0: float = Field(default=0.05, description="Constante de descarga KS0")
    sharpness_m: float = Field(default=2.4, description="Parámetro de agudeza m")
    
    # Nuevas condiciones de alimentación
    feed_solids_tph: float = Field(..., description="Flujo másico total de sólidos (t/h)")
    feed_cp: float = Field(..., description="% de sólidos en peso de la alimentación (Cp)")
    solids_density: float = Field(default=2.7, description="Gravedad específica de sólidos (t/m³)")
    liquid_density: float = Field(default=1.0, description="Densidad del líquido (t/m³)")
    
    # Distribución granulométrica de la nueva alimentación
    feed_psd: List[FeedPsdItem]
    
    # Modo de cálculo de presión: se puede especificar la presión fija (calcula número óptimo) o predecirla
    fixed_pressure_kpa: Optional[float] = Field(default=None, description="Presión fija de manifold en kPa (opcional)")

class CycloneSimulationResponse(BaseModel):
    # Variables operacionales predichas
    predicted_pressure_kpa: float
    total_slurry_flow_m3h: float
    q_per_cyclone_m3h: float
    d50c_microns: float
    water_recovery_underflow_pct: float
    solids_split_underflow_pct: float
    solids_split_overflow_pct: float
    underflow_cv_pct: float
    underflow_cp_pct: float
    overflow_cp_pct: float
    
    # Alerta de régimen de descarga (Spray vs Roping)
    is_roping: bool
    roping_severity: str  # "Normal (Spray)", "Transición Crítica", "Acodalamiento Severo (Roping)"
    
    # Balance de corrientes
    feed: StreamProperties
    overflow: StreamProperties
    underflow: StreamProperties
    
    # Distribuciones granulométricas predichas
    predicted_sieves: List[Dict[str, Any]]
