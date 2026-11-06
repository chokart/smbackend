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

class SieveResidualInfo(BaseModel):
    size_microns: float
    mesh_label: Optional[str] = None
    feed_measured: float
    composite_measured: float
    residual_pct: float            # |f_i - (alpha*o_i + beta*u_i)| en % pasante
    delta_feed_pct: float          # |feed_real - feed_bal|
    delta_over_pct: float          # |over_real - over_bal|
    delta_under_pct: float         # |under_real - under_bal|
    total_delta_pct: float         # suma de correcciones |Δf| + |Δo| + |Δu|
    worst_stream: str              # 'Alimentación' | 'Rebose' | 'Descarga'
    worst_delta_pct: float         # Mayor corrección individual entre las 3 corrientes
    rank: int
    severity: str                  # 'high' (>= 0.5%), 'medium' (0.2% - 0.5%), 'low' (< 0.2%)

class ReconciledSieveRow(BaseModel):
    size_microns: float
    mesh_label: Optional[str] = None
    feed_passing: float
    overflow_passing: float
    underflow_passing: float
    reconciled_feed_passing: float
    reconciled_overflow_passing: float
    reconciled_underflow_passing: float
    composite_passing: Optional[float] = None
    residual_pct: Optional[float] = None     # Desbalance bruto inicial (%)
    delta_feed: Optional[float] = None       # Ajuste reconciliado al Feed (%)
    delta_over: Optional[float] = None       # Ajuste reconciliado al Rebose (%)
    delta_under: Optional[float] = None      # Ajuste reconciliado a la Descarga (%)
    worst_stream: Optional[str] = None       # Corriente con mayor ajuste
    rank: Optional[int] = None               # Ranking de desbalance (1 = mayor)
    severity: Optional[str] = None           # 'high' | 'medium' | 'low'
    partition_raw: Optional[float] = None    # Eficiencia calculada con datos iniciales/crudos (%)
    partition_actual: float                  # Eficiencia calculada con datos balanceados/reconciliados (%)
    partition_reduced: float                 # Eficiencia corregida/reducida balanceada (%)
    partition_fitted: float                  # Curva ajustada por modelo (Lynch-Rao) (%)

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
    top_errors: Optional[List[SieveResidualInfo]] = None
    mean_squared_error: float


# ==============================================================================
# SUBMÓDULO 2: CALIBRACIÓN DEL MODELO DE NAGESWARARAO (PARAMETER ESTIMATION)
# ==============================================================================

class CycloneGeometry(BaseModel):
    dc_mm: float = Field(default=508.0, description="Diámetro interno del hidrociclón Dc (mm) - Krebs GMAX20 (20 pulg)")
    di_mm: float = Field(default=140.0, description="Diámetro equivalente de entrada Di (mm) - Entrada Victaulic 8 pulg")
    do_mm: float = Field(default=171.5, description="Diámetro del vortex finder Do (mm) - Inserto 6.75 pulg")
    du_mm: float = Field(default=101.6, description="Diámetro del apex / tobera Du (mm) - Inserto 4.0 pulg")
    lc_mm: float = Field(default=1086.0, description="Longitud de sección cilíndrica Lc (mm) - Plano 42 3/4 pulg")
    theta_deg: float = Field(default=22.5, description="Ángulo total del cono en grados (°) - Plano 2 x 11.25°")

class CycloneCalibrationRequest(BaseModel):
    geometry: CycloneGeometry
    num_cyclones_operating: int = Field(default=6, description="Número de ciclones operando en batería (Tag 2102-CY-001/006)")
    
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
    sieves: Optional[List[Dict[str, Any]]] = Field(default=None, description="Granulometría balanceada del Submódulo 1 para verificación")

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
    
    # Comparativa detallada: Datos Balanceados vs Modelo Nageswararao
    comparison_kpis: Optional[List[Dict[str, Any]]] = None
    modeled_sieves: Optional[List[Dict[str, Any]]] = None


# ==============================================================================
# SUBMÓDULO 3: SIMULACIÓN / MODELAMIENTO PREDICTIVO (FORWARD SIMULATION)
# ==============================================================================

class FeedPsdItem(BaseModel):
    size_microns: float
    passing_pct: float

class CycloneSimulationRequest(BaseModel):
    geometry: CycloneGeometry
    num_cyclones_operating: int = Field(default=6, description="Número de ciclones operando en batería")
    
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
