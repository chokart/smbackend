from pydantic import BaseModel
from typing import List, Optional

class SieveEntry(BaseModel):
    mesh_size: float
    weight_feed: float
    weight_overflow: float
    weight_underflow: float

class HydrocycloneAnalysisRequest(BaseModel):
    sieves: List[SieveEntry]
    pan_feed: float
    pan_overflow: float
    pan_underflow: float
    
    # Datos de referencia
    pressure: Optional[float] = None
    solid_density: float = 2.65
    # Porcentajes de sólidos (opcionales pero recomendados para balance de agua)
    feed_p_solids: Optional[float] = None
    overflow_p_solids: Optional[float] = None
    underflow_p_solids: Optional[float] = None
    
    # Flujo de alimento (Opcional para balance absoluto)
    feed_flow_rate: Optional[float] = None
    feed_flow_unit: Optional[str] = "tph" # "tph" (masa sólidos) o "m3h" (volumen pulpa)

class PartitionCurvePoint(BaseModel):
    size: float
    actual_recovery: float
    corrected_recovery: float
    adjusted_recovery: Optional[float] = None

class GranulometryPoint(BaseModel):
    size: float
    feed_passing: float
    overflow_passing: float
    underflow_passing: float
    # Versiones ajustadas
    feed_passing_adj: Optional[float] = None
    overflow_passing_adj: Optional[float] = None
    underflow_passing_adj: Optional[float] = None

class BalanceRow(BaseModel):
    size: str
    # Datos experimentales (parciales %)
    feed_pct: float
    overflow_pct: float
    underflow_pct: float
    # Datos ajustados (parciales %)
    feed_pct_adj: Optional[float] = None
    overflow_pct_adj: Optional[float] = None
    underflow_pct_adj: Optional[float] = None
    
    recovery_underflow: float # Eficiencia Ea

class TrompParameters(BaseModel):
    d25c: float
    d50c: float
    d75c: float
    imperfection: float # (d75c - d25c) / (2 * d50c)

class FlowData(BaseModel):
    mass_solids: float
    mass_water: float
    vol_solids: float
    vol_water: float
    vol_pulp: float
    p_solids: float

class GlobalBalance(BaseModel):
    feed: FlowData
    overflow: FlowData
    underflow: FlowData

class WaterBalance(BaseModel):
    solids_recovery_S: float
    water_recovery_Rw: Optional[float] = None
    bypass_Rf: float
    pan_recovery_Epan: Optional[float] = None
    consistency_error: Optional[float] = None
    
    # Flujos relativos (base 100)
    feed_flow: float
    overflow_flow: float
    underflow_flow: float
    
    # Balance global absoluto (si se provee flujo)
    global_balance: Optional[GlobalBalance] = None

class HydrocycloneAnalysisResponse(BaseModel):
    # Resultados (Ajustados si es posible)
    d50c_experimental: float
    d50c_adjusted: Optional[float] = None
    
    # Parámetros de Tromp
    tromp: Optional[TrompParameters] = None
    
    # Diagnóstico y Alertas
    diagnosis_message: Optional[str] = None
    diagnosis_level: Optional[str] = "info"
    
    # Balance de Masa y Agua
    water_balance: WaterBalance
    
    # Gráficos
    partition_curve: List[PartitionCurvePoint]
    granulometry_curve: List[GranulometryPoint]
    
    # Tabla de balance
    balance_table: List[BalanceRow]
    
    # Resumen y mensajes
    summary: dict
