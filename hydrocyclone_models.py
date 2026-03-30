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

class PartitionCurvePoint(BaseModel):
    size: float
    actual_recovery: float
    corrected_recovery: float
    adjusted_recovery: Optional[float] = None # Usando datos reconciliados

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

class WaterBalance(BaseModel):
    solids_recovery_S: float
    water_recovery_Rw: Optional[float] = None
    bypass_Rf: float
    # Flujos calculados (basados en 100 unidades de alimento o pesos totales)
    feed_flow: float
    overflow_flow: float
    underflow_flow: float

class HydrocycloneAnalysisResponse(BaseModel):
    # Resultados (Ajustados si es posible)
    d50c_experimental: float
    d50c_adjusted: Optional[float] = None
    
    # Balance de Masa y Agua
    water_balance: WaterBalance
    
    # Gráficos
    partition_curve: List[PartitionCurvePoint]
    granulometry_curve: List[GranulometryPoint]
    
    # Tabla de balance
    balance_table: List[BalanceRow]
    
    # Resumen y mensajes
    summary: dict
