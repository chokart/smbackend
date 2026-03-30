from pydantic import BaseModel
from typing import List, Optional

class SieveEntry(BaseModel):
    mesh_size: float
    weight_feed: float
    weight_overflow: float
    weight_underflow: float

class HydrocycloneGeometry(BaseModel):
    Dc: float
    Di: float
    Do: float
    Du: float

class RaoLynchRequest(BaseModel):
    sieves: List[SieveEntry]
    pan_feed: float
    pan_overflow: float
    pan_underflow: float
    geometry: HydrocycloneGeometry
    pressure: float
    solid_density: float
    feed_p_solids: float
    # Capacidad experimental medida (opcional, para calibrar K1)
    exp_capacity_Q: Optional[float] = None
    # Constantes manuales (si el usuario ya las conoce)
    k1: Optional[float] = None
    k3: Optional[float] = None

class PartitionCurvePoint(BaseModel):
    size: float
    actual_recovery: float
    corrected_recovery: float

class GranulometryPoint(BaseModel):
    size: float
    feed_passing: float
    overflow_passing: float
    underflow_passing: float

class BalanceRow(BaseModel):
    size: str
    feed_pct: float
    overflow_pct: float
    underflow_pct: float
    recovery_underflow: float

class RaoLynchResponse(BaseModel):
    # Resultados clave
    d50c_experimental: float
    d50c_predicted: float
    water_bypass_Rf: float
    capacity_Q: float
    
    # Constantes de Calibración calculadas
    k1_calculated: float
    k3_calculated: float
    
    # Datos para gráficos
    partition_curve: List[PartitionCurvePoint]
    granulometry_curve: List[GranulometryPoint]
    
    # Tabla de balance
    balance_table: List[BalanceRow]
    
    # Totales
    summary: dict
