from pydantic import BaseModel
from typing import List, Optional

class SieveEntry(BaseModel):
    mesh_size: float
    weight_feed: float
    weight_overflow: float
    weight_underflow: float

class PlittGeometry(BaseModel):
    Dc: float
    Di: float
    Do: float
    Du: float
    h: float
    alpha: Optional[float] = 20.0

class PlittParameters(BaseModel):
    a1: float = 1.0
    a2: float = 1.0
    a3: float = 1.0
    a4: float = 1.0
    l_const: float = 1.0

class HydrocycloneAnalysisRequest(BaseModel):
    mode: Optional[str] = "analysis" 
    sieves: List[SieveEntry]
    pan_feed: float
    pan_overflow: float
    pan_underflow: float
    pressure: Optional[float] = None
    solid_density: float = 2.65
    liquid_density: float = 1.0
    feed_p_solids: Optional[float] = None
    overflow_p_solids: Optional[float] = None
    underflow_p_solids: Optional[float] = None
    feed_flow_rate: Optional[float] = None
    feed_flow_unit: Optional[str] = "tph"
    geometry: Optional[PlittGeometry] = None
    plitt_params: Optional[PlittParameters] = PlittParameters()

class PartitionCurvePoint(BaseModel):
    size: float
    actual_recovery: float
    corrected_recovery: float
    adjusted_recovery: Optional[float] = None
    plitt_recovery: Optional[float] = None

class GranulometryPoint(BaseModel):
    size: float
    feed_passing: float
    overflow_passing: float
    underflow_passing: float
    feed_passing_adj: Optional[float] = None
    overflow_passing_adj: Optional[float] = None
    underflow_passing_adj: Optional[float] = None

class HydrocycloneMetrics(BaseModel):
    d50: float
    d50c: float
    bypass_rf: float
    solids_recovery_s: float
    m_plitt: Optional[float] = None
    s_plitt: Optional[float] = None
    p_plitt: Optional[float] = None

class BalanceRow(BaseModel):
    size: str
    feed_pct_real: float      # % Reconciliado (o experimental)
    overflow_pct_real: float
    underflow_pct_real: float
    feed_pct_sim: float       # % Simulado por Plitt
    overflow_pct_sim: float
    underflow_pct_sim: float
    recovery_real: float      # Ea Reconciliada
    recovery_sim: float       # Ea Simulada (Plitt)

class TrompParameters(BaseModel):
    d25c: float
    d50c: float
    d75c: float

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
    bypass_Rf: float
    feed_flow: float
    overflow_flow: float
    underflow_flow: float
    global_balance: Optional[GlobalBalance] = None

class HydrocycloneAnalysisResponse(BaseModel):
    reconciled_metrics: Optional[HydrocycloneMetrics] = None
    d50c_experimental: float
    d50c_adjusted: Optional[float] = None
    tromp: Optional[TrompParameters] = None
    water_balance: WaterBalance
    partition_curve: List[PartitionCurvePoint]
    granulometry_curve: List[GranulometryPoint]
    comparison_table: List[BalanceRow]
    summary: dict
