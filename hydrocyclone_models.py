from pydantic import BaseModel
from typing import List, Optional

class SieveEntry(BaseModel):
    mesh_size: float
    weight_feed: float
    weight_overflow: float
    weight_underflow: float

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

class PartitionCurvePoint(BaseModel):
    size: float
    adjusted_recovery_mesh: float
    adjusted_recovery_solids: float

class GranulometryPoint(BaseModel):
    size: float
    feed_passing: float
    overflow_passing: float
    underflow_passing: float

class HydrocycloneMetrics(BaseModel):
    d50: float
    d50c: float
    bypass_rf: float
    solids_recovery_s: float

class BalanceRow(BaseModel):
    size: str
    feed_pct_real: float
    overflow_pct_real: float
    underflow_pct_real: float
    
    feed_pct_mesh: float
    overflow_pct_mesh: float
    underflow_pct_mesh: float
    recovery_mesh: float
    
    feed_pct_solids: float
    overflow_pct_solids: float
    underflow_pct_solids: float
    recovery_solids: float

class TrompParameters(BaseModel):
    d25c: float
    d50c: float
    d75c: float

class FlowComparison(BaseModel):
    pulp_mesh: float
    solids_mesh: float
    water_mesh: float
    pulp_solids: float
    solids_solids: float
    water_solids: float

class GlobalFlowBalance(BaseModel):
    feed: FlowComparison
    overflow: FlowComparison
    underflow: FlowComparison

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
    metrics_mesh: Optional[HydrocycloneMetrics] = None
    metrics_solids: Optional[HydrocycloneMetrics] = None
    tromp_mesh: Optional[TrompParameters] = None
    water_balance: WaterBalance
    partition_curve: List[PartitionCurvePoint]
    granulometry_curve: List[GranulometryPoint]
    comparison_table: List[BalanceRow]
    global_flow_balance: Optional[GlobalFlowBalance] = None
    summary: dict