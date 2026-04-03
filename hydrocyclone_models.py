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
    pressure: Optional[float] = None
    solid_density: float = 2.65
    feed_p_solids: Optional[float] = None
    overflow_p_solids: Optional[float] = None
    underflow_p_solids: Optional[float] = None
    feed_flow_rate: Optional[float] = None
    feed_flow_unit: Optional[str] = "tph"

class PartitionCurvePoint(BaseModel):
    size: float
    actual_recovery: float
    corrected_recovery: float
    adjusted_recovery: Optional[float] = None
    solids_recovery: Optional[float] = None

class GranulometryPoint(BaseModel):
    size: float
    feed_passing: float
    overflow_passing: float
    underflow_passing: float
    feed_passing_adj: Optional[float] = None
    overflow_passing_adj: Optional[float] = None
    underflow_passing_adj: Optional[float] = None
    feed_passing_sol: Optional[float] = None
    overflow_passing_sol: Optional[float] = None
    underflow_passing_sol: Optional[float] = None

class BalanceRow(BaseModel):
    size: str
    feed_w: Optional[float] = None
    overflow_w: Optional[float] = None
    underflow_w: Optional[float] = None
    
    # Pesos Re-calculados (Ajustados)
    feed_w_adj: Optional[float] = None
    overflow_w_adj: Optional[float] = None
    underflow_w_adj: Optional[float] = None
    
    # Pesos por Sólidos
    feed_w_sol: Optional[float] = None
    overflow_w_sol: Optional[float] = None
    underflow_w_sol: Optional[float] = None
    
    feed_pct: float
    overflow_pct: float
    underflow_pct: float
    feed_pass: Optional[float] = None
    overflow_pass: Optional[float] = None
    underflow_pass: Optional[float] = None
    
    # Datos Ajustados (para balance reconciliado)
    feed_pct_adj: Optional[float] = None
    overflow_pct_adj: Optional[float] = None
    underflow_pct_adj: Optional[float] = None
    feed_pass_adj: Optional[float] = None
    overflow_pass_adj: Optional[float] = None
    underflow_pass_adj: Optional[float] = None
    
    recovery_underflow: float
    recovery_corrected: Optional[float] = None

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
    water_recovery_Rw: Optional[float] = None
    bypass_Rf: float
    pan_recovery_Epan: Optional[float] = None
    consistency_error: Optional[float] = None
    feed_flow: float
    overflow_flow: float
    underflow_flow: float
    global_balance: Optional[GlobalBalance] = None
    global_balance_solids: Optional[GlobalBalance] = None

class HydrocycloneMetrics(BaseModel):
    d50: float
    d50c: float
    bypass_rf: float
    solids_recovery_s: float

class HydrocycloneAnalysisResponse(BaseModel):
    reconciled_metrics: Optional[HydrocycloneMetrics] = None
    solids_metrics: Optional[HydrocycloneMetrics] = None
    d50c_experimental: float
    d50c_adjusted: Optional[float] = None
    tromp: Optional[TrompParameters] = None
    water_balance: WaterBalance
    partition_curve: List[PartitionCurvePoint]
    granulometry_curve: List[GranulometryPoint]
    balance_table: List[BalanceRow]
    balance_solids_table: Optional[List[BalanceRow]] = None
    summary: dict
