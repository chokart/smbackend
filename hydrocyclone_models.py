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
    actual_recovery: float # Experimental
    corrected_recovery: float # Reconciliado Corregido
    adjusted_recovery: Optional[float] = None # Reconciliado Ajustado
    solids_recovery: Optional[float] = None # Balance por Sólidos

class GranulometryPoint(BaseModel):
    size: float
    # Experimentales
    feed_passing: float
    overflow_passing: float
    underflow_passing: float
    # Ajustados (Reconciliado)
    feed_passing_adj: Optional[float] = None
    overflow_passing_adj: Optional[float] = None
    underflow_passing_adj: Optional[float] = None
    # Por Sólidos
    feed_passing_sol: Optional[float] = None
    overflow_passing_sol: Optional[float] = None
    underflow_passing_sol: Optional[float] = None

class BalanceRow(BaseModel):
    size: str
    # Pesos (solo para tabla experimental)
    feed_w: Optional[float] = None
    overflow_w: Optional[float] = None
    underflow_w: Optional[float] = None
    
    # % Retenidos (Parciales)
    feed_pct: float
    overflow_pct: float
    underflow_pct: float
    
    # % Pasantes (Acumulados)
    feed_pass: Optional[float] = None
    overflow_pass: Optional[float] = None
    underflow_pass: Optional[float] = None
    
    # Datos ajustados (opcionales para balances específicos)
    feed_pct_adj: Optional[float] = None
    overflow_pct_adj: Optional[float] = None
    underflow_pct_adj: Optional[float] = None
    feed_pass_adj: Optional[float] = None
    overflow_pass_adj: Optional[float] = None
    underflow_pass_adj: Optional[float] = None
    
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
    
    # Tablas de balance
    balance_table: List[BalanceRow] # Tabla reconciliada (Ajustada)
    balance_solids_table: Optional[List[BalanceRow]] = None
    
    # Resumen y mensajes
    summary: dict
