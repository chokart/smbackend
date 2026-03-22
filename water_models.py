from pydantic import BaseModel
from typing import List, Optional
from enum import Enum

class NodeType(str, Enum):
    ALIMENTO = "initial"
    PROCESO = "process"
    CONCENTRADO = "concentrate"
    RELAVE = "tailing"
    PRODUCTO = "product"
    AGUA = "water_source" # Opcional: para fuentes de agua pura

class WaterFlowData(BaseModel):
    # En balance de agua, solemos medir Flujo Volumétrico (m3/h) o Másico (tph) de pulpa
    # y el Porcentaje de Sólidos (%S).
    # Asumiremos densidad del agua = 1 para simplificar (1 m3 = 1 ton)
    
    tonelaje_pulpa: float  # Flujo másico total (Agua + Sólido)
    tonelaje_pulpa_error: float # Desviación estándar
    tonelaje_pulpa_fixed: bool = False
    not_measured: bool = False
    
    porcentaje_solidos: float # % de sólidos en peso (Cp)
    porcentaje_solidos_error: float
    porcentaje_solidos_fixed: bool = False
    porcentaje_solidos_not_measured: bool = False

class WaterFlow(BaseModel):
    id: str
    source: str
    target: str
    data: WaterFlowData

class WaterNode(BaseModel):
    id: str
    label: str
    node_type: NodeType

class WaterReconciliationRequest(BaseModel):
    nodes: List[WaterNode]
    flows: List[WaterFlow]

class CorrectedWaterFlowData(BaseModel):
    tonelaje_pulpa_corregido: float
    tonelaje_pulpa_error: float
    
    porcentaje_solidos_corregido: float
    porcentaje_solidos_error: float
    
    # Datos calculados útiles
    tonelaje_agua_calculado: float
    tonelaje_solido_calculado: float

class ReconciledWaterFlow(BaseModel):
    id: str
    corrected_data: CorrectedWaterFlowData

class WaterReconciliationResponse(BaseModel):
    flows: List[ReconciledWaterFlow]
    report: Optional[dict] = None
