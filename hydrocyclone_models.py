from pydantic import BaseModel
from typing import List, Optional

class SieveEntry(BaseModel):
    mesh_size: float  # Tamaño en micrones
    weight_feed: float
    weight_overflow: float
    weight_underflow: float

class HydrocycloneGeometry(BaseModel):
    Dc: float  # Diámetro ciclón (mm)
    Di: float  # Diámetro entrada (mm)
    Do: float  # Diámetro vortex (mm)
    Du: float  # Diámetro apex (mm)

class RaoLynchRequest(BaseModel):
    # Datos experimentales de mallas
    sieves: List[SieveEntry]
    
    # Datos del hidrociclón
    geometry: HydrocycloneGeometry
    
    # Parámetros de operación
    pressure: float      # Presión (kPa)
    solid_density: float # g/cm3 (default 2.65)
    feed_p_solids: float # % sólidos en peso en el alimento
    
    # Parámetros de ajuste (K1, K2, K3, K4) 
    # Si no se envían, usaremos los valores estándar de Rao & Lynch
    k1: Optional[float] = None
    k2: Optional[float] = None
    k3: Optional[float] = None
    k4: Optional[float] = None

class PartitionCurvePoint(BaseModel):
    size: float
    actual_recovery: float    # Recuperación real
    corrected_recovery: float # Recuperación corregida (E_c)

class RaoLynchResponse(BaseModel):
    # Resultados del modelo
    d50c_predicted: float
    d50c_experimental: float
    water_bypass_Rf: float
    capacity_Q: float
    
    # Curva de partición para el gráfico
    partition_curve: List[PartitionCurvePoint]
    
    # Reporte de granulometría calculada (pasantes, retenidos %)
    granulometry_report: dict 
