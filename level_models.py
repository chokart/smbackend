from pydantic import BaseModel, Field
from typing import Optional

class LevelReadingCreate(BaseModel):
    device_id: Optional[str] = Field(default="esp32_tanque_1", description="Identificador del dispositivo ESP32")
    distance_cm: float = Field(..., description="Distancia medida al agua en centímetros")
    level_cm: float = Field(..., description="Nivel estimado del agua en centímetros")
    percentage: float = Field(..., description="Porcentaje de llenado del tanque (0-100%)")

class LevelReadingResponse(BaseModel):
    id: int
    device_id: str
    distance_cm: float
    level_cm: float
    percentage: float
    created_at: str
    is_online: bool = True
    seconds_ago: Optional[float] = None
