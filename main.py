import logging
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from models import ReconciliationRequest, ReconciliationResponse
from reconciliation import reconcile_data_qp, _get_balance_system, format_equations_to_strings
from water_models import WaterReconciliationRequest, WaterReconciliationResponse
from water_reconciliation import reconcile_water_data, _format_physical_equations
from pulp_models import PulpCalculationRequest, PulpCalculationResponse
from pulp_logic import calculate_pulp_properties

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("suiteminerals")

app = FastAPI(
    title="Suite Minerals API",
    version="1.0.0",
    description="Motor de cálculo, balance metalúrgico y reconciliación de datos",
)

# Configuración de CORS estandarizada para consumo web
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/health")
def health_check():
    """Endpoint de salud para balanceadores de carga y Dokploy."""
    return {"status": "ok"}

@app.get("/")
def read_root():
    return {"message": "Servidor de Balance Metalúrgico Suite Minerals funcionando", "version": "1.0.0"}

@app.post("/reconcile", response_model=ReconciliationResponse)
def reconcile(request: ReconciliationRequest):
    try:
        reconciled_flows, report = reconcile_data_qp(request)
        return ReconciliationResponse(flows=reconciled_flows, report=report)
    except ValueError as e:
        logger.warning(f"Error de validación en reconciliación: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Error inesperado durante la optimización:")
        raise HTTPException(status_code=500, detail=f"Error interno en el servidor: {e}")

@app.post("/reconcile/water", response_model=WaterReconciliationResponse)
def reconcile_water(request: WaterReconciliationRequest):
    try:
        reconciled_flows, report = reconcile_water_data(request)
        return WaterReconciliationResponse(flows=reconciled_flows, report=report)
    except ValueError as e:
        logger.warning(f"Error de validación en balance de agua: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Error inesperado en balance de agua:")
        raise HTTPException(status_code=500, detail=f"Error interno en el servidor: {e}")

@app.post("/preview-equations")
def preview_equations(request: ReconciliationRequest):
    try:
        A, b, var_map, _, _, _ = _get_balance_system(request)
        equations = format_equations_to_strings(A, b, var_map)
        return {"equations": equations}
    except Exception as e:
        logger.exception("Error al generar ecuaciones:")
        raise HTTPException(status_code=500, detail=f"Error al generar ecuaciones: {e}")

@app.post("/preview-water-equations")
def preview_water_equations(request: WaterReconciliationRequest):
    try:
        equations = _format_physical_equations(request)
        return {"equations": equations}
    except Exception as e:
        logger.exception("Error al generar ecuaciones de agua:")
        raise HTTPException(status_code=500, detail=f"Error al generar ecuaciones de agua: {e}")

@app.post("/pulp-calculate", response_model=PulpCalculationResponse)
def pulp_calculate(request: PulpCalculationRequest):
    try:
        return calculate_pulp_properties(request)
    except ValueError as e:
        logger.warning(f"Error en datos de pulpa: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Error en el cálculo de pulpa:")
        raise HTTPException(status_code=500, detail=f"Error en el cálculo de pulpa: {e}")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
