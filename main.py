from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from models import ReconciliationRequest, ReconciliationResponse
from reconciliation import reconcile_data_qp, _get_balance_system, format_equations_to_strings
from water_models import WaterReconciliationRequest, WaterReconciliationResponse
from water_reconciliation import reconcile_water_data, _get_water_balance_system, format_water_equations
from hydrocyclone_models import HydrocycloneAnalysisRequest, HydrocycloneAnalysisResponse
from hydrocyclone_logic import analyze_hydrocyclone

app = FastAPI()

# Configuración de CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Permite todas las origenes
    allow_credentials=True,
    allow_methods=["*"],  # Permite todos los métodos
    allow_headers=["*"],  # Permite todas las cabeceras
)

@app.post("/reconcile", response_model=ReconciliationResponse)
def reconcile(request: ReconciliationRequest):
    try:
        reconciled_flows, report = reconcile_data_qp(request)
        return ReconciliationResponse(flows=reconciled_flows, report=report)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        # Captura de errores inesperados durante la optimización
        raise HTTPException(status_code=500, detail=f"Error interno en el servidor: {e}")

@app.post("/reconcile/water", response_model=WaterReconciliationResponse)
def reconcile_water(request: WaterReconciliationRequest):
    try:
        reconciled_flows, report = reconcile_water_data(request)
        return WaterReconciliationResponse(flows=reconciled_flows, report=report)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error interno en el servidor: {e}")

@app.post("/preview-equations")
def preview_equations(request: ReconciliationRequest):
    try:
        A, b, var_map, _, _, _ = _get_balance_system(request)
        equations = format_equations_to_strings(A, b, var_map)
        return {"equations": equations}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error al generar ecuaciones: {e}")

@app.post("/preview-water-equations")
def preview_water_equations(request: WaterReconciliationRequest):
    try:
        A, b, var_map, _, _, _ = _get_water_balance_system(request)
        equations = format_water_equations(A, b, var_map)
        return {"equations": equations}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error al generar ecuaciones de agua: {e}")

@app.post("/model/hydrocyclone/rao-lynch", response_model=HydrocycloneAnalysisResponse)
def hydrocyclone_analysis(request: HydrocycloneAnalysisRequest):
    try:
        return analyze_hydrocyclone(request)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error en el análisis del hidrociclón: {e}")

@app.get("/")
def read_root():
    return {"message": "Servidor de Balance Metalúrgico funcionando"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
