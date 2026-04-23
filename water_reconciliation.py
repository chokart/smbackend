import numpy as np
from scipy.optimize import minimize
from water_models import WaterReconciliationRequest, NodeType
from typing import List

def _get_variable_mapping(request: WaterReconciliationRequest):
    """Mapea cada flujo a índices en el vector de estado y genera valores iniciales."""
    var_map = {}
    initial_guess = []
    bounds = []
    var_idx = 0
    
    for flow in request.flows:
        # Variable: Tonelaje Pulpa (Q)
        var_map[f"q_{flow.id}"] = var_idx
        initial_guess.append(flow.data.tonelaje_pulpa)
        bounds.append((flow.data.tonelaje_pulpa, flow.data.tonelaje_pulpa) if flow.data.tonelaje_pulpa_fixed else (0, None))
        var_idx += 1
        
        # Variable: Porcentaje Sólidos (Cp)
        var_map[f"cp_{flow.id}"] = var_idx
        initial_guess.append(flow.data.porcentaje_solidos)
        bounds.append((flow.data.porcentaje_solidos, flow.data.porcentaje_solidos) if flow.data.porcentaje_solidos_fixed else (0, 100))
        var_idx += 1
        
    return var_map, np.array(initial_guess), bounds

def _get_constraints(y, request: WaterReconciliationRequest, var_map: dict):
    """Calcula el residuo de las restricciones de balance (No Lineales)."""
    process_nodes = [n for n in request.nodes if n.node_type == NodeType.PROCESO]
    residuals = []
    
    for node in process_nodes:
        total_mass_balance = 0.0
        solids_mass_balance = 0.0
        
        for flow in request.flows:
            q_idx = var_map[f"q_{flow.id}"]
            cp_idx = var_map[f"cp_{flow.id}"]
            
            Q = y[q_idx]
            Cp = y[cp_idx]
            
            factor = 0
            if flow.target == node.id: factor = 1
            elif flow.source == node.id: factor = -1
            
            if factor != 0:
                total_mass_balance += factor * Q
                solids_mass_balance += factor * (Q * Cp / 100.0)
        
        residuals.append(total_mass_balance)
        residuals.append(solids_mass_balance)
        
    return np.array(residuals)

def _get_jacobian(y, request: WaterReconciliationRequest, var_map: dict):
    """Calcula la matriz Jacobiana de las restricciones en el punto y."""
    process_nodes = [n for n in request.nodes if n.node_type == NodeType.PROCESO]
    num_vars = len(y)
    num_constraints = len(process_nodes) * 2
    J = np.zeros((num_constraints, num_vars))
    
    row = 0
    for node in process_nodes:
        # Balance Total
        for flow in request.flows:
            q_idx = var_map[f"q_{flow.id}"]
            factor = 1 if flow.target == node.id else (-1 if flow.source == node.id else 0)
            if factor != 0:
                J[row, q_idx] = factor
        row += 1
        
        # Balance Sólidos
        for flow in request.flows:
            q_idx = var_map[f"q_{flow.id}"]
            cp_idx = var_map[f"cp_{flow.id}"]
            Q = y[q_idx]
            Cp = y[cp_idx]
            factor = 1 if flow.target == node.id else (-1 if flow.source == node.id else 0)
            if factor != 0:
                J[row, q_idx] = factor * Cp / 100.0
                J[row, cp_idx] = factor * Q / 100.0
        row += 1
        
    return J

def reconcile_water_data(request: WaterReconciliationRequest):
    var_map, y_init, bounds = _get_variable_mapping(request)
    
    # Construcción de pesos y mediciones
    measured_indices = []
    measured_values = []
    variances = []
    
    for flow in request.flows:
        # Pulpa
        if not flow.data.not_measured and not flow.data.tonelaje_pulpa_fixed:
            idx = var_map[f"q_{flow.id}"]
            measured_indices.append(idx)
            measured_values.append(flow.data.tonelaje_pulpa)
            err = flow.data.tonelaje_pulpa_error if flow.data.tonelaje_pulpa_error > 0 else 1e-4
            variances.append(err**2)
            
        # % Sólidos
        if not flow.data.porcentaje_solidos_not_measured and not flow.data.porcentaje_solidos_fixed:
            idx = var_map[f"cp_{flow.id}"]
            measured_indices.append(idx)
            measured_values.append(flow.data.porcentaje_solidos)
            err = flow.data.porcentaje_solidos_error if flow.data.porcentaje_solidos_error > 0 else 1e-4
            variances.append(err**2)

    # Si no hay mediciones, solo devolvemos los valores iniciales (o error)
    if not measured_indices:
        return _build_response(request, y_init, np.zeros_like(y_init), var_map)

    # Matrices para optimización
    W = np.diag(1.0 / np.array(variances))
    
    def objective(y):
        diff = y[measured_indices] - measured_values
        return diff.T @ W @ diff

    # Restricciones de igualdad: residuals(y) = 0
    constraints = {'type': 'eq', 'fun': _get_constraints, 'args': (request, var_map), 'jac': _get_jacobian}
    
    # Optimización NO lineal directa
    res = minimize(objective, y_init, method='SLSQP', bounds=bounds, constraints=constraints, tol=1e-8, options={'maxiter': 100})
    
    y_hat = res.x
    
    # --- CÁLCULO DE INCERTIDUMBRE (COVARIANZA POSTERIOR) ---
    # Usando el método de multiplicadores de Lagrange / Proyección
    # V_out = V_in - V_in * J^T * (J * V_in * J^T)^-1 * J * V_in
    
    V_in_full = np.zeros((len(y_hat), len(y_hat)))
    for i, idx in enumerate(measured_indices):
        V_in_full[idx, idx] = variances[i]
    
    J = _get_jacobian(y_hat, request, var_map)
    
    try:
        # SVD o pseudo-inversa para estabilidad numérica en sistemas sobre-determinados
        S = J @ V_in_full @ J.T
        V_out = V_in_full - V_in_full @ J.T @ np.linalg.pinv(S) @ J @ V_in_full
        std_devs = np.sqrt(np.maximum(np.diag(V_out), 0))
    except Exception:
        std_devs = np.zeros_like(y_hat)

    return _build_response(request, y_hat, std_devs, var_map, success=res.success)

def _build_response(request, y_hat, std_devs, var_map, success=True):
    reconciled_flows = []
    for flow in request.flows:
        q_idx = var_map[f"q_{flow.id}"]
        cp_idx = var_map[f"cp_{flow.id}"]
        
        q_val = y_hat[q_idx]
        cp_val = y_hat[cp_idx]
        q_err = std_devs[q_idx]
        cp_err = std_devs[cp_idx]
        
        reconciled_flows.append({
            "id": flow.id,
            "corrected_data": {
                "tonelaje_pulpa_corregido": q_val,
                "tonelaje_pulpa_error": q_err,
                "porcentaje_solidos_corregido": cp_val,
                "porcentaje_solidos_error": cp_err,
                "tonelaje_solido_calculado": q_val * (cp_val / 100.0),
                "tonelaje_agua_calculado": q_val * (1.0 - cp_val / 100.0)
            }
        })

    report = {
        "success": success,
        "methodology": [
            "Optimización no lineal directa mediante algoritmo SLSQP.",
            "Cálculo de incertidumbre posterior mediante proyección de matriz de covarianza.",
            "Balance simultáneo de masa total y sólidos secos."
        ],
        "equations": _format_physical_equations(request, var_map)
    }
    return reconciled_flows, report

def _format_physical_equations(request: WaterReconciliationRequest, var_map: dict):
    process_nodes = [n for n in request.nodes if n.node_type == NodeType.PROCESO]
    eqs = []
    for node in process_nodes:
        # Masa Total
        in_flows = [f.id for f in request.flows if f.target == node.id]
        out_flows = [f.id for f in request.flows if f.source == node.id]
        eqs.append(f"Nodo {node.label} (Total): Σ Q[{', '.join(in_flows)}] = Σ Q[{', '.join(out_flows)}]")
        # Masa Sólida
        eqs.append(f"Nodo {node.label} (Sólidos): Σ (Q*Cp/100)_{{in}} = Σ (Q*Cp/100)_{{out}}")
    return eqs
