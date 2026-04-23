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
        
    return var_map, np.array(initial_guess, dtype=float), bounds

def _get_balance_nodes(request: WaterReconciliationRequest):
    """Identifica nodos que deben cumplir balance (tienen entradas y salidas)."""
    # Somos más flexibles: cualquier nodo que no sea explícitamente una fuente pura (initial/water_source) 
    # y tenga flujos de entrada y salida debe balancearse.
    balance_nodes = []
    for node in request.nodes:
        has_in = any(f.target == node.id for f in request.flows)
        has_out = any(f.source == node.id for f in request.flows)
        # Si tiene ambos, es un nodo de proceso/mezcla/división
        if has_in and has_out:
            balance_nodes.append(node)
        # O si es explícitamente de tipo proceso
        elif node.node_type == NodeType.PROCESO:
            balance_nodes.append(node)
    return balance_nodes

def _get_constraints(y, request: WaterReconciliationRequest, var_map: dict):
    """Calcula el residuo de las restricciones de balance (No Lineales)."""
    balance_nodes = _get_balance_nodes(request)
    residuals = []
    
    for node in balance_nodes:
        total_mass_balance = 0.0
        solids_mass_balance = 0.0
        
        for flow in request.flows:
            q_idx = var_map[f"q_{flow.id}"]
            cp_idx = var_map[f"cp_{flow.id}"]
            Q = y[q_idx]
            Cp = y[cp_idx]
            
            factor = 1 if flow.target == node.id else (-1 if flow.source == node.id else 0)
            if factor != 0:
                total_mass_balance += factor * Q
                solids_mass_balance += factor * (Q * Cp / 100.0)
        
        residuals.append(total_mass_balance)
        residuals.append(solids_mass_balance)
        
    return np.array(residuals)

def _get_jacobian(y, request: WaterReconciliationRequest, var_map: dict):
    """Calcula la matriz Jacobiana de las restricciones."""
    balance_nodes = _get_balance_nodes(request)
    num_vars = len(y)
    num_constraints = len(balance_nodes) * 2
    J = np.zeros((num_constraints, num_vars))
    
    row = 0
    for node in balance_nodes:
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
    
    measured_indices, measured_values, variances = [], [], []
    for flow in request.flows:
        # Pulpa
        if not flow.data.not_measured and not flow.data.tonelaje_pulpa_fixed:
            measured_indices.append(var_map[f"q_{flow.id}"])
            measured_values.append(flow.data.tonelaje_pulpa)
            err = flow.data.tonelaje_pulpa_error if flow.data.tonelaje_pulpa_error > 0 else (flow.data.tonelaje_pulpa * 0.02 + 0.1)
            variances.append(err**2)
            
        # % Sólidos
        if not flow.data.porcentaje_solidos_not_measured and not flow.data.porcentaje_solidos_fixed:
            measured_indices.append(var_map[f"cp_{flow.id}"])
            measured_values.append(flow.data.porcentaje_solidos)
            err = flow.data.porcentaje_solidos_error if flow.data.porcentaje_solidos_error > 0 else 0.5
            variances.append(err**2)

    if not measured_indices:
        return _build_response(request, y_init, np.zeros_like(y_init), var_map)

    W = np.diag(1.0 / np.array(variances))
    m_vals = np.array(measured_values)
    
    def objective(y):
        diff = y[measured_indices] - m_vals
        return 0.5 * diff.T @ W @ diff

    def objective_grad(y):
        grad = np.zeros_like(y)
        diff = y[measured_indices] - m_vals
        grad[measured_indices] = W @ diff
        return grad

    constraints = {'type': 'eq', 'fun': _get_constraints, 'jac': _get_jacobian, 'args': (request, var_map)}
    
    res = minimize(
        objective, y_init, jac=objective_grad, method='SLSQP', 
        bounds=bounds, constraints=constraints, tol=1e-8, 
        options={'maxiter': 200, 'ftol': 1e-9}
    )
    
    y_hat = res.x
    
    # Incertidumbre
    V_in = np.zeros((len(y_hat), len(y_hat)))
    for i, idx in enumerate(measured_indices):
        V_in[idx, idx] = variances[i]
    
    J = _get_jacobian(y_hat, request, var_map)
    try:
        if J.size > 0:
            S = J @ V_in @ J.T
            V_out = V_in - V_in @ J.T @ np.linalg.pinv(S) @ J @ V_in
            std_devs = np.sqrt(np.maximum(np.diag(V_out), 1e-9))
        else:
            std_devs = np.sqrt(np.diag(V_in))
    except:
        std_devs = np.zeros_like(y_hat)

    return _build_response(request, y_hat, std_devs, var_map, success=res.success)

def _build_response(request, y_hat, std_devs, var_map, success=True):
    reconciled_flows = []
    for flow in request.flows:
        q_idx, cp_idx = var_map[f"q_{flow.id}"], var_map[f"cp_{flow.id}"]
        q_val, cp_val = y_hat[q_idx], y_hat[cp_idx]
        
        reconciled_flows.append({
            "id": flow.id,
            "corrected_data": {
                "tonelaje_pulpa_corregido": float(q_val),
                "tonelaje_pulpa_error": float(std_devs[q_idx]),
                "porcentaje_solidos_corregido": float(cp_val),
                "porcentaje_solidos_error": float(std_devs[cp_idx]),
                "tonelaje_solido_calculado": float(q_val * cp_val / 100.0),
                "tonelaje_agua_calculado": float(q_val * (1.0 - cp_val / 100.0))
            }
        })

    report = {
        "success": bool(success),
        "methodology": [
            "Optimización no lineal SLSQP con gradientes analíticos.",
            "Detección automática de nodos de balance basada en topología.",
            "Cálculo de varianza posterior mediante proyección de Gauss-Markov."
        ],
        "equations": _format_physical_equations(request)
    }
    return reconciled_flows, report

def _format_physical_equations(request: WaterReconciliationRequest):
    nodes = _get_balance_nodes(request)
    eqs = []
    for node in nodes:
        in_f = [f.id for f in request.flows if f.target == node.id]
        out_f = [f.id for f in request.flows if f.source == node.id]
        eqs.append(f"Nodo {node.label}: Balance de Pulpa (In: {in_f} = Out: {out_f})")
        eqs.append(f"Nodo {node.label}: Balance de Sólidos Secos")
    return eqs
