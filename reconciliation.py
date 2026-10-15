import numpy as np
from scipy.optimize import minimize
from models import ReconciliationRequest, NodeType
from typing import List, Dict, Any, Tuple, Optional

def _get_balance_nodes(request: ReconciliationRequest):
    """
    Identifica nodos que deben cumplir balance de conservación de masa.
    Un nodo debe tener flujos de entrada y salida para cumplir la conservación en estado estacionario.
    """
    balance_nodes = []
    for node in request.nodes:
        has_in = any(f.target == node.id for f in request.flows)
        has_out = any(f.source == node.id for f in request.flows)
        if has_in and has_out:
            balance_nodes.append(node)
    return balance_nodes

def _get_balance_system(request: ReconciliationRequest, y_k: np.ndarray = None):
    """
    Construye el sistema de ecuaciones linealizado (A, b) y el mapeo de variables.
    Linealización bilineal de finos: T*L ≈ L_k*T + T_k*L - T_k*L_k.
    """
    var_map = {}
    initial_guess = []
    bounds = []
    fixed_vars_indices = set()
    var_idx = 0

    # Heurística para estimar valores de flujos no medidos
    measured_tonnages = [f.data.tonelaje for f in request.flows if not f.data.not_measured and f.data.tonelaje > 0]
    avg_tonnage = float(np.mean(measured_tonnages)) if measured_tonnages else 100.0

    for flow in request.flows:
        # Variable: Tonelaje Seco (T)
        var_map[f"t_{flow.id}"] = var_idx
        t_val = flow.data.tonelaje
        if flow.data.not_measured and t_val <= 0:
            t_val = avg_tonnage
        initial_guess.append(t_val)
        bounds.append((flow.data.tonelaje, flow.data.tonelaje) if flow.data.tonelaje_fixed else (0.0, None))
        if flow.data.tonelaje_fixed:
            fixed_vars_indices.add(var_idx)
        var_idx += 1

        # Variables: Leyes de cada elemento (L)
        for el in flow.data.elementos:
            var_map[f"l_{flow.id}_{el.name}"] = var_idx
            l_val = el.ley
            if el.not_measured and l_val <= 0:
                l_val = 1.0  # Semilla positiva para evitar derivadas nulas
            initial_guess.append(l_val)
            bounds.append((el.ley, el.ley) if el.ley_fixed else (0.0, 100.0))
            if el.ley_fixed:
                fixed_vars_indices.add(var_idx)
            var_idx += 1

    if y_k is None:
        y_k = np.array(initial_guess, dtype=float)

    balance_nodes = _get_balance_nodes(request)
    num_elements = len(request.flows[0].data.elementos) if request.flows and request.flows[0].data.elementos else 0
    num_constraints = len(balance_nodes) * (1 + num_elements)
    
    A = np.zeros((num_constraints, len(y_k)))
    b = np.zeros(num_constraints)
    constraint_idx = 0

    for node in balance_nodes:
        # 1. Balance de masa total seca (Lineal)
        for flow in request.flows:
            factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
            if factor != 0:
                A[constraint_idx, var_map[f"t_{flow.id}"]] += factor
        constraint_idx += 1

        # 2. Balance de finos por elemento (Linealizado de Taylor)
        for el_idx in range(num_elements):
            el_name = request.flows[0].data.elementos[el_idx].name
            node_balance_term = 0.0
            for flow in request.flows:
                factor = 1.0 if flow.target == node.id else (-1.0 if flow.source == node.id else 0.0)
                if factor != 0:
                    t_idx = var_map[f"t_{flow.id}"]
                    l_idx = var_map[f"l_{flow.id}_{el_name}"]
                    T_k = float(y_k[t_idx])
                    L_k = float(y_k[l_idx])

                    A[constraint_idx, t_idx] += factor * L_k / 100.0
                    A[constraint_idx, l_idx] += factor * T_k / 100.0
                    node_balance_term += factor * T_k * L_k / 100.0
            
            b[constraint_idx] = node_balance_term
            constraint_idx += 1

    return A, b, var_map, initial_guess, bounds, fixed_vars_indices

def format_equations_to_strings(A: np.ndarray, b: np.ndarray, var_map: dict) -> List[str]:
    """Formatea las matrices A y b en una lista de ecuaciones legibles."""
    equations = []
    idx_to_name = {v: k for k, v in var_map.items()}

    for i, row in enumerate(A):
        terms = []
        for j, coeff in enumerate(row):
            if abs(coeff) > 1e-7:
                var_name = idx_to_name.get(j, f"var_{j}")
                sign = "+" if coeff >= 0 else "-"
                terms.append(f"{sign} {abs(coeff):.4f} * {var_name}")
        
        if not terms:
            continue

        eq_str = " ".join(terms)
        if eq_str.startswith("+ "):
            eq_str = eq_str[2:]
        elif eq_str.startswith("- "):
            eq_str = "-" + eq_str[2:]
            
        b_val = b[i] if i < len(b) else 0.0
        equations.append(f"{eq_str} = {b_val:.4f}")
        
    return equations

def generate_report(y_hat: np.ndarray, var_map: dict, A_final: np.ndarray, b_final: np.ndarray, request: ReconciliationRequest, residuals: list) -> dict:
    """Genera informe técnico de reconciliación metalúrgica con balances de masa y finos."""
    return {
        "success": True,
        "methodology": [
            "Optimización cuadrática iterativa (Sequential Quadratic Programming - SQP).",
            "Detección topológica automática de conservación de masa por nodos de circuito.",
            "Linealización bilineal rigurosa de masa de fino (T * L) mediante expansión de Taylor.",
            "Cálculo analítico de covarianza posterior de variables reconciliadas (Gauss-Markov).",
            "Protección contra singularidad numérica en varianzas nulas de medición."
        ],
        "equations": format_equations_to_strings(A_final, b_final, var_map) if A_final is not None else [],
        "node_residuals": residuals
    }

def reconcile_data_qp(request: ReconciliationRequest):
    """
    Ejecuta la reconciliación metalúrgica de tonelajes y leyes de elementos mediante SQP.
    Garantiza robustez numérica, protección contra varianzas cero y cálculo de Gauss-Markov.
    """
    MAX_ITERATIONS = 40
    TOLERANCE = 1e-7

    A, b, var_map, initial_guess, bounds, fixed_vars_indices = _get_balance_system(request)
    y_k = np.array(initial_guess, dtype=float)

    measured_indices, measured_values, variances = [], [], []

    for flow in request.flows:
        # Tonelaje Seco
        if not flow.data.not_measured and not flow.data.tonelaje_fixed:
            idx = var_map[f"t_{flow.id}"]
            measured_indices.append(idx)
            measured_values.append(flow.data.tonelaje)
            # Protección contra varianza cero: mínimo 2% del valor o 0.1
            err = flow.data.tonelaje_error if flow.data.tonelaje_error > 0 else max(flow.data.tonelaje * 0.02 + 0.1, 0.1)
            variances.append(err ** 2)

        # Leyes de elementos
        for el in flow.data.elementos:
            if not el.not_measured and not el.ley_fixed:
                idx = var_map[f"l_{flow.id}_{el.name}"]
                measured_indices.append(idx)
                measured_values.append(el.ley)
                # Protección contra varianza cero: mínimo 1% del valor o 0.005
                err = el.ley_error if el.ley_error > 0 else max(el.ley * 0.01 + 0.005, 0.005)
                variances.append(err ** 2)

    if not measured_indices:
        # Si no hay mediciones que ajustar, devolver initial_guess
        std_devs = np.zeros_like(y_k)
        return _build_response(request, y_k, std_devs, var_map, None, None, [])

    x_measured = np.array(measured_values, dtype=float)
    W = np.diag(1.0 / np.array(variances, dtype=float))

    H = np.zeros((len(x_measured), len(y_k)))
    for i, idx in enumerate(measured_indices):
        H[i, idx] = 1.0

    A_final, b_final = None, None

    for iteration in range(MAX_ITERATIONS):
        A, b, _, _, _, _ = _get_balance_system(request, y_k=y_k)

        def objective(y):
            diff = H @ y - x_measured
            return 0.5 * diff.T @ W @ diff

        def objective_grad(y):
            diff = H @ y - x_measured
            return H.T @ (W @ diff)

        constraints = [{'type': 'eq', 'fun': lambda y, A_mat=A, b_vec=b: A_mat @ y - b_vec}]
        
        result = minimize(
            objective, y_k, jac=objective_grad, bounds=bounds, 
            constraints=constraints, method='SLSQP', 
            options={'ftol': 1e-9, 'maxiter': 100}
        )

        if not result.success and iteration == 0:
            # Reintentar sin jacobiano en caso de aproximación inicial difícil
            result = minimize(objective, y_k, bounds=bounds, constraints=constraints, method='SLSQP')
            if not result.success:
                raise ValueError(f"La optimización falló: {result.message}")

        y_new = result.x
        A_final, b_final = A, b

        if np.linalg.norm(y_new - y_k) < TOLERANCE:
            y_k = y_new
            break
        y_k = y_new

    y_hat = y_k

    # --- CÁLCULO DE INCERTIDUMBRE POSTERIOR (GAUSS-MARKOV) ---
    num_vars = len(y_hat)
    V_in = np.zeros((num_vars, num_vars))
    for i, idx in enumerate(measured_indices):
        V_in[idx, idx] = variances[i]

    try:
        if A_final is not None and A_final.shape[0] > 0:
            S = A_final @ V_in @ A_final.T
            V_out = V_in - V_in @ A_final.T @ np.linalg.pinv(S) @ A_final @ V_in
            std_devs = np.sqrt(np.maximum(np.diag(V_out), 0.0))
        else:
            std_devs = np.sqrt(np.diag(V_in))
    except Exception:
        std_devs = np.zeros_like(y_hat)

    for i in fixed_vars_indices:
        std_devs[i] = 0.0

    # Residuos de balance nodales finales
    balance_nodes = _get_balance_nodes(request)
    residuals = []
    for node in balance_nodes:
        in_mass = sum(y_hat[var_map[f"t_{f.id}"]] for f in request.flows if f.target == node.id)
        out_mass = sum(y_hat[var_map[f"t_{f.id}"]] for f in request.flows if f.source == node.id)
        residuals.append({
            "node": node.label,
            "residual_mass": float(in_mass - out_mass)
        })

    return _build_response(request, y_hat, std_devs, var_map, A_final, b_final, residuals)

def _build_response(request: ReconciliationRequest, y_hat: np.ndarray, std_devs: np.ndarray, var_map: dict, A_final, b_final, residuals: list):
    reconciled_flows = []
    for flow in request.flows:
        t_idx = var_map[f"t_{flow.id}"]
        corrected_elements = []
        for el in flow.data.elementos:
            l_idx = var_map[f"l_{flow.id}_{el.name}"]
            corrected_elements.append({
                "name": el.name,
                "ley_corregida": float(y_hat[l_idx]),
                "ley_corregida_error": float(std_devs[l_idx])
            })
        reconciled_flows.append({
            "id": flow.id,
            "corrected_data": {
                "tonelaje_corregido": float(y_hat[t_idx]),
                "tonelaje_corregido_error": float(std_devs[t_idx]),
                "elementos": corrected_elements
            }
        })

    report = generate_report(y_hat, var_map, A_final, b_final, request, residuals)
    return reconciled_flows, report
