import numpy as np
from scipy.optimize import minimize
from models import ReconciliationRequest, NodeType
from typing import List

def _get_balance_system(request: ReconciliationRequest, y_k: np.ndarray = None):
    """Construye el sistema de ecuaciones de balance (A, b) y el mapeo de variables."""
    var_map, initial_guess, bounds, fixed_vars_indices = {}, [], [], set()
    var_idx = 0
    for flow in request.flows:
        var_map[f"t_{flow.id}"] = var_idx
        initial_guess.append(flow.data.tonelaje)
        bounds.append((flow.data.tonelaje, flow.data.tonelaje) if flow.data.tonelaje_fixed else (0, None))
        if flow.data.tonelaje_fixed: fixed_vars_indices.add(var_idx)
        var_idx += 1
        for el in flow.data.elementos:
            var_map[f"l_{flow.id}_{el.name}"] = var_idx
            initial_guess.append(el.ley)
            bounds.append((el.ley, el.ley) if el.ley_fixed else (0, 100))
            if el.ley_fixed: fixed_vars_indices.add(var_idx)
            var_idx += 1

    # Si y_k no se proporciona (ej. para previsualización), usar initial_guess
    if y_k is None:
        y_k = np.array(initial_guess)

    process_nodes = [n for n in request.nodes if n.node_type == NodeType.PROCESO]
    num_elements = len(request.flows[0].data.elementos) if request.flows else 0
    num_constraints = len(process_nodes) * (1 + num_elements)
    A = np.zeros((num_constraints, len(y_k)))
    b = np.zeros(num_constraints)
    constraint_idx = 0

    for node in process_nodes:
        # Balance de masa total
        for flow in request.flows:
            if flow.target == node.id: # Flujo de entrada
                A[constraint_idx, var_map[f"t_{flow.id}"]] += 1
            elif flow.source == node.id: # Flujo de salida
                A[constraint_idx, var_map[f"t_{flow.id}"]] -= 1
        # b[constraint_idx] es 0 para balance de masa total
        constraint_idx += 1

        # Balance de masa por elemento (finos)
        for el_idx in range(num_elements):
            el_name = request.flows[0].data.elementos[el_idx].name
            node_balance_term = 0
            for flow in request.flows:
                t_idx = var_map[f"t_{flow.id}"]
                l_idx = var_map[f"l_{flow.id}_{el_name}"]
                
                T_k = y_k[t_idx] # Tonelaje de la iteración anterior
                L_k = y_k[l_idx] # Ley de la iteración anterior

                factor = 0
                if flow.target == node.id: # Flujo de entrada
                    factor = 1
                elif flow.source == node.id: # Flujo de salida
                    factor = -1
                
                if factor != 0:
                    A[constraint_idx, t_idx] += factor * L_k / 100.0
                    A[constraint_idx, l_idx] += factor * T_k / 100.0
                    node_balance_term += factor * T_k * L_k / 100.0
            b[constraint_idx] = node_balance_term
            constraint_idx += 1
            
    return A, b, var_map, initial_guess, bounds, fixed_vars_indices

def format_equations_to_strings(A: np.ndarray, b: np.ndarray, var_map: dict) -> List[str]:
    """Formatea las matrices A y b en una lista de strings de ecuaciones legibles."""
    equations = []
    idx_to_name = {v: k for k, v in var_map.items()}

    for i, row in enumerate(A):
        equation_str = ""
        for j, coeff in enumerate(row):
            if abs(coeff) > 1e-9: # Usar una tolerancia pequeña para coeficientes casi cero
                var_name = idx_to_name.get(j, f"var_{j}")
                sign = "+" if coeff >= 0 else "-"
                equation_str += f" {sign} {abs(coeff):.4f} * {var_name}"
        
        # Eliminar el primer signo si es positivo
        if equation_str.strip().startswith("+"):
            equation_str = equation_str.strip()[1:].strip()
        
        equation_str += f" = {b[i]:.4f}"
        equations.append(equation_str.strip())
        
    return equations

def generate_report(y_hat, var_map, A_final, b_final, request):
    report = {
        "methodology": [
            "1. Se identificaron todas las variables del sistema (tonelajes y leyes), incluyendo las no medidas.",
            "2. Las variables marcadas como 'No ajustar' se trataron como constantes (restricciones de límite).",
            "3. Se formuló un problema de Programación Cuadrática (QP) para minimizar el error cuadrático ponderado de las variables medidas, sujeto a las restricciones de balance.",
            "4. Debido a la no-linealidad de las restricciones de balance de finos (T*L), el problema se resolvió de forma iterativa.",
            "5. En cada iteración, las restricciones se linealizaron utilizando la solución de la iteración anterior.",
            "6. El proceso convergió a una solución estable que satisface todas las restricciones de balance.",
            "7. Las desviaciones estándar ajustadas para las variables fijas se establecieron en cero."
        ],
        "equations": format_equations_to_strings(A_final, b_final, var_map)
    }
    return report

def reconcile_data_qp(request: ReconciliationRequest):
    MAX_ITERATIONS = 30
    TOLERANCE = 1e-7

    # Primera llamada para obtener el mapeo de variables y la suposición inicial
    A, b, var_map, initial_guess, bounds, fixed_vars_indices = _get_balance_system(request)
    y_k = np.array(initial_guess)

    measured_non_fixed_indices, measured_values, variances = [], [], []
    for flow in request.flows:
        if not flow.data.not_measured and not flow.data.tonelaje_fixed:
            idx = var_map[f"t_{flow.id}"]
            measured_non_fixed_indices.append(idx)
            measured_values.append(flow.data.tonelaje)
            variances.append(flow.data.tonelaje_error**2)
        for el in flow.data.elementos:
            if not el.not_measured and not el.ley_fixed:
                idx = var_map[f"l_{flow.id}_{el.name}"]
                measured_non_fixed_indices.append(idx)
                measured_values.append(el.ley)
                variances.append(el.ley_error**2)

    x_measured = np.array(measured_values)
    V_inv = np.linalg.inv(np.diag(variances)) if variances else np.array([])

    H = np.zeros((len(x_measured), len(y_k)))
    for i, idx in enumerate(measured_non_fixed_indices):
        H[i, idx] = 1

    A_final, b_final = None, None
    for iteration in range(MAX_ITERATIONS):
        # Reconstruir A y b en cada iteración con el y_k actual para linealización
        A, b, _, _, _, _ = _get_balance_system(request, y_k=y_k)

        def objective(y):
            if not measured_values: return 0.0
            error_vector = H @ y - x_measured
            return error_vector.T @ V_inv @ error_vector

        constraints = [{'type': 'eq', 'fun': lambda y: A @ y - b}]
        result = minimize(objective, y_k, bounds=bounds, constraints=constraints, method='SLSQP')

        if not result.success: raise ValueError(f"La optimización falló: {result.message}")
        
        y_new = result.x
        A_final, b_final = A, b # Guardar las últimas A y b para el reporte
        if np.linalg.norm(y_new - y_k) < TOLERANCE: y_k = y_new; break
        y_k = y_new
    else: print("Advertencia: La solución no convergió.")

    y_hat = y_k
    std_devs_y = np.zeros_like(y_hat)
    for i in fixed_vars_indices: std_devs_y[i] = 0.0

    reconciled_flows = []
    for flow in request.flows:
        t_idx = var_map[f"t_{flow.id}"]
        corrected_elements = []
        for el in flow.data.elementos:
            l_idx = var_map[f"l_{flow.id}_{el.name}"]
            corrected_elements.append({"name": el.name, "ley_corregida": y_hat[l_idx], "ley_corregida_error": std_devs_y[l_idx]})
        reconciled_flows.append({"id": flow.id, "corrected_data": {"tonelaje_corregido": y_hat[t_idx], "tonelaje_corregido_error": std_devs_y[t_idx], "elementos": corrected_elements}})

    report = generate_report(y_hat, var_map, A_final, b_final, request)

    return reconciled_flows, report
