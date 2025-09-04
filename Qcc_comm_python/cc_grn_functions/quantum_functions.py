"""
Quantum GRN Functions Module

This module provides functions for constructing parameterized quantum circuits
representing gene regulatory networks (GRNs), combining circuits, creating observables
from measurement histograms, and running variational quantum eigensolver (VQE) optimizations
using Qiskit. It is designed for quantum machine learning and quantum biology applications,
especially for modeling cell-type-specific gene interactions.

Main functionalities:
- GRN ansatz circuit construction
- Circuit combination and inter-cell-type interaction modeling
- Observable (Hamiltonian) construction from measurement data
- Parameter management for quantum circuits
- VQE optimization workflow

Dependencies: Qiskit, NumPy, Matplotlib, SciPy
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from qiskit.circuit import QuantumCircuit, Parameter
from qiskit.circuit.library import TwoLocal
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.primitives import StatevectorEstimator, StatevectorSampler
from qiskit.visualization import plot_histogram

def create_grn_ansatz(ng, cell_type):
    """
    Creates a parameterized quantum circuit representing a Gene Regulatory Network (GRN) ansatz.

    Args:
        ng: The number of genes in the network.
        cell_type: A string identifier for the cell type (e.g., "CT1", "CT2").

    Returns:
        A QuantumCircuit object representing the GRN ansatz.
    """
    cell_type = cell_type.lower()
    ansatz_grn = QuantumCircuit(ng, name=f"{cell_type}_GRN_Ansatz")

    # Gene activation probabilities (RY rotations after Hadamard)
    params_act = [Parameter(f'{cell_type}_act_{i}') for i in range(ng)]
    params_post_act = [Parameter(f'{cell_type}_post_acti_{i}') for i in range(ng)]
    params_post_act2 = [Parameter(f'{cell_type}_post_acti2_{i}') for i in range(ng)]
    for i in range(ng):
        ansatz_grn.h(i)
        #ansatz_grn.u(params_act[i], params_post_act[i], params_post_act2[i], i)  
        ansatz_grn.u(params_post_act[i], params_act[i], params_post_act2[i], i)  
        #ansatz_grn.u(params_post_act[i], params_post_act2[i], params_act[i], i) # Not working
        ansatz_grn.x(i)
        #ansatz_grn.u(params_act[i], params_post_act[i],  params_post_act2[i], i)  
        ansatz_grn.u(params_post_act[i], params_act[i], params_post_act2[i], i)  
        #ansatz_grn.u(params_post_act[i], params_post_act2[i], params_act[i], i) # Not working

    # Gene interaction CRX gates
    for i in range(ng):
        for j in range(ng):
            if i != j:
                param_name = f'{cell_type}_grn_{i}_{j}'
                param = Parameter(param_name)
                ansatz_grn.cry(param, i, j)
                #ansatz_grn.crx(param, i, j)

    return ansatz_grn


def create_joint_circuit_from_ansatzes(ansatz_grn_ct1, ansatz_grn_ct2, interactions=None):
    """
    Concatenates two GRN ansatz circuits, preserving their internal structure and parameter names,
    and includes additional inter-cell-type interactions using Controlled-RX rotations with parameterized angles.

    This function directly appends the instructions from the input ansatz circuits,
    avoiding the need to redefine gates and parameters.

    Args:
        ansatz_grn_ct1 (QuantumCircuit): The GRN ansatz circuit for cell type 1.
        ansatz_grn_ct2 (QuantumCircuit): The GRN ansatz circuit for cell type 2.
        interactions (dict, optional): A dictionary where keys are tuples (q1, q2) representing
                                      qubit pairs for inter-cell-type interactions. These interactions
                                      are typically between a qubit from ansatz_grn_ct1 and a qubit
                                      from ansatz_grn_ct2. The values can be arbitrary (e.g., placeholder).
                                      Defaults to None.

    Returns:
        QuantumCircuit: A new QuantumCircuit object representing the combined GRN ansatz
                        with inter-cell-type interactions.
    """
    ng_ct1 = ansatz_grn_ct1.num_qubits
    ng_ct2 = ansatz_grn_ct2.num_qubits
    num_total_qubits = ng_ct1 + ng_ct2
    
    # Create a new quantum circuit for the combined ansatz with the total number of qubits
    ccgrn_circuit = QuantumCircuit(num_total_qubits, name="Combined_GRN_Ansatz")

    # Append the first cell type's ansatz directly.
    # The 'compose' method appends the instructions of ansatz_grn_ct1 to ccgrn_circuit.
    # The 'qubits=range(ng_ct1)' argument ensures these instructions operate on the
    # first 'ng_ct1' qubits (0 to ng_ct1-1) of the combined circuit.
    # 'inplace=True' modifies ccgrn_circuit directly.
    ccgrn_circuit.compose(ansatz_grn_ct1, qubits=range(ng_ct1), inplace=True)

    # Append the second cell type's ansatz, shifting its qubits.
    # The qubits of ansatz_grn_ct2 (which are internally 0 to ng_ct2-1) will now map to
    # qubits from 'ng_ct1' to '(ng_ct1 + ng_ct2 - 1)' in the combined circuit.
    ccgrn_circuit.compose(ansatz_grn_ct2, qubits=range(ng_ct1, num_total_qubits), inplace=True)

    # Add inter-cell-type interactions if provided
    if interactions:
        for (q1, q2) in interactions.keys():
            # Validate qubit indices to ensure they are within the bounds of the combined circuit
            if not (0 <= q1 < num_total_qubits and 0 <= q2 < num_total_qubits):
                raise ValueError(f"Qubit indices ({q1}, {q2}) in interactions are out of range for {num_total_qubits} total qubits.")
            
            # Check if the interaction is indeed between qubits from different cell types.
            # This is a common use case for the 'interactions' parameter.
            is_q1_from_ct1 = q1 < ng_ct1
            is_q2_from_ct1 = q2 < ng_ct1
            
            if is_q1_from_ct1 == is_q2_from_ct1: # If both qubits are from the same cell type
                print(f"Warning: Interaction between qubits ({q1}, {q2}) appears to be within the same cell type. "
                      "Typically, 'interactions' are for inter-cell-type connections. Skipping this interaction.")
                continue # Skip this interaction as it's not an inter-cell interaction

            # Create a new parameter for the ligand-receptor (inter-cell) interaction angle
            param_name = f"lr_{q1}_{q2}" # Parameter name for ligand-receptor interaction
            angle_param = Parameter(param_name)
            
            # Apply a Controlled-RX gate for inter-cell interaction.
            # The control qubit is q1, and the target qubit is q2.
            ccgrn_circuit.crx(angle_param, q1, q2)

    return ccgrn_circuit


from qiskit.quantum_info import SparsePauliOp
from collections import Counter
import itertools # Used for generating all Pauli string combinations
import numpy as np # For numerical precision and calculations

# --- NEW HELPER FUNCTION (placed inside for self-containment, or can be external) ---
def _calculate_pauli_z_eigenvalue_for_basis_state(pauli_string_z_only: str, bit_string: str) -> int:
    """
    Calculates the eigenvalue of a Pauli string (composed ONLY of Z and I)
    for a given computational basis bit string.
    Assumes pauli_string_z_only and bit_string are both ordered from MSB to LSB
    (e.g., "ZIZ" where Z is for q2, I for q1, Z for q0, and "011" where 0 is for q2, 1 for q1, 1 for q0).
    """
    if len(pauli_string_z_only) != len(bit_string):
        raise ValueError("Pauli string length must match bit string length.")

    eigenvalue = 1
    # Iterate from MSB (index 0) to LSB (index num_features-1)
    for i in range(len(pauli_string_z_only)):
        pauli_op = pauli_string_z_only[i].upper() # Pauli op at this qubit position
        bit_val = int(bit_string[i]) # Bit value at this qubit position

        if pauli_op == 'Z':
            eigenvalue *= (1 - 2 * bit_val) # (+1 if bit is 0, -1 if bit is 1)
        # For 'I', eigenvalue *= 1, so no change needed
            
    return eigenvalue


# def create_interaction_observable_from_histogram(
#     joint_counts: Counter,
#     num_features: int,
#     min_ones: int = 0, # Kept as per your original signature
#     unobserved_punishment: float = 1.0,
#     normalization_offset: float = 0.0
# ) -> SparsePauliOp:
#     """Creates a SparsePauliOp representing a diagonal Hamiltonian in the computational basis.
#        The energy of each basis state is derived from its count in the histogram,
#        favoring states that meet `min_ones` and are observed, and punishing others.

#     Args:
#         joint_counts: A Counter object of observed bit string counts.
#         num_features: The total number of qubits.
#         min_ones: The minimum number of '1's required in a bit string for it to be
#                   considered 'favorable' if observed. States not meeting this or unobserved
#                   will be assigned `unobserved_punishment` energy.
#         unobserved_punishment: The positive energy value assigned to unobserved bit strings,
#                                or observed bit strings that don't meet `min_ones` criteria.
#                                These states will be 'punished' (avoided by optimizer).
#         normalization_offset: An optional offset to subtract from counts before calculating energy.
#                               Useful for centering counts, e.g., for 50/50 cases.

#     Returns:
#         A SparsePauliOp observable.
#     """
    
#     # 1. Determine the desired energy (E_b) for *each* computational basis state |b>
#     # This replaces the old `strength` determination logic.
#     state_energies = {} # Maps bit_string -> desired_energy (e.g., "011" -> -500)
    
#     # Iterate through all possible bit strings (2^num_features)
#     for i in range(2**num_features):
#         bit_string = format(i, '0' + str(num_features) + 'b') # Generates MSB...LSB (q_N-1 ... q_0)
#         num_ones = bit_string.count('1') # Count '1's for `min_ones` criteria

#         if num_ones < min_ones:
#             # If it doesn't meet min_ones, treat as undesirable, assign punishment
#             energy_b = unobserved_punishment
#         elif bit_string in joint_counts:
#             # If observed AND meets min_ones: energy proportional to -count
#             energy_b = -float(joint_counts[bit_string] - normalization_offset)
#             #energy_b = -np.log(float(joint_counts[bit_string]) + 1.0)
#         else:
#             # If unobserved (but meets min_ones, or min_ones=0): assign punishment
#             energy_b = unobserved_punishment
        
#         state_energies[bit_string] = energy_b

#     # 2. Generate all possible Z-only Pauli strings (e.g., I, Z, II, IZ, ZI, ZZ for 2 qubits)
#     # in MSB to LSB order (Q_{N-1}...Q_0) for consistency with Qiskit Pauli string order.
#     all_pauli_z_strings = []
#     for pauli_tuple in itertools.product('IZ', repeat=num_features):
#         pauli_string = "".join(pauli_tuple)
#         all_pauli_z_strings.append(pauli_string)
    
#     # 3. Convert these state_energies (E_b) into coefficients (c_P) for Pauli strings
#     # This loop replaces your original `for i in range(2**num_features)` and `pauli_string` construction.
#     pauli_term_coefficients = {}

#     for pauli_str in all_pauli_z_strings: # Iterate through all canonical Pauli-Z strings
#         coeff_p = 0.0
#         # Sum over all basis states 'b' (all 2^N bit strings)
#         for bit_string, energy_b in state_energies.items():
#             # Get eigenvalue of Pauli string P for basis state |b>
#             # <--- CHANGE THIS LINE HERE --- >
#             eigenvalue_p_b = _calculate_pauli_z_eigenvalue_for_basis_state(pauli_str, bit_string) 
#             coeff_p += energy_b * eigenvalue_p_b
        
#         # Normalize by 2^N
#         coeff_p /= (2**num_features)
        
#         if abs(coeff_p) > 1e-9: # Add term only if coefficient is significant
#             pauli_term_coefficients[pauli_str] = coeff_p

#     # Convert to SparsePauliOp format
#     interaction_strength_list = [(pauli_str, coeff) for pauli_str, coeff in pauli_term_coefficients.items()]
#     interaction_observable = SparsePauliOp.from_list(interaction_strength_list)
#     return interaction_observable


# --- Modified create_interaction_observable_from_histogram function ---
def create_interaction_observable_from_histogram(
    joint_counts: Counter,
    num_features: int,
    min_ones: int = 0,
    max_ones: int = None,
    unobserved_punishment: float = 1.0,
    normalization_offset: float = 0.0,
    energy_mapping_strategy: str = 'linear_counts', # New parameter: 'linear_counts', 'log_counts', 'neg_log_prob'
    epsilon: float = 1e-9 # Small value to avoid log(0) for 'neg_log_prob'
) -> SparsePauliOp:
    """Creates a SparsePauliOp representing a diagonal Hamiltonian in the computational basis.
    The energy of each basis state is derived from its count in the histogram,
    favoring states that meet `min_ones` and `max_ones` criteria and are observed,
    and punishing others. The energy mapping strategy can be customized.

    Args:
        joint_counts: A Counter object of observed bit string counts.
        num_features: The total number of qubits.
        min_ones: The minimum number of '1's required in a bit string for it to be
                  considered 'favorable' if observed.
        max_ones: The maximum number of '1's allowed in a bit string for it to be
                  considered 'favorable' if observed. If None, there is no upper limit.
        unobserved_punishment: The positive energy value assigned to unobserved bit strings,
                               or observed bit strings that don't meet `min_ones`/`max_ones` criteria.
        normalization_offset: An optional offset to subtract from counts before calculating energy.
                               Useful for centering counts, primarily for 'linear_counts' strategy.
        energy_mapping_strategy (str): Defines how counts are mapped to energy:
                                       - 'linear_counts': energy = -(count - offset)
                                       - 'log_counts': energy = -log(count + 1.0) (smoother)
                                       - 'neg_log_prob': energy = -log(probability + epsilon) (likelihood-based)
        epsilon (float): A small value added to probabilities/counts to avoid log(0) for log-based strategies.

    Returns:
        A SparsePauliOp observable.
    """
    
    state_energies = {}
    
    # Calculate total counts if 'neg_log_prob' strategy is used
    total_counts = sum(joint_counts.values())
    if total_counts == 0 and energy_mapping_strategy == 'neg_log_prob':
        print("Warning: Total counts are zero. All states will be assigned unobserved_punishment.")

    # Iterate through all possible bit strings (2^num_features)
    for i in range(2**num_features):
        bit_string = format(i, '0' + str(num_features) + 'b')
        num_ones = bit_string.count('1')

        energy_b = unobserved_punishment # Default energy is punishment

        # Check if the number of ones is within the specified range
        is_within_range = (num_ones >= min_ones)
        if max_ones is not None:
            is_within_range = is_within_range and (num_ones <= max_ones)

        # If the state is within the desired range AND it was observed in the histogram
        if is_within_range and bit_string in joint_counts:
            if energy_mapping_strategy == 'linear_counts':
                energy_b = -float(joint_counts[bit_string] - normalization_offset)
            elif energy_mapping_strategy == 'log_counts':
                # Add 1.0 to count to avoid log(0) if count is 0, and to make it always positive
                energy_b = -np.log(float(joint_counts[bit_string]) + 1.0)
            elif energy_mapping_strategy == 'neg_log_prob':
                # Calculate probability, add epsilon to avoid log(0)
                prob_b = joint_counts[bit_string] / total_counts if total_counts > 0 else 0.0
                energy_b = -np.log(prob_b + epsilon)
            else:
                raise ValueError(f"Unknown energy_mapping_strategy: {energy_mapping_strategy}")
        # Else, energy_b remains unobserved_punishment (its default value)
            
        state_energies[bit_string] = energy_b

    # Generate all possible Z-only Pauli strings
    all_pauli_z_strings = []
    for pauli_tuple in itertools.product('IZ', repeat=num_features):
        pauli_string = "".join(pauli_tuple)
        all_pauli_z_strings.append(pauli_string)
    
    # Convert these state_energies (E_b) into coefficients (c_P) for Pauli strings
    pauli_term_coefficients = {}

    for pauli_str in all_pauli_z_strings:
        coeff_p = 0.0
        for bit_string, energy_b in state_energies.items():
            eigenvalue_p_b = _calculate_pauli_z_eigenvalue_for_basis_state(pauli_str, bit_string) 
            coeff_p += energy_b * eigenvalue_p_b
        
        coeff_p /= (2**num_features)
        
        if abs(coeff_p) > 1e-9:
            pauli_term_coefficients[pauli_str] = coeff_p

    # Convert to SparsePauliOp format
    interaction_strength_list = [(pauli_str, coeff) for pauli_str, coeff in pauli_term_coefficients.items()]
    interaction_observable = SparsePauliOp.from_list(interaction_strength_list)
    return interaction_observable


import collections
# The corrected version from previous discussions
def create_interaction_observable_general(interactions, num_features):
    """Creates a SparsePauliOp observable for generalized interactions.
    Args:
        interactions: A dictionary where keys are tuples of node indices 
                      (e.g., (0, 1), (0, 0, 2), (0, 1, 2, 3)) and 
                      values are the corresponding interaction strengths.
        num_features: The total number of qubits.

    Returns:
        A SparsePauliOp observable.
    """
    interaction_strength_list = []
    for nodes, strength in interactions.items():
        strength = -strength # Assuming you want to minimize this energy
        
        pauli_chars = ['I'] * num_features
        
        # This loop correctly maps Qiskit's qubit index (0 for LSB, N-1 for MSB)
        # to the position in an MSB-first string (0 for MSB, N-1 for LSB).
        for node_idx in nodes:
            if not (0 <= node_idx < num_features):
                raise ValueError(
                    f"Node index {node_idx} is out of bounds for {num_features} features."
                )
            pauli_chars[num_features - 1 - node_idx] = "Z" # This is the crucial line
            
        pauli_string = "".join(pauli_chars)
        
        interaction_strength_list.append((pauli_string, strength))

    interaction_observable = SparsePauliOp.from_list(interaction_strength_list)
    return interaction_observable


def evaluate_and_plot_ansatz(ansatz, params, shots=1024, title="VQE Quantum Sampler Results", figsize=(5, 4), filename=None):
    """Evaluates a quantum ansatz, plots results and circuit, and prints counts."""
    try:
        sampler = StatevectorSampler()
        bound_circuit = ansatz.copy()
        bound_circuit.assign_parameters(params, inplace=True)
        bound_circuit.measure_all()

        job = sampler.run([bound_circuit], shots=shots)
        pub_result = job.result()[0]
        data_pub = pub_result.data
        counts = data_pub.meas.get_counts()

        print(f"The counts are: {counts}")

        # Matplotlib customization:
        sorted_counts = dict(sorted(counts.items()))
        x_labels = list(sorted_counts.keys())
        y_values = list(sorted_counts.values())

        plt.figure(figsize=figsize)
        plt.bar(x_labels, y_values)
        plt.xlabel("Measurement Outcomes")
        plt.ylabel("Counts")
        plt.title(title)
        plt.xticks(rotation=45, ha='right')
        #plt.xlabel("Measurement Outcomes", fontsize=16)
        #plt.ylabel("Counts", fontsize=16)
        #plt.title(title, fontsize=18)
        #plt.xticks(rotation=45, ha='right', fontsize=14)
        plt.tight_layout()
        if filename:
            plt.savefig(filename, bbox_inches='tight')
        plt.show()

        return counts, bound_circuit

    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def create_parameter_dictionaries(combined_qc, act_percentages):
    """Creates static and variable parameter dictionaries."""
    # Get Hadamard parameters
    # params_ct = [param for param in combined_qc.parameters if 'ct1_' in param.name or 'ct2_' in param.name
    #                                                         and 'grn' not in param.name 
    #                                                         and 'lr' not in param.name 
    #                                                         and 'act2' not in param.name]
    params_ct = [param for param in combined_qc.parameters if '_act_' in param.name]

    static_params = {}
    for i, val in enumerate(act_percentages):
        static_params[params_ct[i]] = val

    variable_params = [param for param in combined_qc.parameters if param not in static_params]
    x0_interaction = np.ones(len(variable_params))*np.pi*0 # All zeros | Or pi terms
    variable_params = dict(zip(variable_params, x0_interaction))

    return static_params, variable_params


def create_parameter_dictionaries_cust(combined_qc, act_percentages):
    """Creates static and variable parameter dictionaries."""
    params_ct = [param for param in combined_qc.parameters if '_act_' in param.name]

    static_params = {}
    variable_params = [param for param in combined_qc.parameters if param not in static_params]

    # Initialize variable parameters
    #x0_interaction = np.zeros(len(variable_params))  # All zeros
    x0_interaction = np.ones(len(variable_params))*np.pi*0  # All zeros
    variable_params = dict(zip(variable_params, x0_interaction))

    # Now, iterate through the identified 'act' parameters and assign their
    for i, param in enumerate(params_ct):
        variable_params[param] = act_percentages[i]

    return static_params, variable_params

# Create the static and variable parameter dictionaries directly from the circuit.
def create_parameter_dictionaries_from_circuit(circuit):
    """Creates static and variable parameter dictionaries directly from the circuit."""
    static_params = {param: None for param in circuit.parameters if 'lr_' not in param.name}
    variable_params = {param: 0.0 for param in circuit.parameters if 'lr_' in param.name}
    return static_params, variable_params

def cost_func_vqe(params, combined_qc, hamiltonian, estimator):  # combined_qc here
    """Cost function for VQE"""
    bound_qc = combined_qc.assign_parameters(params)  # Assign parameters INSIDE cost_func_vqe
    statevector = Statevector(bound_qc)  # Use bound_qc
    statevector_array = statevector.data
    hamiltonian_matrix = hamiltonian.to_matrix()
    energy = np.real(statevector_array.conjugate() @ hamiltonian_matrix @ statevector_array)
    return energy

def cost_func_wrapper(variable_values, all_params, combined_qc, interaction_observable, estimator, variable_params):
    for i, param in enumerate(variable_params):
        all_params[param] = variable_values[i]
    return cost_func_vqe(all_params, combined_qc, interaction_observable, estimator) # Pass combined_qc

import numpy as np
import matplotlib.pyplot as plt
from qiskit.primitives import StatevectorEstimator
from scipy.optimize import minimize
def vqe_solver(
    histogram_data,
    circuit, # Renamed from 'cirquit' for common convention
    act_percentages, # Renamed from 'act_percentages' for consistency with create_parameter_dictionaries_cust
    cost_func_wrapper, # This function needs to be defined to accept the correct arguments (see comments below)
    min_ones_obs = 0, # Added as an explicit argument for flexibility
    max_ones_obs = None,
    energy_map = 'log_counts', # New parameter: 'linear_counts', 'log_counts', 'neg_log_prob'
    optimizer_method="COBYLA" # "L-BFGS-B | COBYLA
):
    """
    Performs a Variational Quantum Eigensolver (VQE) optimization.

    This function encapsulates the VQE workflow, including:
    1. Creating an interaction observable from histogram data.
    2. Setting up static and variable parameters for the quantum circuit.
    3. Initializing the Qiskit StatevectorEstimator.
    4. Running the optimization using scipy.optimize.minimize.
    5. Collecting cost function values during optimization.
    6. Updating parameters with optimized values.
    7. (Optional) Plotting the energy minimization curve.

    Args:
        histogram_data (dict): Data used to create the interaction observable.
        circuit (QuantumCircuit): The parameterized quantum circuit (ansatz).
        act_percentages (list): Initial percentage values for parameters
                                containing '_act_' in the circuit.
        cost_func_wrapper (callable): The cost function to be minimized.
                                      It MUST accept arguments in the following order:
                                      cost_func_wrapper(current_variable_params_array,
                                                        static_params_dict,
                                                        circuit,
                                                        observable,
                                                        estimator,
                                                        variable_param_objects_list)
                                      Inside this function, you should combine
                                      static_params_dict and map current_variable_params_array
                                      to variable_param_objects_list to form the
                                      full parameter dictionary for circuit binding.
        min_ones_for_observable (int): Minimum number of '1's for filtering
                                       histogram data when creating the observable.
        optimizer_method (str): The optimization method to use for scipy.minimize
                                (e.g., "L-BFGS-B", "COBYLA", "SLSQP").

    Returns:
        tuple: A tuple containing:
            - result_object (OptimizeResult): The result object from scipy.optimize.minimize.
            - optimized_full_params (dict): The dictionary of all parameters
                                            with their final optimized values.
            - cost_history (list): A list of cost function values recorded
                                   during each optimization iteration.
    """

    num_qubits = circuit.num_qubits # Get number of qubits from the circuit

    # 1. Create interaction observable
    #mean_offset = np.mean(np.array(list(histogram_data.values())))*0.75*0
    mean_offset = 0
    interaction_observable = create_interaction_observable_from_histogram(
        joint_counts = histogram_data, num_features = num_qubits, 
        min_ones = min_ones_obs, 
        max_ones = max_ones_obs,
        energy_mapping_strategy = energy_map,
        unobserved_punishment = 1, 
        normalization_offset = mean_offset)

    print("Interaction observable CT from histogram:", interaction_observable)

    # 2. Create static and variable parameter dictionaries
    # This uses the create_parameter_dictionaries_cust function from the immersive
    # which sets static_params to empty and variable_params to all circuit parameters
    # initialized to 0.0.
    # Ensure you are using create_parameter_dictionaries_cust here, not create_parameter_dictionaries
    #static_params, variable_params_dict = create_parameter_dictionaries_cust(circuit, act_percentages)
    static_params, variable_params_dict = create_parameter_dictionaries(circuit, act_percentages) 


    print("Static Parameters:", static_params)
    print("Variable Parameters:", variable_params_dict)

    # 3. Initialize Qiskit Estimator
    estimator = StatevectorEstimator()

    # 4. Prepare initial guess for optimization
    # The 'variable_params_dict' contains parameter objects as keys and 0.0 as values.
    # For scipy.minimize, we need an array of initial values (x0).
    # We extract the initial values from variable_params_dict to form x0.
    x0_interaction = np.array(list(variable_params_dict.values()))

    # To correctly map back optimized values, we also need a list of the parameter objects
    # that correspond to the order of values in x0_interaction.
    variable_param_objects = list(variable_params_dict.keys())


    # 5. Create initial full parameter dictionary (this is mainly for initial print and structure)
    # The actual 'all_params' for binding will be constructed inside cost_func_wrapper.
    initial_all_params_for_display = static_params.copy()
    initial_all_params_for_display.update(dict(zip(variable_param_objects, x0_interaction)))


    # in the context of the vqe_solver function scope.
    iteration_data = {'counter': 0} 
    cost_values = [] # List to store cost values at each iteration
    # Define the callback function for minimize
    def callback_func(xk):
        # Access the counter from the enclosing scope's dictionary
        current_counter = iteration_data['counter']
        if current_counter > 100:
            print_criteria = 100
        else:
            print_criteria = 20

        current_cost = cost_func_wrapper(xk, static_params, circuit, interaction_observable, estimator, variable_param_objects)
        cost_values.append(current_cost)

        # Print the current cost only every 20 iterations
        if current_counter % print_criteria == 0:
            print(f"Iteration {current_counter}: Current cost: {current_cost}")
            
        iteration_data['counter'] += 1 # Increment the counter

    # 6. Call minimize with args
    print(f"Starting optimization with method: {optimizer_method}")
    result_interaction = minimize(
        cost_func_wrapper,
        x0_interaction,
        # IMPORTANT: Pass static_params and variable_param_objects as fixed arguments
        # The cost_func_wrapper will use xk (the first argument) with these to bind parameters.
        args=(static_params, circuit, interaction_observable, estimator, variable_param_objects),
        method=optimizer_method, # Use the passed optimizer_method
        callback=callback_func # Use the defined callback function
        #tol= 1e-5,
    )

    print("\nOptimization Result:")
    print(result_interaction)
    print(f"\nFinal Energy: {result_interaction.fun}")

    # 7. Update the full parameter dictionary with optimized variable parameters
    optimized_variable_parameters = result_interaction.x

    # Construct the final optimized full parameter dictionary
    optimized_full_params = static_params.copy()
    for param_obj, value in zip(variable_param_objects, optimized_variable_parameters):
        optimized_full_params[param_obj] = value

    print("\nOptimized Full Parameters:")
    # Print optimized parameters by name for readability
    for param, value in optimized_full_params.items():
        print(f"  {param.name}: {value}")


    # Return the results
    return result_interaction, optimized_full_params, cost_values



from qiskit.primitives import StatevectorEstimator
from scipy.optimize import minimize
import numpy as np
# from qiskit.circuit import Parameter # Assuming Parameter is imported if needed in other parts
# from qiskit.quantum_info import Statevector # Assuming Statevector is imported if needed in other parts
# from quantum_functions import create_interaction_observable_general # Assuming this is available

def vqe_lr_solver(
    cc_grn_circuit_co,
    optimized_full_params_ct1_co,
    optimized_full_params_ct2_co,
    interactions_lr,
    ng_ct1,
    ng_ct2,
    cost_func_wrapper,
    optimizer_method="COBYLA" # "L-BFGS-B | COBYLA
):
    """
    Performs a VQE-like optimization for LR (Long Range) interactions using the L-BFGS-B method.

    Args:
        cc_grn_circuit_co (QuantumCircuit): The quantum circuit to be optimized, representing the
                                            ansatz or system under study.
        optimized_full_params_ct1_co (dict): A dictionary of pre-optimized parameters for
                                            cell type 1, which will be incorporated as static
                                            parameters in this optimization.
        optimized_full_params_ct2_co (dict): A dictionary of pre-optimized parameters for
                                            cell type 2, also incorporated as static parameters.
        interactions_lr (list): A list defining the structure of the long-range interactions
                                for which the observable is created.
        ng_ct1 (int): The number of qubits associated with cell type 1.
        ng_ct2 (int): The number of qubits associated with cell type 2.
        cost_func_wrapper (callable): The function that calculates the cost (e.g., expectation value)
                                      for a given set of parameters. It should accept arguments in the
                                      order: (variable_parameters_array, full_parameter_dictionary,
                                      circuit, observable, estimator, list_of_variable_parameters).


    Returns:
        tuple: A tuple containing:
            - result_lr_bfgs (OptimizeResult): The result object returned by `scipy.optimize.minimize`,
                                              containing optimization details like the optimized parameters,
                                              function value, number of iterations, etc.
            - all_params_lr_co (dict): The final combined dictionary of all circuit parameters,
                                      including the newly optimized long-range parameters.
            - cost_values (list): A list of cost function values recorded at each iteration
                                  during the optimization process.
    """
    # Initialize the StatevectorEstimator, which will be used to compute expectation values
    estimator = StatevectorEstimator()

    # --- 1. Create Static Parameter Dictionaries ---
    # Extract static and variable parameters from the quantum circuit.
    # 'static_params_lr' will hold parameters that are fixed during this optimization,
    # 'variable_params_lr' will hold parameters that will be optimized.
    static_params_lr, variable_params_lr = create_parameter_dictionaries_from_circuit(cc_grn_circuit_co)

    # Update the 'static_params_lr' dictionary with pre-optimized values from ct1 and ct2.
    # This ensures that parameters already optimized in other stages are fixed for this LR optimization.
    for param in static_params_lr:
        if param.name in [p.name for p in optimized_full_params_ct1_co]:
            # If the parameter's name matches one from ct1 optimized parameters, update its value.
            static_params_lr[param] = optimized_full_params_ct1_co[next(p for p in optimized_full_params_ct1_co if p.name == param.name)]
        elif param.name in [p.name for p in optimized_full_params_ct2_co]:
            # If the parameter's name matches one from ct2 optimized parameters, update its value.
            static_params_lr[param] = optimized_full_params_ct2_co[next(p for p in optimized_full_params_ct2_co if p.name == param.name)]

    # Initialize the array of variable parameters with zeros. This serves as the starting point
    # for the L-BFGS-B optimization algorithm.
    x0_lr = np.zeros(len(variable_params_lr))

    # Create the combined parameter dictionary 'all_params_lr_co'.
    # Start with the static parameters and then add the variable parameters (initially set to zeros).
    # This dictionary will be passed to the cost function.
    all_params_lr_co = static_params_lr.copy()
    all_params_lr_co.update(dict(zip(variable_params_lr, x0_lr)))

    print("Initial combined parameters for LR optimization:", all_params_lr_co)

    # Create the interaction observable for the long-range interactions.
    # This observable is crucial for defining the cost function (e.g., by computing its expectation value).
    # The total number of qubits is the sum of qubits from ct1 and ct2.
    interaction_observable_lr_co = create_interaction_observable_general(interactions_lr, ng_ct1 + ng_ct2)
    print("Interaction observable for LR custom:", interaction_observable_lr_co)

    # Initialize a list to store the cost values at each iteration of the optimization.
    cost_values = []

    # Perform the optimization using scipy's minimize function with the L-BFGS-B method.
    # `cost_func_wrapper`: The function to minimize.
    # `x0_lr`: The initial guess for the variable parameters.
    # `args`: Additional arguments to pass to the `cost_func_wrapper`. These are fixed during optimization.
    # `method`: The optimization algorithm to use.
    # `callback`: A function called after each iteration, used here to record cost values.
    result_lr_bfgs = minimize(
        cost_func_wrapper,
        x0_lr,
        args=(all_params_lr_co, cc_grn_circuit_co, interaction_observable_lr_co, estimator, variable_params_lr),
        method = optimizer_method,
        callback=lambda xk: cost_values.append(
            # Recalculate the cost with the current 'xk' (variable parameters) and append it.
            cost_func_wrapper(xk, all_params_lr_co, cc_grn_circuit_co, interaction_observable_lr_co, estimator, variable_params_lr)
        )
    )

    print("Optimization result for LR VQE:", result_lr_bfgs)

    # --- 6. Results and DataFrame ---
    # Extract the optimized values for the variable parameters from the optimization result.
    optimized_lr_values = result_lr_bfgs.x

    # Update the combined parameter dictionary with the newly optimized long-range parameter values.
    all_params_lr_co.update(dict(zip(variable_params_lr, optimized_lr_values)))

    # Return the optimization result, the final combined parameters, and the list of recorded cost values.
    return result_lr_bfgs, all_params_lr_co, cost_values

import numpy as np
from qiskit.primitives import StatevectorEstimator
from scipy.optimize import minimize
# Ensure other necessary Qiskit imports like Parameter, QuantumCircuit, SparsePauliOp are available
# from qiskit.circuit import Parameter, QuantumCircuit
# from qiskit.quantum_info import SparsePauliOp
# from collections import Counter # Needed for histogram_data_lr type hint
# import itertools # Needed for create_interaction_observable_from_histogram

# Assuming create_parameter_dictionaries_from_circuit is defined as previously:
def create_parameter_dictionaries_from_circuit(circuit):
    """
    Creates static and variable parameter dictionaries directly from the circuit.
    Parameters with 'lr_' in their name are considered variable; others are static.
    Initializes variable parameters to 0.0.
    """
    static_params = {}
    variable_params = {}

    for param in circuit.parameters:
        if 'lr_' in param.name:
            # LR parameters are variable in this VQE step, initialize to 0.0
            variable_params[param] = 0.0
        else:
            # Other parameters are static; get their current value (or 0.0 if not set)
            try:
                param_value = param.value
                if param_value is None:
                    param_value = 0.0
            except AttributeError:
                param_value = 0.0
            static_params[param] = param_value
            
    return static_params, variable_params

# This is the updated vqe_lr_solver2 function with general variable names
def vqe_lr_solver2(
    cc_grn_circuit, # Renamed
    optimized_full_params_ct1, # Renamed
    optimized_full_params_ct2, # Renamed
    histogram_data_lr,
    ng_ct1,
    ng_ct2,
    cost_func_wrapper, # Assumed to be correctly defined elsewhere
    optimizer_method="COBYLA",
    min_ones_obs=0,
    unobserved_punishment_obs=1.0,
    normalization_offset_obs=0.0
):
    """
    Performs a VQE-like optimization for LR (Long Range) interactions.

    Args:
        cc_grn_circuit (QuantumCircuit): The quantum circuit to be optimized.
        optimized_full_params_ct1 (dict): Pre-optimized parameters for cell type 1.
        optimized_full_params_ct2 (dict): Pre-optimized parameters for cell type 2.
        histogram_data_lr (Counter): Histogram data to create the LR interaction observable.
        ng_ct1 (int): Number of qubits associated with cell type 1.
        ng_ct2 (int): Number of qubits associated with cell type 2.
        cost_func_wrapper (callable): The cost function to be minimized.
        optimizer_method (str): The optimization method to use (e.g., "L-BFGS-B", "COBYLA").
        min_ones_obs (int): Minimum number of '1's for filtering histogram data when creating the observable.
        unobserved_punishment_obs (float): Punishment energy for unobserved states in the observable.
        normalization_offset_obs (float): Normalization offset for histogram counts in the observable.

    Returns:
        tuple: A tuple containing:
            - result_lr_bfgs (OptimizeResult): The result object from scipy.optimize.minimize.
            - all_params_lr (dict): The final combined dictionary of all circuit parameters.
            - cost_values (list): A list of cost function values recorded.
    """
    estimator = StatevectorEstimator()

    # --- 1. Create Static and Variable Parameter Dictionaries ---
    # Use create_parameter_dictionaries_from_circuit to ensure 'lr_' parameters are variable
    static_params_lr, variable_params_lr_dict = create_parameter_dictionaries_from_circuit(cc_grn_circuit) # Renamed

    # Update the 'static_params_lr' dictionary with pre-optimized values from ct1 and ct2.
    # This ensures that parameters already optimized in other stages are fixed for this LR optimization.
    for param in static_params_lr:
        if param.name in [p.name for p in optimized_full_params_ct1]: # Renamed
            static_params_lr[param] = optimized_full_params_ct1[next(p for p in optimized_full_params_ct1 if p.name == param.name)] # Renamed
        elif param.name in [p.name for p in optimized_full_params_ct2]: # Renamed
            static_params_lr[param] = optimized_full_params_ct2[next(p for p in optimized_full_params_ct2 if p.name == param.name)] # Renamed

    # Initial guess for optimization (all 'lr_' parameters start at 0.0 as per create_parameter_dictionaries_from_circuit)
    x0_lr = np.array(list(variable_params_lr_dict.values()))
    variable_param_objects = list(variable_params_lr_dict.keys())

    # Create the combined parameter dictionary 'all_params_lr'.
    all_params_lr = static_params_lr.copy() # Renamed
    all_params_lr.update(dict(zip(variable_param_objects, x0_lr)))

    print("Initial combined parameters for LR optimization:", all_params_lr) # Renamed

    # --- 2. Create the LR Interaction Observable from Histogram Data ---
    num_total_qubits = ng_ct1 + ng_ct2
    interaction_observable_lr = create_interaction_observable_from_histogram(
        histogram_data_lr,
        num_total_qubits,
        min_ones=min_ones_obs,
        unobserved_punishment=unobserved_punishment_obs,
        normalization_offset=normalization_offset_obs
    )
    print("Interaction observable for LR from histogram:", interaction_observable_lr) # Fixed missing parenthesis

    # Initialize a list to store the cost values at each iteration
    cost_values = []
    iteration_data = {'counter': 0}

    def callback_func(xk):
        current_counter = iteration_data['counter']
        print_criteria = 100 if current_counter > 100 else 20
        # The cost_func_wrapper will correctly combine xk (variable_values) with static_params_lr
        current_cost = cost_func_wrapper(
            xk, static_params_lr, cc_grn_circuit, interaction_observable_lr, estimator, variable_param_objects # Renamed cc_grn_circuit
        )
        cost_values.append(current_cost)
        if current_counter % print_criteria == 0:
            print(f"Iteration {current_counter}: Current cost: {current_cost}")
        iteration_data['counter'] += 1

    print(f"Starting LR optimization with method: {optimizer_method}")
    result_lr_bfgs = minimize(
        cost_func_wrapper,
        x0_lr,
        args=(static_params_lr, cc_grn_circuit, interaction_observable_lr, estimator, variable_param_objects), # Renamed cc_grn_circuit
        method=optimizer_method,
        callback=callback_func
    )

    print("\nOptimization result for LR VQE:", result_lr_bfgs)
    print(f"\nFinal LR Energy: {result_lr_bfgs.fun}")

    # Update the full parameter dictionary with optimized variable parameters
    optimized_lr_values = result_lr_bfgs.x
    all_params_lr.update(dict(zip(variable_param_objects, optimized_lr_values))) # Renamed

    print("\nOptimized Full Parameters (including LR):")
    for param, value in all_params_lr.items(): # Renamed
        print(f"  {param.name}: {value}")

    return result_lr_bfgs, all_params_lr, cost_values # Renamed