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


def create_interaction_observable_from_histogram(
    joint_counts: Counter,
    num_features: int,
    min_ones: int = 0, # Kept as per your original signature
    # Added new parameters for fine-grained control over energy assignment
    unobserved_punishment: float = 10.0,
    normalization_offset: float = 0.0
) -> SparsePauliOp:
    """Creates a SparsePauliOp representing a diagonal Hamiltonian in the computational basis.
       The energy of each basis state is derived from its count in the histogram,
       favoring states that meet `min_ones` and are observed, and punishing others.

    Args:
        joint_counts: A Counter object of observed bit string counts.
        num_features: The total number of qubits.
        min_ones: The minimum number of '1's required in a bit string for it to be
                  considered 'favorable' if observed. States not meeting this or unobserved
                  will be assigned `unobserved_punishment` energy.
        unobserved_punishment: The positive energy value assigned to unobserved bit strings,
                               or observed bit strings that don't meet `min_ones` criteria.
                               These states will be 'punished' (avoided by optimizer).
        normalization_offset: An optional offset to subtract from counts before calculating energy.
                              Useful for centering counts, e.g., for 50/50 cases.

    Returns:
        A SparsePauliOp observable.
    """
    
    # 1. Determine the desired energy (E_b) for *each* computational basis state |b>
    # This replaces the old `strength` determination logic.
    state_energies = {} # Maps bit_string -> desired_energy (e.g., "011" -> -500)
    
    # Iterate through all possible bit strings (2^num_features)
    for i in range(2**num_features):
        bit_string = format(i, '0' + str(num_features) + 'b') # Generates MSB...LSB (q_N-1 ... q_0)
        num_ones = bit_string.count('1') # Count '1's for `min_ones` criteria

        if num_ones < min_ones:
            # If it doesn't meet min_ones, treat as undesirable, assign punishment
            energy_b = unobserved_punishment
        elif bit_string in joint_counts:
            # If observed AND meets min_ones: energy proportional to -count
            energy_b = -float(joint_counts[bit_string] - normalization_offset)
        else:
            # If unobserved (but meets min_ones, or min_ones=0): assign punishment
            energy_b = unobserved_punishment
        
        state_energies[bit_string] = energy_b

    # 2. Generate all possible Z-only Pauli strings (e.g., I, Z, II, IZ, ZI, ZZ for 2 qubits)
    # in MSB to LSB order (Q_{N-1}...Q_0) for consistency with Qiskit Pauli string order.
    all_pauli_z_strings = []
    for pauli_tuple in itertools.product('IZ', repeat=num_features):
        pauli_string = "".join(pauli_tuple)
        all_pauli_z_strings.append(pauli_string)
    
    # 3. Convert these state_energies (E_b) into coefficients (c_P) for Pauli strings
    # This loop replaces your original `for i in range(2**num_features)` and `pauli_string` construction.
    pauli_term_coefficients = {}

    for pauli_str in all_pauli_z_strings: # Iterate through all canonical Pauli-Z strings
        coeff_p = 0.0
        # Sum over all basis states 'b' (all 2^N bit strings)
        for bit_string, energy_b in state_energies.items():
            # Get eigenvalue of Pauli string P for basis state |b>
            # <--- CHANGE THIS LINE HERE --- >
            eigenvalue_p_b = _calculate_pauli_z_eigenvalue_for_basis_state(pauli_str, bit_string) 
            coeff_p += energy_b * eigenvalue_p_b
        
        # Normalize by 2^N
        coeff_p /= (2**num_features)
        
        if abs(coeff_p) > 1e-9: # Add term only if coefficient is significant
            pauli_term_coefficients[pauli_str] = coeff_p

    # Convert to SparsePauliOp format
    interaction_strength_list = [(pauli_str, coeff) for pauli_str, coeff in pauli_term_coefficients.items()]
    interaction_observable = SparsePauliOp.from_list(interaction_strength_list)
    return interaction_observable

import collections
# Import SparsePauliOp for returning the result in the desired format
from qiskit.quantum_info import SparsePauliOp

def create_interaction_observable_from_histogram2(joint_counts, min_ones, num_features):
    """
    Analyzes a histogram of bit strings, accumulating strengths for specific Pauli patterns.
    The function automatically generates Pauli patterns for:
    1. All 'I's (e.g., 'IIII' for 4 features). Its strength is set to the negative
       of the count of the all-zero bit string ('0'*num_features) and does not
       accumulate from other bit strings.
    2. Exactly one 'Z' and the rest 'I's (e.g., 'ZIII', 'IZII', 'IIZI, IIIZ' for 4 features).
       These patterns accumulate strengths from all matching bit strings.

    A bit string contributes to a Pauli pattern's strength if it matches the pattern based on:
    - 'I' in the pattern matches both '0' and '1' in the bit string at that position (wildcard).
    - 'Z' in the pattern matches only '1' in the bit string at that position.

    Args:
        joint_counts (dict): A dictionary where keys are bit strings (e.g., '00', '01')
                             and values are their counts (e.g., {'00': 50, '11': 120}).
        min_ones (int): The minimum number of '1's required in a bit string
                        for it to be considered for analysis.
        num_features (int): The expected total length of the bit strings and Pauli patterns.

    Returns:
        SparsePauliOp: A Qiskit SparsePauliOp object where each Pauli term corresponds
                       to a generated Pauli pattern and its coefficient is the cumulative strength.
    """
    # Initialize a dictionary to store cumulative strengths for each target Pauli pattern
    cumulative_pauli_strengths = collections.defaultdict(float)

    # --- Generate the target Pauli patterns based on num_features ---
    generated_pauli_patterns = []

    # 1. Add the all 'I's pattern
    all_i_pattern = 'I' * num_features
    generated_pauli_patterns.append(all_i_pattern)
    cumulative_pauli_strengths[all_i_pattern] = 0.0 # Initialize its strength

    # 2. Add patterns with exactly one 'Z'
    # These will be processed for cumulative sums
    single_z_patterns = []
    for i in range(num_features):
        single_z_pattern_list = ['I'] * num_features
        single_z_pattern_list[i] = 'Z'
        single_z_pattern = "".join(single_z_pattern_list)
        generated_pauli_patterns.append(single_z_pattern)
        single_z_patterns.append(single_z_pattern) # Store separately for easier iteration
        cumulative_pauli_strengths[single_z_pattern] = 0.0 # Initialize its strength
    # --- End of pattern generation ---


    # Helper function to determine if a bit string matches a given Pauli pattern
    def does_bitstring_match_pauli_pattern(bit_string, pauli_pattern):
        # Length check for safety
        if len(bit_string) != len(pauli_pattern):
            return False

        for i in range(len(bit_string)):
            pattern_char = pauli_pattern[i]
            bit_char = bit_string[i]

            if pattern_char == 'I':
                # 'I' (Identity) in the pattern acts as a wildcard, matching both '0' and '1'
                continue
            elif pattern_char == 'Z':
                # 'Z' (Pauli-Z) in the pattern requires the bit string to have a '1' at this position
                if bit_char != '1':
                    return False # Mismatch if bit is '0' when 'Z' is required
            else:
                # Handle unexpected characters in the Pauli pattern (e.g., 'X', 'Y' not supported here)
                return False
        return True # All characters matched the pattern rules

    # Determine the all-zero bit string for special handling of 'IIII'
    all_zero_bit_string = '0' * num_features

    # Iterate through each bit string and its count in the histogram
    for bit_string, count in joint_counts.items():
        # Ensure the bit string's length matches the expected number of features
        if len(bit_string) != num_features:
            # Skip bit strings that do not conform to the expected length
            continue

        # Count the number of '1's in the current bit string
        num_ones = bit_string.count('1')

        # Only process bit strings that have at least `min_ones` '1's
        if num_ones >= min_ones:
            # Special handling for the 'all I's' pattern (e.g., 'IIII')
            if bit_string == all_zero_bit_string:
                # If the bit string is all zeros, set the strength of the 'all I's' pattern
                # to the negative of its count. This pattern will not accumulate further.
                cumulative_pauli_strengths[all_i_pattern] = -float(count)
                continue # Move to the next bit string, as '0'*N only affects 'I'*N in this specific way

            # For all other bit strings (not all zeros) and other Pauli patterns (single 'Z's)
            for target_pattern in single_z_patterns: # Only iterate through patterns with single 'Z's
                if does_bitstring_match_pauli_pattern(bit_string, target_pattern):
                    # If it matches, cumulatively add its count to the pattern's strength
                    cumulative_pauli_strengths[target_pattern] += float(count)

    # Convert defaultdict to a list of (label, complex(value)) tuples
    # to ensure compatibility with SparsePauliOp.from_list
    final_pauli_terms_list = [(label, complex(value)) for label, value in cumulative_pauli_strengths.items()]
    interaction_observable = SparsePauliOp.from_list(final_pauli_terms_list)
    return interaction_observable


def create_interaction_observable_from_histogram_simple(joint_counts, num_features, min_ones=0, standardize=False, rm_all_ones=False):
    """Creates a SparsePauliOp from joint histogram counts, 
       considering interactions with a minimum number of '1's.
    Args:
        joint_counts: A Counter object from create_joint_histogram.
        num_features: The total number of qubits.
        min_ones: The minimum number of '1's required in a bit string 
                  for the interaction to be included.
        standardize: If True, standardizes the counts before calculating strengths.
    Returns:
        A SparsePauliOp observable.
    """
    interaction_strength_list = []

    # Prepare counts for standardization if needed
    counts_array = np.array(list(joint_counts.values()))
    mean_count = np.mean(counts_array)
    if standardize:
        mean_count = np.mean(counts_array)
        std_count = np.std(counts_array)
        if std_count == 0:  # Handle case where all counts are the same
            standardized_counts = np.zeros_like(counts_array)
        else:
            standardized_counts = (counts_array - mean_count) / std_count
        
        # Create a dictionary to map bitstrings to standardized counts
        bitstring_to_std_count = dict(zip(joint_counts.keys(), standardized_counts))

    for bit_string, count in joint_counts.items():
        num_ones = bit_string.count('1')  # Count the number of '1's

        if num_ones == num_features and rm_all_ones:
            continue
        
        if num_ones >= min_ones:  # Consider only if at least min_ones '1's are present
            nodes = tuple(i for i, bit in enumerate(bit_string) if bit == '1')
            
            if standardize:
                strength = -bitstring_to_std_count[bit_string]  # Use standardized count
            else:
                #strength = -(count - mean_count)
                strength = (count - mean_count)

            pauli_string = ""
            for i in range(num_features):
                if i in nodes:
                    pauli_string += "Z"
                else:
                    pauli_string += "I"
            interaction_strength_list.append((pauli_string, strength))

    interaction_observable = SparsePauliOp.from_list(interaction_strength_list)
    return interaction_observable



# You also provided the create_interaction_observable_general function, which is separate
# and was not part of the problem. Keeping it as is.
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
        pauli_string = ""
        for i in range(num_features):
            if i in nodes:  # Check if the current qubit is in the interaction
                pauli_string += "Z"
            else:
                pauli_string += "I"
        interaction_strength_list.append((pauli_string, strength))

    interaction_observable = SparsePauliOp.from_list(interaction_strength_list)
    return interaction_observable

def create_interaction_observable_general2(interactions, num_features):
    """Creates a SparsePauliOp observable for generalized interactions.

    Args:
        interactions: A dictionary where keys are tuples of node indices
                      (e.g., (0, 1), (0, 2), (1,)) and
                      values are the corresponding interaction strengths.
                      For each tuple, a multi-qubit Pauli Z operator is constructed
                      where 'Z' acts on the specified qubits and 'I' (identity)
                      acts on all other qubits.
        num_features: The total number of qubits.

    Returns:
        A SparsePauliOp observable.
    """
    interaction_strength_list = []
    for nodes, strength in interactions.items():
        # Negate the strength: a positive input strength will result in a
        # negative coefficient in the Hamiltonian, encouraging minimization
        # of that interaction term's energy.
        strength = -strength 
        
        # Initialize the Pauli string with 'I' for all qubits
        # The string is built from MSB (highest qubit index) to LSB (qubit 0)
        pauli_chars = ['I'] * num_features 
        
        # For each qubit index specified in 'nodes', change 'I' to 'Z'
        for node_idx in nodes:
            if not (0 <= node_idx < num_features):
                raise ValueError(
                    f"Node index {node_idx} is out of bounds for {num_features} features."
                )
            # Qiskit Pauli string order is typically MSB...LSB (q_N-1 ... q_0)
            # So, to place 'Z' on qubit 'node_idx', we need to calculate its position
            # from the left (MSB side).
            # Example: for num_features=3, node_idx=0 (q0) -> position 2 (pauli_chars[2])
            #          node_idx=1 (q1) -> position 1 (pauli_chars[1])
            #          node_idx=2 (q2) -> position 0 (pauli_chars[0])
            pauli_chars[num_features - 1 - node_idx] = "Z"
            
        pauli_string = "".join(pauli_chars)
        
        interaction_strength_list.append((pauli_string, strength))

    # Create the SparsePauliOp from the list of (Pauli string, coefficient) pairs
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

        # Plot histogram:
        plot_histogram(counts, bar_labels=True, title=title).show()

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
    min_ones_obs=0, # Added as an explicit argument for flexibility
    optimizer_method="L-BFGS-B"
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
        histogram_data, num_qubits, min_ones=min_ones_obs, unobserved_punishment = 10, 
        normalization_offset = mean_offset)
    
    # interaction_observable = create_interaction_observable_from_histogram2(
    #     histogram_data, min_ones = min_ones_obs, num_features = num_qubits)
    
    # interaction_observable = create_interaction_observable_from_histogram_simple(
    #     histogram_data, num_qubits, min_ones_obs)
    
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
    create_parameter_dictionaries_from_circuit,
    create_interaction_observable_general # Explicitly pass this function too
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
        create_parameter_dictionaries_from_circuit (callable): A utility function that extracts
                                                                and separates static and variable
                                                                parameters from the given quantum circuit.
                                                                Expected signature:
                                                                `static_params, variable_params = func(circuit)`.
        create_interaction_observable_general (callable): A utility function to construct the
                                                          interaction observable based on the
                                                          defined interactions and total number of qubits.
                                                          Expected signature:
                                                          `observable = func(interactions_list, total_qubits)`.

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
        method="L-BFGS-B",
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
from qiskit.primitives import Sampler
from qiskit.circuit import QuantumCircuit, Parameter # Ensure these are imported if not already
from qiskit.quantum_info import Statevector # For StatevectorEstimator context


# --- MODIFIED cost_func_probability_matching ---
def cost_func_probability_matching(
    params_array,
    all_circuit_params_objects,
    target_probabilities: dict,
    circuit: QuantumCircuit,
    sampler: Sampler,
    shots: int = 1024,
    epsilon: float = 1e-9
) -> float:
    """
    Cost function for minimizing the distance between the quantum circuit's
    measurement probabilities and a target probability distribution (histogram).

    It calculates the Mean Squared Error (MSE) between the two distributions.

    Args:
        params_array (np.ndarray): The current numerical values of the *variable* parameters
                                   being optimized by scipy.optimize.minimize.
        all_circuit_params_objects (list): A list of all Parameter objects in the `circuit`,
                                           in the order they appear when `circuit.parameters` is called.
                                           This is used to correctly bind `params_array` to the circuit.
        target_probabilities (dict): A dictionary mapping bit strings (e.g., "011") to
                                     their target probability values. This is your normalized
                                     experimental histogram data.
        circuit (QuantumCircuit): The parameterized quantum circuit (ansatz) to be optimized.
        sampler (Sampler): A Qiskit Sampler primitive for running the circuit and getting probabilities.
        shots (int): The number of shots to run the quantum circuit for measurement.
        epsilon (float): A small value for numerical stability, especially if using KL divergence.

    Returns:
        float: The calculated Mean Squared Error (MSE) between the two probability distributions.
                Returns infinity if any target probability is non-zero but the ansatz yields zero,
                or if the sum of probabilities is not close to 1.
    """
    # 1. Bind parameters to the circuit
    bound_circuit = circuit.assign_parameters(dict(zip(all_circuit_params_objects, params_array)))
    
    # 2. Add measurements to the bound circuit for sampling
    # This creates a new circuit with measurements.
    # It's important that the original 'circuit' passed in has classical bits defined.
    circuit_with_measurements = bound_circuit.measure_all(inplace=False) # Use inplace=False to not modify bound_circuit directly

    # 3. Run the circuit on the sampler to get measurement probabilities
    try:
        job = sampler.run(circuit_with_measurements, shots=shots)
        result = job.result()
        ansatz_probabilities = result.quasi_dists[0].binary_probabilities()
    except Exception as e:
        print(f"Error during quantum circuit execution: {e}")
        return np.inf # Return a very high cost on error

    # 4. Calculate the Mean Squared Error (MSE)
    mse_cost = 0.0
    num_qubits = circuit.num_qubits # Number of qubits for iterating all states

    # Iterate through all possible bit strings (2^num_qubits) to cover all states
    for i in range(2**num_qubits):
        bit_string = format(i, '0' + str(num_qubits) + 'b')

        p_ansatz = ansatz_probabilities.get(bit_string, 0.0)
        p_target = target_probabilities.get(bit_string, 0.0)

        mse_cost += (p_ansatz - p_target)**2

    return mse_cost

# --- MODIFIED cost_func_wrapper_for_prob_matching ---
def cost_func_wrapper_for_prob_matching(
    variable_values,
    static_params: dict,
    circuit: QuantumCircuit,
    sampler: Sampler,
    variable_param_objects: list,
    target_probabilities: dict,
    shots: int = 1024
) -> float:
    """
    Wrapper function to adapt the probability matching cost function for scipy.optimize.minimize.
    It combines static and variable parameters, binds them to the circuit, and calls the
    `cost_func_probability_matching`.
    """
    all_params_for_binding = static_params.copy()
    for i, param_obj in enumerate(variable_param_objects):
        all_params_for_binding[param_obj] = variable_values[i]

    # Pass the full parameter dictionary to cost_func_probability_matching
    # The `all_circuit_params_objects` for `cost_func_probability_matching`
    # should be `list(circuit.parameters)` to ensure correct order for binding.
    return cost_func_probability_matching(
        params_array=list(all_params_for_binding.values()), # Pass values in circuit.parameters order
        all_circuit_params_objects=list(circuit.parameters), # Pass parameter objects in circuit.parameters order
        target_probabilities=target_probabilities,
        circuit=circuit,
        sampler=sampler,
        shots=shots
    )

# --- MODIFIED prob_dist_matching_solver ---
# (No changes needed here for the measurement fix, but including for completeness)
def prob_dist_matching_solver(
    target_histogram_data: Counter,
    circuit: QuantumCircuit,
    act_percentages: list,
    min_ones_obs: int = 0,
    optimizer_method: str = "L-BFGS-B",
    shots: int = 1024,
):
    num_qubits = circuit.num_qubits

    total_counts = sum(target_histogram_data.values())
    if total_counts == 0:
        raise ValueError("Target histogram data is empty or all counts are zero.")

    filtered_counts = Counter()
    for bit_string, count in target_histogram_data.items():
        if bit_string.count('1') >= min_ones_obs:
            filtered_counts[bit_string] = count

    if sum(filtered_counts.values()) == 0:
        print("Warning: Filtering resulted in an empty histogram. Using original counts for normalization.")
        target_probabilities = {bs: count / total_counts for bs, count in target_histogram_data.items()}
    else:
        target_probabilities = {bs: count / sum(filtered_counts.values()) for bs, count in filtered_counts.items()}

    print("\nTarget Probabilities (P_obs):", target_probabilities)

    static_params, variable_params_dict = create_parameter_dictionaries(circuit, act_percentages)

    initial_full_params = static_params.copy()
    initial_full_params.update(variable_params_dict)

    all_circuit_param_objects = list(circuit.parameters)
    x0_initial = np.array([initial_full_params.get(p, 0.0) for p in all_circuit_param_objects])

    print("\nInitial All Parameters (for optimizer x0):", x0_initial)
    print("Parameter Objects (ordered):", [p.name for p in all_circuit_param_objects])


    cost_history = []
    iteration_data = {'counter': 0}
    def callback_func(xk):
        current_cost = cost_func_wrapper_for_prob_matching(
            xk,
            static_params,
            circuit,
            all_circuit_param_objects,
            target_probabilities,
            shots
        )
        cost_history.append(current_cost)

        current_counter = iteration_data['counter']
        print_criteria = 20 if current_counter <= 100 else 100
        if current_counter % print_criteria == 0:
            print(f"Iteration {current_counter}: Current cost: {current_cost}")
        iteration_data['counter'] += 1

    print(f"\nStarting probability distribution matching optimization with method: {optimizer_method}")
    result = minimize(
        fun=lambda xk: cost_func_wrapper_for_prob_matching(
            xk,
            static_params,
            circuit,
            all_circuit_param_objects,
            target_probabilities,
            shots
        ),
        x0=x0_initial,
        method=optimizer_method,
        callback=callback_func
    )

    print("\nOptimization Result:")
    print(result)
    print(f"\nFinal Cost (MSE): {result.fun}")

    optimized_params_array = result.x
    optimized_full_params = {}
    for i, param_obj in enumerate(all_circuit_param_objects):
        optimized_full_params[param_obj] = optimized_params_array[i]

    print("\nOptimized Full Parameters:")
    for param, value in optimized_full_params.items():
        print(f"  {param.name}: {value}")

    final_bound_circuit = circuit.assign_parameters(optimized_full_params)
    
    # Ensure this circuit also has measurements for final sampling/plotting
    final_circuit_with_measurements = final_bound_circuit.measure_all(inplace=False)

    job_final = sampler.run(final_circuit_with_measurements, shots=shots)
    final_ansatz_probabilities = job_final.result().quasi_dists[0].binary_probabilities()
    print("\nFinal Optimized Ansatz Probabilities:", final_ansatz_probabilities)

    return result, optimized_full_params, cost_history, final_ansatz_probabilities

