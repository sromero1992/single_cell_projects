def create_interaction_observable_from_histogram(
    joint_counts: Counter,
    num_features: int,
    min_ones: int = 0,
    max_ones: int = None, # New parameter: maximum number of '1's
    unobserved_punishment: float = 1.0,
    normalization_offset: float = 0.0
) -> SparsePauliOp:
    """Creates a SparsePauliOp representing a diagonal Hamiltonian in the computational basis.
    The energy of each basis state is derived from its count in the histogram,
    favoring states that meet `min_ones` and `max_ones` criteria and are observed,
    and punishing others.

    Args:
        joint_counts: A Counter object of observed bit string counts.
        num_features: The total number of qubits.
        min_ones: The minimum number of '1's required in a bit string for it to be
                  considered 'favorable' if observed.
        max_ones: The maximum number of '1's allowed in a bit string for it to be
                  considered 'favorable' if observed. If None, there is no upper limit.
        unobserved_punishment: The positive energy value assigned to unobserved bit strings,
                               or observed bit strings that don't meet `min_ones`/`max_ones` criteria.
                               These states will be 'punished' (avoided by optimizer).
        normalization_offset: An optional offset to subtract from counts before calculating energy.
                               Useful for centering counts, e.g., for 50/50 cases.

    Returns:
        A SparsePauliOp observable.
    """
    
    state_energies = {} # Maps bit_string -> desired_energy (e.g., "011" -> -500)
    
    # Iterate through all possible bit strings (2^num_features)
    for i in range(2**num_features):
        bit_string = format(i, '0' + str(num_features) + 'b') # Generates MSB...LSB (q_N-1 ... q_0)
        num_ones = bit_string.count('1') # Count '1's for criteria

        energy_b = unobserved_punishment # Default energy is punishment

        # Check if the number of ones is within the specified range
        is_within_range = (num_ones >= min_ones)
        if max_ones is not None:
            is_within_range = is_within_range and (num_ones <= max_ones)

        # If the state is within the desired range AND it was observed in the histogram
        if is_within_range and bit_string in joint_counts:
            # Assign energy proportional to -count (or -log(count) if that's preferred)
            energy_b = -float(joint_counts[bit_string] - normalization_offset)
            # If you want to use log smoothing, uncomment the line below and comment the one above:
            # energy_b = -np.log(float(joint_counts[bit_string]) + 1.0) # Added 1.0 to avoid log(0) if counts can be 0

        state_energies[bit_string] = energy_b

    # 2. Generate all possible Z-only Pauli strings
    all_pauli_z_strings = []
    for pauli_tuple in itertools.product('IZ', repeat=num_features):
        pauli_string = "".join(pauli_tuple)
        all_pauli_z_strings.append(pauli_string)
    
    # 3. Convert these state_energies (E_b) into coefficients (c_P) for Pauli strings
    pauli_term_coefficients = {}

    for pauli_str in all_pauli_z_strings: # Iterate through all canonical Pauli-Z strings
        coeff_p = 0.0
        # Sum over all basis states 'b' (all 2^N bit strings)
        for bit_string, energy_b in state_energies.items():
            # Get eigenvalue of Pauli string P for basis state |b>
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
