def create_interaction_observable_from_histogram(
    joint_counts: Counter,
    num_features: int,
    min_ones: int = 0, # Kept as per your original signature
    # Added new parameters for fine-grained control over energy assignment
    unobserved_punishment: float = 1.0,
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
            #energy_b = -np.log(float(joint_counts[bit_string]) + 1.0)
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
