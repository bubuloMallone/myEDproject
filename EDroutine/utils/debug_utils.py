import numpy as np
from utils.hamiltonian_utils import apply_hamiltonian_serial, apply_hamiltonian_parallel
from utils.input_utils import parse_lattice_file, parse_matrix_file, parse_input_string, parse_tags_factors_matrices
from utils.observable_utils import translation_op, hadamard_op

def check_duality_symmetry_fly(preprocessed_hamiltonian_dict, sites, movelist):

    hilbert_size = 2 ** len(sites)
    state = np.random.randn(hilbert_size) + 1j * np.random.randn(hilbert_size)
    norm_squared = state.conj().T @ state
    state /= np.sqrt(norm_squared)  # Normalize the state
    # state1 = apply_hamiltonian(hadamard_op(state, len(sites)), preprocessed_hamiltonian_dict)  # Test the Hamiltonian application
    # state2 = hadamard_op(apply_hamiltonian(state, preprocessed_hamiltonian_dict), len(sites))  # Test the Hamiltonian application
    # print("[H,S] = 0 --> ", np.allclose(state1-state2, np.zeros_like(state)))  # Check if the two states are close
    state1 = apply_hamiltonian_serial(hadamard_op(translation_op(state, len(sites), movelist), len(sites)), preprocessed_hamiltonian_dict)  # Test the Hamiltonian application
    state2 = hadamard_op(translation_op(apply_hamiltonian_serial(state, preprocessed_hamiltonian_dict), len(sites), movelist), len(sites))  # Test the Hamiltonian application
    # state1 = apply_hamiltonian_serial(hadamard_op(state, len(sites)), preprocessed_hamiltonian_dict)  # Test the Hamiltonian application
    # state2 = hadamard_op(apply_hamiltonian_serial(state, preprocessed_hamiltonian_dict), len(sites))  # Test the Hamiltonian application
    print("[H,S] = 0 --> ", np.allclose(state1, state2))  # Check if the two states are close
    state3 = hadamard_op(translation_op(hadamard_op(translation_op(state, len(sites), movelist), len(sites)), len(sites), movelist), len(sites))
    print("(S^2 = Id) --> ", np.allclose(state, state3))  # Check if the two states are close
    state4 = translation_op(translation_op(state, len(sites), movelist), len(sites), movelist)
    print("(T^2 = Id) --> ", np.allclose(state, state4))  # Check if the two states are close
    state5 = hadamard_op(hadamard_op(state, len(sites)), len(sites))
    print("(Had^2 = Id) --> ", np.allclose(state, state5))  # Check if the two states are close


def check_duality_symmetry_sparse(hamiltonian_sparse, sites, movelist):

    hilbert_size = 2 ** len(sites)
    state = np.random.randn(hilbert_size) + 1j * np.random.randn(hilbert_size)
    norm_squared = state.conj().T @ state
    state /= np.sqrt(norm_squared)  # Normalize the state

    # Test the Hamiltonian application
    state1 = hadamard_op(translation_op(state, len(sites), movelist), len(sites))  
    state1 = hamiltonian_sparse @ state1
    # Test the Hamiltonian application
    state2 = hamiltonian_sparse @ state
    state2 = hadamard_op(translation_op(state2, len(sites), movelist), len(sites))  
    
    print("[H,S] = 0 --> ", np.allclose(state1, state2))  # Check if the two states are close
    state3 = hadamard_op(translation_op(hadamard_op(translation_op(state, len(sites), movelist), len(sites)), len(sites), movelist), len(sites))
    print("(S^2 = Id) --> ", np.allclose(state, state3))  # Check if the two states are close
    state4 = translation_op(translation_op(state, len(sites), movelist), len(sites), movelist)
    print("(T^2 = Id) --> ", np.allclose(state, state4))  # Check if the two states are close
    state5 = hadamard_op(hadamard_op(state, len(sites)), len(sites))
    print("(Had^2 = Id) --> ", np.allclose(state, state5))  # Check if the two states are close



def check_single_bit_hadamard():
    # check hadamard op is correct
    dim=2
    simple_state=np.zeros(dim, dtype=np.complex128)
    basis_state = 1
    simple_state[basis_state] = 1.0 + 0.0j
    out_state = hadamard_op(simple_state, 1)
    print(f"H|{basis_state}>: ", out_state)


def check_hamiltonian_action_on_basis(preprocessed_hamiltonian_dict, sites):
    hilbert_size = 2 ** len(sites)
    simple_state = np.zeros(hilbert_size, dtype=np.complex128)
    basis_state = 7
    simple_state[7] = 1.0 + 0.0j
    # out_state = hadamard_op(translation_op(simple_state, len(sites), movelist), len(sites))
    out_state = apply_hamiltonian_serial(simple_state, preprocessed_hamiltonian_dict)
    print("non-vanishing coefficients in S|0...0>: ", np.nonzero(out_state))
    print("coefficients in S|0...0>: ", out_state)
    print("true coefficient: ", 1/np.sqrt(2**len(sites)))
    print(f"coefficients in H|{basis_state}>: ", out_state)


def print_hamiltonian_terms(preprocessed_hamiltonian_dict):
    """
    Print the Hamiltonian terms for debugging purposes.
    """
    for arity, terms_list in preprocessed_hamiltonian_dict.items():
        print(f"type {arity}-BODY: {len(terms_list)} terms \n")
        for term in terms_list:
            sites_tuple, factor, matrix_elements_array, matrix_values_array = term
            print(f"Sites: {sites_tuple} \nFactor: {factor}")
            print("Matrix elements (input -> output : value):")
            for i in range(len(matrix_elements_array)):
                input_local_state = matrix_elements_array[i, 0]
                output_local_state = matrix_elements_array[i, 1]
                matrix_value = matrix_values_array[i]
                print(f"{input_local_state} -> {output_local_state} : {matrix_value}")
            print("")
        print("--------------------")