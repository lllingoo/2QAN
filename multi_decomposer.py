# We can use three different decomposers
# Import cirq decomposition modules
import cirq
from CirqToQiskit import qiskit_to_cirq, cirq_to_qiskit
from decompose_to_fsim import decompose_circuit_to_fsim, decompose_into_fsim

# Import qiskit decomposition modules
from qiskit import transpile

# Import NuOp decomposition modules
from parallel_two_qubit_gate_decomposition_v2 import *
from gates_numpy import cnot_gate, fsim_gate, cphase_gate, xy_gate, get_gate_unitary_qiskit


def cirq_decompose(circ, bgate='syc'):
    # Decomposition into other two-qubit gates using cirq ['sycamore', 'xmon', 'xmon_partial_cz', 'sqrt_iswap', cirq.FSimGate(np.pi/2, 0)]
    cirq_circ = qiskit_to_cirq(circ.copy())
    cirq_circ = decompose_circuit_to_fsim(cirq_circ, fsim_gates=['sycamore'])
    qns = [q.index for q in circ.qubits]
    qiskit_circ = cirq_to_qiskit(cirq_circ, max(qns)+1)
    return qiskit_circ

def nuop_decompose(circ, params=[np.pi/2, np.pi/6], gname='syc', num_threads=1):
    """ Using the numerical decomposition NuOp to decompse two-qubit gates into any two-qubit gate set"""
    pgrp_fsim = ParallelGateReplacementPass([fsim_gate], [params], [gname], 1, 1, 1e-7)
    nuop_circ, nuop_fid, g_counts, all_gatecircs = pgrp_fsim.run(circ.copy(), num_threads=num_threads, exact_decom=True, max_num_layers=4, trials=5)
    return nuop_circ

def run_decompose(circ, bgate='syc', params=[np.pi/2, np.pi/6], num_threads=4):
    # syc [np.pi/2, np.pi/6], sqrtiswap [np.pi/4, 0], iswap [np.pi/2, 0]
    if bgate == 'cx':
        # Using qiskit to decompose circuits into CNOT gate set
        basis_gates = ['id', 'rz', 'u3', 'u2', 'cx', 'reset']
        circ_d = transpile(circ.copy(), basis_gates=basis_gates, optimization_level=3)
        cxs = circ_d.count_ops()['cx']
    elif bgate == 'syc':
        # Using cirq to decompose circuits into SYC/sqrt_iSWAP gate set
        circ_d = cirq_decompose(circ.copy(), bgate)
        cxs = circ_d.count_ops()['unitary']
    else:
        # Using NuOp to decompose circuits into any two-qubit gate set
        # bgate == 'iswap'
        # params = [np.pi/2, 0]
        circ_d = nuop_decompose(circ.copy(), params, bgate, num_threads)
        cxs = circ_d.count_ops()['unitary']
    return circ_d, cxs
