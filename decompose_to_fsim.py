from typing import (
    Sequence,
    Union,
    Any,
    List,
    Iterator,
    TypeVar,
    Iterable,
    Optional,
)

import numpy as np

import cirq
from cirq import optimizers
# from cirq.optimizers import decompose_cphase_into_two_fsim

def decompose_into_fsim(
        interaction: Union[cirq.Operation, cirq.Gate, np.ndarray, Any],
        fsim_gate: Union[str, cirq.FSimGate]='sycamore',
        qubits: Sequence[cirq.Qid] = None) -> cirq.Circuit:
    """Decompose operations into one of the FSimGates and single-qubit rotations.
    Args:
        interaction: The two qubit operation to synthesize. This can either be
            a cirq object (such as a gate, operation, or circuit) or a raw numpy
            array specifying the 4x4 unitary matrix.
        fsim_gate: The type of FSimGates decompositions and corresponding two-qubit gates are:
            'sycamore': SYC = fsim(theta=pi/2, phi=pi/6))
            'xmon_partial_cz': czpow(t) = fsim(theta=0, phi=t*pi))
            'xmon': cz = fsim(theta=0, phi=pi)
            'sqrt_iswap': sqrt_iswap = fsim(theta=pi/2, phi=0)
            cirq.FSimGate(theta,phi) 
        qubits: The qubits that the resulting operations should apply the
            desired interaction to. If not set then defaults to either the
            qubits of the given interaction (if it is a `cirq.Operation`) or
            else to `cirq.LineQubit.range(2)`.

    Returns:
        A list of operations implementing the desired two qubit unitary. The
        list will include four operations of the given fsim gate, various single
        qubit operations, and a global phase operation.
    """
    if isinstance(interaction, cirq.Operation):
        qubits = interaction.qubits
    elif isinstance(interaction, cirq.Circuit):
        for op in interaction.all_operations():
            qubits = op.qubits
    else:
        # input interaction is a gate
        if qubits is None:
            qubits = cirq.devices.LineQubit.range(2)
        interaction = interaction.on(*qubits)

    if len(qubits) != 2:
        raise ValueError(f'Expected a pair of qubits, but got {qubits!r}.')

    if isinstance(interaction, cirq.Circuit):
        input_c = interaction
    else:
        input_c = cirq.Circuit(interaction)

    if fsim_gate == 'sycamore' or fsim_gate == cirq.google.SYC:
        decomposed_c = cirq.google.optimized_for_sycamore(circuit=input_c, 
                                                          optimizer_type='sycamore')
    elif fsim_gate == 'xmon':
        decomposed_c = cirq.google.optimized_for_sycamore(circuit=input_c, 
                                                          optimizer_type='xmon')
    elif fsim_gate == 'xmon_partial_cz':
        decomposed_c = cirq.google.optimized_for_sycamore(circuit=input_c, 
                                                          optimizer_type='xmon_partial_cz')
    elif fsim_gate == 'sqrt_iswap':
        decomposed_c = cirq.google.optimized_for_sycamore(circuit=input_c, 
                                                          optimizer_type='sqrt_iswap')
    elif isinstance(fsim_gate, cirq.FSimGate):
        if isinstance(interaction.gate, (cirq.XXPowGate, cirq.ZZPowGate, cirq.YYPowGate, cirq.CZPowGate)):
            try:
                # try to decompose a xx/yy/zz interaction into two fsim gates
                ops = cirq.two_qubit_matrix_to_operations(qubits[0], qubits[1], cirq.unitary(interaction), True)
                decomposed_c = cirq.Circuit()
                for op in ops:
                    new_ops = op
                    if len(op.qubits) == 2:
                        new_ops = cirq.decompose_cphase_into_two_fsim(op.gate, fsim_gate=fsim_gate, qubits=op.qubits)
                    decomposed_c.append(new_ops)
                decomposed_c = optimise_mix_decomposed_circuit(decomposed_c)
                return decomposed_c
            except:
                pass
        # try to decompose a two-qubit interaction into four fsim gates
        if not 3 / 8 * np.pi <= abs(fsim_gate.theta) <= 5 / 8 * np.pi:
            raise ValueError('Must have 3π/8 ≤ |fsim_gate.theta| ≤ 5π/8')
        if abs(fsim_gate.phi) > np.pi / 4:
            raise ValueError('Must have abs(fsim_gate.phi) ≤ π/4')
        decomposed_c = cirq.decompose_two_qubit_interaction_into_four_fsim_gates_via_b(
                        interaction, fsim_gate=fsim_gate)
        decomposed_c = optimise_mix_decomposed_circuit(decomposed_c)
    else:
        raise ValueError(f'Unexpected native two-qubit gate {fsim_gate}.')

    return decomposed_c


def optimise_mix_decomposed_circuit(circuit: cirq.Circuit) -> cirq.Circuit:
    """Last-step optimisation: merging single-qubit rotations if a decomposed circuit consists of different two-qubit natives. 
    """
    optimizers.merge_single_qubit_gates_into_phxz(circuit)
    optimizers.EjectPhasedPaulis().optimize_circuit(circuit)
    optimizers.EjectZ().optimize_circuit(circuit)
    try:
        optimizers.DropNegligible().optimize_circuit(circuit)
    except:
        pass

    return circuit

def decompose_circuit_to_fsim(circuit: cirq.Circuit, fsim_gates: List[Union[str, cirq.FSimGate, Any]]=['sycamore', 'xmon', 'xmon_partial_cz','sqrt_iswap']) -> cirq.Circuit:
    """Decompose a given circuit in cirq into one that only has fsim gates. 
    If multiple fsim gates are available, for each two-qubit gate, choose the decomposition
    that has least number of fsim gates and least number of gates.
    TODO, add fidelity check to choose the decomposition that has highest fidelity
    """
    def find_gatecount(circuit, metric='count2qgate'):
        if metric == 'count2qgate':
            count = 0
            for op in circuit.all_operations():
                if len(op.qubits) == 2:
                    count += 1
            return count
        elif metric == 'countallgate':
            count = 0
            for op in circuit.all_operations():
                count += 1
            return count
        else:

            raise ValueError('Currently only find the decomposition with minimum number of all gates or two-qubit gates')

    new_circuit = cirq.Circuit()
    for op in circuit.all_operations():
        if len(op.qubits) == 1:
            new_circuit.append(op)
        elif len(op.qubits) == 2:
            decomposed_gates = []
            gate_counts = []
            for fsim_gate in fsim_gates:
                try:
                    new_decom = decompose_into_fsim(interaction=op, fsim_gate=fsim_gate)
                    decomposed_gates.append(new_decom)
                except:
                    pass
            if decomposed_gates:
                for cir in decomposed_gates:
                    gate_counts.append(find_gatecount(cir, metric='count2qgate'))
                num_min = gate_counts.count(min(gate_counts))
                if num_min != 1:
                    gate_counts2 = []
                    for i, item in enumerate(gate_counts):
                        if item == min(gate_counts):
                            gate_counts2.append(find_gatecount(decomposed_gates[i], metric='countallgate'))
                        else: 
                            gate_counts2.append(float(np.inf))
                else:
                    gate_counts2 = gate_counts
                best_cir = decomposed_gates[gate_counts2.index(min(gate_counts2))]
                new_circuit.append(best_cir)
            else:
                raise ValueError(f'Could not find any decompositions for gate {op} based on the given native gates {fsim_gates}')
    new_circuit = optimise_mix_decomposed_circuit(new_circuit)
    return new_circuit
    

if __name__ == '__main__':
    # testing
    q0, q1 = cirq.LineQubit.range(2)
    interaction = cirq.XXPowGate(exponent=1/16)(q0, q1)
    # from scipy import stats
    # mat = stats.unitary_group.rvs(4, random_state=np.random.default_rng())
    # print(mat)
    # interaction = cirq.MatrixGate(mat).on(cirq.LineQubit(0), cirq.LineQubit(1))
    circuit = cirq.Circuit(interaction)
    fsim_gate = cirq.FSimGate(-np.pi/2, 0)
    # fsim_gate = 'xmon'
    new_circuit = decompose_into_fsim(interaction, fsim_gate=fsim_gate)
    # new_circuit = decompose_circuit_to_fsim(circuit)
    print(new_circuit)

    # circuits = []
    # for op in [cirq.ZZPowGate(exponent=1/4)(q0, q1)]:
    #     for fsim in ['sycamore', 'xmon', 'xmon_partial_cz', 'sqrt_iswap']:
    #         circuit = decompose_into_fsim(interaction=op, fsim_gate=fsim)
    #         circuits.append(circuit)
    #         print(cirq.equal_up_to_global_phase(cirq.unitary(op), cirq.unitary(circuit)))
    # circuits[0].append(circuits[4])
    # print('The number of operations after merge is {0}\n'
    #       .format(len([str(op) for op in circuits[0].all_operations()])))
    # optimise_mix_decomposed_circuit(circuits[0])
    # print('The number of operations after optimisation is {0}'
    #       .format(len([str(op) for op in circuits[0].all_operations()])))