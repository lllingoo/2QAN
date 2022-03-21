# Circuit conversion between Cirq and Qiskit
import numpy as np
import cirq # version 0.8.0
from cirq import LineQubit
from cirq.circuits.qasm_output import QasmUGate, QasmTwoQubitGate
from cirq import protocols
from cirq.google.optimizers.convert_to_sycamore_gates import swap_rzz, rzz

from qiskit import QuantumCircuit
# from fsim import fsim
from qiskit.quantum_info import Operator
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.dagcircuit import DAGCircuit
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
import qiskit
from qiskit import Aer, execute

from qiskit.transpiler.passes.optimization.optimize_1q_gates import Optimize1qGates
optimise1qgates = Optimize1qGates()

def fsim(theta, phi):
    matrix = np.array([[1,             0,                         0, 0],
                       [0,       np.cos(theta), -1j * np.sin(theta), 0],
                       [0, -1j * np.sin(theta),       np.cos(theta), 0],
                       [0, 0, 0,         np.cos(phi) - 1j*np.sin(phi)]], 
                      dtype=complex)
    return Operator(matrix)

def rzz_unitary(theta):
    return np.array([[np.exp(-1j*theta/2), 0, 0, 0],
                     [0, np.exp(1j*theta/2), 0, 0],
                     [0, 0, np.exp(1j*theta/2), 0],
                     [0, 0, 0, np.exp(-1j*theta/2)]], dtype=complex)


def extract_fsim_params(mat):
    v00 = mat[0,0] # in case there is a global phase in mat
    v12 = mat[1,2]/v00 #theta term
    v33 = mat[3,3]/v00 #phi term
    theta = np.arcsin(np.real((1j)*v12))
    phi = np.real((1j)*np.log(v33))
    return theta, phi

def extract_zz_params(mat):
    v11 = mat[1,1] 
    v33 = mat[3,3] # in case there is a global phase in mat
    t = np.real((-1j)*np.log(v11/v33)/np.pi)
    global_shift = np.real((-1j)*np.log(v33)/(t*np.pi))
    return t, global_shift

from cirq.circuits.qasm_output import QasmUGate, QasmTwoQubitGate
def qiskit_to_cirq(circuit):
    n_qubits = len(circuit.qubits)
    n_cbits = len(circuit.clbits)
    dag = circuit_to_dag(circuit)
    output = cirq.Circuit()
    for gate in dag.topological_op_nodes():
        op = gate.op
        qargs = gate.qargs
        cargs = gate.cargs
        p = gate.op.params
        # print(p)
        qidx = list(reversed([q.index for q in qargs]))
        #Single-qubit gate
        if len(qargs) == 1 and len(cargs) == 0:            
            if isinstance(op, qiskit.circuit.library.standard_gates.rx.RXGate):
                output.append(cirq.rx(p[0])(LineQubit(qidx[0])))
            elif isinstance(op, qiskit.circuit.library.standard_gates.ry.RYGate):
                output.append(cirq.ry(p[0])(LineQubit(qidx[0])))
            elif isinstance(op, qiskit.circuit.library.standard_gates.rz.RZGate):
                output.append(cirq.rz(p[0])(LineQubit(qidx[0])))
            elif isinstance(op, qiskit.circuit.library.standard_gates.h.HGate):
                output.append(cirq.H(LineQubit(qidx[0])))  
            elif isinstance(op, qiskit.circuit.library.standard_gates.x.XGate):
                output.append(cirq.X(LineQubit(qidx[0])))        
            elif isinstance(op, qiskit.circuit.library.standard_gates.y.YGate):
                output.append(cirq.Y(LineQubit(qidx[0])))           
            elif isinstance(op, qiskit.circuit.library.standard_gates.z.ZGate):
                output.append(cirq.Z(LineQubit(qidx[0])))                  
            elif isinstance(op, qiskit.circuit.library.standard_gates.SGate):
                output.append(cirq.S(LineQubit(qidx[0])))  
            elif isinstance(op, qiskit.circuit.library.standard_gates.SdgGate):
                output.append((cirq.S**-1)(LineQubit(qidx[0])))  
            elif isinstance(op, qiskit.circuit.library.standard_gates.TGate):
                output.append(cirq.T(LineQubit(qidx[0])))  
            elif isinstance(op, qiskit.circuit.library.standard_gates.TdgGate):
                output.append((cirq.T**-1)(LineQubit(qidx[0])))                  
            elif isinstance(op, qiskit.circuit.library.standard_gates.u1.U1Gate):
                output.append(QasmUGate(0,0,p[0]/np.pi)(LineQubit(qidx[0])))                  
            elif isinstance(op, qiskit.circuit.library.standard_gates.u2.U2Gate):
                output.append(QasmUGate(1/2, p[0]/np.pi,p[1]/np.pi)(LineQubit(qidx[0])))      
            elif isinstance(op, qiskit.circuit.library.standard_gates.u3.U3Gate):
                output.append(QasmUGate(theta=p[0]/np.pi,phi=p[1]/np.pi,lmda=p[2]/np.pi)(LineQubit(qidx[0])))   
            else:
            #     print("Warning: Unknown gate type", op, op.label)
                output.append(cirq.MatrixGate(op.to_matrix()).on(LineQubit(qidx[0])))    
        #Two-qubit gate             
        elif len(qargs) == 2 and len(cargs) == 0:
            if isinstance(op, qiskit.circuit.library.standard_gates.rzz.RZZGate):
                # mat = rzz_unitary(p[0])
                output.append(cirq.ZZPowGate(exponent=p[0], global_shift=-0.5).on(LineQubit(qidx[0]), LineQubit(qidx[1]))) 
                # output.append(cirq.MatrixGate(mat).on(LineQubit(qidx[0]), LineQubit(qidx[1])))
            elif isinstance(op, qiskit.circuit.library.standard_gates.rxx.RXXGate):
                output.append(cirq.XXPowGate(exponent=p[0], global_shift=-0.5).on(LineQubit(qidx[0]), LineQubit(qidx[1]))) 
            elif isinstance(op, qiskit.circuit.library.standard_gates.ryy.RYYGate):
                output.append(cirq.YYPowGate(exponent=p[0], global_shift=-0.5).on(LineQubit(qidx[0]), LineQubit(qidx[1]))) 
            elif isinstance(op, qiskit.circuit.library.standard_gates.SwapGate):
                output.append(cirq.SWAP(LineQubit(qidx[0]), LineQubit(qidx[1])))                
            elif isinstance(op, qiskit.circuit.library.standard_gates.CXGate):
                output.append(cirq.CNOT(LineQubit(qidx[1]), LineQubit(qidx[0])))                
            elif isinstance(op, qiskit.extensions.unitary.UnitaryGate):
                if op.label == 'SYC':
                    output.append(cirq.google.SYC(LineQubit(qidx[0]), LineQubit(qidx[1])))
                elif op.label == 'CZ':
                    output.append(cirq.CZ(LineQubit(qidx[0]), LineQubit(qidx[1])))  
                elif op.label == 'CZPowGate':
                    theta, phi = extract_fsim_params(op.to_matrix())
                    output.append(cirq.CZPowGate(exponent=-phi/np.pi).on(LineQubit(qidx[0]), LineQubit(qidx[1])))
                elif op.label[:3] == 'dZZ':
                    # decompose a ZZswap into syc and rzz
                    t = -np.pi/(float(op.label[3:])/2-np.pi/24+np.pi/4)
                    output.append(cirq.google.SYC(LineQubit(qidx[0]), LineQubit(qidx[1])))
                    output.append(rzz(t, LineQubit(qidx[0]), LineQubit(qidx[1]))) 
                    # output.append(cirq.CZPowGate(exponent=-t).on(LineQubit(qidx[0]), LineQubit(qidx[1]))) 
                elif op.label in ['sqrtISWAP', 'ISWAP','ISwapPowGate']:
                    theta, phi = extract_fsim_params(op.to_matrix())
                    output.append(cirq.ISwapPowGate(exponent=-2*theta/(np.pi)).on(LineQubit(qidx[0]), LineQubit(qidx[1])))                    
                elif op.label == 'ZZPowGate':
                    t, shift = extract_zz_params(op.to_matrix())
                    output.append(cirq.ZZPowGate(exponent=t).on(LineQubit(qidx[0]), LineQubit(qidx[1])))   
                elif op.label in ['FsimGate', 'XY', 'fSWAP']:
                    theta, phi = extract_fsim_params(op.to_matrix())
                    output.append(cirq.FSimGate(theta, phi).on(LineQubit(qidx[0]), LineQubit(qidx[1])))                 
                elif op.label in ['YYPowGate', 'XXPowGate']:
                    output.append(cirq.MatrixGate(op.to_matrix()).on(LineQubit(qidx[0]), LineQubit(qidx[1])))    
                else:
                    # print("Warning: Unknown gate type", op, op.label)
                    output.append(cirq.MatrixGate(op.to_matrix()).on(LineQubit(qidx[0]), LineQubit(qidx[1])))
            else:
                # print("Warning: Unknown gate type", op, op.label)
                output.append(cirq.MatrixGate(op.to_matrix()).on(LineQubit(qidx[0]), LineQubit(qidx[1])))                  
    return output

def u3_decompose(op):
    mat = cirq.unitary(op)
    lam = QasmUGate.from_matrix(mat).lmda
    theta = QasmUGate.from_matrix(mat).theta
    phi = QasmUGate.from_matrix(mat).phi
    return theta, phi, lam

def cirq_to_qiskit(circuit, qubit_count=0, cbit_count=0, measure_qubit_labels=[], measure_cbit_labels=[], unitary_label=""):
    assert len(measure_qubit_labels) == len(measure_cbit_labels)
    if qubit_count != 0 and cbit_count != 0:
        output_circuit = QuantumCircuit(qubit_count, cbit_count)
    if cbit_count == 0:
        output_circuit = QuantumCircuit(qubit_count)
    for op in circuit.all_operations():
        qubits = list(reversed([int(str(q)) for q in op.qubits]))        
        if op.gate == cirq.google.SYC:
            output_circuit.unitary(fsim(np.pi/2,np.pi/6), qubits, label='SYC')  
        elif isinstance(op.gate, cirq.ops.common_gates.CZPowGate):
            e = op.gate.exponent
            if e in [1, -1]:
                output_circuit.unitary(fsim(0,-e*np.pi), qubits, label='CZ')
            else:
                output_circuit.unitary(fsim(0,-e*np.pi), qubits, label='CZPowGate')
        elif isinstance(op.gate, cirq.ops.swap_gates.ISwapPowGate):
            e = op.gate.exponent
            if e in [1, -1]:
                output_circuit.unitary(fsim(-e*np.pi/2, 0), qubits, label='ISWAP')
            elif e in [0.5, -0.5]:
                output_circuit.unitary(fsim(-e*np.pi/2, 0), qubits, label='sqrtISWAP')
            else:
                output_circuit.unitary(fsim(-e*np.pi/2, 0), qubits, label='ISwapPowGate')
        elif isinstance(op.gate, cirq.ops.ZZPowGate):
            this_unitary = cirq.unitary(op)
            output_circuit.unitary(Operator(this_unitary), qubits, label='ZZPowGate')
        elif isinstance(op.gate, cirq.ops.YYPowGate):
            this_unitary = cirq.unitary(op)
            output_circuit.unitary(Operator(this_unitary), qubits, label='YYPowGate')
        elif isinstance(op.gate, cirq.ops.XXPowGate):
            this_unitary = cirq.unitary(op)
            output_circuit.unitary(Operator(this_unitary), qubits, label='XXPowGate')
        elif isinstance(op.gate, cirq.FSimGate):
            theta = op.gate.theta
            phi = op.gate.phi
            if phi == 0:
                label = 'XY'
            elif theta == 0.5*np.pi and phi == np.pi:
                label = 'fSWAP'
            else:
                label = 'FSimGate'
            output_circuit.unitary(fsim(theta, phi), qubits, label=label)
        elif op.gate == cirq.SWAP:
            output_circuit.swap(*qubits)
        elif op.gate == cirq.CNOT:
            output_circuit.cx(*list(reversed(qubits)))            
        elif isinstance(op.gate,cirq.PhasedXZGate) or isinstance(op.gate, cirq.circuits.qasm_output.QasmUGate):
            theta, phi, lam = u3_decompose(op)
            output_circuit.u3(np.pi*theta, np.pi*phi, np.pi*lam, qubits)
            # output_circuit.u3(2*np.pi*theta, 2*np.pi*phi, 2*np.pi*lam, qubits)
        elif op.gate == cirq.H:
            output_circuit.h(qubits)
        elif isinstance(op.gate, cirq.ops.common_gates.XPowGate):
            rads = op.gate._exponent
            if op.gate._global_shift ==-0.5:
                output_circuit.rx(np.pi*rads, qubits)
            else:
                print(f'Cannot convert XPowGate with global_shift {op.gate._global_shift}')
        elif isinstance(op.gate, cirq.ops.common_gates.YPowGate):
            rads = op.gate._exponent
            if op.gate._global_shift ==-0.5:
                output_circuit.ry(np.pi*rads, qubits)
                print(f'Cannot convert YPowGate with global_shift {op.gate._global_shift}')
        elif isinstance(op.gate, cirq.ops.common_gates.ZPowGate):
            rads = op.gate._exponent
            output_circuit.u1(np.pi*rads, qubits)
        else:
            if qubits != []:
                if len(qubits) == 1:
                    unitary_label = 'su4'
                this_unitary = cirq.unitary(op)
                output_circuit.unitary(Operator(this_unitary), qubits, label=unitary_label)
    if measure_qubit_labels:
        output_circuit.measure(measure_qubit_labels, measure_cbit_labels)
    optimized_dag = optimise1qgates.run(circuit_to_dag(output_circuit))
    output_circuit = dag_to_circuit(optimized_dag)
    return output_circuit

def print_unitary(circuit):
    simulator = Aer.get_backend('unitary_simulator')
    result = execute(circuit, simulator).result()
    unitary = result.get_unitary(circuit)
    return unitary
