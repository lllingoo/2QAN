import numpy as np
import qiskit
from qiskit.extensions.unitary import UnitaryGate
from qiskit.circuit.library.standard_gates import RZZGate
def cphase_gate(theta):
    return np.matrix([
        [
            1,
            0,
            0,
            0
        ],
        [
            0,
            1,
            0,
            0
        ],
        [
            0,
            0,
            1,
            0
        ],
        [
            0,
            0,
            0,
            np.cos(theta) + 1j*np.sin(theta)
        ]])
        
def cnot_gate():
     return np.matrix([
        [
            1,
            0,
            0,
            0
        ],
        [
            0,
            1,
            0,
            0
        ],
        [
            0,
            0,
            0,
            1
        ],
        [
            0,
            0,
            1,
            0
        ]])   
    
def fsim_gate(theta, phi):
    return np.matrix([
        [
            1,
            0,
            0,
            0
        ],
        [
            0,
            np.cos(theta),
            -1j * np.sin(theta),
            0
        ],
        [
            0,
            -1j * np.sin(theta),
            np.cos(theta),
            0
        ],
        [
            0,
            0,
            0,
            np.cos(phi) - 1j*np.sin(phi)
        ]])

def xy_gate(theta):
    return np.matrix([
        [
            1,
            0,
            0,
            0
        ],
        [
            0,
            np.cos(theta/2),
            1j * np.sin(theta/2),
            0
        ],
        [
            0,
            1j * np.sin(theta/2),
            np.cos(theta/2),
            0
        ],
        [
            0,
            0,
            0,
            1
        ]
    ])

def cz_gate():
    return np.matrix([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, -1]])

def rzz_unitary(theta):
    return np.array([[np.exp(-1j*theta/2), 0, 0, 0],
                     [0, np.exp(1j*theta/2), 0, 0],
                     [0, 0, np.exp(1j*theta/2), 0],
                     [0, 0, 0, np.exp(-1j*theta/2)]], dtype=complex)
    
def get_gate_unitary_qiskit(gate_op):
    # Let's assume all the default unitary matrices in qiskit, which has different endianness from our 
    # convention, so we will need to reverse the qubit order when we apply our decomposition pass.
    # if isinstance(gate_op, qiskit.circuit.library.standard_gates.x.CXGate):
    #     return cnot_gate()
    # elif isinstance(gate_op, qiskit.circuit.library.standard_gates.z.CZGate):
    #     return cz_gate()
    if isinstance(gate_op, UnitaryGate):
        return gate_op.to_matrix()
    elif isinstance(gate_op, RZZGate):
        return rzz_unitary(gate_op.params[0])
    else:
        return gate_op.to_matrix()
        
    