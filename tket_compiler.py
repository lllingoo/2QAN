import numpy as np
import pytket
from pytket import *
from pytket.extensions.qiskit import qiskit_to_tk, tk_to_qiskit
from pytket.passes import (FullMappingPass,DefaultMappingPass, DecomposeMultiQubitsIBM, SequencePass, OptimisePhaseGadgets, BasePass,
                            DecomposeBoxes, SynthesiseIBM, FullPeepholeOptimise, RebaseIBM, PauliSimp)
from pytket.predicates import CompilationUnit
from pytket.routing import Architecture, SquareGrid
from pytket.device import Device

# We compare 2QAN with tket
class TketCompile(object):
    """Compile circuits using tket [1]
    [1] tket compiler: 
    Paper: Seyon Sivarajah, Silas Dilkes, Alexander Cowtan, Will Simmons, Alec Edgington, Ross Duncan. 
    t|ket⟩: a retargetable compiler for NISQ devices, QST21.
    Source code: https://github.com/CQCL/tket
    """  
    def __init__(self, benchmark, lattice_xy=None, decom=False, coupling_map=None):
        """
        Args:
            benchmark: currently in qiskit circuit format
            lattice_xy: (x,y) if coupling_map is not given and the topology is grid
            coupling_map: connectivity between qubits such as 
            [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (0,7), (0,8), (1,9), (0,10)]
        """
        qc = benchmark.copy()
        self.tk_circuit = qiskit_to_tk(qc)
        
        self.n_qbts = len(benchmark.qubits)
        if coupling_map:
            self.topology = Architecture(coupling_map)
        else:
            self.topology = SquareGrid(lattice_xy[1], lattice_xy[0])
        self.device = Device(self.topology)

    def placement(self):
        """Find qubit initial map"""
        circ = self.tk_circuit.copy()
        if self.n_qbts < 24:
            graph_device = pytket.routing.GraphPlacement(self.device)
        else:
            graph_device = pytket.routing.LinePlacement(self.device)
        vpmap = graph_device.get_placement_map(circ)
        init_map = {}
        for key in vpmap.keys():
            init_map[key.index[0]]=vpmap[key].index[0]

        final_map = self.fill_partial_mapping(init_map)
        return final_map

    def fill_partial_mapping(self, partial_map):
        final_map = partial_map.copy()
        device_qubits = [node.index[0] for node in self.topology.nodes]
        unused_dqs = [q for q in device_qubits if q not in partial_map.values()]
        for i in range(len(unused_dqs)):
            final_map[self.n_qbts+i]=unused_dqs[i]
        return final_map

    def route_circuit(self):
        # Do not perform gate decomposition
        circ = self.tk_circuit.copy()
        # Apply smart initial placement/mapping
        if self.n_qbts < 16:
            graph_device = pytket.routing.GraphPlacement(self.device)
        else:
            graph_device = pytket.routing.LinePlacement(self.device)
        graph_device.place(circ)
        
        basic_parameters = dict(decompose_swaps=False)
        routed_tk_circuit = pytket.routing.route(circ, self.device, **basic_parameters)
        qiskit_circuit = tk_to_qiskit(routed_tk_circuit)
        return qiskit_circuit

    def default_tket(self):
        # Decompose to cx gate set
        circ = self.tk_circuit.copy()

        if self.n_qbts < 20:
            dft = [DecomposeBoxes(), FullPeepholeOptimise(), DefaultMappingPass(self.device), SynthesiseIBM()]
        else:
            graph_device = pytket.routing.LinePlacement(self.device)
            dft = [DecomposeBoxes(), FullPeepholeOptimise(), FullMappingPass(self.device, graph_device), SynthesiseIBM()]
        dpass = SequencePass(dft)
        cu = CompilationUnit(circ)
        dpass.apply(cu)

        fnl_tk_circuit = cu.circuit
        qiskit_circuit = tk_to_qiskit(fnl_tk_circuit)
        return qiskit_circuit
