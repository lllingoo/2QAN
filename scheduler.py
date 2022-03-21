# from qiskit.dagcircuit import DAGCircuit
from qiskit.circuit import Barrier, Measure

class Scheduler(object):
    """Scheduling utilities"""
    def __init__(self, dag):
        """
        Parse the dag of a qiskit circuit
        Args:
            dag: the directed acyclic graph of the input circuit
        """
        self.dag = dag
        self.all_tp_gates()
        self.asap_moments = []
        self.alap_moments = []

    def find_moment(self, op, layers):
        k = 0
        for l in range(len(layers)):
            if op in layers[l]:
                k = l
        return k
      
    def cx_tuple(self, gate):
        """
        Representation for two-qubit gate
        """
        physical_q_0 = gate.qargs[0].index
        physical_q_1 = gate.qargs[1].index
        r_0 = min(physical_q_0, physical_q_1)
        r_1 = max(physical_q_0, physical_q_1)
        return (r_0, r_1)

    def singleq_tuple(self, gate):
        """
        Representation for single-qubit gate
        """
        physical_q_0 = gate.qargs[0].index
        tup = (physical_q_0,)
        return tup

    def gate_tuple(self, gate):
        """
        Representation for gate
        """
        if len(gate.qargs) == 2:
            return self.cx_tuple(gate)
        else:
            return self.singleq_tuple(gate)
            
    def all_tp_gates(self):
        all_qbts = []
        all_gates = []
        self.msmt = []
        for gate in list(self.dag.topological_nodes()):
            if not gate.type == 'op':
                continue
            if isinstance(gate.op, Barrier):
                continue
            if isinstance(gate.op, Measure):
                self.msmt.append(gate)
                continue
            all_gates.append(gate)
            for q in list(self.gate_tuple(gate)):
                if q not in all_qbts:
                    all_qbts.append(q)
        self.all_gates = tuple(all_gates)
        self.all_qbts = tuple(all_qbts)
        all_gates_r = []
        for i in range(len(self.all_gates)):
            gate = self.all_gates[len(self.all_gates)-1-i]
            all_gates_r.append(gate)
        self.all_gates_r = tuple(all_gates_r)

    def find_predecessors(self, ops):
        predecessors = []
        for gate in list(self.dag.predecessors(ops)):
            if not gate.type == 'op':
                continue
            if isinstance(gate.op, Barrier):
                continue
            predecessors.append(gate)
        return predecessors

    def find_successors(self, ops):
        successors = []
        for gate in list(self.dag.successors(ops)):
            if not gate.type == 'op':
                continue
            if isinstance(gate.op, Barrier):
                continue
            successors.append(gate)
        return successors

    def find_asap_moments(self):
        asap_moments = [[]]
        for gate in self.all_gates:
            predecessors = self.find_predecessors(gate)
            if len(predecessors):
                prd_k = []
                for op in predecessors:
                    prd_k.append(self.find_moment(op, asap_moments))
                k = max(prd_k)
                if len(asap_moments) == (k+1):
                    asap_moments.append([])
                    asap_moments[k+1] = [gate]
                else:
                    asap_moments[k+1].append(gate)
            else:
                asap_moments[0].append(gate)
        for moment in asap_moments:
            self.asap_moments.append(tuple(moment))

        self.asap_moments = tuple(self.asap_moments)

    def find_alap_moments(self):
        alap_moments = [[]]
        for gate in self.all_gates_r:
            successors = self.find_successors(gate)
            if len(successors):
                prd_k = []
                for op in successors:
                    prd_k.append(self.find_moment(op, alap_moments))
                k = max(prd_k)
                if len(alap_moments) == (k+1):
                    alap_moments.append([])
                    alap_moments[k+1] = [gate]
                else:
                    alap_moments[k+1].append(gate)
            else:
                alap_moments[0].append(gate)
        for md in range(len(alap_moments)):
            moment = alap_moments[len(alap_moments)-1-md]
            self.alap_moments.append(tuple(moment))
        self.alap_moments = tuple(self.alap_moments)
        return self.alap_moments