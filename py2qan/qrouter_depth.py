import numpy as np
import networkx as nx
import copy

from bench_arch import BenchArch
from heuristic_mapper import HeuristicMapper
from scheduler import Scheduler

from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator
from qiskit.circuit.library.standard_gates import SwapGate, RZZGate
from qiskit.extensions.unitary import UnitaryGate
from qiskit.converters import *


class QuRouter(BenchArch):
    """
    The permutation-aware routing and scheduling heuristics developed in 2QAN [1]
    **References:**
    [1] Lingling Lao, Dan E. Browne. 2QAN: A quantum compiler for 2-local Hamiltonian simulation algorithms, 
    ISCA'22, arXiv:2108.02099.
    """
    def __init__(self, qasm, init_map, coupling_map=None, verbose=False):
        """
        Args:
            qasm: circuit in OpenQASM format
            lattice_xy(tuple): (x,y) if coupling_map is not given and the topology is grid
            init_map(dict): the initial qubit map {circuit qubit index: device qubit index}
            coupling_map(list): connectivity between qubits such as 
            [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (0,7), (0,8), (1,9), (0,10)]
            verbose(bool): True if print out intermediate mapping info
        """
        super().__init__(qasm=qasm, coupling_map=coupling_map)

        self.init_map = init_map
        self.init_instrs, self.unroute_instrs = self.find_r_ur_instrs(self.pairs, self.init_map)
        self.route_instrs = copy.deepcopy(self.init_instrs)

        self.maps = {}
        self.map_instrs = {}
        self.map_moves = {}
        self.swap_count = 0
        self.dswap_count = 0
        self.verbose = verbose
        if verbose:
            print('Initial vp-map ', self.init_map)
            print('NN gate (qubit pairs) ', self.init_instrs)
            print('Un-routed gate (qubit pairs) ', self.unroute_instrs)

    def find_r_ur_instrs(self, unrouted_pairs, vpmap):
        """
        Find routed and unrouted gates in a qubit map
        """
        route_instrs = []
        unroute_instrs = []
        for qpair in unrouted_pairs:
            q0, q1 = self.find_qlocs(qpair, vpmap)
            if self.distmat[q0][q1] == 1:
                route_instrs.append(qpair)
            else:
                unroute_instrs.append(qpair)
        return route_instrs, unroute_instrs

    def find_qlocs(self, qpair, vpmap):
        """
        Find the device locations of a qubit pair
        """
        q0 = vpmap[qpair[0]]
        q1 = vpmap[qpair[1]]
        ql0 = min(q0, q1)
        ql1 = max(q0, q1)
        return ql0, ql1

    def check_nn(self, qpair, vpmap):
        q0, q1 = self.find_qlocs(qpair, vpmap)
        d = self.distmat[q0][q1]
        return d-1

    def dist_instrs(self, instrs, vpmap):
        """
        Calculate the distances of a set of instructions in a given qubit map
        """
        dist_instrs = {}
        instrs_dists = {}
        for qpair in instrs:
            q0, q1 = self.find_qlocs(qpair, vpmap)
            dist = int(self.distmat[q0][q1])
            instrs_dists[qpair] = dist
            if dist in dist_instrs.keys():
                dist_instrs[dist].append(qpair)
            else:
                dist_instrs[dist] = [qpair]
        dmin = min(dist_instrs.keys())
        dsum = sum(instrs_dists.values())
        return dist_instrs[dmin], dsum

    def reverse_dict(self, vpmap):
        rev_map = dict([reversed(i) for i in vpmap.items()])
        return rev_map

    def swap_q_map(self, q0, q1, map):
        cmap = copy.deepcopy(map)
        p0 = cmap[q0]
        p1 = cmap[q1]
        cmap[q0] = p1
        cmap[q1] = p0
        return cmap

    def qbt_overlap(self, new_move, old_moves):
        overlap = 0
        if old_moves:
            for mov in old_moves:
                if set(mov) & set(new_move):
                    overlap += 1
        return overlap

    def router(self):
        """The routing heuristic in 2QAN:
        For each qubit map, all NN gates will be routed regardless their order in the initial circuit
        One finds best movement set for each un-NN gate
        """
        route_instrs = copy.deepcopy(self.init_instrs)
        unroute_instrs = copy.deepcopy(self.unroute_instrs)
        instrs_canswap = copy.deepcopy(self.pairs)
        maps = {}
        current_map = copy.deepcopy(self.init_map)
        instr_sets = {}
        movements = {}
        map_id = 0
        new_cycles = [[]]
        def schedule_moves_insts(inst, cycles):
            new_cycles = copy.deepcopy(cycles)
            for i in range(len(new_cycles)):
                j = len(new_cycles) - 1 - i
                if self.qbt_overlap(inst, new_cycles[j]):
                    break
            if j+1 < len(new_cycles):
                new_cycles[j+1].append(inst)
            else:
                new_cycles.append([inst])
            return new_cycles, j

        while unroute_instrs:
            # Find the gates which are closest
            closest_unroute_ins, dsum = self.dist_instrs(unroute_instrs, current_map)
            pv_map = self.reverse_dict(current_map) 
            instr = closest_unroute_ins[0]

            # Find the shortest paths for this pair;
            # find all possible movement sets for each path;
            # find all possible first movement for each path, either from 
            # evaluate which movement leads to minimum overhead for all the unrouted gates.
            p0, p1 = self.find_qlocs(instr, current_map)
            paths = [p for p in nx.all_shortest_paths(self.topology, source=p0, target=p1)]
            tmaps = []
            moves = []
            dsums = []
            dps = []
            for path in paths:
                p0nn = path[1]
                p1nn = path[-2]
                move01 = [self.zz_tuple((pv_map[p0], pv_map[p0nn])), 
                            self.zz_tuple((pv_map[p1], pv_map[p1nn]))]
                tmap01 = [self.swap_q_map(pv_map[p0], pv_map[p0nn], current_map), 
                          self.swap_q_map(pv_map[p1], pv_map[p1nn], current_map)]
                dsum01 = [self.dist_instrs(unroute_instrs, tmap01[0])[1], 
                          self.dist_instrs(unroute_instrs, tmap01[1])[1]]
                if dsum01[0] == dsum01[1]:
                    cycles, dp0 = schedule_moves_insts(move01[0], new_cycles)
                    cycles, dp1 = schedule_moves_insts(move01[1], new_cycles)
                    if dp0 == dp1:
                        if move01[0] in instrs_canswap and move01[1] not in instrs_canswap:
                            iid = 0
                        elif move01[0] not in instrs_canswap and move01[1] in instrs_canswap:
                            iid = 1
                        else:
                            iid, ids = locate_min(dsum01)
                    else:
                        iid = np.argmin([dp0, dp1])
                    dp = [dp0,dp1][iid]
                else:
                    iid = np.argmin(dsum01)
                    cycles, dp = schedule_moves_insts(move01[iid], new_cycles)
                tmap = tmap01[iid]
                move = move01[iid]
                d = dsum01[iid]
                tmaps.append(tmap)
                moves.append(self.zz_tuple(move))
                dsums.append(d)
                dps.append(dp)
            
            dmin_id, ids = locate_min(dsums)
            if len(ids) >1:
                new_idps = {}
                for i in ids:
                    new_idps[i] = dps[i]
                dp_minid, idps = locate_min(list(new_idps.values()))
                if len(idps) > 1:
                    for i in idps:
                        if moves[list(new_idps.keys())[i]] in instrs_canswap:
                            dmin_id = i
                            break

            # TODO: if there are multiple minimal sets, 
            # select the one which leads to shortest depth, 
            # or select the one which has more movements on existing gate pairs.
            fmove = moves[dmin_id]
            if fmove in instrs_canswap:
                instrs_canswap.remove(fmove)
            current_map = self.swap_q_map(fmove[0], fmove[1], current_map)
            new_routed, unroute_instrs = self.find_r_ur_instrs(unroute_instrs, current_map)
            if map_id in movements.keys():
                movements[map_id].append(fmove)
            else:
                movements[map_id] = [fmove]

            new_cycles, dp = schedule_moves_insts(fmove, new_cycles)
            for inst in new_routed:
                new_cycles, dp = schedule_moves_insts(inst, new_cycles)
            if new_routed:
                maps[map_id] = current_map
                instr_sets[map_id] = new_routed
                map_id += 1

        if self.verbose:
            print('All sets of SWAPs ', movements)
            # print('All vp-maps after each SWAP set ', maps)
            print('New NN gate (qubit pairs) for each vp-map', instr_sets)
        self.maps = maps
        self.map_instrs = instr_sets
        self.map_moves = movements

    def critical_gates(self, unroute_instrs):
        """Find the gates in critical path"""
        qbt_gates = {}
        for qpair in unroute_instrs:
            for q in qpair:
                if q not in qbt_gates:
                    qbt_gates[q] = [qpair]
                else:
                    qbt_gates[q].append(qpair)
        values = list(qbt_gates.values())
        lens = [len(value) for value in values]
        keys = list(qbt_gates.keys())
        iid = np.argmax(lens)
        busy_q = keys[iid]
        cgates = values[iid]
        return cgates

    def graph_color(self, routed_instrs):
        # Schedule init_routed_instrs using graph coloring, all intructions could be exchangeable.
        G0 = self.instruction_graph(routed_instrs)
        LG_G0 = nx.line_graph(G0)
        Cl = nx.coloring.greedy_color(LG_G0, strategy="largest_first")
        cmax = max(Cl.values())
        order_g0 = {}
        for c in range(cmax+1): 
            order_g0[c] = []
            for ed in Cl:
                if Cl[ed] == c:
                    order_g0[c].append(self.zz_tuple(ed))
        new_init_routed = []
        for c in order_g0:
            new_init_routed += order_g0[c]
        return order_g0, new_init_routed

    def schedule_cycle(self, inst, cycle_insts, map_insts, vpmap):
        schedule_at = False
        nn = self.check_nn(inst, vpmap)
        if not nn:
            if not self.qbt_overlap(inst, cycle_insts):
                schedule_at = True
                map_insts.remove(inst)
                cycle_insts.append(inst)
        return map_insts, cycle_insts, schedule_at

    def scheduler(self):
        """
        Hybrid scheduling: using graph coloring algorithms to schedule the commutable gates and 
        using the normal dag-based approach to schedule uncommutable gates.
        """
        init_routed_instrs = copy.deepcopy(self.route_instrs)
        all_routed_instrs = copy.deepcopy(self.route_instrs)

        # Combine movement with circuit gate, and then specify movement types
        # In g_types, 1 represents a dressed movement (move+circuit gate), 
        # 0 represents move gate
        mv_types = {}
        g_types = {}
        swap_count = 0
        dswap_count = 0
        n_maps = len(self.map_moves)
        for i in range(len(self.map_moves)):
            mv_types[i] = {}
            swap_count += len(self.map_moves[i])
            for j in range(len(self.map_moves[i])):
                mv = self.map_moves[i][j]
                mv_types[i][mv] = 0
                if mv in self.pairs:
                    if mv not in g_types:
                        mv_types[i][mv] = 1
                        dswap_count += 1
                        g_types[mv] = 1
                        if mv in init_routed_instrs:
                            init_routed_instrs.remove(mv)
                        else:
                            for mn in range(n_maps):
                                if mv in self.map_instrs[mn]:
                                    self.map_instrs[mn].remove(mv)

        self.swap_count = swap_count
        self.dswap_count = dswap_count
        if self.verbose:
            print(f'Total number of SWAPs is {swap_count} and the dressed SWAP count is {dswap_count}')
            print(f'Total number of two-qubit gates is {len(self.instrs)+(swap_count-dswap_count)}')
            print(f'Total number of CX gates is {len(self.instrs)*2+3*(swap_count-dswap_count)+dswap_count}')

        # We first schedule the initially routed gates using graph coloring, 
        # this will help to reduce the circuit depth
        if init_routed_instrs:
            G0 = self.instruction_graph(init_routed_instrs)
            LG_G0 = nx.line_graph(G0)
            Cl = nx.coloring.greedy_color(LG_G0, strategy="largest_first")
            cmax = max(Cl.values())
            order_g0 = {}
            for c in range(cmax+1): 
                order_g0[c] = []
                for ed in Cl:
                    if Cl[ed] == c:
                        order_g0[c].append(self.zz_tuple(ed))
            new_init_routed = []
            for c in order_g0:
                new_init_routed += order_g0[c]
            init_routed_instrs =list(reversed(new_init_routed))
            init_routed_instrs =new_init_routed
        #End of graph coloring schedule

        if not n_maps:
            # no routing is required for this circuit
            return init_routed_instrs

        # Cycles, each cycle consists of the gates that are scheduled at this cycle
        cycles = [[]] 
        # Corresponding cycles, gtypes indicates the gate types for each instruction at each cycle
        gtypes = [[]] 
        # 2 represents an original circuit gate, 
        # 1 represents a move+circuit gate, 0 represents move gate
        cycle_maps = [self.maps[n_maps-1]]
        phy_cycles = [[]]
        dp = 0
        gn0 = len(self.instrs)+(swap_count-dswap_count) # The number of gates need to be scheduled
        gn1 = 0 # The number of gates have been scheduled
        while gn0 != gn1:
            # Check whether any gate that were routed undirectly can be scheduled at each cycle
            for i in range(n_maps):
                n = n_maps-1-i
                map_instrs = copy.deepcopy(self.map_instrs[n])
                for inst in map_instrs:
                    self.map_instrs[n], cycles[dp], scheduled = self.schedule_cycle(inst, cycles[dp], 
                                                                self.map_instrs[n], cycle_maps[dp])
                    if scheduled:
                        gtypes[dp].append(2)
                        phy_cycles[dp].append(self.find_qlocs(inst, cycle_maps[dp]))

            # Check whether any movement can be scheduled at each cycle
            for i in range(n_maps):
                n = n_maps-1-i
                moves = copy.deepcopy(self.map_moves[n])
                unscheduled = []
                for j in range(n, n_maps-1):
                    unscheduled += self.map_instrs[j]
                    if j < n_maps-1:
                        unscheduled += self.map_moves[j+1]
                for move in reversed(moves):
                    # Movement can be scheduled only if 
                    # the circuit gates that depend on it have been scheduled
                    if not self.qbt_overlap(move, unscheduled):
                        self.map_moves[n], cycles[dp], scheduled = self.schedule_cycle(move, cycles[dp], 
                                                                    self.map_moves[n], cycle_maps[dp])
                        if scheduled:
                            gtypes[dp].append(mv_types[n][move])
                            cycle_maps[dp] = self.swap_q_map(move[0], move[1], cycle_maps[dp])
                            phy_cycles[dp].append(self.find_qlocs(move, cycle_maps[dp]))

            # Check whether any initially routed gate can be scheduled at each cycle
            routeds = copy.deepcopy(init_routed_instrs)
            for inst in routeds:
                init_routed_instrs, cycles[dp], scheduled = self.schedule_cycle(inst, cycles[dp], 
                                                                init_routed_instrs, cycle_maps[dp])
                if scheduled:
                    gtypes[dp].append(2)
                    phy_cycles[dp].append(self.find_qlocs(inst, cycle_maps[dp]))

            current_map = cycle_maps[-1]
            gn1 = sum([len(gates) for gates in cycles])
            if gn0 != gn1:
                cycle_maps.append(current_map)
                cycles.append([])
                gtypes.append([])
                phy_cycles.append([])
                dp += 1
    
        if gn0 != gn1:
            print('Scheduled cycles', cycles)
            raise ValueError('Scheduled circuit doesnt cover all gates')
        ordered_all_instrs = []
        dressed_instrs = []
        ordered_all_instrs_phy = []
        for c in range(len(cycles)):
            i = len(cycles) - 1 - c
            ordered_all_instrs += cycles[i]
            dressed_instrs += gtypes[i]
            ordered_all_instrs_phy += phy_cycles[i]

        return ordered_all_instrs, dressed_instrs, ordered_all_instrs_phy
        
    def no_routing(self, bench):
        """
        Assuming all-to-all connectivity, only perform scheduling.
        """
        init_routed_instrs = self.pairs
        G0 = self.instruction_graph(init_routed_instrs)
        LG_G0 = nx.line_graph(G0)
        Cl = nx.coloring.greedy_color(LG_G0, strategy="largest_first")
        cmax = max(Cl.values())
        order_g0 = {}
        for c in range(cmax+1): 
            order_g0[c] = []
            for ed in Cl:
                if Cl[ed] == c:
                    order_g0[c].append(self.zz_tuple(ed))
        new_init_routed = []
        for c in order_g0:
            new_init_routed += order_g0[c]

        qc = QuantumCircuit(self.b_qbts, self.b_qbts)
        if self.q1_instrs[0]:
            for gate in self.q1_instrs[0]:
                qc.append(gate[0], [gate[1]])

        for qpair in new_init_routed:
            qc.append(self.instrs[self.zz_tuple(qpair)], qpair)
        
        if self.q1_instrs[1]:
            for gate in self.q1_instrs[1]:
                qc.append(gate[0], [gate[1]])

        return qc

    def construct_circ(self, ordered_all_instrs, dressed_instrs, 
                        ordered_all_instrs_phy, layers=1, msmt=False):
        """Construct the qiskit circuit after compilation, 
        for multiple trotter steps, we only perform compilation for the first layer,
        for odd layers, we use the same circuit as first layer,
        for even layers, we reverse the circuit.
        """
        swap_mat = SwapGate().to_matrix()
        new_circ = QuantumCircuit(self.n_qbits, self.b_qbts)

        if self.q1_instrs[0]:
            for gate in self.q1_instrs[0]:
                new_circ.append(gate[0], [self.init_map[gate[1]]])
        
        zz_id = 0
        dzz = 0
        for l in range(layers):
            for gid in range(len(ordered_all_instrs)):
                if (l%2):
                    gid = len(ordered_all_instrs) - 1 - gid
                qpair = ordered_all_instrs[gid]
                phqs = list(ordered_all_instrs_phy[gid])
                # apply operations on the hardware qubits 
                if dressed_instrs[gid] == 2:
                    new_circ.append(self.instrs[self.zz_tuple(qpair)], phqs)
                    zz_id += 1
                elif dressed_instrs[gid] == 1:
                    u_mat = np.matmul(self.instrs[qpair].to_matrix(), swap_mat)
                    if isinstance(self.instrs[qpair], RZZGate):
                        angle = self.instrs[qpair].params[0]
                        new_circ.unitary(Operator(u_mat), phqs, label='dZZ'+str(angle))
                    else:
                        new_circ.unitary(Operator(u_mat), phqs, label='d'+self.instrs[qpair].name)
                    dzz += 1
                elif dressed_instrs[gid] == 0:
                    new_circ.swap(phqs[0], phqs[1])
                else:
                    raise ValueError('unsupported instruction types')

            if self.maps:
                if (l%2):
                    fnl_vpmap = self.init_map
                else:
                    fnl_vpmap = self.maps[list(self.maps.keys())[-1]]
            else:
                fnl_vpmap = self.init_map

            for gate in self.q1_instrs[1]:
                new_circ.append(gate[0], [fnl_vpmap[gate[1]]])
        if msmt:
            # The measurement outcome goes to classical bits 
            # which are indexed by their virtual qubit indices
            for cv in range(self.b_qbts):
                new_circ.measure(fnl_vpmap[cv], cv)
        return new_circ

    def construct_qaoa(self, ordered_all_instrs, dressed_instrs, ordered_all_instrs_phy, 
                        layers=1, gammas=None, betas=None, msmt=False, init_map=None, maps=None):
        swap_mat = SwapGate().to_matrix()
        new_circ = QuantumCircuit(self.n_qbits, self.b_qbts)

        if self.q1_instrs[0]:
            for gate in self.q1_instrs[0]:
                new_circ.append(gate[0], [init_map[gate[1]]])
        
        zz_id = 0
        dzz = 0
        for l in range(layers):
            for gid in range(len(ordered_all_instrs)):
                if (l%2):
                    gid = len(ordered_all_instrs) - 1 - gid
                qpair = ordered_all_instrs[gid]
                phqs = list(ordered_all_instrs_phy[gid])
                # apply operations on the hardware qubits 
                if dressed_instrs[gid] == 2:
                    # new_circ.append(self.instrs[self.zz_tuple(qpair)], phqs)
                    new_circ.rzz(2*gammas[l], phqs[0], phqs[1])
                    zz_id += 1
                elif dressed_instrs[gid] == 1:
                    rzz_mat = RZZGate(2*gammas[l]).to_matrix()
                    u_mat = np.matmul(rzz_mat, swap_mat)
                    new_circ.unitary(Operator(u_mat), phqs, label='dZZ'+str(2*gammas[l]))
                    dzz += 1
                elif dressed_instrs[gid] == 0:
                    new_circ.swap(phqs[0], phqs[1])
                else:
                    raise ValueError('unsupported instruction types')

            if maps:
                if (l%2):
                    fnl_vpmap = init_map
                else:
                    fnl_vpmap = maps[list(maps.keys())[-1]]
            else:
                fnl_vpmap = init_map

            for gate in self.q1_instrs[1]:
                new_circ.rx(2*betas[l], [fnl_vpmap[gate[1]]])
        if msmt:
            # The measurement outcome goes to classical bits which are indexed 
            # by their virtual qubit indices
            for cv in range(self.b_qbts):
                new_circ.measure(fnl_vpmap[cv], cv)
        return new_circ

    def run_qaoa(self, layers=1, gammas=None, betas=None, msmt='False'):
        self.router()
        ordered_all_instrs, dressed_instrs, ordered_all_instrs_phy = self.scheduler()
        new_circ = self.construct_qaoa(ordered_all_instrs, dressed_instrs, 
                            ordered_all_instrs_phy, layers, gammas, betas, msmt, self.init_map, self.maps)
        return new_circ, (layers*self.swap_count, layers*self.dswap_count)

    def run(self, layers=1, msmt='False'):
        self.router()
        ordered_all_instrs, dressed_instrs, ordered_all_instrs_phy = self.scheduler()
        new_circ = self.construct_circ(ordered_all_instrs, dressed_instrs, 
                                        ordered_all_instrs_phy, layers, msmt)

        return new_circ, (layers*self.swap_count, layers*self.dswap_count)


def locate_min(a):
    smallest = min(a)
    ids = [index for index, element in enumerate(a) if smallest == element]
    return np.random.choice(ids), ids
