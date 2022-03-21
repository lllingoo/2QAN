import numpy as np
import random
import copy

from bench_arch import BenchArch

class HeuristicMapper(BenchArch):
    """
    Formualte the qubit placement problem as a quadratic assignment problem (QAP) [1]
    and use the Tabu search heuristic to solve the QAP problem [2].
    **References:**
    [1] Lingling Lao, Dan E Browne. 2QAN: A quantum compiler for 2-local Hamiltonian simulation algorithms, 
    ISCA'22, arXiv:2108.02099.
    [2] https://github.com/zeman412/Tabu_Search_QAP_20
    """
    def __init__(self, benchmark, lattice_xy, coupling_map):
        """
        Args:
            benchmark: currently in qiskit circuit format
            lattice_xy: (x,y) if coupling_map is not given and the topology is grid
            coupling_map: connectivity between qubits such as 
            [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (0,7), (0,8), (1,9), (0,10)]
        """
        super().__init__(benchmark, lattice_xy, coupling_map)
        self.flow = np.array(copy.deepcopy(self.adjmat))
        self.dist = np.array(copy.deepcopy(self.distmat))
        self.idx = -1
        self.N = int(self.n_qbits*(self.n_qbits-1)/2)
        self.neighbors = np.zeros((self.N, self.n_qbits +2), dtype=int)

    def run(self, output=None, num_iter=200, lst_len=20):
        
        vp_map, cost = self.place_tabu(num_iter, lst_len)
        if output:
            output_placement = './maps/' + str(output) + '_vpmap_heuristic.txt'
            target = open(output_placement, "w")
            target.write("# name = %s \n" %output)
            target.write("# map_type = Heuristic Mapper \n")
            target.write("# Z = %d \n" % cost)
            target.write("VP MAP %s \n" % vp_map)
            target.close()
        self.reset()
        return vp_map, cost

    def assignmt_cost(self, sol):
        """Calculate the placement cost"""
        cost=0
        for i in range(self.n_qbits):
            for j in range(self.n_qbits):
                cost+=self.dist[i][j] *self.flow[sol[i]][sol[j]]
        return cost

    def swap_move(self, sol_n):
        for i in range (self.n_qbits):
            j=i+1
            for j in range(self.n_qbits):
                if i<j:
                    self.idx=self.idx+1
                    sol_n[j], sol_n[i] = sol_n[i], sol_n[j]
                    self.neighbors[self.idx, :-2] = sol_n
                    self.neighbors[self.idx, -2:] = [sol_n[i], sol_n[j]]
                    sol_n[i], sol_n[j] = sol_n[j], sol_n[i]

    def not_in_tabu(self, solution, tabu):
        not_found = False
        if not solution.tolist() in tabu:
            solution[0], solution[1] = solution[1], solution[0]
            if not solution.tolist() in tabu:
                not_found = True
        return not_found

    def place_tabu(self, num_iter, lst_len):
        """The tabu search heuristic"""
        curnt_sol = random.sample(range(self.n_qbits), self.n_qbits)
        best_soln = curnt_sol
        Tabu = []
        frequency = {}
        # print("Initial: %s cost %s " % (curnt_sol, assignmt_cost(curnt_sol)))
        while num_iter > 0:
            self.idx = -1
            self.swap_move(curnt_sol)  # make a move to self.neighbors
            cost = np.zeros((len(self.neighbors)))  # holds the cost of the self.neighbors
            for index in range(len(self.neighbors)):
                # evaluate the cost of the candidate self.neighbors
                cost[index] = self.assignmt_cost(self.neighbors[index, :-2])  
            rank = np.argsort(cost)  # sorted index based on  cost
            self.neighbors = self.neighbors[rank]

            for j in range(self.N):
                not_tabu = self.not_in_tabu(self.neighbors[j, -2:], Tabu)
                if (not_tabu):
                    curnt_sol = self.neighbors[j, :-2].tolist()
                    Tabu.append(self.neighbors[j, -2:].tolist())
                    if len(Tabu) > lst_len-1:
                        Tabu = Tabu[1:]
                    #frequency based
                    if not tuple(curnt_sol) in frequency.keys():
                        frequency[tuple(curnt_sol)] = 1 # set key->penality -> to One
                        if self.assignmt_cost(curnt_sol) < self.assignmt_cost(best_soln):
                            best_soln = curnt_sol
                    else:
                        cur_cost= self.assignmt_cost(curnt_sol) + frequency[tuple(curnt_sol)] # penalize by frequency
                        frequency[tuple(curnt_sol)] += 1   # increament the frequency for the current visit
                        if cur_cost < self.assignmt_cost(best_soln):
                            best_soln = curnt_sol
                    break
                #Aspiration
                elif self.assignmt_cost(self.neighbors[j, :-2]) < self.assignmt_cost(best_soln):
                    curnt_sol = self.neighbors[j, :-2].tolist()
                    Tabu.insert(0, Tabu.pop(Tabu.index(self.neighbors[j, -2:].tolist())))
                    if len(Tabu) > lst_len - 1:
                        Tabu = Tabu[1:]
                        # frequency based
                    if not tuple(curnt_sol) in frequency.keys():
                        frequency[tuple(curnt_sol)] = 1  # set key->penality -> to One
                        best_soln = curnt_sol
                    else:
                        cur_cost= self.assignmt_cost(curnt_sol) + frequency[tuple(curnt_sol)] # penalize by frequency
                        frequency[tuple(curnt_sol)] += 1   # increament the frequency for the current visit
                        if cur_cost < self.assignmt_cost(best_soln):
                            best_soln = curnt_sol
            num_iter -= 1
        fnl_cost = self.assignmt_cost(best_soln)-sum(sum(self.flow))
        rvp_map = {}
        for i in range(len(best_soln)):
            rvp_map[best_soln[i]] = i
        vp_map = {}
        for i in range(len(best_soln)):
            vp_map[i] = rvp_map[i]
        # print("Best placement %s cost: %s " % (vp_map, fnl_cost/2))
        return vp_map, fnl_cost/2
        
    def reset(self):
        self.benchmark = None
        self.n_qbits = 0
        self.pairs = []
        self.G_circ = None
        self.adjmat = None
        self.topology = None
        self.distmat = None
        self.pqbits_locs = {}
        self.locs_pqbits = {}
        self.flow = None
        self.dist = None
