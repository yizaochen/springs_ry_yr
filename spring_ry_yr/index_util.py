from enmspring.graphs_bigtraj import StackMeanModeAgent
import numpy as np

class TableAgent:
    d_atomname_i_lst = {'A': ['C5', 'C4', 'C2', 'N1', 'C6'],
                        'G': ['N9', 'C4', 'N3', 'C2', 'N1'],
                        'T': ['N1', 'C2', 'N3', 'C4'],
                        'C': ['N1', 'C2', 'N3', 'C4', 'O2']}
    d_atomname_j_lst = {'T': ['C5', 'C4', 'C2', 'N3'],
                        'C': ['C5', 'C4', 'N3', 'C2'],
                        'A': ['N7', 'C5', 'C6', 'N6'],
                        'G': ['N7', 'C5', 'C6', 'O6', 'C4']}

    def __init__(self, host):
        self.host = host
        self.table = None

        self.s_agent = None

        self.idx_i_array = None
        self.idx_j_array = None

        self.lambda_array = None # eigenvalues
        self.n_node = None
        self.eigvector_mat = None
        self.n_mode = None

        self.d_idx = None
        self.d_idx_inverse = None
        self.resid_map = None
        self.strandid_map = None
        self.atomname_map = None

        self.d_seq = None

    def make_table_to_npy(self):
        d_table = dict()
        for idx_i, idx_j in zip(self.idx_i_array, self.idx_j_array):
            d_table[(idx_i, idx_j)] = self.get_spring(idx_i, idx_j)

    def get_spring(self, idx_i, idx_j):
        spring = Spring(idx_i, idx_j)
        spring.set_cgname(self.d_idx_inverse)
        spring.set_strandid(self.strandid_map)
        if spring.check_not_st:
            return spring
        spring.set_resid(self.resid_map)
        spring.check_resid()
        if spring.check_not_st:
            return spring
        spring.set_resname(self.d_seq)
        spring.set_RY_or_YR()
        spring.set_atomname(self.atomname_map)
        if spring.RY_or_YR == 'RY':
            spring.set_bondtype_RY()
        else:
            spring.set_bondtype_YR()
        return spring

    def read_table_from_npy(self):
        pass

    def initialize_s_agent(self):
        rootfolder = '/home/ytcdata/bigtraj_fluctmatch/500ns'
        interval_time = 500
        self.s_agent = StackMeanModeAgent(self.host, rootfolder, interval_time)
        self.s_agent.load_mean_mode_laplacian_from_npy()
        self.s_agent.eigen_decompose()
        self.s_agent.initialize_nodes_information()
        self.s_agent.set_benchmark_array()
        self.s_agent.set_strand_array()

        self.lambda_array = self.s_agent.w # eigenvalues
        self.n_node = self.s_agent.n_node
        self.eigvector_mat = self.s_agent.v
        self.n_mode = self.lambda_array.shape[0]
        self.idx_i_array, self.idx_j_array = np.triu_indices(self.n_node, k=1)

        self.d_idx = self.s_agent.d_idx
        self.d_idx_inverse = {y:x for x,y in self.d_idx.items()}
        self.resid_map = self.s_agent.resid_map
        self.strandid_map = self.s_agent.strandid_map


class Spring:
    def __init__(self, idx_i, idx_j):
        self.idx_i = idx_i
        self.idx_j = idx_j

        self.bond_type = None

        self.cgname_i = None
        self.cgname_j = None

        self.strandid_i = None
        self.strandid_j = None

        self.resid_i = None
        self.resid_j = None # resid_i < resid_j

        self.resname_i = None
        self.resname_j = None
        self.RY_or_YR = None

        self.atomname_i = None
        self.atomname_j = None
        self.bondname = None

    def set_cgname(self, d_idx_inverse):
        self.cgname_i = d_idx_inverse[self.idx_i]
        self.cgname_j = d_idx_inverse[self.idx_j]

    def set_strandid(self, strandid_map):
        self.strandid_i = strandid_map[self.cgname_i]
        self.strandid_j = strandid_map[self.cgname_j]
        if self.strandid_i != self.strandid_j:
            self.bond_type = 'not-st' # Cross-strand, Not belong to Stack

    def set_resid(self, resid_map):
        self.resid_i = resid_map[self.cgname_i]
        self.resid_j = resid_map[self.cgname_j]

    def set_resname(self, d_seq):
        self.resname_i = d_seq[self.strandid_i][self.resid_i-1]
        self.resname_j = d_seq[self.strandid_j][self.resid_j-1]

    def set_RY_or_YR(self):
        if self.resname_i in ['A', 'G']:
            self.RY_or_YR = 'RY'
        else:
            self.RY_or_YR = 'YR'

    def set_atomname(self, atomname_map):
        self.atomname_i = atomname_map[self.cgname_i]
        self.atomname_j = atomname_map[self.cgname_j]
        self.bondname = f'{self.atomname_i}-{self.atomname_j}'

    def set_bondtype_RY(self): 
        if self.resname_i == 'A':
            if self.bondname not in ['C5-C4', 'C4-C5', 'C2-C2', 'N1-N3', 'C6-C4']:
                self.bond_type = 'others'
            else:
                self.bond_type = self.bondname
        else:
            if self.bondname not in ['N9-C5', 'C4-C4', 'N3-C2', 'C2-C2', 'N1-N3']:
                self.bond_type = 'others'
            else:
                self.bond_type = self.bondname

    def set_bondtype_YR(self): 
        if self.resname_i == 'T':
            if self.bondname not in ['N1-N7', 'C2-C5', 'C2-N7', 'N3-C6', 'C4-N6']:
                self.bond_type = 'others'
            else:
                self.bond_type = self.bondname
        else:
            if self.bondname not in ['N1-N7', 'C2-C5', 'N3-C6', 'C4-O6', 'O2-C4']:
                self.bond_type = 'others'
            else:
                self.bond_type = self.bondname

    def check_resid(self):
        resid_diff = self.resid_j - self.resid_i
        if np.abs(resid_diff) != 1:
            self.bond_type = 'not-st' # Not belong to Stack
        elif resid_diff == -1: # make resid_i < resid_j
            self.swap_all()

    def check_not_st(self):
        if self.bond_type == 'not-st':
            return True
        else:
            return False

    def swap_all(self):
        self.idx_i, self.idx_j = self.idx_j, self.idx_i
        self.cgname_i, self.cgname_j = self.cgname_j, self.cgname_i
        self.strandid_i, self.strandid_j = self.strandid_j, self.strandid_i
        self.resid_i, self.resid_j = self.resid_j, self.resid_i 

    

    