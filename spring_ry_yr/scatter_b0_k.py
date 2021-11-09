from os import path
import numpy as np
import matplotlib.pyplot as plt
from enmspring.graphs_bigtraj import StackMeanModeAgent

class DataAgent:

    def __init__(self, host):
        self.host = host
        self.table = self.read_table_from_npy()

        self.s_agent = None
        self.n_node = None
        self.idx_i_array = None
        self.idx_j_array = None
        self.eigvector_mat = None

    def make_k_b0_data(self):
        data_obj = Data(self.host)
        data_obj.process_spring(self.idx_i_array, self.idx_j_array, self.s_agent.laplacian_mat, self.s_agent.b0_mean_mat, self.table)
        k_b0_data_npy = self.get_k_b0_data_npy()
        np.save(k_b0_data_npy, data_obj)
        print(f'Save data_obj into {k_b0_data_npy}')

    def read_k_b0_data(self):
        k_b0_data_npy = self.get_k_b0_data_npy()
        return np.load(k_b0_data_npy, allow_pickle='TRUE').item()

    def get_k_b0_data_npy(self):
        datafolder = '/home/yizaochen/codes/dna2021paper/notebooks/SI/prominet_modes_st/chord_diagram_st/heterogeneous/data'
        return path.join(datafolder, f'{self.host}_k_b0_data.npy')

    def read_table_from_npy(self):
        table_npy = self.get_table_npy()
        return np.load(table_npy, allow_pickle='TRUE').item()

    def get_table_npy(self):
        datafolder = '/home/yizaochen/codes/dna2021paper/notebooks/SI/prominet_modes_st/chord_diagram_st/heterogeneous/data'
        return path.join(datafolder, f'{self.host}_table.npy')

    def initialize_s_agent(self):
        rootfolder = '/home/ytcdata/bigtraj_fluctmatch/500ns'
        interval_time = 500
        self.s_agent = StackMeanModeAgent(self.host, rootfolder, interval_time)
        self.s_agent.load_mean_mode_laplacian_from_npy()
        self.s_agent.load_b0_mean_std_from_npy()
        self.s_agent.eigen_decompose()
        self.s_agent.initialize_nodes_information()
        self.s_agent.set_benchmark_array()
        self.s_agent.set_strand_array()

        self.n_node = self.s_agent.n_node
        self.idx_i_array, self.idx_j_array = np.triu_indices(self.n_node, k=1)

class Data:
    k_criteria = 0.5

    def __init__(self, host):
        self.host = host
        self.k_array_RY = list()
        self.b0_array_RY = list()
        self.k_array_YR = list()
        self.b0_array_YR = list()

    def process_spring(self, idx_i_array, idx_j_array, K_mat, b0_mat, table):
        for idx_i, idx_j in zip(idx_i_array, idx_j_array):
            k_value = K_mat[idx_i, idx_j]
            b0_value = b0_mat[idx_i, idx_j]
            if k_value < self.k_criteria:
                continue
            bondtype = table[(idx_i, idx_j)].bond_type
            if bondtype == 'not-st':
                continue
            if table[(idx_i, idx_j)].RY_or_YR == 'RY':
                self.k_array_RY.append(k_value)
                self.b0_array_RY.append(b0_value)
            elif table[(idx_i, idx_j)].RY_or_YR == 'YR':
                self.k_array_YR.append(k_value)
                self.b0_array_YR.append(b0_value)

    def __repr__(self):
        return f'{self.host}: k-b0 Data'

class Scatter:
    # Fontsize
    lbfz = 6 # xy labels
    tickfz = 6

    def __init__(self, host):
        self.host = host
        self.k_b0_data = self.read_k_b0_data()

    def plot_main(self, figsize):
        fig, ax = plt.subplots(figsize=figsize)
        self.scatter_RY(ax)
        self.scatter_YR(ax)
        xlabel, ylabel = self.get_xlabel_ylabel()
        ax.set_xlabel(xlabel, fontsize=self.lbfz)
        ax.set_ylabel(ylabel, fontsize=self.lbfz)
        ax.tick_params(axis='both', labelsize=self.tickfz, length=1.5, pad=1)
        return fig, ax

    def get_xlabel_ylabel(self):
        xlabel = r'$b^{0}$ ($\mathrm{\AA}$)'
        ylabel = r'$k$ (kcal/mol/$\mathrm{\AA^{2}}$)'
        return xlabel, ylabel

    def scatter_RY(self, ax):
        x_array = self.k_b0_data.b0_array_RY
        y_array = self.k_b0_data.k_array_RY
        ax.scatter(x_array, y_array, s=2, color='blue')

    def scatter_YR(self, ax):
        x_array = self.k_b0_data.b0_array_YR
        y_array = self.k_b0_data.k_array_YR
        ax.scatter(x_array, y_array, s=2, color='red')

    def read_k_b0_data(self):
        k_b0_data_npy = self.get_k_b0_data_npy()
        return np.load(k_b0_data_npy, allow_pickle='TRUE').item()

    def get_k_b0_data_npy(self):
        datafolder = '/home/yizaochen/codes/dna2021paper/notebooks/SI/prominet_modes_st/chord_diagram_st/heterogeneous/data'
        return path.join(datafolder, f'{self.host}_k_b0_data.npy')