import matplotlib.pyplot as plt
import numpy as np

class BarBoxPlot:
    width = 1
    tickfz = 6

    d_bondtype_lst_RY = {
        'atat_21mer': [('C5', 'C4'), ('C4', 'C5'), ('C2', 'C2'), ('N1', 'N3'), ('C6', 'C4')],
        'gcgc_21mer': [('N9', 'C5'), ('C4', 'C4'), ('N3', 'C2'), ('C2', 'C2'), ('N1', 'N3')]
    }
    d_bondtype_lst_YR = {
        'atat_21mer': [('N1', 'N7'), ('C2', 'C5'), ('C2', 'N7'), ('N3', 'C6'), ('C4', 'N6')],
        'gcgc_21mer': [('N1', 'N7'), ('C2', 'C5'), ('N3', 'C6'), ('C4', 'O6'), ('O2', 'C4')]
    }

    abbr_host = {'atat_21mer': 'TATA', 'gcgc_21mer': 'CpG'}

    def __init__(self, host, mode_id, k_mat):
        self.host = host
        self.mode_id = mode_id
        self.k_mat = k_mat

    def plot_main(self, figsize):
        fig, axes = plt.subplots(nrows=2, figsize=figsize, facecolor='white')
        d_axes = self.get_d_axes(axes)
        self.bar_RY(d_axes['counts'])
        self.bar_YR(d_axes['counts'])
        self.boxplot_RY(d_axes['boxplot'])
        self.boxplot_YR(d_axes['boxplot'])
        self.set_xtick(d_axes)
        self.set_xlim(d_axes)
        self.set_ylabel(d_axes)
        self.set_title(d_axes)
        return fig, d_axes

    def bar_RY(self, ax):
        x_array = self.get_xarray_RY()
        height_array = np.random.randint(1, 10, size=6)
        ax.bar(x_array, height_array, self.width, color='blue', edgecolor='white', linewidth=0.5)

    def bar_YR(self, ax):
        x_array = self.get_xarray_YR()
        height_array = np.random.randint(1, 10, size=6)
        ax.bar(x_array, height_array, self.width, color='red', edgecolor='white', linewidth=0.5)

    def boxplot_RY(self, ax):
        x_array = self.get_xarray_RY()
        data = [np.random.rand(10) for i in range(6)]
        ax.boxplot(data, positions=x_array, whis=1., widths=0.6, boxprops={'facecolor': 'steelblue'}, medianprops={'color': 'black'}, patch_artist=True)

    def boxplot_YR(self, ax):
        x_array = self.get_xarray_YR()
        data = [np.random.rand(10) for i in range(6)]
        ax.boxplot(data, positions=x_array, whis=1., widths=0.6, boxprops={'facecolor': 'steelblue'}, medianprops={'color': 'black'}, patch_artist=True)

    def set_ylabel(self, d_axes):
        d_axes['counts'].set_ylabel('Number of springs', fontsize=6)
        d_axes['boxplot'].set_ylabel('k (kcal/mol/Å$^2$)', fontsize=6)

    def get_xarray_RY(self):
        return range(1, 7)

    def get_xarray_YR(self):
        return range(8, 14)

    def set_title(self, d_axes):
        abbr = self.abbr_host[self.host]
        d_axes['counts'].set_title(f'{abbr}: Mode {self.mode_id}', fontsize=8)

    def set_xlim(self, d_axes):
        d_axes['counts'].set_xlim(0, 14)
        d_axes['boxplot'].set_xlim(0, 14)

    def set_xtick(self, d_axes):
        d_axes['counts'].set_xticks(list(self.get_xarray_RY())+list(self.get_xarray_YR()))
        d_axes['boxplot'].set_xticks(list(self.get_xarray_RY())+list(self.get_xarray_YR()))
        xticklabels = self.get_xticklabels()
        d_axes['counts'].set_xticklabels(xticklabels)
        d_axes['boxplot'].set_xticklabels(xticklabels)
        d_axes['counts'].tick_params(axis='both', labelsize=self.tickfz, length=1.5, pad=1)
        d_axes['boxplot'].tick_params(axis='both', labelsize=self.tickfz, length=1.5, pad=1)

    def get_xticklabels(self):
        RY_xticklabels = [f'{name_i}-{name_j}' for name_i, name_j in self.d_bondtype_lst_RY[self.host]] + ['others']
        YR_xticklabels = [f'{name_i}-{name_j}' for name_i, name_j in self.d_bondtype_lst_YR[self.host]] + ['others']
        return RY_xticklabels + YR_xticklabels

    def get_d_axes(self, axes):
        return {'counts': axes[0], 'boxplot': axes[1]}


class RY:
    d_bondtype_lst = {
        'atat_21mer': [('C5', 'C4'), ('C4', 'C5'), ('C2', 'C2'), ('N1', 'N3'), ('C6', 'C4')],
        'gcgc_21mer': [('N9', 'C5'), ('C4', 'C4'), ('N3', 'C2'), ('C2', 'C2'), ('N1', 'N3')]
    }

    def __init__(self, host):
        self.host = host
        self.bondtype_lst = self.d_bondtype_lst[host]

        self.n_spring_container = self.ini_n_spring_container()
        self.k_container = self.ini_k_container()

    def ini_n_spring_container(self):
        n_spring_container = {'others': 0}
        for bondtype in self.bondtype_lst:
            n_spring_container[bondtype] = 0
        return n_spring_container

    def set_n_spring_container(self):
        pass

    def ini_k_container(self):
        k_container = {'others': list()}
        for bondtype in self.bondtype_lst:
            k_container[bondtype] = list()
        return k_container


class YR:
    d_bondtype_lst = {
        'atat_21mer': [('N1', 'N7'), ('C2', 'C5'), ('C2', 'N7'), ('N3', 'C6'), ('C4', 'N6')],
        'gcgc_21mer': [('N1', 'N7'), ('C2', 'C5'), ('N3', 'C6'), ('C4', 'O6'), ('O2', 'C4')]
    }

    def __init__(self, host):
        self.host = host
        self.bondtype_lst = self.d_bondtype_lst[host]



class Spring:
    def __init__(self, resid_i, resid_j, resname_i, resname_j, atomname_i, atomname_j, k):
        self.resid_i = resid_i
        self.resid_j = resid_j # resid_i > resid_j
        self.resname_i = resname_i
        self.resname_j = resname_j
        self.atomname_i = atomname_i
        self.atomname_j = atomname_j
        self.k = k

    def get_RY_or_YR(self):
        if self.resid_i in ['A', 'G']:
            return 'RY'
        else:
            return 'YR'