from os import path
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from prominent_modes.prominent import AAProminentModes, TTProminentModes, GGProminentModes, CCProminentModes, TATAProminentModes, CpGProminentModes

class TwoPlot:
    group_lst = ['boxplot', 'scatter']
    data_folder = '/home/yizaochen/codes/dna2021paper/notebooks/SI/prominet_modes_st/data'
    d_prominent = {'a_tract_21mer': {'STRAND1': AAProminentModes, 'STRAND2': TTProminentModes},
                   'g_tract_21mer': {'STRAND1': GGProminentModes, 'STRAND2': CCProminentModes}}

    # Fontsize
    lbfz = 6 # xy labels
    tickfz = 6

    def __init__(self, host, strand_id):
        self.host = host
        self.strand_id = strand_id
        self.f_df = path.join(self.data_folder, f'{self.host}_{self.strand_id}_st.csv')
        self.df = pd.read_csv(self.f_df)

        self.prominent_agent = self.d_prominent[self.host][self.strand_id]()

    def plot_main(self, figsize, xlim=None, xticks=None, ylim=None, yticks=None, label_txt=False):
        fig = plt.figure(figsize=figsize, facecolor='white')
        d_axes = self.get_d_axes(fig)
        x_array, y_array = self.get_x_y_array()
        x_criteria = self.get_x_outlier_criteria(x_array)
        self.do_boxplot(d_axes['boxplot'], x_array, x_criteria)
        self.do_scatter_plot(d_axes['scatter'], x_array, y_array, x_criteria)
        self.remove_ticks(d_axes)
        if xlim is not None:
            self.set_xlim_xticks(d_axes, xlim, xticks)
        if ylim is not None:
            self.set_ylim_yticks(d_axes, ylim, yticks)
        if label_txt:
            self.label_text(d_axes['scatter'])
        return fig, d_axes

    def get_x_outlier_criteria(self, x_array):
        whis = 1
        q1 = np.quantile(x_array, .25)
        q3 = np.quantile(x_array, .75)
        return q3 + whis * (q3 - q1)

    def get_x_y_array(self):
        x_array = self.df['lambda'] # λ
        y_array = self.df['<r>'] # <r>
        return x_array, y_array

    def set_xlim_xticks(self, d_axes, xlim, xticks):
        d_axes['scatter'].set_xticks(xticks)
        d_axes['boxplot'].set_xlim(xlim)
        d_axes['scatter'].set_xlim(xlim)

    def set_ylim_yticks(self, d_axes, ylim, yticks):
        d_axes['scatter'].set_yticks(yticks)
        d_axes['scatter'].set_ylim(ylim)

    def do_scatter_plot(self, ax, x_array, y_array, x_criteria):
        x_array_others, y_array_others = self.get_other_xy_array(x_array, y_array, x_criteria)
        #ax.scatter(x_array_prominent, y_array_prominent, s=2, color='tab:orange')
        self.scatter_prominent(ax)
        ax.scatter(x_array_others, y_array_others, s=2, color='dimgray')
        xlabel = r'$\lambda_{\alpha' + "'}" + r'$'
        ylabel = r'$\left< r_{\alpha' + "'}" + r'\right>$' 
        ax.set_xlabel(xlabel, fontsize=self.lbfz)
        ax.set_ylabel(ylabel, fontsize=self.lbfz)
        ax.axvline(x_criteria, linestyle='--', linewidth=1, color='red', alpha=0.8)
        ax.tick_params(axis='both', labelsize=self.tickfz, length=1.5, pad=1)

    def scatter_prominent(self, ax):
        mode_container = self.prominent_agent.get_mode_container()
        for mode in mode_container:
            ax.scatter(mode.x, mode.y, s=2, color=mode.color)

    def label_text(self, ax):
        mode_container = self.prominent_agent.get_mode_container()
        for mode in mode_container:
            if mode.label_ornot:
                txt = f'{mode.mode_id}'
                ax.text(mode.x, mode.y, txt, fontsize=4)

    def get_prominent_xy_array(self, x_array, y_array, x_criteria):
        mask_prominent = (x_array >= x_criteria)
        x_array_prominent = x_array[mask_prominent]
        y_array_prominent = y_array[mask_prominent]
        return x_array_prominent, y_array_prominent

    def get_other_xy_array(self, x_array, y_array, x_criteria):
        mask_prominent = (x_array >= x_criteria)
        x_array_others = x_array[~mask_prominent]
        y_array_others = y_array[~mask_prominent]
        return x_array_others, y_array_others

    def do_boxplot(self, ax, x_array, x_criteria):
        ax.boxplot([x_array], whis=1., widths=0.6, vert=False, boxprops={'facecolor': 'steelblue'}, medianprops={'color': 'black'}, patch_artist=True, sym='d', flierprops={'markerfacecolor': 'tab:blue', 'markersize': 3, 'markeredgewidth': 0.8})
        ax.axvline(x_criteria, linestyle='--', linewidth=1, color='red', alpha=0.8)

    def get_d_axes(self, fig):
        d_axes = {group_name: dict() for group_name in self.group_lst}
        outer_grid = gridspec.GridSpec(265, 1, wspace=0, hspace=0) 
        d_axes['boxplot'] = fig.add_subplot(outer_grid[:38])
        d_axes['scatter'] = fig.add_subplot(outer_grid[45:])
        return d_axes

    def remove_ticks(self, d_axes):
        d_axes['boxplot'].tick_params(axis='both', bottom=False, left=False, labelbottom=False, labelleft=False)


class TwoPlotHetero(TwoPlot):
    d_prominent = {'atat_21mer': TATAProminentModes,
                   'gcgc_21mer': CpGProminentModes}

    def __init__(self, host):
        self.host = host

        self.f_df = path.join(self.data_folder, f'{self.host}_st.csv')
        self.df = pd.read_csv(self.f_df)

        self.prominent_agent = self.d_prominent[self.host]()

    def plot_main(self, figsize, xlim=None, xticks=None, ylim=None, yticks=None, label_txt=False):
        fig = plt.figure(figsize=figsize, facecolor='white')
        d_axes = self.get_d_axes(fig)
        x_array, y_array = self.get_x_y_array()
        x_criteria = self.get_x_outlier_criteria(x_array)
        self.do_boxplot(d_axes['boxplot'], x_array, x_criteria)
        self.do_scatter_plot(d_axes['scatter'], x_array, y_array, x_criteria)
        self.remove_ticks(d_axes)
        if xlim is not None:
            self.set_xlim_xticks(d_axes, xlim, xticks)
        if ylim is not None:
            self.set_ylim_yticks(d_axes, ylim, yticks)
        if label_txt:
            self.label_text(d_axes['scatter'])
        return fig, d_axes

    def do_scatter_plot(self, ax, x_array, y_array, x_criteria):
        #x_array_prominent, y_array_prominent = self.get_prominent_xy_array(x_array, y_array, x_criteria)
        x_array_others, y_array_others = self.get_other_xy_array(x_array, y_array, x_criteria)
        self.scatter_prominent(ax)
        #ax.scatter(x_array_prominent, y_array_prominent, s=2, color='tab:orange')
        ax.scatter(x_array_others, y_array_others, s=2, color='dimgray')
        xlabel = r'$\lambda_{\alpha' + "'}" + r'$'
        ylabel = r'$\left< r_{\alpha' + "'}" + r'\right>$' 
        ax.set_xlabel(xlabel, fontsize=self.lbfz)
        ax.set_ylabel(ylabel, fontsize=self.lbfz)
        ax.axvline(x_criteria, linestyle='--', linewidth=1, color='red', alpha=0.8)
        ax.tick_params(axis='both', labelsize=self.tickfz, length=1.5, pad=1)

    def get_mode_array(self):
        return self.df['mode-id']

    def get_x_y_array(self):
        x_array = self.df['lambda']
        y_array = self.df['<r>']
        return x_array, y_array

class TwoPlotHeteroV1(TwoPlot):
    group_lst = ['n(RY)>n(YR)', 'n(RY)<n(YR)', 'n(RY)=n(YR)', 'all k<0.5']
    d_color = {'n(RY)>n(YR)': 'blue', 'n(RY)<n(YR)': 'red', 'n(RY)=n(YR)': 'magenta', 'all k<0.5': 'dimgray'}

    def __init__(self, host):
        self.host = host

        self.f_df = path.join(self.data_folder, f'{self.host}_st.csv')
        self.df = pd.read_csv(self.f_df)

        self.df_prob = self.get_df_prob()
        self.df_concat = self.get_df_concat()

        #self.prominent_agent = self.d_prominent[self.host][self.strand_id]()

    def plot_main(self, figsize, xlim=None, xticks=None, ylim=None, yticks=None, label_txt=False):
        fig = plt.figure(figsize=figsize, facecolor='white')
        d_axes = self.get_d_axes(fig)
        x_array, y_array = self.get_x_y_array()
        x_criteria = self.get_x_outlier_criteria(x_array)
        self.do_boxplot(d_axes['boxplot'], x_array, x_criteria)
        self.do_scatter_plot_v1(d_axes['scatter'])
        self.remove_ticks(d_axes)
        if xlim is not None:
            self.set_xlim_xticks(d_axes, xlim, xticks)
        if ylim is not None:
            self.set_ylim_yticks(d_axes, ylim, yticks)
        if label_txt:
            self.label_text(d_axes['scatter'])
        return fig, d_axes

    def do_scatter_plot_v1(self, ax):
        d_xy_array = self.get_d_x_y_array_by_group()
        for group_id in self.group_lst:
            x_array, y_array = d_xy_array[group_id]
            color = self.d_color[group_id]
            ax.scatter(x_array, y_array, s=2, color=color)
        xlabel = r'$\lambda_{\alpha' + "'}" + r'$'
        ylabel = r'$\left< r_{\alpha' + "'}" + r'\right>$' 
        ax.set_xlabel(xlabel, fontsize=self.lbfz)
        ax.set_ylabel(ylabel, fontsize=self.lbfz)
        #ax.axvline(x_criteria, linestyle='--', linewidth=1, color='red', alpha=0.8)
        ax.tick_params(axis='both', labelsize=self.tickfz, length=1.5, pad=1)

    def get_d_x_y_array_by_group(self):
        d_xy_array = {group_id: None for group_id in self.group_lst}

        small_k_mask = self.get_small_k_mask()
        df_small_k = self.df_concat[small_k_mask]
        d_xy_array['all k<0.5'] = (df_small_k['lambda'], df_small_k['<r>'])

        df_no_small_k = self.df_concat[~small_k_mask]
        df_RY = df_no_small_k[df_no_small_k['Count-RY'] > df_no_small_k['Count-YR']]
        df_YR = df_no_small_k[df_no_small_k['Count-RY'] < df_no_small_k['Count-YR']]
        d_xy_array['n(RY)>n(YR)'] = (df_RY['lambda'], df_RY['<r>'])
        d_xy_array['n(RY)<n(YR)'] = (df_YR['lambda'], df_YR['<r>'])

        df_equal = df_no_small_k[df_no_small_k['Count-RY'] == df_no_small_k['Count-YR']]
        d_xy_array['n(RY)=n(YR)'] = (df_equal['lambda'], df_equal['<r>'])
        return d_xy_array

    def get_d_df_by_group(self):
        d_df = {group_id: None for group_id in self.group_lst}

        small_k_mask = self.get_small_k_mask()
        df_small_k = self.df_concat[small_k_mask]
        d_df['all k<0.5'] = df_small_k

        df_no_small_k = self.df_concat[~small_k_mask]
        df_RY = df_no_small_k[df_no_small_k['Count-RY'] > df_no_small_k['Count-YR']]
        df_YR = df_no_small_k[df_no_small_k['Count-RY'] < df_no_small_k['Count-YR']]
        d_df['n(RY)>n(YR)'] = df_RY
        d_df['n(RY)<n(YR)'] = df_YR

        df_equal = df_no_small_k[df_no_small_k['Count-RY'] == df_no_small_k['Count-YR']]
        d_df['n(RY)=n(YR)'] = df_equal
        return d_df

    def get_small_k_mask(self):
        return (self.df_concat['Count-RY'] == 0) & (self.df_concat['Count-YR'] == 0)

    def do_scatter_plot(self, ax, x_array, y_array, x_criteria):
        x_array_prominent, y_array_prominent = self.get_prominent_xy_array(x_array, y_array, x_criteria)
        x_array_others, y_array_others = self.get_other_xy_array(x_array, y_array, x_criteria)
        #self.scatter_prominent(ax)
        ax.scatter(x_array_prominent, y_array_prominent, s=2, color='tab:orange')
        ax.scatter(x_array_others, y_array_others, s=2, color='dimgray')
        xlabel = r'$\lambda_{\alpha' + "'}" + r'$'
        ylabel = r'$\left< r_{\alpha' + "'}" + r'\right>$' 
        ax.set_xlabel(xlabel, fontsize=self.lbfz)
        ax.set_ylabel(ylabel, fontsize=self.lbfz)
        ax.axvline(x_criteria, linestyle='--', linewidth=1, color='red', alpha=0.8)
        ax.tick_params(axis='both', labelsize=self.tickfz, length=1.5, pad=1)

    def get_mode_array(self):
        return self.df['mode-id']

    def get_x_y_array(self):
        x_array = self.df['lambda']
        y_array = self.df['<r>']
        return x_array, y_array

    def get_df_prob(self):
        datafolder = '/home/yizaochen/codes/dna2021paper/notebooks/SI/prominet_modes_st/chord_diagram_st/heterogeneous/data'
        f_df = path.join(datafolder, f'{self.host}.prob.RY_YR.csv')
        return pd.read_csv(f_df)

    def get_df_concat(self):
        return pd.concat([self.df, self.df_prob], axis=1)