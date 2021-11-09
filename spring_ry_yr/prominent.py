from os import path
import pandas as pd

class AAProminentModes:
    host = 'a_tract_21mer'
    strand_id = 'STRAND1'

    data_folder = '/home/yizaochen/codes/dna2021paper/notebooks/SI/prominet_modes_st/data'

    d_colors = {1: 'tab:orange', 2: 'tab:green'}
    mode_lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13]
    d_category = {1: [1, 2, 3, 4, 5, 6, 9], 2: [7, 8, 11, 13]}
    labelornot_lst = [1, 7]

    def __init__(self):
        self.f_df = path.join(self.data_folder, f'{self.host}_{self.strand_id}_st.csv')
        self.df = pd.read_csv(self.f_df)

    def get_mode_container(self):
        x_array, y_array = self.get_x_y_array()
        mode_container = list()
        for idx, mode_id in enumerate(self.mode_lst):
            group_id = self.determine_group_id(mode_id)
            color = self.d_colors[group_id]
            label_ornot = self.determine_labelornot(mode_id)
            x = x_array[idx]
            y = y_array[idx]
            mode_container.append(ProminentMode(x, y, mode_id, color, label_ornot))
        return mode_container

    def get_x_y_array(self):
        x_array = self.df['lambda'] # λ
        y_array = self.df['<r>'] # <r>
        return x_array, y_array

    def determine_group_id(self, mode_id):
        if mode_id in self.d_category[1]:
            return 1
        else:
            return 2

    def determine_labelornot(self, mode_id):
        if mode_id in self.labelornot_lst:
            return True
        else:
            return False

class TTProminentModes(AAProminentModes):
    host = 'a_tract_21mer'
    strand_id = 'STRAND2'

    d_colors = {1: 'magenta'}
    mode_lst = [10, 12, 14, 17, 18, 21, 22]
    d_category = {1: [10, 12, 14, 17, 18, 21, 22]}
    labelornot_lst = [10]

class GGProminentModes(AAProminentModes):
    host = 'g_tract_21mer'
    strand_id = 'STRAND1'

    d_colors = {1: 'tab:orange', 2: 'tab:green', 3: 'tab:cyan'}
    mode_lst = [1, 2, 3, 4, 5, 7, 8, 10, 11, 13, 14]
    d_category = {1: [1, 2, 4, 14], 2: [3, 5, 7, 11, 13], 3: [8, 10, 14]}
    labelornot_lst = [1, 3, 8]

    def determine_group_id(self, mode_id):
        if mode_id in self.d_category[1]:
            return 1
        elif mode_id in self.d_category[2]:
            return 2
        else:
            return 3

class CCProminentModes(AAProminentModes):
    host = 'g_tract_21mer'
    strand_id = 'STRAND2'

    d_colors = {1: 'magenta'}
    mode_lst = [6, 9, 12, 16, 18, 21, 24, 27]
    d_category = {1: [6, 9, 12, 16, 18, 21, 24, 27]}
    labelornot_lst = [6]

class TATAProminentModes(AAProminentModes):
    host = 'atat_21mer'
    strand_id = None

    d_colors = {1: 'tab:orange', 2: 'tab:green'}
    mode_lst = list(range(1, 21))
    d_category = {1: list(range(1, 15)), 2: list(range(15, 21))}
    labelornot_lst = [1, 15]

    def __init__(self):
        self.f_df = path.join(self.data_folder, f'{self.host}_st.csv')
        self.df = pd.read_csv(self.f_df)

class CpGProminentModes(TATAProminentModes):
    host = 'gcgc_21mer'
    strand_id = None

    d_colors = {1: 'tab:orange'}
    mode_lst = list(range(1, 41))
    d_category = {1: list(range(1, 41))}
    labelornot_lst = [1]

class ProminentMode:
    def __init__(self, x, y, mode_id, color, label_ornot):
        self.x = x
        self.y = y
        self.mode_id = mode_id
        self.color = color
        self.label_ornot = label_ornot # True or False