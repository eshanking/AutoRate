import pandas as pd
import os

class Experiment():
    def __init__(self,folder_path,moat=False):
        self.moat = moat
        self.folder_path = folder_path
        self.plate_data_paths = self.load_plate_data()
        self.plates = []

        for pdp in self.plate_data_paths:
            self.plates.append(Plate(pdp))

    def load_plate_data(self):
        plate_files = sim_files = os.listdir(path=self.folder_path)
        plate_files.sort()

        plate_data_paths = []

        for pf in plate_files:
            if pf != '.DS_Store':
                plate_path = self.folder_path + os.sep + pf
                plate_data_paths.append(plate_path)

        return plate_data_paths


class Plate():
    def __init__(self,data_path,moat=False):

        self.moat = moat
        # print(data_path)
        self.data = pd.read_csv(data_path)

        self.background_keys = self.get_background_keys()
        self.data_keys = self.get_data_keys()

    def get_background_keys(self):
        """Gets the dataframe keys for the background assuming a 1-well moat

        Returns:
            list: list of background keys
        """
        # row A, row H, col 1, and col 12

        if self.moat:
            k = self.data.keys()

            k = k[2:]
            bg_keys = [y for y in k if int(y[1:]) == 1] # col 1
            bg_keys = bg_keys + [y for y in k if (int(y[1:]) == 12 and y not in bg_keys)]
            bg_keys = bg_keys + [y for y in k if (y[0] == 'A' and y not in bg_keys)]
            bg_keys = bg_keys + [y for y in k if (y[0] == 'H' and y not in bg_keys)]
        
        else:
            bg_keys = None

        return bg_keys

    def get_data_keys(self):
        """Gets the dataframe keys for the data assuming a 1-well moat

        Args:
            df (pandas dataframe): datafram containing raw OD data

        Returns:
            list: list of keys
        """
        bg_keys = self.get_background_keys()
        if self.background_keys is None:
            data_keys = self.data.keys()
        else:
            data_keys = [k for k in self.data.keys() if k not in self.background_keys]
        
        data_keys = data_keys[2:]

        return data_keys
