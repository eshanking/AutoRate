import pandas as pd

class Plate():
    def __init__(self,data_path,moat=False):

        self.moat = moat
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

        data_keys = [k for k in self.data.keys() if k not in bg_keys]
        data_keys = data_keys[2:]

        return data_keys
