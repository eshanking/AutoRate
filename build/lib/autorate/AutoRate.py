import pandas as pd
import os
import scipy.optimize as sciopt
import scipy.interpolate as sciinter
import matplotlib.pyplot as plt
import numpy as np

class Experiment():
    """Experiment class for a given plate reader experiment
    """
    def __init__(self,
                folder_path,
                moat=False,
                replicate_arrangement='rows',
                drug_conc = None,
                units = 'ug/mL',
                debug=False):
        """Initializer

        Args:
            folder_path (str): path of plate reader data
            moat (bool, optional): If true, assumes the outer row of the plate is a moat. Defaults to False.
            replicate_arrangement (str, optional): Determines if replicates are arranged in rows or columns. Defaults to rows.
            drug_conc (list, optional): Drug concentrations corresponding to the drug diluation scheme. If none, defaults to 
            [0,0.003,0.0179,0.1072,0.643,3.858,23.1481,138.8889,833.3333,5000].
            units (str, optional): Drug concentration units. Defaults to ug/mL
        """
        self.moat = moat
        self.folder_path = folder_path
        self.replicate_arrangement = replicate_arrangement
        self.plate_data_paths = self.load_plate_data()
        self.plates = []
        self.units = units
        self.debug = debug

        if drug_conc is None:
            self.drug_conc = [0,0.003,0.0179,0.1072,0.643,3.858,23.1481,138.8889,833.3333,5000]
        else:
            self.drug_conc = drug_conc

        for pdp in self.plate_data_paths:
            self.plates.append(Plate(pdp,self.drug_conc,debug=debug,moat=moat))
        
        self.growth_rate_lib = self.gen_growth_rate_lib()
        self.seascape_lib = self.gen_seascape_lib()
        


    def load_plate_data(self):
        """Gets plate data paths

        Returns:
            list: list of plate data paths
        """
        plate_files = os.listdir(path=self.folder_path)

        #Need to make sure we are only attempting to load .csv data
        plate_files = [i for i in plate_files if '.csv' in i]

        plate_files.sort()

        plate_data_paths = []

        for pf in plate_files:
            if pf != '.DS_Store':
                plate_path = self.folder_path + os.sep + pf
                plate_data_paths.append(plate_path)

        return plate_data_paths
    
    def gen_growth_rate_lib(self):
        """Generates the complete growth rate library by compining the growth rate libraries of each plate

        Returns:
            dict: growth rate library
        """
        gl = {}
        offset = 0
        for p in self.plates:
            
            for k in p.growth_rate_lib.keys():

                rep_num = int(k)
                gl[str(rep_num+offset)] = p.growth_rate_lib[k]
                #print(str(rep_num+offset))
            keys = [int(k) for k in p.growth_rate_lib.keys()]
            offset += max(keys)

        return gl
    
    def gen_seascape_lib(self):
        """Fits raw estimated growth rate values to a Hill dose-response curve

        Args:
            pop (population class object, optional): population class object. Defaults to None.
            debug (bool, optional): generates plots useful for debugging if true. Defaults to False.

        Raises:
            ValueError: raises error if there is no growth rate library in the experiment object

        Returns:
            dict: seascape library
        """


        if not 'growth_rate_lib' in self.__dict__:
            raise ValueError('No growth rate library in population.')
        else:
            gl = self.growth_rate_lib
        
        sl = {}
        dc = self.drug_conc

        replicates = [int(k) for k in self.growth_rate_lib.keys()]
        
        for r in replicates:

            popt = self.fit_hill_curve(dc,gl[str(r)])
            ic50 = popt[0]
            g_drugless = popt[1]
            hill_coeff = popt[2]
            d_t = {'ic50':ic50,
                'g_drugless':g_drugless,
                'hill_coeff':hill_coeff}

            sl[str(r)] = d_t

        return sl

    def fit_hill_curve(self,xdata,ydata):
        """Fits dose-response curve to growth rate data

        Args:
            xdata (list or numpy array): drug concentration curve from plate experiment
            ydata (list or numpy array): growth rate versus drug concetration for a given replicate

        Returns:
            list: List of optimized paramters: IC50, drugless growth rate, and Hill coefficient
        """
        # interpolate data
        xd_t = xdata
        yd_t = ydata
        f = sciinter.interp1d(xdata,ydata)

        if min(xdata) == 0:
            xmin = np.log10(xdata[1])
        else:
            xmin = np.log10(min(xdata))
        xmax = np.log10(max(xdata))

        xdata = np.logspace(xmin,xmax)
        if not xdata[0] == 0:
            xdata = np.insert(xdata,0,0)

        ydata = f(xdata)

        
        p0 = [0,ydata[0],-0.08]

        if ydata[0] == 0:
            g_drugless_bound = [0,1]
        else:
            # want the estimated drugless growth rate to be very close to the value given in ydata
            g_drugless_bound = [ydata[0]-0.0001*ydata[0],ydata[0]+0.0001*ydata[0]]

        bounds = ([-5,g_drugless_bound[0],-1],[4,g_drugless_bound[1],-0.001])
        popt, pcov = sciopt.curve_fit(self.logistic_pharm_curve_vectorized,
                                            xdata,ydata,p0=p0,bounds=bounds)

        if self.debug:
            est = [self.logistic_pharm_curve(x,popt[0],popt[1],popt[2]) for x in xdata]
            fig,ax = plt.subplots()
            ax.plot(xd_t,yd_t)
            ax.plot(xdata,ydata)
            ax.plot(xdata,est)
            ax.set_xscale('log')
            ax.set_title('IC50 = ' + str(popt[0]))

        return popt
    
    def logistic_pharm_curve_vectorized(self,x,IC50,g_drugless,hill_coeff):
        """Defines the logistic dose-response curve. Use if the input is a vector of drug concentration curves

        Args:
            x (numpy array): drug concentration vector
            IC50 (float)): IC50
            g_drugless (float): drugless growth rate
            hill_coeff (float): Hill coefficient

        Returns:
            numpy array: array of growth rates
        """
        g = []

        for x_t in x:
            if x_t == 0:
                g.append(g_drugless)
            else:
                g.append(g_drugless/(1+np.exp((IC50-np.log10(x_t))/hill_coeff)))

        return g

    def save_results(self):
        seascape_df = pd.DataFrame.from_dict(self.seascape_lib)
        growth_rate_df = pd.DataFrame.from_dict(self.growth_rate_lib)

        seascape_df.to_csv(path_or_buf='sl.csv')
        growth_rate_df.to_csv(path_or_buf='gr.csv')


class Plate():
    """96-well plate object
    """
    def __init__(self,data_path,drug_conc,replicate_arrangement='rows',moat=False,debug=False):
        """Initializer

        Args:
            data_path (str): csv file path
            drug_conc (list of floats): drug concentration gradient
            replicates (str, optional): Determines if replicates are arranged in rows or columns. Defaults to rows.
            moat (bool, optional): If true, assumes the outer row of the plate is a moat. Defaults to False.
            debug (bool, optional): If true, plots growth curve and estimated growth curve. Defaults to False.
        """
        self.moat = moat
        self.data = pd.read_csv(data_path)

        self.background_keys = self.get_background_keys()
        self.data_keys = self.get_data_keys()
        self.replicate_arrangement = replicate_arrangement
        self.drug_conc = drug_conc
        self.debug = debug
        self.growth_rate_lib = self.gen_growth_rate_lib()

    def get_background_keys(self):
        """Gets the dataframe keys for the background (aka moat)

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
        """Gets the dataframe keys for the data

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
    
    def est_growth_rate(self,growth_curve,t=None):
        """Estimates growth rate from OD growth curve

        Args:
            growth_curve (list or numpy array): vector of OD data
            t (list or numpy array, optional): Time vector. If None, algorithm assumes each time step is 1 s. Defaults to None.

        Returns:
            float: Growth rate in units of 1/s
        """
        
        if t is None:
            t = np.arange(len(growth_curve))
        
        p0 = [10**-6,0.05,1] # starting parameters

        popt, pcov = sciopt.curve_fit(self.logistic_growth_curve,
                                            t,growth_curve,p0=p0,
                                            bounds=(0,1))

        r = popt[0]
        if r < 0:
            r = 0
        if popt[2] < popt[1]: # if the carrying capacity is less than the initial population size
            r = 0
        if popt[2] < 0.1:
            r = 0
        
        if self.debug:
            fig,ax = plt.subplots()

            ax.plot(t,growth_curve)

            est = self.logistic_growth_curve(t,popt[0],popt[1],popt[2])
            
            ax.plot(t,est)
            # print(popt[0])
            p0 = round(popt[1]*10**5)/10**5
            k = round(popt[2]*10**5)/10**5
            r = round(r*10**5)/10**5
            title = 'rate = ' + str(r*(60**2)) + ' cc = ' + str(k)
            ax.set_title(title)        

        return r

    def get_growth_rates_from_df(self):
        
        """Estimates the growth rates from timeseries growth data in a dataframe

        Returns:
            growth_rates: dict
                dictionary of growth rates for each experimental condition
        """

        growth_rates = {}
        df = self.data

        data_keys = self.get_data_keys()
        time = df['Time [s]']

        for k in data_keys:
            growth_rates[k] = self.est_growth_rate(df[k],t=time)

        return growth_rates

    def gen_growth_rate_lib(self):
        """Generates growth rate library from OD data

        Returns:
            dict: Dict of dose-response curves indexed by replicate
        """
        growth_rates = self.get_growth_rates_from_df()
        replicate_num = 0
        growth_rate_lib = {}

        if self.replicate_arrangement == 'rows': # letters represent individual replicates
            # get all the letters in the data keys
            replicates = []
            concentrations = []
            for key in self.data_keys:
                replicates.append(key[0])
                concentrations.append(key[1:])

            replicates = list(set(replicates))
            concentrations = list(set(concentrations))


            for r in replicates:
                gr_vect = []
                i = 0 # concentration index
                for c in self.drug_conc:
                    key = r + concentrations[i]
                    gr_vect.append(growth_rates[key])
                    i += 1

                growth_rate_lib[str(replicate_num)] = gr_vect
                replicate_num += 1

        else:
            replicates = []
            concentrations = []
            for key in self.data_keys:
                replicates.append(key[1:])
                concentrations.append(key[0])

            for r in replicates:
                gr_vect = []
                i = 0 # concentration index
                for c in self.drug_conc:
                    key = concentrations[i] + r
                    gr_vect.append(growth_rates[key])
                    i += 1

                growth_rate_lib[str(replicate_num)] = gr_vect
                replicate_num += 1

        return growth_rate_lib

    def logistic_growth_curve(self,t,r,p0,k):
        """Logistic growth equation

        Args:
            t (float): time
            r (float): growth rate
            p0 (float): starting population size
            k (float): carrying capacity

        Returns:
            float: population size at time t
        """
        p = k/(1+((k-p0)/p0)*np.exp(-r*t))

        return p