from autorate import Experiment
from importlib_resources import files

source = files("autorate.test.data")
exp = Experiment(folder_path,debug=False,moat=True)

p = exp.plates[0]
gl = p.gen_growth_rate_lib()