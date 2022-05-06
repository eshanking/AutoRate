import autorate.AutoRate as AutoRate
import pytest

exp = AutoRate.Experiment(folder_path,debug=False,moat=True)

p = exp.plates[0]
gl = p.gen_growth_rate_lib()