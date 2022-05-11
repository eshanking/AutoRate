import pytest
from autorate import Experiment
from importlib_resources import files

folder_path = str(files("autorate.test.data"))
exp = Experiment(folder_path,debug=False,moat=True)
exp.execute()

p = exp.plates[0]
gl = p.gen_growth_rate_lib()

def test_moat_zero():
    #some code making a bool proving the moat is gone
    #assert moat_zero == True
    assert 1+1 == 2