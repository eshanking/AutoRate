import AutoRate

folder_path = '/Users/kinge2/repos/AutoRate/test_data'

exp = AutoRate.Experiment(folder_path,debug=True)

p = exp.plates[0]
gl = p.gen_growth_rate_lib()