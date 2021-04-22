from utilities import read_gravity_field
from utilities import degree_varience_sq
from utilities import degree_varience_disturbing_pot
from utilities import degree_variance_un
from utilities import degree_varience_GD
from utilities import degree_varience_GG
from utilities import read_LLN_internalformat, read_LLN_from_LOADEV
from utilities import sigma_EWH

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Read love numbers (k)
path_to_LLN = "/home/fine/PROJECTS/GFZ/GRACE_postprocessing/LoadLove_PG_CF_oct.dat"
path_to_LOADDEV = "/home/fine/PROJECTS/GFZ/LoadDef/output/Love_Numbers/LLN/lln_.txt"



f_static = "/home/fine/PROJECTS/GFZ/GRACE_postprocessing/static_models/GGM05C.gfc"
path_to_GSM = "/home/fine/PROJECTS/GFZ/GRACE_postprocessing/GRACE/GFZ/RL06_60"
list_of_files = os.listdir(path_to_GSM)
LLN_internal = read_LLN_internalformat(path_to_LLN, 60)
LLN_loaddev = read_LLN_from_LOADEV(path_to_LOADDEV, 60)

# for fn in list_of_files[0:1]:
# 	fp = os.path.join(path_to_GSM, fn)
# 	C,S, gm, r = read_gravity_field(fp, f_static)
# 	sigma_sq = degree_varience_sq(C, S, 60)
# 	sigma = np.sqrt(sigma_sq)
# 	disturbingPot = degree_varience_disturbing_pot(gm,r, sigma)
# 	undulations = degree_variance_un(r,sigma)
# 	degreeDist = degree_varience_GD(gm,r,sigma, 60)
# 	degreeGrad = degree_varience_GG(gm, r, sigma, 60)
# 	EWH = sigma_EWH(5500, 1000, r, LLN_loaddev, 60, sigma)
#
# 	df = pd.DataFrame()
# 	df['Disturbing Pot'] = disturbingPot
# 	df['Undulations'] = undulations
# 	df['degree Disturbances'] = degreeDist
# 	df['Degree Gradient'] = degreeGrad
# 	df['EWH'] = EWH
# 	# ax = df['EWH'][3::].plot()
# 	# ax.set_yscale('log')
# 	# plt.show()

df = pd.DataFrame([[22,33,44]])
df2 = pd.DataFrame([[2,3,4]])
print([[22,33,44]])

df3 = pd.concat([df2, df])

ind = [0,1]
df3.index = ind
print(df3)
print('******************')
print(df3[2])
df = pd.DataFrame()
DATE = []
firstread=True
df_list = []

for fn in list_of_files:
	fp = os.path.join(path_to_GSM, fn)
	C,S, gm, r, date_start, date_end = read_gravity_field(fp, f_static)

	DATE.append(date_start.strftime(format = '%Y-%m-%d %H:%M:%S')[0:10])


	sigma_sq = degree_varience_sq(C, S, 60)
	sigma = np.sqrt(sigma_sq)
	EWH = [list(sigma_EWH(5500, 1000, r, LLN_loaddev, 60, sigma).flatten())]

	df = pd.DataFrame(EWH)

	df_list.append(df)
DF = pd.concat(df_list)
DF.index = pd.to_datetime(DATE)
DF_sorted = DF.sort_index()
DF_sorted_transp = DF_sorted.T
DF_sorted.to_csv('EWH.csv')
DF_sorted_transp.to_csv('EWH_transposed.csv')





