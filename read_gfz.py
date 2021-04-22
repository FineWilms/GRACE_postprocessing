import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pyshtools as pysh
from pyshtools import constants

from pyshtools import shtools

def read_gravity_field(fp):
	# path_to_GSM = os.path.abspath(os.path.join(os.path.dirname('__file__'), 'GRACE', 'GFZ', 'RL06_60'))
	# list_of_GSM_files = os.listdir(path_to_GSM)
	# for fn in list_of_GSM_files[10:11]:
	fn = os.path.basename(fp)
	string_date = fn.strip().split('_')[1]
	start = string_date.split('-')[0]
	end = string_date.split('-')[1]
	year = np.int(start[0:4])
	day = np.int(start[4::])
	date_start = datetime.datetime.strptime('{} {}'.format(day, year), '%j %Y')
	year = np.int(end[0:4])
	day = np.int(end[4::])
	date_end = datetime.datetime.strptime('{} {}'.format(day, year), '%j %Y')
	print(date_start, date_end)
	string_eoY = '# End of YAML header'
	string_gm = 'earth_gravity_param'
	string_r = 'mean_equator_radius'
	linecount = 0
	GM_r = []
	with open(fp, 'r') as f:
		for i, line in enumerate(f):
			if string_gm in line:
				GM_r.append(i + 3)
			if string_r in line:
				GM_r.append(i + 3)

	with open(fp, 'r') as f:
		gm_value = np.float(f.readlines()[GM_r[0]].strip().split(':')[1])
	with open(fp, 'r') as f:
		r_value = np.float(f.readlines()[GM_r[1]].strip().split(':')[1])

	with open(fp, 'r') as f:
		for line in f:
			if string_eoY in line:
				header2skip = linecount
			linecount += 1
	df = pd.read_csv(fp, skiprows=header2skip, delim_whitespace=True, usecols=[1, 2, 3, 4], index_col=False,
					 names=['degree_l', 'order_m', 'C_lm', 'S_lm'], header=0)
	L = df['degree_l'].values
	M = df['order_m'].values
	LM = list(zip(L, M))

	arr_C = np.empty((61, 61))
	arr_C[:] = np.NaN

	arr_S = np.empty((61, 61))
	arr_S[:] = np.NaN
	for v in LM:
		C = df.loc[(df['degree_l'] == v[0]) & (df['order_m'] == v[1])]['C_lm'].values[0]
		arr_C[v[0], v[1]] = C

		S = df.loc[(df['degree_l'] == v[0]) & (df['order_m'] == v[1])]['S_lm'].values[0]
		arr_S[v[0], v[1]] = S
	clm = np.empty((2, 61, 61))
	clm[:] = np.NaN

	clm[0, :, :] = arr_C
	clm[1, :, :] = arr_S
	# Get the change only (anomaly not total gravity field)



	fp = os.path.abspath(os.path.join(os.path.dirname('__file__'), 'static_models', 'GGM05C.gfc'))
	rows2skip = list(range(45)) + [46]
	df_static = pd.read_csv(fp, header=0, skiprows=rows2skip, delim_whitespace=True)

	L = df_static['L'].values
	M = df_static['M'].values

	LM = list(zip(L, M))
	LMFiltered = [(x, y) for x, y in LM if (x <= 60)]
	LMFiltered = [(x, y) for x, y in LMFiltered if (y <= 60)]

	arr_C_static = np.empty((61, 61))
	arr_C_static[:] = np.NaN

	arr_S_static = np.empty((61, 61))
	arr_S_static[:] = np.NaN

	for v in LMFiltered:
		C = df_static.loc[(df_static['L'] == v[0]) & (df_static['M'] == v[1])]['C'].values[0].replace('D', 'E')
		arr_C_static[v[0], v[1]] = C
		S = df_static.loc[(df_static['L'] == v[0]) & (df_static['M'] == v[1])]['S'].values[0].replace('D', 'E')
		arr_S_static[v[0], v[1]] = S

	arr_C_delta = arr_C - arr_C_static

	arr_S_delta = arr_S - arr_S_static
	clm_delta = np.empty((2, 61, 61))

	clm_delta[:] = np.NaN

	clm_delta[0, :, :] = arr_C_delta
	clm_delta[1, :, :] = arr_S_delta
