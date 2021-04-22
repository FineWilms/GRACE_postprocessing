import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pyshtools as pysh
from pyshtools import constants

from pyshtools import shtools

def read_gravity_field(fp_full, fp_static):
	fn = os.path.basename(fp_full)

	string_date = fn.strip().split('_')[1]
	start = string_date.split('-')[0]
	end = string_date.split('-')[1]
	year = np.int(start[0:4])
	day = np.int(start[4::])
	date_start = datetime.datetime.strptime('{} {}'.format(day, year), '%j %Y')
	year = np.int(end[0:4])
	day = np.int(end[4::])
	date_end = datetime.datetime.strptime('{} {}'.format(day, year), '%j %Y')
	print('reading gravity coefficients for period %s to %s' %(date_start, date_end))
	string_eoY = '# End of YAML header'
	string_gm = 'earth_gravity_param'
	string_r = 'mean_equator_radius'
	linecount = 0
	GM_r = []
	with open(fp_full, 'r') as f:
		for i, line in enumerate(f):
			if string_gm in line:
				GM_r.append(i + 3)
			if string_r in line:
				GM_r.append(i + 3)

	with open(fp_full, 'r') as f:
		gm_value = np.float(f.readlines()[GM_r[0]].strip().split(':')[1])
	with open(fp_full, 'r') as f:
		r_value = np.float(f.readlines()[GM_r[1]].strip().split(':')[1])

	with open(fp_full, 'r') as f:
		for line in f:
			if string_eoY in line:
				header2skip = linecount
			linecount += 1
	df = pd.read_csv(fp_full, skiprows=header2skip, delim_whitespace=True, usecols=[1, 2, 3, 4], index_col=False,
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



	fp = fp_static
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

	return arr_C_delta, arr_S_delta, gm_value, r_value, date_start, date_end

# squared degree variances
def degree_varience_sq(arr_C_delta, arr_S_delta, maxdeg):
	sigma_l_squared = []
	[sigma_l_squared.append(np.nansum(np.power(arr_C_delta[ll, :], 2) + np.power(arr_S_delta[ll, :], 2))) for ll in
	 range(0, maxdeg)]
	arr_sigma_squared = np.array(sigma_l_squared)
	return arr_sigma_squared

# Degree variances for disturbing potential
def degree_varience_disturbing_pot(gm,r, sigma):
	sigma_DP = gm/ r * sigma
	return sigma_DP

# Degree variances for undulations
def degree_variance_un(r, sigma):
	sigma_Un = r * sigma
	return sigma_Un

# Degree variances for Gravity disturbances
def degree_varience_GD(gm, r, sigma, maxdeg):
	GD = []
	for l in range(0, maxdeg):
		GD.append(gm / r ** 2 * sigma[l] * (l + 1))
	sigma_GD = np.array(GD)
	return sigma_GD

# Degree variances for Gravity gradient
def degree_varience_GG(gm, r, sigma, maxdeg):
	GG = []
	for l in range(0, maxdeg):
		GG.append(gm / r ** 3 * sigma[l] * (l + 1) * (l + 2))
	sigma_GG = np.array(GG)
	return sigma_GG

def read_LLN_internalformat(fp, maxdeg):
	df = pd.read_csv(fp, skiprows=14, delim_whitespace=True,
					 names=['n', 'h', 'nl', 'nk', 'h_asymp', 'nl_asymp', 'nk_asymp'])
	df.drop(['n', 'h', 'nl', 'h_asymp', 'nl_asymp', 'nk_asymp'], axis=1, inplace=True)
	df = df[0:maxdeg]
	return df.to_numpy()

def read_LLN_from_LOADEV(fp, maxdeg):
	df = pd.read_csv(fp, delim_whitespace=True, skiprows=14,
					  names=['n', 'h', 'nl', 'nk', 'h_asymptotic', 'nl_asymptotic', 'nk_asymptotic'])
	df['kndivn'] = df['nk'] / df['n']
	df_out = pd.DataFrame()
	df_out['kn'] = df['kndivn'][0:maxdeg]
	return df_out.to_numpy()

def sigma_EWH(rho_e, rho_w, r, LLN, maxdeg, sigma):
	EWH = []
	C = r*rho_e/(3*rho_w)
	for l in range(0, maxdeg):
		EWH.append(C*(2*l+1)/(1+LLN[l])*sigma[l])
	sigma_EWH = np.array(EWH)
	return sigma_EWH










