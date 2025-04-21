import math

import pandas as pd
import numpy as np

def choice_a2f():
	"""
	输出a2f谱函数
	:return:
	"""
	column_names = ['E(THz)', '0.005', '0.010', '0.015', '0.020', '0.025', '0.030', '0.035', '0.040', '0.045', '0.050']
	a2f_total = pd.read_csv('alpha2F.dat', delim_whitespace=True, skiprows=1, names=column_names)
	data_fre = a2f_total["E(THz)"]

	num1 = input("选择相应的展宽: ")
	num1 = '%.3f' % float(num1)
	data_total = a2f_total[num1]
	all_data = pd.concat([data_fre, data_total], axis=1)
	"""
	修改列表名
	"""
	all_data.columns.values[1] = 'total-a2f'
	all_data.to_csv('a2f_total.dat', sep=" ", index=False, float_format='%.5f')


def cal_Allen_Dynes_equation():
	"""
    计算Allen_Dynes_equation的Tc，包括带f1和f2的
    :return:
    """
	"""
    读取a2f_total.dat文件，先计算a2f*2/w。
    """
	# choice_a2f()
	u = float(input("请输入u*的值："))
	data = pd.read_csv('a2f_total.dat', sep=' ')
	a2f_freq = data['E(THz)']
	new_freq = []
	for i in a2f_freq:
		if float(i) == 0:
			i = 1e-2
			new_freq.append(i)
		else:
			new_freq.append(i)
	a2f_freq = new_freq
	a2f_total = data['total-a2f']
	data['total_cal'] = (a2f_total * 2) / a2f_freq
	data.replace([np.nan, np.inf, -np.inf], 0, inplace=True)
	# data.to_csv("output.dat", sep=' ', index=False, float_format='%.5f')
	data_values = data['total_cal']
	# data_values = pd.read_csv("output.dat",sep=' ')['total_cal']
	""" 计算对应展宽下的λ和ωlog的值"""
	# for i in range(len(a2f_total)):
	# 	integ = np.trapz(a2f_total[:i + 1], x=a2f_freq[:i + 1])
	# 	integ = round(integ, 5)
	total_lambda = np.trapz(data_values, x=a2f_freq)
	total_lambda = round(total_lambda, 5)

	e = math.e
	array1 = np.array(np.log(a2f_freq))
	array2 = np.array(a2f_total)
	array3 = np.array(a2f_freq)
	values = array1 * array2 * 2 / (array3 * total_lambda)
	va = np.trapz(values, a2f_freq)
	THz_to_K = 33.356 * 1.439
	omega_log = math.e ** va * THz_to_K
	omega_log = round(omega_log, 5)

	"""计算ω2"""
	array_new = (2 / total_lambda) * array2 * array3
	omega_2 = np.sqrt(np.trapz(array_new, x=array3)) * THz_to_K
	omega_2 = round(omega_2, 5)
	"""计算f1和f2"""
	f1 = (1 + (total_lambda / (2.46 * (1 + 3.8 * u))) ** (3 / 2)) ** (1 / 3)
	f2 = 1 + (((omega_2 / omega_log) - 1) * (total_lambda ** 2)) / (
			(total_lambda ** 2) + ((1.82 * (1 + 6.3 * u)) * (omega_2 / omega_log)) ** 2)
	f1 = round(f1, 5)
	f2 = round(f2, 5)
	"""计算Allen-Dynes-equation计算Tc"""
	Tc = (omega_log / 1.2) * (
			e ** ((-1.04 * (1 + total_lambda)) / (total_lambda - u * (1 + 0.62 * total_lambda))))
	Tc_with_f1f2 = (f1 * f2 * omega_log / 1.2) * (
			e ** ((-1.04 * (1 + total_lambda)) / (total_lambda - u * (1 + 0.62 * total_lambda))))
	Tc = round(Tc, 5)
	Tc_with_f1f2 = round(Tc_with_f1f2, 5)
	message = {'λ:': total_lambda, 'ωlog': f'{omega_log}K', 'ω2': f'{omega_2}K', 'f1': f1, 'f2': f2, 'Tc': f'{Tc}K',
			   'Tc_with_f1f2': f'{Tc_with_f1f2}K'}
	for k, v in message.items():
		print(f'{k} : {v}')

choice_a2f()
cal_Allen_Dynes_equation()
