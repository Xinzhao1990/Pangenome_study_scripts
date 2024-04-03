import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt


# Custom nonlinear regression function
def custom_function_pan(params, x, y1):
    a, b, c = params
    y1_fit = a * np.power(x, b) + c
    residual_pan = y1_fit - y1
    return residual_pan

def custom_function_core(params, x, y2):
    A, B, C = params
    y2_fit = A * np.exp(-B * x) + C
    residual_core = y2_fit - y2
    return residual_core

    
# input data 
file_path = '/Users/jiangwengyao/Desktop/application/SURF23/pangenome_for_panGP_plot.txt'
data = np.loadtxt(file_path, skiprows=1) 
data_pan = data[:, 1]
data_core = data[:, 2]

grouped_data_pan = [data_pan[i:i+10] for i in range(0, len(data_pan), 10)]
mean_values_pan = [sum(group) / len(group) for group in grouped_data_pan]
#std_values_pan = [np.std(group) for group in grouped_data_pan]

grouped_data_core = [data_core[i:i+10] for i in range(0, len(data_core), 10)]
mean_values_core = [sum(group) / len(group) for group in grouped_data_core]
#std_values_core = [np.std(group) for group in grouped_data_core]

x = np.arange(1, 62, 1)
y1 = mean_values_pan
y2 = mean_values_core

# Use least_squares for fitting
initial_guess = [1000, 0.4, 200]  # Initial parameter guess value

# Carry out nonlinear least square fitting
params_fit_pan = least_squares(custom_function_pan, initial_guess, args=(x, y1)).x
a_fit, b_fit, c_fit = params_fit_pan

params_fit_core = least_squares(custom_function_core, initial_guess, args=(x, y2)).x
A_fit, B_fit, C_fit = params_fit_core

# Calculate the fitting result
y_pan = a_fit * np.power(x, b_fit) + c_fit
y_core = A_fit * np.exp(-B_fit * x) + C_fit


# Print fitting parameters
print("parameters A_pan = ", params_fit_pan[0])
print("parameters B_pan = ", params_fit_pan[1])
print("parameters C_pan = ", params_fit_pan[2])

print("parameters A_core = ", params_fit_core[0])
print("parameters B_core = ", params_fit_core[1])
print("parameters C_core = ", params_fit_core[2])

# Draw a fitting curve
plt.figure(figsize=(16, 13))
# Draw a vertical box diagram
positions = range(1, 62)
plt.rcParams['font.family'] = 'Arial'

#plt.boxplot(grouped_data_pan, vert=True, showfliers=False, positions=positions)
plt.boxplot(grouped_data_pan, vert=True, showfliers=False, positions=positions, medianprops={'color': 'black'})

plt.boxplot(grouped_data_core, vert=True, showfliers=False, positions=positions, medianprops={'color': 'black'})

plt.plot(x, y_pan, label="Pan genome size", color='blue', linewidth=3)
plt.plot(x, y_core, label="Core genome size", color='red', linewidth=3)

plt.xlabel("Number of genome", fontsize=20, labelpad=15)
plt.ylabel("Number of gene cluster", fontsize=20, labelpad=15)

plt.legend(fontsize=18, frameon=False)

plt.yticks(fontsize=16)
#plt.xticks(ticks=np.arange(0,62,2), fontsize=10)
plt.xticks(ticks=np.arange(0,62,2), labels=np.arange(0,62,2), fontsize=16)

plt.savefig('aaa.png', dpi=300)
plt.show()