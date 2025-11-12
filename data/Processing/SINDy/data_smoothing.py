from integrate_velocity import displacement
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from pysindy.differentiation import SmoothedFiniteDifference

# Time base for displacement dataset
rec_time = 0.1
sampling_rate = 1250000

# # Time base for voltage dataset
# rec_time = 10
# sampling_rate = 96000

start_time = 0
t_array = np.arange(start_time, rec_time, 1/sampling_rate)
t_array = np.array(t_array)
t_array= t_array[:,np.newaxis]

# # Loading voltage dataset
# data_file = loadmat("10filtered_700_13000.mat")
# read_out_raw = data_file["low_filtered"]

# # Creating voltage training dataset 
# x_rec = read_out_raw[196001:200001]
# t_train = t_array[196001:200001]

# Creating displacement training dataset 
x_rec = displacement[46001:52001]
t_train = t_array[46001:52001]

#---------------------------------------------#
# Smoothing is being performed via Savgol Golay filter, which performs a moving 3rd order polynomial approximation on the dataset.
# Trying out different differentiators

#x_rec =np.abs(signal_filtered)
#fd = FiniteDifference(drop_endpoints=True)
sfd = SmoothedFiniteDifference(smoother_kws={'window_length': 5} )
#sfd = SmoothedFiniteDifference()
#diff = ps.SINDyDerivative(kind="trend_filtered", order=1, alpha=1e-3)

#xdot = fd._differentiate(x_rec, t_train)
xdot = sfd._differentiate(x_rec, t_train[:, 0])
#xdot = diff._differentiate(x_rec,t_train)

# Here dataset is being prepared for feeding to SINDy, with which SINDy knows that we are loooking for a 2nd order ODE. 
data_x_train = np.hstack((x_rec, xdot, t_train))

print(data_x_train.shape)
print(t_train.shape)
print(xdot.shape)
# plt.plot(x_rec [:500, 0])