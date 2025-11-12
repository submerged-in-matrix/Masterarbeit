from data_smoothing import sfd
import numpy as np
from integrate_velocity import displacement, t_array 
## validation 100001:103001. test_set varies (102k-105k)
x_test = displacement[100001:103001]
t_test = t_array[100001:103001]
# t_test = np.array(t_test)
# t_test = t_test[:,np.newaxis]

xdot_test = sfd._differentiate(x_test, t_test[:, 0])
data_x_train_test = np.hstack((x_test, xdot_test, t_test))
print(data_x_train_test.shape)