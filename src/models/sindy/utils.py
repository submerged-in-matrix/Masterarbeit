# Finding best dataset fragment for training and validation. 
# Only shown for the displacement training dataset which is also the same approach for the voltage dataset.
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from .train import generalized_library, sfd, stlsq
from data.integrate_velocity import displacement, sampling_rate, t_array

window_size = 6000
start_indices = np.arange(40001, 100001, window_size)  # Start indices for each portion
t_array_new = t_array
errors = []
errors_p = []
for start_index in start_indices:
    end_index = start_index + window_size
    portion_sig =displacement[start_index:end_index]
    portion_time = t_array_new[start_index:end_index]
    portion_xdot = sfd._differentiate(portion_sig, portion_time[:, 0] )
    data_x_train_new = np.hstack((portion_sig, portion_xdot, portion_time))
    model_sine=ps.SINDy(
    optimizer=stlsq,
    #optimizer=ensemble_optimizer,
    #optimizer=sr3,
    differentiation_method= sfd,
    # differentiation_method=diff,
    feature_names=["x", "xdot", "t"],
    discrete_time=False,
    feature_library=generalized_library 
    )
    model_sine.fit(data_x_train_new, 1/sampling_rate)
    
    sim_new = model_sine.simulate(np.array(data_x_train_new[0]), portion_time[:, 0], integrator="odeint") 
    # Calculate and store the error metrics of each portion

    error = r2_score(data_x_train_new[:, 0], sim_new[:, 0])
    error_p = pearsonr(data_x_train_new[:, 0], sim_new[:, 0])
    errors.append(error)
    errors_p.append(error_p)
    
    print(f"Start index: {start_index}, Error: {error, error_p}")
