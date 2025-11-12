#  Optimizing the hyperparameters. 
# The SINDy models have 4 hyperparameters: Threshold , regularization strength (alpha),  window length,  
# and the total number of epochs (has almost no effect due to the library settings).
# The optimization scheme is shown for the threshold value which is to be repeated for all the hyperparameters.

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import root_mean_squared_error, mean_squared_error
from data.data_smoothing import sampling_rate, sfd
from data.validation_test_set import data_x_train, data_x_train_test, t_test
from models.sindy.train import poly_order, max_iter, stlsq, lib_xdot, lib_id, funct, names, generalized_library

# Make coefficient plot for threshold scan
def plot_pareto(coefs, opt, model, threshold_scan, x_test, t_test):
    dt = 1/sampling_rate
    mse = np.zeros(len(threshold_scan))
    mse_sim = np.zeros(len(threshold_scan))
    for i in range(len(threshold_scan)):
        opt.coef_ = coefs[i]
        mse[i] = model.score(x_test, t=dt, metric=mean_squared_error)
        x_test_sim = model.simulate(x_test[0], t_test[:,0], integrator="odeint")
        # if np.any(x_test_sim > 1e4):
        #     x_test_sim = 1e4
        mse_sim[i] = np.sum((x_test[:,0] - x_test_sim[:, 0]) ** 2)
    
    plt.figure()
    plt.semilogy(threshold_scan, mse, "ro")
    plt.semilogy(threshold_scan, mse, "g")
    plt.title(" Effect of threshold value on xdot")
    plt.ylabel(r"$\dot{X}$ RMSE", fontsize=20)
    plt.xlabel(r"threshold values", fontsize=20)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.grid(True)
    
    plt.figure()
    plt.semilogy(threshold_scan, mse_sim, "yo")
    plt.semilogy(threshold_scan, mse_sim, "b")
    plt.title(" Effect of threshold value on simulated xdot")
    plt.ylabel(r"$\dot{X}$ RMSE", fontsize=20)
    plt.xlabel(r"threshold values", fontsize=20)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.grid(True)

#--------------------------------------------------------------#
threshold_scan = np.arange(0.0, 1.05, 0.05)
coefs = []
rmse = root_mean_squared_error(data_x_train, np.zeros(data_x_train.shape), squared=False)

#x_train_added_noise = x_train + np.random.normal(0, rmse / 10.0,  x_train.shape)

for i, threshold in enumerate(threshold_scan):
    sparse_regression_optimizer = ps.STLSQ(threshold=threshold)
    model_sine = ps.SINDy(
    optimizer=sparse_regression_optimizer,
    #optimizer=sr3,
    differentiation_method= sfd,
    feature_names=["x", "xdot", "t"],
    discrete_time=False,
    feature_library=generalized_library
    )
    model_sine.fit(data_x_train, 1/sampling_rate, quiet=True)
    #model_sine.fit(data_x_train, 1/sampling_rate, ensemble=True, replace=False, n_subset=1000, n_models=3, quiet=True, unbias=False, ensemble_aggregator=np.median)
    coefs.append(model_sine.coefficients())

plot_pareto(coefs, sparse_regression_optimizer, model_sine,
            threshold_scan, data_x_train_test, t_test)