import numpy as np
import pandas as pd
from train import model_sine
from validation_test_set import data_x_train_test, t_test
from sklearn.metrics import r2_score
from scipy.stats import pearsonr

# Training evaluation 
# sim = model_sine.simulate(np.array(data_x_train[0]), t_train[:,0]) 

# Model stability evaluation
sim = model_sine.simulate(np.array(data_x_train_test[0]), t_test[:,0]) 

sim.shape

#------------------------------------------------------------#
# # R2_score_train = r2_score(data_x_train[:, 1], sim[:, 1])
# # R2_score_train = "{:.3f}".format(R2_score_train)
# Pearson_correlation_train = pearsonr(data_x_train[:, 0], sim[:, 0])
# Pearson_correlation_train_noisy = pearsonr(data_x_train_noisy[:, 0], sim_noisy[:, 0])
# Pearson_correlation_train = "{:.3f}".format(Pearson_correlation_train[0])

# #print(f"The observed R2_score for the training_set is:",R2_score_train)
# print(f"The Pearson_correlation for the training_set is :", Pearson_correlation_train)
# print(f"The Pearson_correlation for the noisy training_set is :", Pearson_correlation_train_noisy)

R2_score_test = r2_score(data_x_train_test[:, 0], sim[:, 0])
R2_score_test = "{:.3f}".format(R2_score_test)
Pearson_correlation_test = pearsonr(data_x_train_test[0:, 0], sim[:, 0])
Pearson_correlation_test = "{:.3f}".format(Pearson_correlation_test[0])

print(f"The observed R2_score for the validation set is:",R2_score_test)
print(f"The Pearson_correlation for the validation set is :", Pearson_correlation_test )