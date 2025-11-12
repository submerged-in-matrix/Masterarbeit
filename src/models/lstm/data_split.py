from data.Processing.ANNs.seq_2_samples import split_sequence
from data.Processing.ANNs.load_data import x_rec

# choose a number of time steps
n_steps = 3
# split into samples
X, y = split_sequence(x_rec[40000:42004, 0], n_steps)

# reshape from [samples, timesteps] into [samples, timesteps, features]
n_features = 1
X = X.reshape((X.shape[0], X.shape[1], n_features))
print(X.shape)
print(y.shape)

n_features = 1
# validation set 100k-112k
X_test, y_test = split_sequence(x_rec[112000:113004, 0], n_steps)
X_test = X_test.reshape((X_test.shape[0], X_test.shape[1], n_features))
print(X_test.shape)
print(y_test.shape)

