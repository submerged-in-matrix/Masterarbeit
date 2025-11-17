from data.Processing.ANNs.seq_2_samples import split_sequence
from data.Processing.ANNs.load_data import x_rec
from data.Processing.ANNs.seq_2_samples import n_steps

n_features = 1
# Validation data set 60k-63k.Test set varies
X_test, y_test = split_sequence(x_rec[112000:113004, 0], n_steps)
print(X_test.shape)
print(y_test.shape)