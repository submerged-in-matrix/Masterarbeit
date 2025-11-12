from src.models.mlp.dependencies_mlp import *
from tensorflow.keras.models import load_model
from data.Processing.ANNs.seq_2_samples import n_steps
from data.Processing.ANNs.test_data import X_test, y_test
from data.Processing.ANNs.seq_2_samples import split_sequence, X, y
from src.models.mlp.train_mlp import history
model = load_model('MLP_sine_final.hdf5')

#input = X_test
input = X
prediction_list = []
for row in input:
    x_input = array(row)
    x_input = x_input.reshape((1, n_steps))
    
    # making predictions
    yhat = model.predict(x_input, verbose=0)
    prediction_list.append(yhat)
predicted_values = np.array(prediction_list)

R2_score_train = r2_score(X[4:, 0], predicted_values[0:-4, 0], force_finite=False)
R2_score_train = "{:.4f}".format(R2_score_train)
print(f"The observed R2_score for the training_set is:",R2_score_train)

# R2_score_test = r2_score(X_test[4:, 0], predicted_values[0:-4, 0], force_finite=False)
# R2_score_test = "{:.4f}".format(R2_score_test)
# print(f"The observed R2_score for the test_set is:",R2_score_test)

# Evaluation Plots
from matplotlib import pyplot
import matplotlib.pyplot as plt
fig,ax = plt.subplots(1,1, figsize=(8,6))
ax.set_ylabel('volatage (mV)', fontsize=12)
ax.set_xlabel('Time (sec)', fontsize=12)
print(predicted_values.shape)

# For training set
# plt.plot (X[4:1004, 0], 'r', label="signal response")
# plt.plot(predicted_values[0:1000, 0], 'g', label= "MLP prediction")

# # For Laser data 
# ax.set_title(f'R2_score for displacement training_set is : {R2_score_train} ')

# # For voltage data
# ax.set_title(f'R2_score for voltage training_set with MLP is : {R2_score_train}')
 
# For test_set
plt.plot (X_test[4:1004, 0], 'r', label="signal response")
plt.plot(predicted_values[0:1000, 0], 'g', label= "MLP prediction")

# # For Laser data 
# ax.set_title(f'R2_score for displacement test set with MLP is : {R2_score_test}')

#For voltage data
ax.set_title(f'R2_score for voltage test_set with MLP is : {R2_score_test}')

ax.legend(loc='upper left', bbox_to_anchor=(1.0005, 1))

## History plotting

# plot history of training and validation/test loss
pyplot.plot(history.history['loss'], label='train')
pyplot.plot(history.history['val_loss'], label='test')
#pyplot.legend()

pyplot.title('model loss',size=15)
pyplot.ylabel('loss',size=15)
pyplot.xlabel('epochs',size=15)
pyplot.legend(loc='upper right',fontsize=15)

pyplot.show()