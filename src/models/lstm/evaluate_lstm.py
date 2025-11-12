from  src.models.lstm.lstm_dependencies import *
from src.models.lstm.data_split import X, y, X_test, y_test, n_steps, n_features
from src.models.lstm.train_lstm import model, history
from keras.models import load_model


# demonstrate prediction
#input = X_test
input = X
prediction_list = []
for row in input:
    x_input = array(row)
    x_input = x_input.reshape((1, n_steps, n_features))
    
    # making predictions
    yhat = model.predict(x_input, verbose=0)
    prediction_list.append(yhat)


predicted_values = np.array(prediction_list)

R2_score_train = r2_score(y, predicted_values[:, 0])
R2_score_train = "{:.4f}".format(R2_score_train)
print(f"The observed R2_score for the voltage training_set is:",R2_score_train)


# R2_score_test = r2_score(X_test[4:, 0], predicted_values[0:-4, 0], force_finite=False)
# R2_score_test = "{:.4f}".format(R2_score_test)
# print(f"The observed R2_score for the test_set is:",R2_score_test)

fig,ax = plt.subplots(1,1, figsize=(8,6))
ax.set_ylabel('Voltage (mV)', fontsize=12)
ax.set_xlabel('Time (sec)', fontsize=12)
print(predicted_values.shape)

#For training set
plt.plot (X[4:1004, 0], 'r', label="signal response")
plt.plot(predicted_values[0:1000, 0], 'g', label= "LSTM prediction")

# For Laser data 
ax.set_title(f'R2_score for displacement training_set is : {R2_score_train} ')

# # For voltage data
# ax.set_title(f'R2_score for voltage training_set with LSTM is : {R2_score_train}')
 
#For test_set
# plt.plot (X_test[4:1004, 0], 'r', label="signal response")
# plt.plot(predicted_values[0:1000, 0], 'g', label= "LSTM prediction")

# # For Laser data 
# ax.set_title(f'R2_score for displacement validation set with LSTM is : {R2_score_test} ')

#For voltage data
# ax.set_title(f'R2_score for voltage test_set with LSTM is : {R2_score_test}')

ax.legend(loc='upper left', bbox_to_anchor=(1.0005, 1))

## Plot History
# plot history
from matplotlib import pyplot
pyplot.plot(history.history['loss'], label='train')
pyplot.plot(history.history['val_loss'], label='test')
#pyplot.legend()

pyplot.title('model loss',size=15)
pyplot.ylabel('loss',size=15)
pyplot.xlabel('epochs',size=15)
pyplot.legend(loc='upper right',fontsize=15)

pyplot.show()
print(predicted_values.shape)