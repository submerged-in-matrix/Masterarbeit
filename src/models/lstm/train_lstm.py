from src.models.lstm.lstm_dependencies import *
from src.models.lstm.data_split import X, y, X_test, y_test, n_steps, n_features

model = Sequential()
model.add(LSTM(150, activation='relu', return_sequences= True, input_shape=(n_steps, n_features)))
model.add(LSTM(150, activation='relu'))
model.add(Dense(1))
model.compile(optimizer='adam', loss='mse')
# fit model
history = model.fit(X, y, epochs=2000, validation_data=(X_test, y_test), verbose=0)
model.save('LSTM_sine_final.hdf5')
model.summary()