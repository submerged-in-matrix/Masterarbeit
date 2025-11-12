from data.Processing.ANNs.seq_2_samples import n_steps, split_sequence, X, y
from data.Processing.ANNs.load_data import x_rec
from data.Processing.ANNs.test_data import X_test, y_test
from src.models.mlp.dependencies_mlp import *

model = Sequential()
model.add(Dense(200, activation='relu', input_dim=n_steps))
model.add(Dense(200, activation='relu', input_dim=n_steps))
model.add(Dense(1))
model.compile(optimizer='adam', loss='mse')
# fit model
history = model.fit(X, y, epochs=3000, validation_data=(X_test, y_test), verbose=0)
#history = model.fit(trainX, trainY, epochs=50, batch_size=100, validation_data=(testX, testY), verbose=1, shuffle=False)
model.save('MLP_sine_final.hdf5')
model.summary()