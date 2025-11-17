# Masterarbeit
Mathematical modeling of source-nonlinearity in Electro-Magnetic Sensors.
Main steps and results are found under summary_thesis. The focus of the thesis was leveraging SINDy algorithm to predict the source non-linearity (i.e., harmonics) arising due experimental setup (e.g., Amplifiers, interfaces etc.) even the actuator is being in its linear regime operated. 

Sparse Identification of Nonlinear Dynamics (SINDY)
We have the following equation for SINDY modelling;
<img width="979" height="198" alt="image" src="https://github.com/user-attachments/assets/08256044-f28e-4249-950b-be84228d7c5e" /> 
To apply SINDy we need a set of measurement data that collected at different times (matrix X). Then we create a time derivative of the X matrix numerically. 
The matrix Θ(X) could be chosen from the library for different set of basis functions that we want to apply on the data.
Then we use a kind of inverse problem optimization to find the set of sparse coefficient vectors as follow and fit the model.
 <img width="331" height="146" alt="image" src="https://github.com/user-attachments/assets/5247f14e-eea0-4c8c-b312-eb28d547cb22" />

And at the final step we have a differential equation that could predict the effect of input parameters on the output parameters.

 

