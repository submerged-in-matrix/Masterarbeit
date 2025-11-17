from scipy.io import loadmat
from data.Processing.ANNs.Vel_2_disp import displacement
# # For voltage data
# data_file =loadmat("10filtered_700_13000.mat") 
# x_rec = data_file["low_filtered"]
# print(x_rec[1:5, 0])
# print(x_rec.shape)

# for displacement data
#data_file =loadmat("integrated displacement sine.mat") 
x_rec = displacement
print(x_rec[1:5, 0])
print(x_rec.shape)