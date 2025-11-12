from src.models.mlp.dependencies_mlp import *
from src.models.mlp.grid_search_mlp import *

### Effect of layer_size on the model performance
scores = np.array(scores)
plt.figure()
plt.semilogy(scores[:, 0], scores[:, 1], "ro")
plt.semilogy(scores[:, 0], scores[:, 1], "g")
plt.title(" Effect of layer size on the model")
plt.ylabel(r" RMSE", fontsize=15)
plt.xlabel(r"Configurations", fontsize=15)
plt.xticks(fontsize=10, rotation=30)
plt.yticks(fontsize=10)
plt.grid(True)

### Effect of input_size size (series to supervised) on the model performance
# scores = np.array(scores)
# # Plot the data
# plt.figure()
# plt.semilogy(scores[:, 0], scores[:, 1], "ro")
# plt.semilogy(scores[:, 0], scores[:, 1], "g")
# plt.title(" Effect of input_size size on the model")
# plt.ylabel(r" RMSE", fontsize=15)
# plt.xlabel(r"Configurations", fontsize=15)
# plt.xticks(fontsize=10, rotation=30)
# plt.yticks(fontsize=10)
# plt.grid(True)

### Effect of training_epochs on the model performance
scores = np.array(scores)
plt.figure()
plt.semilogy(scores[:, 0], scores[:, 1], "ro")
plt.semilogy(scores[:, 0], scores[:, 1], "g")
plt.title(" Effect of number of epochs on the model")
plt.ylabel(r" RMSE", fontsize=15)
plt.xlabel(r"Configurations", fontsize=15)
plt.xticks(fontsize=10, rotation=30)
plt.yticks(fontsize=10)
plt.grid(True)