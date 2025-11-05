import numpy as np
import pysindy as ps
from data.data_smoothing import sampling_rate, data_x_train, t_train

poly_order = 3
max_iter = 10

# for STLSQ
regularization_on_weight = 0.05
threshold = 0.2

#variable thresholding for SR3 optimizers
# thresholds = 0.05* np.ones((11, 3))
# thresholds[2 :, :] = 0.5
# nu = 0.3
# tol= 0.00001

#ensemble_optimizer = ps.EnsembleOptimizer(ps.STLSQ(threshold=threshold, alpha=regularization_on_weight, normalize_columns=False, max_iter=max_iter), bagging=True, replace=False, n_subset=1000, n_models=6, ensemble_aggregator=np.median)
stlsq = ps.STLSQ(threshold=threshold, max_iter=max_iter, alpha=regularization_on_weight)
#sr3 = ps.SR3(thresholder="weighted_l0", thresholds=thresholds, max_iter=max_iter, nu=nu, tol=tol)

#Initialize Libraries

#lib_poly = ps.feature_library.PolynomialLibrary(degree = poly_order, include_interaction=False, include_bias=False)

#lib_fourier = ps.feature_library.FourierLibrary(n_frequencies = 1, include_cos=False, include_sin=True)

funct = [lambda x : x,
         lambda x : x ** 3 ,
         #lambda x,y : (x* y)
        ]

names = [lambda x : '(' + x + ')',
         lambda x : '(' + x + '^' '3' ')',
         #lambda x, y : '(' + x + '*' + y + ')'
         ]
lib_xdot = ps.CustomLibrary(library_functions=funct, function_names=names)

lib_id = ps.IdentityLibrary()

# For displacement dataset
functions = [lambda x :  np.cos(x * 7420),
            lambda x :  np.cos(x * 14840),
            lambda x :  np.cos(x * 22270),
            lambda x :  np.cos(x * 29680),
            lambda x :  np.cos(x * 37090),
            # lambda x :  np.cos(x * 9648),
            # lambda x :  np.cos(x * 10389),
            # lambda x :  np.cos(x * 11880),
            # lambda x :  np.cos(x * 13362),
            # lambda x :  np.cos(x * 14095),
            # lambda x :  np.cos(x * 14837),
            ]

function_names = [lambda x : 'cos(' + x + '*' '7420' ')',
                  lambda x : 'cos(' + x + '*' '14840' ')',
                  lambda x : 'cos(' + x + '*' '22270' ')',
                  lambda x : 'cos(' + x + '*' '29680' ')',
                  lambda x : 'cos(' + x + '*' '37090' ')',
                #   lambda x : 'cos(' '12590', '*', + x + ')',
                #   lambda x : 'cos(' '15100', '*', + x+ ')'
                    ]
# # For voltage dataset
# functions = [lambda x :  np.cos(x * 2517),
#              lambda x :  np.cos(x * 5034),
#              lambda x :  np.cos(x * 7551),
#              lambda x :  np.cos(x * 10068),
#              lambda x :  np.cos(x * 12585)
#             ]
# function_names = [lambda x : 'cos(' + x + '*' '2517' ')',
#                   lambda x : 'cos(' + x + '*' '5034' ')',
#                   lambda x : 'cos(' + x + '*' '7551' ')',
#                   lambda x : 'cos(' + x + '*' '10068' ')',
#                   lambda x : 'cos(' + x + '*' '12585' ')'
#                  ]  

lib_custom = ps.CustomLibrary(library_functions=functions, function_names=function_names)
 
inputs_per_library = [[0, 0, 1],  [2, 2, 2]]

inputs_per_library = np.reshape(inputs_per_library, (2, 3))

# Tensor all the polynomial and Fourier library terms together
tensor_array = [[0, 1, 1]]

generalized_library = ps.GeneralizedLibrary([lib_xdot, lib_custom], 
                                            #tensor_array=tensor_array,
                                            #exclude_libraries=[1],
                                            inputs_per_library=inputs_per_library)
# Lets define model architecture 
model_sine=ps.SINDy(
    optimizer=stlsq,
    #optimizer=ensemble_optimizer,
    #optimizer=sr3,
    differentiation_method= sfd,
    feature_names=["x", "xdot", "t"],
    discrete_time=False,
    feature_library=generalized_library
    
    )
# Lets train the model
model_sine.fit(data_x_train, 1/sampling_rate)

# See the discovered equation and the used features, i.e., the candidate terms
model_sine.print(lhs=['dx/dt', 'd(xdot)/dt', 'dt/dt'])
model_sine.get_feature_names()