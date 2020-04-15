from root_numpy import root2array, tree2array, list_structures
from root_numpy import testdata
import numpy as np
import pandas
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression


def get_index(name_string, df_names):

    for i in range(0,len(df_names.values[0][0])):
        if(df_names.values[0][0][i] == name_string):
            index = i
            break

    return index

filename = 'GenLevelOutput.root'

# Convert a TTree in a ROOT file into a NumPy structured array + pandas Dataframe, with which you can create an array to be given to Keras
branches = ["EFTWeights"]
branches_names = ["EFTWeightsNames"]
# one needs probably a flat tree, to access information

weight_dataframe = pandas.DataFrame(root2array(filename, treename='boosted/events',branches=branches))
weight_dataframe_names = pandas.DataFrame(root2array(filename, treename='boosted/events',branches=branches_names))

parameter_to_scan = "ctq8" #"ctG"

weight_SM = weight_dataframe.applymap(lambda x: x[get_index("rwgt_SM",weight_dataframe_names)])
weight_EFT_min5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_"+str(parameter_to_scan)+"_min5",weight_dataframe_names)])
weight_EFT_min4 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_"+str(parameter_to_scan)+"_min4",weight_dataframe_names)])
weight_EFT_min3 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_"+str(parameter_to_scan)+"_min3",weight_dataframe_names)])
weight_EFT_min2 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_"+str(parameter_to_scan)+"_min2p0",weight_dataframe_names)])
weight_EFT_min1 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_"+str(parameter_to_scan)+"_min1p0",weight_dataframe_names)])
weight_EFT_1 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_"+str(parameter_to_scan)+"_1p0",weight_dataframe_names)])
weight_EFT_2 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_"+str(parameter_to_scan)+"_2p0",weight_dataframe_names)])
weight_EFT_3 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_"+str(parameter_to_scan)+"_3p0",weight_dataframe_names)])
weight_EFT_4 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_"+str(parameter_to_scan)+"_4p0",weight_dataframe_names)])
weight_EFT_5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_"+str(parameter_to_scan)+"_5p0",weight_dataframe_names)])

merged_df = weight_SM.merge(weight_EFT_1, left_index=True, right_index=True)
merged_df = merged_df.merge(weight_EFT_2, left_index=True, right_index=True)
print(merged_df.corr())

fig_corr = plt.figure(figsize=(12,8))
plt.matshow(merged_df.corr())
plt.xticks(range(merged_df.shape[1]), merged_df.columns, fontsize=14, rotation=45)
plt.yticks(range(merged_df.shape[1]), merged_df.columns, fontsize=14)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=14)
plt.title('Correlation Matrix', fontsize=16);
plt.savefig('corr.png')

x_values = [-5, -4,-3, -2, -1, 0, 1, 2, 3, 4, 5]
#['rwgt_ctG_1p0', 'rwgt_ctG_4p0', 'rwgt_ctG_3p0', 'rwgt_ctG_min3','rwgt_ctG_min4', 'rwgt_ctG_min2p0', 'rwgt_ctG_min1p0', 'rwgt_SM','rwgt_ctG_5p0', 'rwgt_ctG_min5', 'rwgt_ctG_2p0']
#'rwgt_cQq8_0p0_ctq8_0p0_ctG_min5p0_ctGI_0p0', 'rwgt_ctG_1p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_5p0_ctGI_0p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_min2p5_ctGI_5p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_min5p0_ctGI_5p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_0p0_ctGI_2p5', 'rwgt_ctG_4p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_min2p5_ctGI_0p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_5p0_ctGI_2p5',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_0p0_ctGI_5p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_2p5_ctGI_2p5',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_2p5_ctGI_2p5', 'rwgt_ctG_3p0',
#       'rwgt_ctG_min3', 'rwgt_ctG_min4', 'rwgt_ctG_min2p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_5p0_ctGI_min5p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_2p5_ctGI_min5p0', 'rwgt_ctG_min1p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_min5p0_ctGI_2p5',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_min5p0_ctGI_min5p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_0p0_ctGI_min5p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_min2p5_ctGI_2p5',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_min2p5_ctGI_min2p5',
#      'rwgt_cQq8_0p0_ctq8_0p0_ctG_2p5_ctGI_min2p5',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_min5p0_ctGI_min2p5',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_2p5_ctGI_0p0', 'rwgt_SM',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_min2p5_ctGI_min5p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_0p0_ctGI_min2p5', 'rwgt_ctG_5p0',
#       'rwgt_cQq8_0p0_ctq8_0p0_ctG_5p0_ctGI_5p0', 'rwgt_ctG_min5',
#       'rwgt_ctG_2p0', 'rwgt_cQq8_0p0_ctq8_0p0_ctG_5p0_ctGI_min2p5']

y_values = [weight_EFT_min5.mean(), weight_EFT_min4.mean(), weight_EFT_min3.mean(), weight_EFT_min2.mean(), weight_EFT_min1.mean(), weight_SM.mean(), weight_EFT_1.mean(), weight_EFT_2.mean(), weight_EFT_3.mean(), weight_EFT_4.mean(), weight_EFT_5.mean()]

y_values_err = [weight_EFT_min5.std(), weight_EFT_min4.std(), weight_EFT_min3.std(), weight_EFT_min2.std(), weight_EFT_min1.std(), weight_SM.std(), weight_EFT_1.std(), weight_EFT_2.std(), weight_EFT_3.std(), weight_EFT_4.std(), weight_EFT_5.std()]

y_values = y_values/weight_SM.mean().values[0] 
y_values_err = y_values_err/weight_SM.mean().values[0] 

print 'Variables setup'

#trying a NN implementation

print('Plotting stuff')
 
# summarize history for accuracy                                                                                                                                                                           
fig_acc = plt.figure(figsize=(12,8))
plt.hist(weight_SM.values,color="blue",bins=10,range=(0,5))
plt.hist(weight_EFT_1.values,color="red",bins=10,range=(0,5))
plt.hist(weight_EFT_1.values,color="orange",bins=10,range=(0,5))
plt.hist(weight_EFT_1.values,color="yellow",bins=10,range=(0,5))
plt.hist(weight_EFT_1.values,color="brown",bins=10,range=(0,5))
plt.title('EFT weights', fontsize=20)
plt.ylabel('Event weight', fontsize=20)
plt.xlabel('Values', fontsize=20)
plt.legend(['SM', '-2','-1','+1','+2'], loc='upper left')
fig_acc.savefig('plot.png')
#plt.show()

x_values = np.asarray(x_values)
y_values = np.asarray(y_values)

# fit the transformed features to Linear Regression
poly_model = LinearRegression()
poly_model.fit(x_values.reshape(-1, 1), y_values.reshape(-1, 1))

print('Score 1-degree '+str(poly_model.score(x_values.reshape(-1, 1), y_values.reshape(-1, 1))))

regression_line = poly_model.coef_ * x_values.reshape(-1, 1) + poly_model.intercept_

print('Simple linear regression')
print(str(poly_model.coef_)+ ' ' + str(poly_model.intercept_))

from sklearn.preprocessing import PolynomialFeatures
poly_deg_2 = PolynomialFeatures(2)
x_degree_2 = poly_deg_2.fit_transform(x_values.reshape(-1, 1))

print('Simple degree-2 polynomial regression')
poly_model.fit(x_degree_2, y_values.reshape(-1, 1))

print('Score 2-degree '+str(poly_model.score(x_degree_2, y_values.reshape(-1, 1))))

print(str(poly_model.coef_)+ ' ' + str(poly_model.intercept_))

regression_line_2 = poly_model.coef_[0,1] * x_degree_2[:,1] + poly_model.coef_[0,2] * x_degree_2[:,2] + poly_model.intercept_

fig_fit = plt.figure(figsize=(12,8))
plt.title('EFT weights', fontsize=20)
plt.ylabel('Event weights', fontsize=20)
plt.xlabel('Values', fontsize=20)
axes = plt.gca()
axes.set_xlim([-5.5,5.5])
axes.set_ylim([0.1,4.5])
plt.errorbar(x_values, y_values, y_values_err, linestyle='None', marker='^')
plt.plot(x_values, regression_line, color="red", label='regression line')
plt.plot(x_values, regression_line_2, color="orange", label='regression line degree 2')
plt.legend(loc="upper right")
fig_fit.savefig('fit_1D_variations.png')

###### 2D-variations and fits

weight_EFT_min5_0 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_min5p0_ctGI_0p0",weight_dataframe_names)])
weight_EFT_min2p5_0 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_min2p5_ctGI_0p0",weight_dataframe_names)])

weight_EFT_2p5_0 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_2p5_ctGI_0p0",weight_dataframe_names)])
weight_EFT_5_0 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_5p0_ctGI_0p0",weight_dataframe_names)])

weight_EFT_min5_min5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_min5p0_ctGI_min5p0",weight_dataframe_names)])
weight_EFT_min2p5_min5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_min2p5_ctGI_min5p0",weight_dataframe_names)])
weight_EFT_0_min5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_0p0_ctGI_min5p0",weight_dataframe_names)])

weight_EFT_2p5_min5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_2p5_ctGI_min5p0",weight_dataframe_names)])
weight_EFT_5_min5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_5p0_ctGI_min5p0",weight_dataframe_names)])

weight_EFT_min5_min2p5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_min5p0_ctGI_min2p5",weight_dataframe_names)])
weight_EFT_min2p5_min2p5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_min2p5_ctGI_min2p5",weight_dataframe_names)])
weight_EFT_0_min2p5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_0p0_ctGI_min2p5",weight_dataframe_names)])

weight_EFT_2p5_min2p5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_2p5_ctGI_min2p5",weight_dataframe_names)])
weight_EFT_5_min2p5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_5p0_ctGI_min2p5",weight_dataframe_names)])

weight_EFT_min5_5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_min5p0_ctGI_5p0",weight_dataframe_names)])
weight_EFT_min2p5_5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_min2p5_ctGI_5p0",weight_dataframe_names)])
weight_EFT_0_5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_0p0_ctGI_5p0",weight_dataframe_names)])

weight_EFT_2p5_5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_2p5_ctGI_2p5",weight_dataframe_names)])
weight_EFT_5_5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_5p0_ctGI_5p0",weight_dataframe_names)])

weight_EFT_min5_2p5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_min5p0_ctGI_2p5",weight_dataframe_names)])
weight_EFT_min2p5_2p5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_min2p5_ctGI_2p5",weight_dataframe_names)])
weight_EFT_0_2p5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_0p0_ctGI_2p5",weight_dataframe_names)])

weight_EFT_2p5_2p5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_2p5_ctGI_2p5",weight_dataframe_names)])
weight_EFT_5_2p5 = weight_dataframe.applymap(lambda x: x[get_index("rwgt_cQq8_0p0_ctq8_0p0_ctG_5p0_ctGI_2p5",weight_dataframe_names)])

x_values = [-5, -2.5, 0, 2.5, 5]
y_values = [-5, -2.5, 0, 2.5, 5]

X, Y = np.meshgrid(x_values, y_values)

z_values = [[weight_EFT_min5_min5.mean(), weight_EFT_min2p5_min5.mean(), weight_EFT_0_min5.mean(), weight_EFT_2p5_min5.mean(), weight_EFT_5_min5.mean()],
            [weight_EFT_min5_min2p5.mean(), weight_EFT_min2p5_min2p5.mean(), weight_EFT_0_min2p5.mean(), weight_EFT_2p5_min2p5.mean(), weight_EFT_5_min2p5.mean()],
            [weight_EFT_min5_0.mean(), weight_EFT_min2p5_0.mean(), weight_SM.mean(), weight_EFT_2p5_0.mean(), weight_EFT_5_0.mean()],
            [weight_EFT_min5_2p5.mean(), weight_EFT_min2p5_2p5.mean(), weight_EFT_0_2p5.mean(), weight_EFT_2p5_2p5.mean(), weight_EFT_5_2p5.mean()],
            [weight_EFT_min5_5.mean(), weight_EFT_min2p5_5.mean(), weight_EFT_0_5.mean(), weight_EFT_2p5_5.mean(), weight_EFT_5_5.mean()]]

z_values_err = [[weight_EFT_min5_min5.std(), weight_EFT_min2p5_min5.std(), weight_EFT_0_min5.std(), weight_EFT_2p5_min5.std(), weight_EFT_5_min5.std()],
            [weight_EFT_min5_min2p5.std(), weight_EFT_min2p5_min2p5.std(), weight_EFT_0_min2p5.std(), weight_EFT_2p5_min2p5.std(), weight_EFT_5_min2p5.std()],
            [weight_EFT_min5_0.std(), weight_EFT_min2p5_0.std(), weight_SM.std(), weight_EFT_2p5_0.std(), weight_EFT_5_0.std()],
            [weight_EFT_min5_2p5.std(), weight_EFT_min2p5_2p5.std(), weight_EFT_0_2p5.std(), weight_EFT_2p5_2p5.std(), weight_EFT_5_2p5.std()],
            [weight_EFT_min5_5.std(), weight_EFT_min2p5_5.std(), weight_EFT_0_5.std(), weight_EFT_2p5_5.std(), weight_EFT_5_5.std()]]

z_values = z_values/weight_SM.mean().values[0]
z_values_err = z_values_err/weight_SM.mean().values[0]

z_values_plot = z_values.reshape(5,5)
z_values_fit = z_values.reshape(-1,1)

print(z_values_plot)

# fit the transformed features to Linear Regression                                                                                                                                                                                           
cartesian_product = [[x0, y0] for x0 in x_values for y0 in y_values]
cartesian_product = np.asarray(cartesian_product)

poly_model_2D = LinearRegression()
poly_model_2D.fit(cartesian_product, z_values_fit)

print('Score 1-degree 2D '+str(poly_model_2D.score(cartesian_product, z_values_fit)))

regression_line_2D = poly_model_2D.coef_ * cartesian_product + poly_model_2D.intercept_

print('Simple linear regression for 2D chains')
print(str(poly_model_2D.coef_)+ ' ' + str(poly_model_2D.intercept_))

poly_deg_2_2D = PolynomialFeatures(2)
x_degree_2_2D = poly_deg_2_2D.fit_transform(cartesian_product)

print('Simple degree-2 polynomial regression')
poly_model_2D.fit(x_degree_2_2D, z_values_fit)

print('Score 2-degree 2D '+str(poly_model_2D.score(x_degree_2_2D, z_values_fit)))

print(str(poly_model_2D.coef_)+ ' ' + str(poly_model_2D.intercept_))

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig_fit_2D = plt.figure(figsize=(12,8))
plt.title('EFT weights', fontsize=20)
plt.ylabel('Event weights', fontsize=20)
plt.xlabel('Values', fontsize=20)

ax = fig_fit_2D.gca(projection='3d')
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Plot the surface.
surf = ax.plot_surface(X, Y, z_values_plot, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#surf = ax.hist2d(X, Y, z_values_plot)#, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# Make data.

# Customize the z axis.
ax.set_zlim(0.5, 3.0)
ax.set_xlim([-5.5,5.5])
ax.set_ylim([-5.5,5.5])

fig_fit_2D.colorbar(surf)#, shrink=0.5, aspect=5)

fig_fit_2D.savefig('fit_2D_variations.png')







