from root_numpy import root2array, tree2array, list_structures
from root_numpy import testdata
import numpy as np
import pandas
import matplotlib.pyplot as plt

filename = 'GenLevelOutput.root'

# Convert a TTree in a ROOT file into a NumPy structured array + pandas Dataframe, with which you can create an array to be given to Keras
branches = ["EFTWeights"]
# one needs probably a flat tree, to access information

weight_dataframe = pandas.DataFrame(root2array(filename, treename='boosted/events',branches=branches))
print(weight_dataframe)

weight_SM = weight_dataframe.applymap(lambda x: x[0])
weight_EFT_1 = weight_dataframe.applymap(lambda x: x[1])
weight_EFT_2 = weight_dataframe.applymap(lambda x: x[2])
weight_EFT_3 = weight_dataframe.applymap(lambda x: x[3])
weight_EFT_4 = weight_dataframe.applymap(lambda x: x[4])

x_values = [0, -2, -1, 1, 2]

y_values = [weight_SM.mean(), weight_EFT_1.mean(), weight_EFT_2.mean(), weight_EFT_3.mean(), weight_EFT_4.mean()]

y_values_err = [weight_SM.std(), weight_EFT_1.std(), weight_EFT_2.std(), weight_EFT_3.std(), weight_EFT_4.std()]

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

fig_fit = plt.figure(figsize=(12,8))
plt.title('EFT weights', fontsize=20)
plt.ylabel('Event weights', fontsize=20)
plt.xlabel('Values', fontsize=20)
axes = plt.gca()
axes.set_xlim([-2.5,2.5])
axes.set_ylim([0.5,1.5])
plt.errorbar(x_values, y_values, y_values_err, linestyle='None', marker='^')
fig_fit.savefig('fit.png')



