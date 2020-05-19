from root_numpy import root2array, tree2array, list_structures
from root_numpy import testdata
import numpy as np
import pandas
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import math

def get_index(name_string, df_names):

    for i in range(0,len(df_names.values[0][0])):
        if(df_names.values[0][0][i] == name_string):
            index = i
            break

    return index

filename = 'GenLevelOutput.root'

# Convert a TTree in a ROOT file into a NumPy structured array + pandas Dataframe, with which you can create an array to be given to Keras
branches = ["EFTWeights"]
branches_obs = ["DiTopMass"]
branches_names = ["EFTWeightsNames"]
# one needs probably a flat tree, to access information

weight_dataframe = pandas.DataFrame(root2array(filename, treename='boosted/events',branches=branches))
weight_dataframe_obs = pandas.DataFrame(root2array(filename, treename='boosted/events',branches=branches_obs))
weight_dataframe_names = pandas.DataFrame(root2array(filename, treename='boosted/events',branches=branches_names))

parameter_to_scan = "ctGI"#"ctG"

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

weight_obs = weight_dataframe_obs.applymap(lambda x: x[0])

print('Plotting stuff')

# summarize history for accuracy                                                                                                                                                                           
#fig_acc = plt.figure(figsize=(12,8))
#plt.hist(weight_dataframe_obs["DiTopMass"],weights=weight_SM.values,color="blue",bins=10,range=(0,5))
fig_acc, (ax1, ax2) = plt.subplots(nrows=2)

#bins = [0,100,200,300,400,500,600,700,800,900,1000,1200,1500,2000,2500,3000]

#ns_1 = plt.bar(bincenters, y, width=width, color='r', yerr=menStd)
#print(ns_1)

ns_1, bins_1, patches_1 = ax1.hist(weight_obs.values,color="red",bins=10, alpha = 0.1, log=True, label='SM', normed=True)
ns_2, bins_2, patches_2 = ax1.hist(weight_obs.values,weights=(1+weight_EFT_5.values),color="orange",bins=10, alpha = 0.1, log=True, label='+5', normed=True)




ratio = []
ratio_err = []
for i in range(0,len(ns_2)):
    if ns_2[i] !=0 :
        ratio.append(ns_2[i] / ns_1[i])
        ratio_err.append(math.sqrt((ns_2[i])/(ns_1[i]*len(weight_obs.values)*len(weight_obs.values)*ns_1[i])+ns_2[i]*ns_2[i]/(ns_1[i]*len(weight_obs.values)*len(weight_obs.values)*ns_1[i]*ns_1[i])))
    else:
        ratio.append(1)
        ratio_err.append(1)

ax2.bar(bins_1[:-1],     # this is what makes it comparable
        ratio,
        alpha=0.4,
        width = 200,
        yerr = ratio_err
)

ax1.set_ylabel('A.U.')
ax2.set_ylabel('Ratio (+5/SM)')

plt.xlabel('Ditop quark mass (GeV)', fontsize=20)
ax1.legend(loc='upper right')
fig_acc.savefig('ditopmass-'+str(parameter_to_scan)+'.png')


def chi2_square(new_predictions, SM, SM_error):

    chi_2 = 0

    for i in range(0,len(new_predictions)):
        
        chi_2 += (new_predictions[i]-SM[i])*(new_predictions[i]-SM[i]) / SM_error[i]

    return chi_2

#plot the chis for different predictions and fit as a function of the values

def plot_chi2_square(chi2_vector, EFT_vector):

    fig_acc = plt.figure(figsize=(12,8))                                                                                                                                                                                                   
    plt.plot(EFT_vector, chi2_vector,color="blue",bins=10,range=(0,5))     

    #fit and plot
    fit = LinearRegression()
    


    return fit


#extract the values of the parameters for which the chi2 variation is +/- 6.63 for 95%, +/- 2.71 for 90% and +/- 1 for 68% 
#https://www.astro.ubc.ca/people/jvw/ASTROSTATS/lectures/PL_10_2013_www.pdf


