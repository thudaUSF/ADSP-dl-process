import traceback
import pandas as pd
import numpy as np
import os
import sys
#from lifelines.statistics import logrank_test
#from lifelines.utils import median_survival_times #I couldn't get this to work
#from lifelines import KaplanMeierFitter
#kmf = KaplanMeierFitter()
#from lifelines import CoxPHFitter
import scipy.stats as stats
#import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')

##### recoveries data ##### 
path_recoveries = "/Users/thuda/Desktop/ADphyschem/AD_pivot_pchem_with_clinical_data.csv"
df = pd.read_csv(path_recoveries)

###### clinical data ####### (NOTE: FOR BOTTOM50/TOP50, WE CAN COMPARE BRAAK WITHIN AD ONLY, THEN COMPARE AD VS CONTROL)
survival = "/Users/thuda/Desktop/ADphyschem/AD_pivot_pchem_with_clinical_data.csv"
data_surv = pd.read_csv(survival, sep=",")
data_clinical = data_surv[["SUBJID","BODY_SITE",
                           'APOE',
                           "Race",'AD']] 

###### SURVIVAL DATA #######
path_survival = "/Users/thuda/Desktop/ADphyschem/AD_pivot_pchem_with_clinical_data.csv"
data_surv = pd.read_csv(path_survival, sep=",")
data_survival = data_surv[["SUBJID","Braak"]] #REMOVED AGE_BASELINE BECAUSE CAUSES AROUND 2000 BLANKS, then not included in processing

##### physchem data ##### 
path_pchem = "/Users/thuda/Desktop/ADphyschem/AD_pchem.csv"
df_pchem = pd.read_csv(path_pchem)

AD_patients = data_clinical.loc[data_clinical["AD"] == 1]["SUBJID"].drop_duplicates().tolist() #NOTE: TO FILTER ALL PATIENTS WITH VALUE OF 1;CHECK IF WORKING

def sort_and_get_hi_lo(table,sampletype,receptor,proper):
    try:
        patients = table.loc[(sampletype,receptor), proper].sort_values(ascending=True)
        length = len(patients)
        middle_index = length//2
        bottom_half = patients[:middle_index].index.tolist()
        top_half = patients[middle_index:].index.tolist()
    except:
        bottom_half = []
        top_half = []
    
    return bottom_half, top_half

#NOTE: REPLACE WITH T-TEST
def get_logrank(df1,df2):
    results = stats.ttest_ind(df1['Braak'], df2['Braak'])
    return results.pvalue


###### COMPARE BETTER SURVIVOR FN #### 
def get_better_survivor(mean_braak_group1,mean_braak_group2):
    if mean_braak_group1 > mean_braak_group2:
        better_survivor = "Top50%"
    elif mean_braak_group1 < mean_braak_group2:
        better_survivor = "Bottom50%"
    elif mean_braak_group1 == mean_braak_group2:
        better_survivor = "Equal Survival"
    else:
        better_survivor = ""
    return better_survivor


###### GENERATE RESULTS OF ANALYSIS IN LIST FORMAT ##########    
def get_data(sampletype,receptor,proper,df1,df2):
    
    N1 = len(df1["SUBJID"].to_list())
    N2 = len(df2["SUBJID"].to_list())
    
    if N1 and N2 >= 2:
        
        Median_braak_group1 = df1["Braak"].median()
        Median_braak_group2 = df2["Braak"].median()
        Mean_braak_group1 = df1["Braak"].mean()
        Mean_braak_group2 = df2["Braak"].mean()
        pvalue = get_logrank(df1,df2)
        
        higher_survival = get_better_survivor(Mean_braak_group1,Mean_braak_group2)
        
        data = [sampletype,receptor,proper,
                N1, N2, Median_braak_group1, Median_braak_group2,
                Mean_braak_group1,Mean_braak_group2,pvalue,higher_survival]
        
    else:
        data = [sampletype,receptor,proper,N1, N2,"", "","",""]
        
    return data

# get patients ids for Bcell ALL
sampletypes = df_pchem["BODY_SITE"].drop_duplicates().tolist()
receptors = ["TRA","TRB","TRD","TRG","IGH","IGK","IGL"]
properties = df_pchem.columns[3:] ###NOTE: DOES THIS DO ALL COLUMNS?
#filetype = "wxs"

results = []

for sampletype in sampletypes:
    for receptor in receptors:

        #### filter physchem df based on whether they have AD or not #### 
        df_filtered_1 = df_pchem[df_pchem["SUBJID"].isin(AD_patients)]
        
        #### filter reads based on sampletype and receptor type #### 
        df_filtered_2 = df_filtered_1[(df_filtered_1["BODY_SITE"] == sampletype) & (df_filtered_1["Receptor"] == receptor)]

        #### generate corresponding pivot #####
        ind = ["BODY_SITE","Receptor","SUBJID"]
        table = pd.pivot_table(df_filtered_2, index = ind, values = df_filtered_2.columns, aggfunc = np.mean)
        for proper in properties:
            bottom_half, top_half = sort_and_get_hi_lo(table,sampletype,receptor,proper)
            #print(bottom_half[:5])
            df_survival_top50 = data_survival[data_survival.SUBJID.isin(top_half)].dropna().drop_duplicates()
            df_survival_bottom50 = data_survival[data_survival.SUBJID.isin(bottom_half)].dropna().drop_duplicates()
            if proper == "fraction_negative" and sampletype == "Brain" and receptor== "TRA":
                #patients = patients.merge(df_pchem[['SUBJID','fraction_negative']],on='SUBJID', how='left') #cc = control, case
                df_survival_top50.to_csv("/Users/thuda/Desktop/ADphyschem/trabrainfraction_negativepchemt50.csv",index=False)
                df_survival_bottom50.to_csv("/Users/thuda/Desktop/ADphyschem/trabrainfraction_negativepchemb50.csv",index=False)
            
            
            
            
            data = get_data(sampletype,receptor,proper,df_survival_top50,df_survival_bottom50)
            results.append(data)



cols_results = ["Sampletype","Receptor","Property",
                "N - top50%","N - bottom50%",
                "Median braak - top50%", "Median braak - bottom50%",
                "Mean braak - top50%", "Mean braak - bottom50%",
                "P-value","Worse Braak"]
df_results = pd.DataFrame(results)
df_results.columns = cols_results 
df_results = df_results.dropna().sort_values(by=["P-value"]) #NOT WORKING
results_path = "/Users/thuda/Desktop/ADphyschem/results.xlsx" #GET DATE AND TIME AS STRING
df_results.to_excel(results_path, index = False)
