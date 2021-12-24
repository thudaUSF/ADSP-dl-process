import traceback
import pandas as pd
import numpy as np
import csv
import os
import sys
#from lifelines.statistics import logrank_test
#from lifelines.utils import median_survival_times #I couldn't get this to work
#from lifelines import KaplanMeierFitter
#kmf = KaplanMeierFitter()
#from lifelines import CoxPHFitter
import scipy.stats as stats
import rpy2.robjects as R
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
#import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')

path = ""
##### recoveries data ##### 
path_recoveries = path + "AD_pivot_pchem_with_clinical_data.csv"
df = pd.read_csv(path_recoveries)

###### clinical data ####### (NOTE: FOR BOTTOM50/TOP50, WE CAN COMPARE BRAAK WITHIN AD ONLY, THEN COMPARE AD VS CONTROL)
survival = path + "AD_pivot_pchem_with_clinical_data.csv"
data_surv = pd.read_csv(survival, sep=",")
data_clinical = data_surv[["SUBJID","BODY_SITE",
                           'APOE',
                           "Race",'AD']] 
data_clinical = data_clinical.drop_duplicates(subset="SUBJID") 

###### SURVIVAL DATA #######
path_survival = path + "AD_pivot_pchem_with_clinical_data.csv"
data_surv = pd.read_csv(path_survival, sep=",")
data_surv = data_surv.drop_duplicates(subset="SUBJID") 
#data_survival = data_surv[["SUBJID","Braak","BODY_SITE"]] #REMOVED AGE_BASELINE BECAUSE CAUSES AROUND 2000 BLANKS, then not included in processing

##### physchem data ##### 
path_pchem = path + "AD_pchem.csv"
df_pchem = pd.read_csv(path_pchem)
AD_patients = data_clinical.loc[data_clinical["AD"] == 1]["SUBJID"].drop_duplicates().tolist() #NOTE: TO FILTER ALL PATIENTS WITH VALUE OF 1;CHECK IF WORKING

def runKW(pdf, proper):
    pdf = pdf.rename(columns={proper: 'characteristic'})
    groupsizes = pdf.groupby(["Braak"]).size().reset_index(name='Number of people')
    finallist = groupsizes["Number of people"].tolist()
    for i,braakn in groupsizes.iterrows():
        if braakn[1] < 15:
            pdf = pdf[pdf.Braak != braakn[0]]
    with localconverter(R.default_converter + pandas2ri.converter):
        df = R.conversion.py2rpy(pdf)
    if proper=="fraction_negative":
        pdf.to_csv(path + "something.csv")
    Kruskall = R.r(r'''
        function(df) {
            options(warn=-1)
            df$Braak <- as.factor(df$Braak)
            levels(df$Braak)
            df$Braak <- ordered(df$Braak, levels = c("0","3","5","6"))

            # run kruskal-wallis for property
            pvalue <- kruskal.test(characteristic~Braak, data=df)$p.value
            }
    ''')
    posthoc = R.r(r'''
        function(df) {
            options(warn=-1)
            df$Braak <- as.factor(df$Braak)
            levels(df$Braak)
            df$Braak <- ordered(df$Braak, levels = c("0","3","5","6"))

            library(FSA)
            dunn <- dunnTest(characteristic ~ Braak,
                    data=df,
                    method="bonferroni")$res
            return(dunn)
            }
    ''')
    rresults = Kruskall(df)
    results=np.asarray(rresults)
    if(results<1):
        rposthoc = posthoc(df)
    with localconverter(R.default_converter + pandas2ri.converter):
        posthoc = R.conversion.rpy2py(rposthoc)
    posthoclist = posthoc["P.adj"]
    finallist.extend(results)
    finallist.extend(posthoclist)
    return finallist

###### GENERATE RESULTS OF ANALYSIS IN LIST FORMAT ##########    
def get_data(sampletype,receptor,proper,df): 

    list = runKW(df,proper)         
    data = [sampletype,receptor,proper]
    data.extend(list)
    return data

# get patients ids for Bcell ALL
sampletypes = df_pchem["BODY_SITE"].drop_duplicates().tolist()
receptors = ["TRA","TRB","TRD","TRG","IGH","IGK","IGL"]
properties = df_pchem.columns[3:] ###NOTE: DOES THIS DO ALL COLUMNS?
#filetype = "wxs"

bloodresults = []
brainresults = []

for sampletype in sampletypes:
    for receptor in receptors:

        #### filter physchem df based on whether they have AD or not #### 
        df_filtered_1 = data_surv[data_surv["SUBJID"].isin(AD_patients)]
        
        #### filter reads based on sampletype and receptor type #### 
        df_filtered_2 = df_filtered_1[(df_filtered_1["BODY_SITE"] == sampletype) & (df_filtered_1["Receptor"] == receptor)]
        if (df_filtered_2['SUBJID'].count() < 50):
            continue
        #### generate corresponding pivot #####
        #ind = ["BODY_SITE","Receptor","SUBJID"]
        #table = pd.pivot_table(df_filtered_2, index = ind, values = df_filtered_2.columns, aggfunc = np.mean)
        for proper in properties:
            table = df_filtered_2[["SUBJID","Braak",proper]] 
            table.to_csv(path + "someting.csv")
            table.loc[(table['Braak'] == 4), 'Braak'] = 3
            if sampletype=='Blood':
                table.loc[(table['Braak'] < 3), 'Braak'] = 0                
           
            data = get_data(sampletype,receptor,proper,table)
            if sampletype=='Blood':
                bloodresults.append(data)
            if sampletype=='Brain':
                brainresults.append(data)

Detailsblood = ["Sampletype","Receptor","Property",
                "N - Braak 0, I, II","N - Braak III, IV",
                "N - Braak V", "N - Braak VI", "KW p-value",
                "Post-hoc Comparison 0,I,II/III,IV", "Post-hoc Comparison 0,I,II/V",
                "Post-hoc Comparison 0,I,II/VI","Post-hoc Comparison III,IV/V",
                "Post-hoc Comparison III,IV/VI","Post-hoc Comparison V/VI"]
Detailsbrain = ["Sampletype","Receptor","Property",
                "N - Braak 0, I, II","N - Braak III, IV",
                "N - Braak V", "N - Braak VI", "KW p-value",
                "Post-hoc Comparison III,IV/V",
                "Post-hoc Comparison III,IV/VI","Post-hoc Comparison V/VI"]

df_bloodresults = pd.DataFrame(bloodresults)
df_brainresults = pd.DataFrame(brainresults)
df_bloodresults.columns = Detailsblood 
df_brainresults.columns = Detailsbrain
df_bloodresults = df_bloodresults.dropna().sort_values(by=["KW p-value"])
df_brainresults = df_brainresults.dropna().sort_values(by=["KW p-value"])
df_bloodresults.to_csv(path + "bloodresults.csv", index = False)
df_brainresults.to_csv(path + "brainresults.csv", index = False)
