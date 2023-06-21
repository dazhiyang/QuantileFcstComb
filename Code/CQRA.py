# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 18:15:33 2022
@email: yangguoming1995@gmail.com
@author: Dazhi Yang, Guoming Yang
"""

import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import xlsxwriter


##-----------data process----------------------------------
station = 'fpk'   ###------------modify here, if you want to get the result of the other location
file_path = 'D:\Doctor\paper\supervisor\yang_combiningforecast\github\Originial_data'+ '\\'+ station +'_PP_2019_2020.txt'  #data path
df = pd.read_table(file_path, sep='\t',)
df['Time'] = pd.to_datetime(df['Time'])   #time data type
df = df.set_index('Time')                 #set index
df2019_train = df['2019']                 #get training data
df2020_predict = df['2020']               #get forecast data

###----------------------creat model------------------------
Probabilistic=gp.Model("Probabilistic")
# define variables
N = 10    #model number
quantile = [0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99]  #quantile
Q = len(quantile)   #quantile number
Wnq = Probabilistic.addVars(N, Q, vtype=GRB.CONTINUOUS, name='weights')
Theta_positive = Probabilistic.addVars(len(df2019_train), Q, vtype=GRB.CONTINUOUS, name='positive_residual')
Theta_negative = Probabilistic.addVars(len(df2019_train), Q, vtype=GRB.CONTINUOUS, name='negative_residual')

# set objective function
Probabilistic.setObjective(sum(quantile[q] * Theta_positive[t, q] + (1 - quantile[q]) * Theta_negative[t, q] for t in range(len(df2019_train)) for q in range(Q)), GRB.MINIMIZE)
#    sum(Vtq[t, q] for t in range(len(df2019_train)) for q in range(Q)), GRB.MINIMIZE)


# add constraints
#the sum of the weights of the ten model for each quantile is equal to one.
Probabilistic.addConstrs((sum(Wnq[n, q] for n in range(N)) == 1 for q in range(Q)), name= 'weights_sum_one')
#equality constraint
for t in range(len(df2019_train)):
    for q in range(Q):
        Probabilistic.addConstr((df2019_train.iloc[t,0] == sum(Wnq[n,q] * df2019_train.iloc[t,2+q+n*Q] for n in range(N)) + Theta_positive[t,q] - Theta_negative[t,q]),'equal_con1')

# Optimize model
Probabilistic.optimize()


##--------------------get the value of Wnq, see table 3---------------
Wnq_c=np.zeros((N,Q))
for n in range(N):
    for q in range(Q):
        Wnq_c[n,q] = Wnq[n,q].x
#保留一位小数
Wnq_c0=np.zeros((N,Q))
for n in range(N):
    for q in range(Q):
        Wnq_c0[n,q] = round(Wnq_c[n,q],2)
#列名称
names = []
for q in range(Q):
    names.append('qua.'+str(int(round(quantile[q]*100, 0))))
Wnq_c0 = pd.DataFrame(Wnq_c0, columns=names)
##save data
new_file_path1 = 'D:\Doctor\paper\supervisor\yang_combiningforecast\github\Result'+ '\\'+ station + '_CQRA_Wnq' +'.txt'  #file path
Wnq_c0.to_csv(new_file_path1, header=True, sep='\t', index=False)
##-------------------get all forecasts-------------------------------
#in-sample forcasts for the year of 2019
Y_inpre = np.zeros((len(df2019_train), Q))
for q in range(Q):
    for t in range(len(df2019_train)):
        Y_inpre[t, q] = sum(Wnq_c[n,q] * df2019_train.iloc[t,2+q+n*Q] for n in range(N))
        Y_inpre[t, q] = round(Y_inpre[t, q],1)


#out-of-sample combined forecasts for the year of 2020
Y_pre = np.zeros((len(df2020_predict), Q))
for q in range(Q):
    for t in range(len(df2020_predict)):
        Y_pre[t, q] = sum(Wnq_c[n,q] * df2020_predict.iloc[t,2+q+n*Q] for n in range(N))
        Y_pre[t, q] = round(Y_pre[t, q], 1)
#simple quartile averaging
Y_sqa_2019 = np.zeros((len(df2019_train), Q))
for q in range(Q):
    for t in range(len(df2019_train)):
        Y_sqa_2019[t,q] =  round((sum(df2019_train.iloc[t,2+q+n*Q]for n in range(N)))/N,1)


for q in range(Q):
    df2019_train['SQA.'+str(int(round(quantile[q]*100, 0)))] = Y_sqa_2019[:,q]


Y_sqa_2020 = np.zeros((len(df2020_predict), Q))
for q in range(Q):
    for t in range(len(df2020_predict)):
        Y_sqa_2020[t,q] =  round((sum(df2020_predict.iloc[t,2+q+n*Q]for n in range(N)))/N,1)

for q in range(Q):
    df2020_predict['SQA.'+str(int(round(quantile[q]*100, 0)))] = Y_sqa_2020[:,q]

#incorporate Y_inpre into the dataframe df2019_train
for q in range(Q):
    df2019_train['CQRA.'+str(int(round(quantile[q]*100, 0)))] = Y_inpre[:,q]

#incorprate Y_pre into the dataframe df2020_predict
for q in range(Q):
    df2020_predict['CQRA.'+str(int(round(quantile[q]*100, 0)))] = Y_pre[:,q]

#merge data together
Final_data = pd.concat([df2019_train,df2020_predict], axis=0)

#save data
new_file_path2 = 'D:\Doctor\paper\supervisor\yang_combiningforecast\github\Result'+ '\\'+ station + '_CQRA_2019_2020' +'.txt'  #file path
Final_data.to_csv(new_file_path2, sep='\t')
