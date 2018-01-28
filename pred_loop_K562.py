#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
import numpy as np
from sklearn.datasets import load_iris
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_blobs
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
from numpy import *

import math

file="/home/zhanghui/project1/zhanghui_K562/pred_loop_with_sampled_reads/K562/m1h1/rep0/new_np/pn_new_exp.bed.withdistance"


fr=open(file,'r')
n=len(fr.readlines())
returnMat1=zeros((n,106))
fr=open(file,'r')
index = 0
for line in fr.readlines():
    line = line.strip()
    listFromLine = line.split('\t')
    returnMat1[index,:] = listFromLine[6:112] 
    index += 1

lables=zeros(10000)
lables[0:4999]=1


clf = GradientBoostingClassifier(n_estimators=500)
#clf = AdaBoostClassifier(n_estimators=500)
#clf = RandomForestClassifier(n_estimators=300,max_depth=80)
#clf = ExtraTreesClassifier(n_estimators=300,max_depth=50)

clf.fit(returnMat1,lables)



import sys
file2=sys.argv[1]  
fileout=sys.argv[2]  
    



fr=open(file2,'r')
n=len(fr.readlines())
returnMat2=zeros((n,106))
position=[]
fr=open(file2,'r')
index = 0
for line in fr.readlines():
    line = line.strip()
    listFromLine = line.split('\t')
    returnMat2[index,:] = listFromLine[6:112] 
    position.append(listFromLine[0:6])
    index += 1
    
result=clf.predict(returnMat2)
fw=open(fileout,"w")
for i in range(n):
    fw.write("\t".join(position[i])+"\t"+str(result[i])+"\n")
    
    
  
    