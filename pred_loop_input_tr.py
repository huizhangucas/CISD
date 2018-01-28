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


import sys
trfile=sys.argv[1] #training set
file2=sys.argv[2]  #candidate with hic,fft,distance
fileout=sys.argv[3]  #output predicted loop, 7 fields, the 7th field: 1 represent loop, and 0 represent not loop


fr=open(trfile,'r')
n=len(fr.readlines())
returnMat1=zeros((n,106))
lables=zeros(n)
fr=open(trfile,'r')
index = 0
for line in fr.readlines():
    line = line.strip()
    listFromLine = line.split('\t')
    returnMat1[index,:] = listFromLine[6:112]
    lables[index]=listFromLine[112]
    index += 1




#clf = GradientBoostingClassifier(n_estimators=500)
#clf = AdaBoostClassifier(n_estimators=500)
#clf = RandomForestClassifier(n_estimators=300,max_depth=80)
#clf = ExtraTreesClassifier(n_estimators=300,max_depth=50)

#clf.fit(returnMat1,lables)






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


clf = GradientBoostingClassifier(n_estimators=500)
clf.fit(returnMat1,lables)    
result=clf.predict(returnMat2)
fileout1=fileout+".gbdt"
fw=open(fileout1,"w")
for i in range(n):
    fw.write("\t".join(position[i])+"\t"+str(result[i])+"\n")
    
clf = AdaBoostClassifier(n_estimators=500)
clf.fit(returnMat1,lables)    
result=clf.predict(returnMat2)
fileout1=fileout+".ada"
fw=open(fileout1,"w")
for i in range(n):
    fw.write("\t".join(position[i])+"\t"+str(result[i])+"\n")

clf = RandomForestClassifier(n_estimators=300,max_depth=80)
clf.fit(returnMat1,lables)    
result=clf.predict(returnMat2)
fileout1=fileout+".rf"
fw=open(fileout1,"w")
for i in range(n):
    fw.write("\t".join(position[i])+"\t"+str(result[i])+"\n")

clf = ExtraTreesClassifier(n_estimators=300,max_depth=50)
clf.fit(returnMat1,lables)    
result=clf.predict(returnMat2)
fileout1=fileout+".ext"
fw=open(fileout1,"w")
for i in range(n):
    fw.write("\t".join(position[i])+"\t"+str(result[i])+"\n")

   
  
    