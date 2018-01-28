#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
file_RAWkr=sys.argv[1]  # input eg. "chr1_5kb.KRnorm.nonan"
file_RAWexpected=sys.argv[2] # input eg. chr1_5kb.RAWexpected  
file_RAWobserved=sys.argv[3]  # input eg. chr1_5kb.RAWobserved  
file_output_normalized=sys.argv[4]  # output normalized in dict format

                               
                               
fr0=open(file_RAWkr,'r')
kr=[float(i.strip()) for i in fr0]
fr1=open(file_RAWexpected,'r')
expected=[float(i.strip()) for i in fr1]

fr2=open(file_RAWobserved,'r')
fw=open(file_output_normalized,"w")

for line in fr2:
    (anchor1,anchor2,value0)=line.strip().split("\t")
    idx1=int(anchor1)/5000
    idx2=int(anchor2)/5000
    value1=float(value0)/(kr[idx1]*kr[idx2])
    expected0=float(expected[abs(idx2-idx1)])
    value2 = value1/expected0
    fw.write(anchor1+"_"+anchor2+"\t"+str(round(value2,2))+"\n")
#    fw.write(anchor1+"_"+anchor2+"\t"+line.strip()+"\t"+str(idx1)+"\t"+str(idx2)+"\t"+str(value1)+"\t"+str(expected0)+"\t"+str(round(value2,2))+"\n")

