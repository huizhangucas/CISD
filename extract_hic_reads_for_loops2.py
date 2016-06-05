#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
file_RAWexpected=sys.argv[1]  # input eg. "chr1_5kb.RAWexpected"
file_RWWobserved_dict=sys.argv[2]  # input eg. chr1_5kb.RAWobserved.dict  must be observed reads file in dict format such as: anchor1_anchor2 13
file_loop=sys.argv[3]  # input
file_loop_with_hic_reads=sys.argv[4] # output

fr1=open(file_RAWexpected,'r')
expected=[float(i.strip()) for i in fr1]
fr2=open(file_RWWobserved_dict,'r')
observed=[i.strip() for i in fr2]

#dict0=dict([[i[0]+'_'+i[1],i[2]] for i  in [i.split('\t') for i in observed]])
dict0=dict([[i[0],i[1]] for i  in [i.split('\t') for i in observed]])
import itertools

def to5k(loop_str):
        line=loop_str.split('\t')
        center1=(int(line[1])+int(line[2]))//10000*5000
        center2=(int(line[4])+int(line[5]))//10000*5000
        return([key0 for key0 in itertools.product([center1-5000,center1,center1+5000],[center2-5000,center2,center2+5000])])

def getvalues(keys):
	values=[]
	for i in range(len(keys)):
		key0=str(keys[i][0])+'_'+str(keys[i][1])
		if dict0.has_key(key0):
			value0=float(dict0[key0])
		else:
			value0=0
		expected0=float(expected[(keys[i][1]-keys[i][0])/5000])
		value_normalized= value0/expected0
		values.append(value_normalized)
	value_mean=sum(values)/len(values)
	return(value_mean)

frloop=open(file_loop,'r')
loops=[loop.strip() for loop in frloop]
fw=open(file_loop_with_hic_reads,'w')
for loop in loops:
        keys=to5k(loop)
        value_mean=getvalues(keys)
        fw.write(loop+"\t"+str(value_mean)+'\n')
fw.close()


