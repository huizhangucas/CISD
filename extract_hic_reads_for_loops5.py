#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
file_RWWobserved_dict=sys.argv[1]  # input eg. chr1_5kb.RAWobserved.dict  must be observed reads file in dict format such as: anchor1_anchor2 13
file_loop=sys.argv[2]  # input
file_loop_with_hic_reads=sys.argv[3] # output

                                 

fr2=open(file_RWWobserved_dict,'r')
observed=[i.strip() for i in fr2]

#dict0=dict([[i[0]+'_'+i[1],i[2]] for i  in [i.split('\t') for i in observed]])
dict0=dict([[i[0],i[1]] for i  in [i.split('\t') for i in observed]])
import itertools

def to5k(loop_str):
        line=loop_str.split('\t')
        center1=(int(line[1])+int(line[2]))//10000*5000
        center2=(int(line[4])+int(line[5]))//10000*5000
        return([key0 for key0 in itertools.product([center1-10000, center1-5000, center1, center1+5000, center1+10000],[center2-10000, center2-5000, center2, center2+5000, center2+10000])])

def getvalues(keys):
	values=[]
	for i in range(len(keys)):
		key0=str(keys[i][0])+'_'+str(keys[i][1])
		if dict0.has_key(key0):
			value0=float(dict0[key0])
		else:
			value0=0
		values.append(value0)
	idx1=[12]
	idx2=[6,7,8,11,12,13,16,17,18]
	idx3=range(25)
	value_level1=sum([values[i] for i in idx1])
	value_level2=sum([values[i] for i in idx2])
	value_level3=sum([values[i] for i in idx3])
	value_21 = value_level2 - value_level1
	value_32 = value_level3 - value_level2
	value_str="\t".join([str(i) for i in [value_level1,value_level2,value_level3,value_21,value_32]])
	return(value_str)

frloop=open(file_loop,'r')
loops=[loop.strip() for loop in frloop]
fw=open(file_loop_with_hic_reads,'w')
for loop in loops:
        keys=to5k(loop)
        value_mean=getvalues(keys)
        fw.write(loop+"\t"+str(value_mean)+'\n')
fw.close()


