#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import sys

pred_sites = sys.argv[1] # must be 3 columns
output_randomloop = sys.argv[2] # output is 6 columns
fr=open(pred_sites,'r')
sites=[line.strip() for line in fr]

fw=open(output_randomloop,'w')
random_pairs=["\t".join(i) for i in itertools.combinations(sites,2)]
fw.writelines([loop+'\n' for loop in random_pairs])

