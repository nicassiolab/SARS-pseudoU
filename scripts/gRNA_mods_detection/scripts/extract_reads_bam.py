#!/usr/bin/python

import pysam
import sys, getopt
from collections import Counter
import re
import numpy as np
import glob
import os


inputfile = sys.argv[1];
outputfile = sys.argv[2]
grepping_list = sys.argv[3]

print(inputfile)
print(outputfile)
print(grepping_list)
pos_cov = [ 0 for _ in range(0, 30000)]

# extract reads from a list

listnames=glob.glob(grepping_list)

for q in listnames:
       qnames = np.genfromtxt(q,dtype='str')
       bamfile =pysam.AlignmentFile(inputfile, "rb")
       out = pysam.AlignmentFile(outputfile, "wb",template=bamfile)
       for x in bamfile.fetch(until_eof=True):
               if x.query_name in qnames:
                       out.write(x)
       bamfile.close()
       out.close()



