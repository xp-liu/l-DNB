import os
import scipy.stats as stat
import math
import numpy as np
import random
import time




begin=time.strftime("%Y-%m-%d:%H-%M-%S",time.localtime(time.time()))

normal={}
f=open("Reference_samples.txt")
flag=0
for p in f:
    flag+=1
    if flag==1:
        continue
    t=p.split()
    normal[t[0]]=[float(t[i]) for i in range(1,len(t))]
f.close()


network={}
keys=list(normal.keys())
n=len(keys)

fw=open("reference_network.txt","w")
for i in range(n-1):
    print(i)
    for j in range(i+1,n):
        r=stat.pearsonr(normal[keys[i]],normal[keys[j]])
        # The threshold of P-value need be set in here for Pearson Correlation Coefficient
        if r[1] < 0.01 / (20501*20501): 
            fw.write(keys[i]+"\t"+keys[j]+"\t"+str(r[0])+"\n")
            
fw.close()

        
print("Begin time is "+begin)
print("End time is "+time.strftime("%Y-%m-%d:%H-%M-%S",time.localtime(time.time())))

