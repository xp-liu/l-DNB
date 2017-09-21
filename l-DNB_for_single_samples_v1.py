import os
import scipy.stats as stat
import math
import numpy as np
import random
import time
import multiprocessing


def ssn_score(deta,pcc,nn):
    if pcc==1:
        pcc=0.99999999
    if pcc==-1:
        pcc=-0.99999999
    z=deta/((1-pcc*pcc)/(nn-1))
    return z
    

def parallel_procedure(stage,normal,disease,title,ref,sd_mean,j,refnum,pvalue):
    begin=time.strftime("%Y-%m-%d:%H-%M-%S",time.localtime(time.time()))


    print("Stage: ",stage," Sample: ",j+1)

    network={}
    ssn={}
    for p in ref.keys():
        tmp=[]
        t=p.split()
        r1=ref[p]
        r2=stat.pearsonr(normal[t[0]]+[disease[t[0]][j]],normal[t[1]]+[disease[t[1]][j]])[0]
        r=r2-r1
        z = ssn_score(r,r1,refnum)
        p_value=1-stat.norm.cdf(abs(z))
        if p_value < pvalue:
            r=r if r>0 else -r
            ssn[p]=r
            ssn[t[1]+"\t"+t[0]]=r
            
            if t[0] not in network.keys():
                network[t[0]]=[]
            network[t[0]].append(t[1])

            if t[1] not in network.keys():
                network[t[1]]=[]
            network[t[1]].append(t[0])


    ci={}
    for p in network.keys():
        if len(network[p])<3:
            continue
        
        sd=abs(disease[p][j]-sd_mean[p][1])/sd_mean[p][0]
        pcc_in=0
        pcc_out=0
        count=0
        for q in network[p]:
            sd+=abs(disease[q][j]-sd_mean[q][1])/sd_mean[q][0]
            pcc_in+=ssn[p+"\t"+q]

            for m in network[q]:
                if m!=p:
                    pcc_out+=ssn[q+"\t"+m]
                    count+=1
        sd/=len(network[p])+1
        pcc_in/=len(network[p])
        if count==0:
            continue
        pcc_out/=count
        if pcc_out==0:
            continue
        ci[p]=[sd*pcc_in/pcc_out,sd,pcc_in,pcc_out]

    ci=sorted(ci.items(),key=lambda d:d[1][0],reverse=True)


    fw=open("Max_score_module in %s for sample %d" % (stage,j+1),"w")
    for k in range(len(ci)):
        fw.write(ci[k][0]+"\t"+str(ci[k][1][0])+"\t"+str(ci[k][1][1])+"\t"+str(ci[k][1][2])+"\t"+str(ci[k][1][3])+"\n")
    fw.close()


    print("Begin time is "+begin)
    print("End time is "+time.strftime("%Y-%m-%d:%H-%M-%S",time.localtime(time.time())))


if __name__=="__main__":
   
    pvalue=0.05  #p-value threshold is set
    
    refnum=0  
    normal={}
    f=open("Reference_samples.txt")
    flag=0
    for p in f:
        flag+=1
        if flag==1:
            t=p.split()
            refnum=len(t)-1
            continue
        t=p.split()
        normal[t[0]]=[float(t[i]) for i in range(1,len(t))]
    f.close()

    sd_mean={}
    for key in normal.keys():
        sd_mean[key]=[np.std(normal[key]),np.mean(normal[key])]

    f=open("reference_network.txt")
    network={}
    ref={}
    for p in f:
        t=p.split()
        ref[t[0]+"\t"+t[1]]=float(t[2])

    f.close()
            
    stages=["Stage IA.txt","Stage IB.txt","Stage IIA.txt","Stage IIB.txt","Stage IIIA.txt","Stage IIIB.txt","Stage IV.txt"]

    pool=multiprocessing.Pool(3)

    for stage in stages:
        file=stage
        f=open(file)
        disease={}
        title=[]
        flag=0
        for p in f:
            flag+=1
            t=p.split()
            if flag==1:
                title=[str(k) for k in range(1,len(t))]
                continue
            disease[t[0]]=[float(t[k]) for k in range(1,len(t))]

        f.close()

        for j in range(len(title)):
            pool.apply_async(parallel_procedure,(stage,normal,disease,title,ref,sd_mean,j,refnum,pvalue,))
    pool.close()
    pool.join()

        
