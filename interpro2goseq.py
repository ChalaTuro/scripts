#!/usr/bin/env python
# coding: utf-8

#this is used to extract GO ID and Protein ID from interproscan annotation
import os
import re
from collections import OrderedDict
import shutil, sys 
import itertools
from collections import defaultdict
import pandas as pd
import collections



data_dict = defaultdict(list)  # that is most important 
data_list = []
outid =  open("C:/Users/254357H/Documents/IPR/ptt_17FRG026_id_leng.txt", 'w')
outfileFinal = open("C:/Users/254357H/Documents/IPR/ptt_17FRG026_GoMAppedFinalDefaultDict_out.txt", 'w')
gseq = open("R:/CCDM-FGR-SE00542/Chala/RNASEqAnalysis2022/FinalData/ptt_17FRG026_goseq_out_17052022.csv", 'w')
with open('C:/Users/254357H/Documents/IPR/17FRG026_augustus.hints2.aa.tsv', 'r') as paf:
    for line in paf:
        if not re.match("#",line):
            line = line.replace('\n', '')
            if (re.search("GO:\d+\|GO:\d+", line)):
                data_list = list(line.split("\t",20))
                #line = list(line.split("\t",20))
                #print(line)
                protid = data_list[0] # protein id  Key
                start = data_list[6]
                stop = data_list[7]
                length = data_list[2]
                goes = data_list[13]
                goCopy = goes
                #print(protid, '\t', goes)
                goCopy = goCopy.replace('|', ',') # comented after dict written to file
                goCopy2 = goCopy.split('|') # make a copy to create a dict 
                #goCopy2 = goCopy.replace('|','\n') # make a copy to create a dict 
                #for i in goCopy2:
                    #print(protid, '\t', i)
                for i in goCopy2:
                    got = i.split(',')
                    for x in got:
                        #print(protid + "\t" + x)
                        gseq.write(protid + "\t" + x + "\n")
                    #for x in got:
                        #print(x)
                
                    #gseq.write(protid + '\t' + i + '\n')  # write the go data to  and id to foile
               # outid.write(protid + '\t'+ length + '\n')
                data_dict[protid].append(goCopy)
                #goCopy = goCopy.replace('|','\n')
               # goCopyString ='\t' .join([str(item) for item in goCopy])
                #print(goCopy, end=",")  # print on single liens 
               # go_out.write(goCopy)              
#print (data_dict,"n")
for key, value in data_dict.items():
    #print(key , str(value))
    #print(key , value)
    print (key, '\t', ','.join(value))
    outfileFinal.write(key + '\t' + ','.join(value) + "\n")  # write to file 


goseqUniq = open("C:/Users/254357H/Documents/IPR/ptt_17FRG026_GoseqUniq.txt", 'w')
goUnq = []
with open("C:/Users/254357H/Documents/IPR/ptt_17FRG026_goseq_out.txt", 'r') as f:
    for line in f:
        if line not in goUnq:
            goUnq.append(line)  # append uniq data only
for i in goUnq:
    print(i)
    goseqUniq.write(i)  # writr go seq format to a file
    

lenUniq = open("C:/Users/254357H/Documents/IPR/ptt_17FRG026_GoseqLength.txt", 'w')
LengUnq = []
with open("C:/Users/254357H/Documents/IPR/ptt_17FRG026_id_leng.txt", 'r') as f:
    for line in f:
        if line not in LengUnq:
            LengUnq.append(line)  # append uniq data only
for i in LengUnq:
    #print(i) 
    lenUniq.write(i)  # writr protien id and length