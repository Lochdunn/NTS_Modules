#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PMOD_26 import *
import subprocess

subprocess.call("rm mas14.srt",shell=True)
subprocess.call("rm max.srt",shell=True)
subprocess.call("rm caus.srt",shell=True)
subprocess.call("rm eosout.don",shell=True)
subprocess.call("rm outnts.srt",shell=True)

reset = 1
counter = 1

tstart = time.time()

#subprocess.call("f90 $F90FLAGS -o run -s -w nts_os.f nts_routines.f $LINK_FNL",shell=True)
 
pareos_grab = list_file_grab("parameos.don",[],False,True)
pareos_grab = pareos_grab[0]
par_grab = list_file_grab("par.don",[],False,True)
par_grab = par_grab[1]

pareos_str = ""
for i in range(len(pareos_grab)-2): pareos_str = pareos_str+str(pareos_grab[i])+" "
if (reset == 1):
    reset_str = list_file_grab('parameos.don',[],False,False)
    reset_str = [str(reset_str[0])]  

nos = int(par_grab[0])
n   = int(par_grab[1])
m   = int(par_grab[2])
gam1_o = float(par_grab[3])
gam1_inc = float(par_grab[4])
gam2_o = float(par_grab[5])
gam2_inc = float(par_grab[6])
nprint = int(par_grab[7])  
ntsprint = int(par_grab[8])
p = n*m

neos = int(pareos_grab[1])
ncd = 30

gam1s=[]
gam2s=[]

for i in range(n):
    g1 = gam1_o + i*gam1_inc    
    gam1s.append(g1)

for i in range(m): 
    g2 = gam2_o + i*gam2_inc
    gam2s.append(g2)

pareos_inc = '' 
mas14_group, max_group, caus_group, eos_group, nts_group = empty_lists(5)
mas14_col, max_col, caus_col, eos_col, nts_col = empty_lists(5)

if(nos != 1):
   n = 0
   m = 0
   subprocess.call("./run",shell=True)

for i in range(n):        
    eos_file_name = "EoS_G_"+str(gam1s[i])+".txt"
    nts_file_name = "NTS_G_"+str(gam1s[i])+".txt" 
    remove_cmd_eos = "rm "+eos_file_name
    remove_cmd_nts = "rm "+nts_file_name
    subprocess.call(remove_cmd_eos,shell=True)
    subprocess.call(remove_cmd_nts,shell=True)

# Generating Loop

for i in range(n):
    for j in range(m):
        pareos_inc = [pareos_str + str(gam1s[i]) + " " + str(gam2s[j])]
        text_file_replace("parameos.don","parameos.don",[1],pareos_inc)
        subprocess.call("./run",shell=True)
    
        mas14_grab = list_file_grab("mas14.srt",[],False,True)[0] 
        max_grab = list_file_grab("max.srt",[],False,True)[0]
        caus_grab = list_file_grab("caus.srt",[],False,True)[0] 
        eos_grab = list_file_grab("eosout.don",[],False,True) 
        nts_grab = list_file_grab("outnts.srt",[],False,True)

        mas14_group.append(mas14_grab)
        max_group.append(max_grab)
        caus_group.append(caus_grab)
        eos_group.append(eos_grab)
        nts_group.append(nts_grab)

        counter+=1
     
    mas14_col.append(mas14_group)
    max_col.append(max_group)
    caus_col.append(caus_group)
    eos_col.append(eos_group)
    nts_col.append(nts_group)   
    
    mas14_group, max_group, caus_group, eos_group, nts_group = empty_lists(5) 

if(nos == 1):
    subprocess.call("rm mas14.srt",shell=True)
    subprocess.call("rm max.srt",shell=True)
    subprocess.call("rm caus.srt",shell=True)
    subprocess.call("rm eosout.don",shell=True)
    subprocess.call("rm outnts.srt",shell=True)
    
# Printing Loops


mas14_title = "Γ(1) Γ(2) Radius    Cen_Den   Cen_Eng   SoS\n"
with open('mas14.srt','a') as mas14_file:    
    mas14_file.write(mas14_title)
    for i in range(n):
        mas14_file.write("\n")
        for j in range(m):
            head_str = str(gam1s[i]) + "  "+str(gam1s[j])+"  "
            mas14_str = head_str + float_list_to_string(mas14_col[i][j])+str("\n")
            mas14_file.write(mas14_str)
            
max_title = "Γ(1) Γ(2) Max_Mass  Radius    Cen_Den   Cen_Eng   Cen_Prs   SoS\n"
with open('max.srt','a') as max_file:    
    max_file.write(max_title)
    for i in range(n):
        max_file.write("\n")
        for j in range(m):
            head_str = str(gam1s[i]) + "  "+str(gam1s[j])+"  "
            max_str = head_str + float_list_to_string(max_col[i][j])+str("\n")
            max_file.write(max_str)  

caus_title = "Γ(1) Γ(2) CausMass  Radius    Cen_Den\n"
with open('caus.srt','a') as caus_file:    
    caus_file.write(caus_title)
    for i in range(n):
        caus_file.write("\n")
        for j in range(m):
            head_str = str(gam1s[i]) + "  "+str(gam1s[j])+"  "
            caus_str = head_str + float_list_to_string(caus_col[i][j])+str("\n")
            caus_file.write(caus_str)  
            
eos_count_run = 1
eos_header_string ='Den       '
eos_out_str = ''

for i in range(m):
    if(i<10):
        eos_header_string = eos_header_string + "Prs_G" + str(gam2s[i])+"  "+ "Eng_G" +str(gam2s[i])+"  " 
    else: 
        eos_header_string = eos_header_string + "Prs_G" + str(gam2s[i])+"  "+ "Eng_G" +str(gam2s[i])+"  " 

if(nprint == 1):
    for i in range(n):        
        eos_file_name = "EoS_G_"+str(gam1s[i])+".txt"
        with open(eos_file_name,'a') as eos_file_write:
            eos_file_write.write(eos_header_string+"\n")
            for k in range(neos):
                for j in range(m):                          
                    if(eos_count_run == 1):
                        eos_add_1 = float(eos_col[i][j][k][0])
                        eos_add_2 = float(eos_col[i][j][k][1])
                        eos_add_3 = float(eos_col[i][j][k][2])         
                        eos_out_str = eos_out_str + numeric_uni_format(eos_add_3) + \
                        "  " + numeric_uni_format(eos_add_1) + "  "+ numeric_uni_format(eos_add_2) + "  "
                    else:
                        eos_add_1 = float(eos_col[i][j][k][0])
                        eos_add_2 = float(eos_col[i][j][k][1])
                        eos_out_str = eos_out_str + numeric_uni_format(eos_add_1) + \
                        "  " + numeric_uni_format(eos_add_2) + "  "
                    eos_count_run+=1
                eos_out_str = eos_out_str+"\n"
                eos_file_write.write(eos_out_str)
                eos_out_str = ''  
                eos_count_run = 1
            eos_count_run = 1
                
            
nts_count_run = 1
nts_header_string ='Den       '
nts_out_str = ''

for i in range(m):
    if(i<10):
        nts_header_string = nts_header_string + "Rad_G" + str(gam2s[i])+"  "+ "Mas_G" + str(gam2s[i])+ \
                            "  " + "SOS_G" + str(gam2s[i]) + "  "
    else: 
        nts_header_string = nts_header_string + "Rad_G" + str(gam2s[i])+"  "+ "Mas_G" + str(gam2s[i])+ \
                            "  " + "SOS_G" + str(gam2s[i]) + "  "

#Γ
if(ntsprint == 1):
    for i in range(n):
        nts_file_name = "NTS_G_"+str(gam1s[i])+".txt"
        with open(nts_file_name,'a') as nts_file_write:
            nts_file_write.write(nts_header_string+"\n")
            for k in range(ncd):
                for j in range(m):                    
                    if(nts_count_run == 1):
                        nts_add = nts_col[i][j][k]
                        nts_add_1 = float(nts_col[i][j][k][0])
                        nts_add_2 = float(nts_col[i][j][k][1])
                        nts_add_3 = float(nts_col[i][j][k][2])         
                        nts_add_4 = float(nts_col[i][j][k][3])
                        nts_out_str = nts_out_str + numeric_uni_format(nts_add_1) + \
                        "  " + numeric_uni_format(nts_add_3) + "  " + numeric_uni_format(nts_add_2) + \
                        "  " + numeric_uni_format(nts_add_4) + "  "
                    else:
                        nts_add_1 = float(nts_col[i][j][k][1])
                        nts_add_2 = float(nts_col[i][j][k][2])
                        nts_add_3 = float(nts_col[i][j][k][3])
                        nts_out_str = nts_out_str + numeric_uni_format(nts_add_2) + \
                        "  " + numeric_uni_format(nts_add_1) + "  " + numeric_uni_format(nts_add_3) + "  "
                    nts_count_run+=1                 
                nts_out_str = nts_out_str+"\n"
                nts_file_write.write(nts_out_str) 
                nts_out_str = ''
                nts_count_run=1
            nts_count_run=1

if(reset == 1):
    text_file_replace('parameos.don','parameos.don',[1],reset_str)
        
tend = time.time()
deltat = tend - tstart
time_took = convert_time(deltat,False)
print(" ")
print("Success!")
print('The process took '+time_took+' to finish.'+"\n")


