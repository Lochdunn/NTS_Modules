from PMOD_26 import *
import subprocess

opti_pars = list_file_grab("par.don",[],False,True)
pars_reset = list_file_grab("par.don",[],False,False)
opti_pars = opti_pars[0]

n_first = opti_pars[0]
assert int(opti_pars[-1]) == 1 , "Error: 'match_type' in 'par.don' must be '1' "

par_front = str(opti_pars[0])+" "+str(opti_pars[1])+" "
par_end = " " + str(opti_pars[-2])+" "+str(opti_pars[-1])

pars = list_file_grab("mme_par.don",[],False,True)
rho_match_vals = pars[1] 
pars = pars[0]

beos = list_file_grab("beos.don",[],False,False)

ro = 0.16
r2o = 0.32


match_number = int(pars[0])

match_val_list = []
for i in range(match_number):
    match_val_list.append(float(rho_match_vals[i]))

count = 0

subprocess.call("make",shell=True)
subprocess.call("f90 $F90FLAGS -o run -s -w eos_match.f $LINK_FNL",shell=True)


vals_list = []

for i in range(match_number):
    replace_str = par_front + str(match_val_list[i]) + par_end
    text_file_replace('par.don',[1],[replace_str])
    
    subprocess.call("python opti.py",shell=True)
    subprocess.call("./run",shell=True)
    subprocess.call("rm gam.srt",shell=True)
    
    if(count == 0):
       par_front = str(45)+" "+str(opti_pars[1])+" " 
 
    vals = list_file_grab("val_out.don",[],False,False)
    vals = vals[0]      
    vals_list.append(vals)
    
    new_beos = list_file_grab("eos_beta.don",[],False,False)
    subprocess.call("rm beos.don",shell=True)
    text_file_print('beos.don',new_beos,False)
    
text_file_print('tot_vals.don',vals_list,False)
subprocess.call("rm val_out.don",shell=True)

subprocess.call("rm beos.don",shell=True)
n_reset = len(beos)
text_file_print('beos.don',beos,False)

subprocess.call("rm par.don",shell=True)
text_file_print('par.don',pars_reset,False)


