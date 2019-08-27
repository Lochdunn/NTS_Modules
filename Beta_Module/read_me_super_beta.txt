
Read me for "Beta Module"

1.) Input the nuclear asym and sym EoS into 'ex_nxlo.don':

den  e_sym  e_asym
.    .      . 
.    .      .
.    .      .

2.) Input parameters into 'par.don':

n  m  nprint

n: number of EoS points
m: number of output points if nprint = 2,3
nprint: 1,2,3 (If 1, the output EoS at the same density points as the input)
              (If 2, n # of output points ending at the density point 1.8 fm^-3)
              (If 3, n # of output points ending at the density point 0.64 fm^-3)
        else  no output EoS is printed

2.) Run with "./sb"

Output files: 

'frac_nxlo.don' : The fraction of particle by species

'beta_eos.don': Output EoS in beta-equalibrium

