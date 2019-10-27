        
c      Version 2.5.0
        
       program nts_generator
       implicit real*8 (a-h,o-z)
       real*8 ndens,last_density 
       parameter (number_central_densities=30)                
       parameter (neos=97)                
       dimension cen_dens(number_central_densities)
       dimension radius(10000)
       dimension smass(number_central_densities)
       dimension ndens(neos)
       dimension eps(neos)
       dimension pr(neos)
       dimension y(2), ynew(2), ders(2)
       dimension gamout(100)            
       
       integer iii
       logical format_val
       logical ntsprint
             
       dimension gam(100), alf(100), c(100), alf2(100), c2(100)
       dimension eng(100), prs(100), den(100)
       dimension rho(100), drucke(100), dichte(100)
       dimension drucke1(100), dichte1(100)
       
       dimension cnd14(100), rad14(100), smm14(100)
       dimension dehors14(100), facteur14(4,100)
       dimension dehors14red(100), facteur14red(4,100)
       dimension dehors14prs(100), facteur14prs(4,100)
       dimension dehors14eng(100), facteur14eng(4,100)
       
       dimension cnd(100), rad(100), smm(100)
       dimension cnds(100), rads(100), smms(100)
       dimension ced(100), ceds(100)
        
       dimension dehors(100), facteur(4,100)
       dimension dehors1(100), facteur1(4,100)
       dimension dehors2(100), facteur2(4,100)
       dimension dehors3(100), facteur3(4,100)
       
       dimension outptd1(100), outptd2(100)
       dimension coeffd1(100), coeffd2(100)
       dimension dpde(100)
       
       dimension set1d(100)
       dimension set1r(100),set2r(100),set3r(100),set4r(100)
       dimension set5r(100),set6r(100),set7r(100),set0r(100)
       
       dimension set1m(100),set2m(100),set3m(100),set4m(100)
       dimension set5m(100),set6m(100),set7m(100),set0m(100)
       
       dimension den66(100), prs66(100), eng66(100)
        
       dimension outpt4(100), coeff4(4,100), css(100), css_comp(100)
       dimension outpt5(100), coeff5(4,100)
       dimension denex(100), engex(100), prsex(100)
       
       dimension vsl(100), denvsl(100), outpt6(100), coeff6(4,100)
      
       character(100) :: par_title
       
       common/eosval1/outpt1(neos), coeff1(4,neos)
       common/eosval2/outpt2(neos), coeff2(4,neos)
       common/eosval3/outpt3(neos), coeff3(4,neos)
       common/block1/g,fourpi,n2,sm_solar
       common/block2/number_differential_eqs,max_rk4_steps,diff_eq_step
       common/block3/first_density,last_density, xkg, density_step
       common/block4/central_energy_density,const_1,const_2
       common/eos_total/den_all(neos),prs_all(neos),eng_all(neos)
       common/printing/nprint  
        
       ncd=number_central_densities
                       
       open(626,file='max.srt')
       open(727,file='mas14.srt')
       open(888,file='caus.srt')    
       open(666,file='outnts.srt')

       open(555,file='par.don')
       read(555,'(A)') par_title
       read(555,*) nos, ng1, ng2, g1_int, g1_inc, g2_int, g2_inc,
     1    nprint, mprint, opNAN
       
c       write(*,*) 'opNAN', opNAN
       flag = 0
       g=6.67259D-45*197.327
       pi = 3.14597d0
       fourpi=12.56637061     
       n2= neos-1    
       
       n=97
       
       if(nos .eq. 1) then
          call eos_gen              
          do i=1,n
             eng(i)=eng_all(i)
             prs(i)=prs_all(i)
             den(i)=den_all(i)
          end do
       else
          open(300,file='eosin.don')
          do i=1,n
             read(300,*) eng(i), prs(i), den(i)
          end do
       end if
       
c number of differntial equations             
       number_differential_eqs=2   
c number of Runge-Kutta iterations 
       max_rk4_steps=100000        
c dimensionless step in RK4 procedure 
       diff_eq_step=0.001     
c chosen first and last number densities 
       first_density=0.01      
       last_density=1.8
c conversion factor
       xkg=1.E+30/1.78266270
c solar mass in MeV/c^2
       sm_solar=1.989E+30*xkg
       density_step=(last_density - first_density)/
     1       number_central_densities
       
c---------------------------------------------------------------------------------------------------
c                                                                                                  |
c             Splunky Spliner Group                                                                |
c                                                                                                  | 
c---------------------------------------------------------------------------------------------------
       
       call dcsakm(neos,prs,eng,outpt1,coeff1)
       call dcsakm(neos,den,prs,outpt2,coeff2)
       call dcsakm(neos,den,eng,outpt3,coeff3)
       
c---------------------------------------------------------------------------------------------------
c                                                                                                  |
c             Just Cause Group: Causality Limit                                                    |
c                                                                                                  |
c---------------------------------------------------------------------------------------------------

       bvr=1.d0

       n73=24 !24
       n273=n73-1
       n37=73 !73

       do i=1,n73
          den66(i)=den(i+n37)
          eng66(i)=eng(i+n37)
          prs66(i)=prs(i+n37)
       end do

       nst=25 !24
       nst2=nst-1
       ninv=60 !73

       do i=1,nst
          denex(i)=den(i+ninv)
          engex(i)=eng(i+ninv)
          prsex(i)=prs(i+ninv)
       end do
                
       call dcsakm(nst,engex,prsex,outpt4,coeff4)
       call dcsakm(nst,denex,engex,outpt5,coeff5)

       call dcsakm(n73,eng66,prs66,outptd1,coeffd1)
      
       do i=1,n73
          dpde(i)=dcsder(1,eng66(i),n273,outptd1,coeffd1)
       end do
        
       call dcsakm(n73,dpde,den66,outptd2,coeffd2)
       denlim=dcsval(bvr,n273,outptd2,coeffd2)
              
c---------------------------------------------------------------------------------------------------
                                                                                                   |
c             Nunts Crumper Group: Neutron Star Mass, Radius, CenDen                               |
                                                                                                   |
c---------------------------------------------------------------------------------------------------           

       DO 200 i=1, number_central_densities
          ync=first_density+i*density_step
          central_energy_density=dcsval(ync,n2,outpt3,coeff3) 
          ced = central_energy_density                    
          pressure=dcsval(ync,n2,outpt2,coeff2) 
          const_1=1.E-19/sqrt(fourpi*g*central_energy_density)
          const_2=fourpi*central_energy_density/   
     1       (sqrt(fourpi*g*central_energy_density)**3)
          y(1)=pressure/central_energy_density
          star_radius=diff_eq_step/10.
          y(2)=((star_radius)**3)/3.
c    Start Runge-Kutta solution for fixed central density 
          k=0
          do while((pressure > 0.).and.(k<= max_rk4_steps))
             k=k+1 
             star_radius=star_radius+diff_eq_step       
             call derivatives(star_radius,y,ders)
             call runge_kutta_4(star_radius,y,ynew,ders)
             y(1)=ynew(1)
             y(2)=ynew(2)
             pressure=y(1)*central_energy_density
          end do 
          tp = 10.d0
          cen_dens(i)=ync
          radius(i)=star_radius*const_1
          smass(i)=y(2)*const_2/sm_solar

          cnd(i)=cen_dens(i)
          rad(i)=radius(i)*tp  
          smm(i)=smass(i)  

          neos1=neos-1
          css(i)= dsqrt(dcsder(1,dcsval(ync,nst2,outpt5,coeff5)
     1                  ,nst2,outpt4,coeff4))         
           
c    Printing subgroup
          if (flag .EQ. 1 .and. opNAN .EQ. 1) then 
             WRITE(666,500) cen_dens(i), radius(i)*tp, 
     1                   smass(i), 1.9999099
          else
             WRITE(666,500) cen_dens(i), radius(i)*tp, 
     1                   smass(i),css(i)
          end if
          if (css(i) .GT. 1.d0) then 
             flag = 1
          end if 

c    Just Cause subgroup: 

          if(cen_dens(i) .LT. denlim) then
             causden=cen_dens(i) 
             causrad=radius(i)*tp
             causmass=smass(i)
          end if         
200    continue 
500    format(2x,E10.4,2x,E10.4,2x,F6.4,2x,F10.4)     
c    Just Cause subgroup

       write(888,*) causmass,causrad,causden
	   
       ivsl = 2
       nedy = ncd-2
       nedy1 = nedy-1
       do i=1,nedy
          denvsl(i) = cnd(i+ivsl)
          vsl(i) = css(i+ivsl)
       end do
       call dcsakm(nedy,denvsl,vsl,outpt6,coeff6) 

c---------------------------------------------------------------------------------------------------
c             Splunky Spliners Mad Max and 14 Group                                                |
c---------------------------------------------------------------------------------------------------
                                    
       smmmax=3.0d0
       smmmin=0.9d0	   
       n2=n-1
       ncd2=ncd-1	   
       smmtemp=0.d0
       mm_c=0

       do i=1,ncd2
          mm_c=mm_c+1
          if(smm(i) .LE. smmmax .AND. smm(i) .GE. smmmin) then
             if(smm(i) .GT. smmtemp) then
                xmaxim=smm(i)
                nmax=mm_c 
             end if
          end if
          smmtemp=smm(i) 
       end do   	   
       h_mms=smm(nmax) 
         
       call dcsakm(ncd,rad,cnd,dehors1,facteur1)
       cen_num_den=dcsval(rad(nmax),ncd2,dehors1,facteur1)                         
       call dcsakm(n,den,eng,dehors2,facteur2)
       cen_eng_den=dcsval(cen_num_den,n2,dehors2,facteur2)          
       call dcsakm(n,den,prs,dehors3,facteur3)
       cen_prs=dcsval(cen_num_den,n2,dehors3,facteur3)	   
       cen_vsl = dcsval(cen_num_den,nedy1,outpt6,coeff6)
	   
       write(626,8934) h_mms,rad(nmax),
     1                 cen_num_den,cen_eng_den,cen_prs,cen_vsl
       
8934   FORMAT(2x,F9.4,2x,F9.4,2x,
     $        F9.4,2x,F9.4,2x,F9.4,2x,F9.4)           

c      Mas-14 Subgroup
       
       nmax1=nmax-1 
       nmax2=nmax1-1                   
            
       if(smm(nmax) .LT. 1.4d0) then
          cen_den_14=0.d0
          rad_14=0.d0
          go to 223
       end if
                                  
       rad_cutoff=50.d0

       do i=1,nmax1
          cnd14(i)=cnd(i)
          rad14(i)=rad(i)
          smm14(i)=smm(i)
          if(rad14(i) .GT. rad_cutoff) then
             smm14(i)=0.05d0
             rad14(i)=30.d0
          end if
       end do
          
       
       call dcsakm(nmax1,smm14,cnd14,dehors14,facteur14)
       cen_den_14 = dcsval(1.4d0,nmax2,dehors14,facteur14)
       
       call dcsakm(nmax1,cnd14,rad14,
     1             dehors14red,facteur14red)
       rad_14 = dcsval(cen_den_14,nmax2,
     1                 dehors14red,facteur14red)
              
       eng_14 = dcsval(cen_den_14,n2,dehors2,facteur2) 
          
       cen_vsl_14=dcsval(cen_den_14,nedy1,outpt6,coeff6)
       
223    continue 
           
       write(727,8933) rad_14, cen_den_14, eng_14, cen_vsl_14          
8933   FORMAT(2x,F9.4,2x,F9.4,2x,F9.4,2x,F9.4)

       end