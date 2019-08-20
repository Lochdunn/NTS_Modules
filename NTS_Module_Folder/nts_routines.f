c---------------------------------------------------------------------------------------------------
c                                                                                                  | 
c      Function\Subroutine                                                                         | 
c                                                                                                  |  
c---------------------------------------------------------------------------------------------------

       subroutine eos_gen
       implicit real*8 (a-h,o-z)
       parameter (neos=97) 

       common/eos_total/den_all(neos),prs_all(neos),eng_all(neos)
       common/printing/nprint

       dimension eng(100), prs(100), den(100)
       dimension alf(100), c(100), alf2(100), c2(100)
       dimension rho(100), drucke(100), dichte(100)
       dimension drucke1(100), dichte1(100)
       
c      open(000,file='dump.don')
       open(323,file='parameos.don')
       open(300,file='eosin.don')
       open(333,file='eosout.don')
       read(323,*) m, n, rho_i, rho_f, rho_fin, 
     1             match, gam, gam2

       if (match .EQ. 1) then
          valle = rho_i
       else if (match .NE. 1) then
          valle = rho_f
       end if
        
       rho_inc = (rho_f - rho_i)/4.d0  
       rho_pas = (rho_fin - rho_f)/m
       rho_com = rho_f
        
       inac = 0
       do i=1,n
          read(300,*) eng(i), prs(i), den(i)
          if (den(i) .LT. valle) then
             rhobc = den(i)
             pbc = prs(i)
             ebc = eng(i) 
             inac = inac +1
c          write(000,*) ebc, pbc, rhobc
          end if     
       end do
c       write(*,*) inac       
       alf = pbc/(rhobc**gam)
       c = (1/rhobc)*(ebc-((alf/(gam-1))*(rhobc**gam)))
       
       if (match .EQ. 1) then
          do i=1,4
             den(i+inac)=rho_i+rho_inc*(i-1.d0)
             prs(i+inac)=druck(den(i+inac),gam,alf)
             eng(i+inac)=avg(den(i+inac),gam,alf,c)
          end do
             alf2 = prs(inac+4)/(den(inac+4)**gam2)
             c2 = (1/den(inac+4))*(eng(inac+4)-((alf2/
     1               (gam2-1))*(den(inac+4)**gam2)))
          do i=1,23
             den(i+4+inac)=rho_com+rho_pas*(i-1)
             prs(i+4+inac)=druck(den(i+4+inac),gam2,alf2)
             eng(i+4+inac)=avg(den(i+4+inac),gam2,alf2,c2)
          end do
       else if (match .NE. 1) then
          do i=1,27
             den(i+inac)=rho_i+rho_inc*(i-1.d0)
             prs(i+inac)=druck(den(i+inac),gam,alf)
             eng(i+inac)=avg(den(i+inac),gam,alf,c)
          end do
       end if
       
       if(nprint.EQ.1) then
          do i=1,neos
             write(333,*) eng(i), prs(i), den(i)
          end do
       end if

       do i=1,neos
          den_all(i) = den(i)
          prs_all(i) = prs(i)
          eng_all(i) = eng(i) 
       end do

       return
       end
                  
       function avg(r,gam,alfa,cc)
          implicit real*8(a-h,o-z)
          avg=(alfa/(gam-1.d0))*r**gam+(cc*r) 
          return
       end
       
       function druck(r,gam,alfa)
          implicit real*8(a-h,o-z) 
          druck=(alfa*(r**gam))
          return 
       end       

       SUBROUTINE runge_kutta_4(x,y,yout,dydx)
          implicit real*8(a-h,o-z) 
          real*8 last_density 
          dimension yt(2), dyt(2), dym(2) 
          dimension y(2), dydx(2), yout(2) 
          common/block1/g,fourpi,n2,sm_solar 
          common/block2/number_differential_eqs,max_rk4_steps,
     1                  diff_eq_step 
          common/block3/first_density,last_density,xkg,density_step 
          common/block4/central_energy_density,const_1,const_2  

          hh=diff_eq_step*0.5                                  
          h6=diff_eq_step/6.                 
          xh=x+hh 
          yt(1)=y(1)+hh*dydx(1)
          yt(2)=y(2)+hh*dydx(2) 
          call derivatives(xh,yt,dyt)
          yt(1)=y(1)+hh*dyt(1) 
          yt(2)=y(2)+hh*dyt(2) 
          call derivatives(xh,yt,dym)
          yt(1)=y(1)+diff_eq_step*dym(1)              
          yt(2)=y(2)+diff_eq_step*dym(2)              
          dym(1)=dyt(1)+dym(1) 
          dym(2)=dyt(2)+dym(2) 
          call derivatives(x+diff_eq_step,yt,dyt)
          yout(1)=y(1) + h6*(dydx(1)+dyt(1)+2.*dym(1))
          yout(2)=y(2) + h6*(dydx(2)+dyt(2)+2.*dym(2))
       return 
       end                              

       subroutine derivatives(r,y,ders)                 
          implicit real*8(a-h,o-z) 
          real*8 last_density 
          parameter (neos=97)                
          dimension y(2), ders(2)                      
          common/eosval1/outpt1(neos), coeff1(4,neos)                  
          common/eosval2/outpt2(neos), coeff2(4,neos)                  
          common/eosval3/outpt3(neos), coeff3(4,neos)                  
          common/block1/g,fourpi,n2,sm_solar 
          common/block2/number_differential_eqs,max_rk4_steps,
     1                  diff_eq_step 
          common/block3/first_density,last_density, xkg,density_step
          common/block4/central_energy_density,const_1,const_2   
          
          if(y(1).gt.0.) then 
             p1=y(1)*central_energy_density
             e_rho=dcsval(p1,n2,outpt1,coeff1)
             e_rho=e_rho/central_energy_density
c         
             ders(1)=tov(r,e_rho,y)
             ders(2)=(r**2)*e_rho
          end if
          return                      
       end 

       DOUBLE PRECISION FUNCTION tov(r,e_rho,y)
          implicit real*8(a-h,o-z)
          dimension y(2)   
          tov=-(e_rho+y(1))*                 
     1       (y(2) + (r**3)*y(1))/(r*r-2*r*y(2))
          return            
       end     
