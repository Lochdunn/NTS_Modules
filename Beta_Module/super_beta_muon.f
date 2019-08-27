
      program super_beta                                             
      implicit real*8 (a-h,o-z)           
      common/paspoi/pas(200),poi(200),xfs(200),wfs(200)
      common/const/hc,pi,pi2,fmu,fmu2,fmn,fmp,c1,c2   
      common/fracs/oh,ot,ht,qt,ft,deltm     
      common/ys/Y_u(100),Y_e(100),Y_p(100),Y_n(100),ierr(100)     
        
      parameter(nn=100)
   
      dimension rho(nn), ratio(nn), xrho(nn)                     
      dimension enm(nn), e0(nn),xkf(nn)                                   
      dimension e_n(nn),e_0(nn),e_s(nn)                
      dimension press(nn), break1(nn), coeff(4,nn)              
      dimension e_el(nn),e_u(nn),e_B(nn),e_tot(nn),eps(nn)       

      dimension coef2(4,nn),coef3(4,nn)
      dimension break2(nn), break3(nn)         
          
      open(200,file='par.don')
      open(20,file='ex_nxlo.don')   
      open(32,file='frac_nxlo.don')                              
      open(33,file='beta_eos.don') 

      read(200,*) n, m, nprint
      n2=n-1      
       
      hc=197.327d0
      pi=3.1415927d0
      pi2=pi**2
      
      oh=1.d0/2.d0
      ot=1.d0/3.d0
      ht=2.d0/3.d0
      qt=4.d0/3.d0
      ft=5.d0/3.d0
          
      fmu=105.658d0
      fmu2=fmu*fmu
      fmn=939.57d0
      fmp=938.28d0
      deltm=fmn-fmp
          
      c1=hc*(3.d0*pi2)**(ot)
      c2=hc**2/(6.d0*pi2*fmu)*(3.d0*pi2)**ft
          
c-----Read-in-data-----c

      do i=1,n
         read(20,*) xkf(i), e0(i), enm(i)
         rho(i)=2.d0*xkf(i)**3/(3.d0*pi2)
         e_s(i)=enm(i)-e0(i)
      end do

c calculate the fractions

      call frac_eval(n,rho,e_s)

c fractions output
      
      do i=1,n
         add_frac=Y_p(i)+Y_n(i)
         write(32,104) rho(i),Y_u(i),Y_e(i),Y_n(i),
     1                 Y_p(i),add_frac,ierr(i)
      end do
  
c energy calculatations    
  
c calculate e_el

      do i=1,n 
         e_el(i)=e_el_func(rho(i),Y_e(i))
      end do
      
c calculate e_u

c--------Approx. Muon expression: Ultra-Relativistic Approx.

c      do i=1,n
c         e_u(i) = e_u_ura(rho(i),Y_u(i))
c      end do

c--------Analytical Muon Energy per Particle Expression

c      do i=1,n
c         e_u(i) = e_u_anal(rho(i),Y_u(i))
c      end do

c--------Numeric Muon Energy per Particle Expression 

       do i=1,n
          e_u(i) = e_u_num(rho(i),Y_u(i))
       end do

c--------Baryon Energy per Particle Expression
                 
      do i=1,n
         alp=Y_n(i)-Y_p(i) 
         easy=e0(i)+e_s(i)*alp**2
         e_B(i) = easy + e_el(i) + e_u(i)
      end do

c calculate e_tot=total energy/baryon 

      do i=1,n
         e_tot(i) = e_B(i) + fmp*Y_p(i) + fmn*Y_n(i)         
      end do

c calculate energy_density 
     
      do i=1,n
         eps(i) = rho(i)*e_tot(i)
      end do

c output                            
      
      call dcsakm(n,rho,eps,break1,coeff)
      do i=1,n
         der=dcsder(1,rho(i),n2,break1,coeff)
         press(i)=rho(i)*der - eps(i)
         if(nprint .EQ. 1) then
            write(33,500) eps(i), press(i), rho(i), e_tot(i)
         end if
      end do
      call dcsakm(n,rho,press,break2,coef2)
      
      if(nprint .EQ. 2) then
         rho0 = rho(1)
         do i=1,m
            rho_ext = rho0+(i-1)*(1.8d0/(m-1))         
            prs_ext = dcsval(rho_ext,n2,break2,coef2)
            eng_ext = dcsval(rho_ext,n2,break1,coeff)
            write(33,560) eng_ext, prs_ext, rho_ext 
         end do
      end if
      
      if(nprint .EQ. 3) then
         rho0 = rho(1)
         do i=1,m
            rho_ext = rho0+(i-1)*(0.64d0/(m-1))         
            prs_ext = dcsval(rho_ext,n2,break2,coef2)
            eng_ext = dcsval(rho_ext,n2,break1,coeff)
            write(33,560) eng_ext, prs_ext, rho_ext 
         end do
      end if      
         
  104 format(F11.5,2x,F11.5,2x,F11.5,2x,F11.5,2x,F11.5,2x,F11.5,2x,i3)  
  500 format(6x,e14.8,4x,e14.8,4x,e14.8,4x,e14.8)
  560 format(6x,e14.8,4x,e14.8,4x,e14.8)

      end         


c---------------------------------------------------------------------------------------------------|

      function e_u_num(rhoi,y_u)
      implicit real*8 (a-h,o-z)
      common/paspoi/pas(200),poi(200),xfs(200),wfs(200)  
      common/const/hc,pi,pi2,fmu,fmu2,fmn,fmp,c1,c2   
      common/fracs/oh,ot,ht,qt,ft,deltm        
      
      tfac=1.d0/(pi2*hc**3)
      sum=0.d0
      x1=0.d0
      xkfmu=(3.d0*pi2*y_u*rhoi)**ot
      x2=xkfmu*hc
      xinf=0.d0
      nps=90
      xn=1.d0
      call lgauss(nps)
      call paspoid(x1,x2,1,nps,xinf)
      sum=0.d0
      do i=1,nps
         r=pas(i)
         ww=poi(i)
         funct=eps_mu(r)*ww
         sum=sum+funct
      end do
      EV_u=tfac*sum
      e_u_num= EV_u/rhoi      
      end  

c------------------------------------------------------------------|

      function e_el_func(rhoi,y_e)
      implicit real*8(a-h,o-z)
      common/const/hc,pi,pi2,fmu,fmu2,fmn,fmp,c1,c2   
      common/fracs/oh,ot,ht,qt,ft,deltm        
      eel=(hc/(4.d0*pi2))*(3.d0*pi2*y_e)**(qt)
      e_el_func = eel*rhoi**ot
      end 

c------------------------------------------------------------------|

      function e_u_ura(rhoi,y_u)
      implicit real*8(a-h,o-z)
      common/const/hc,pi,pi2,fmu,fmu2,fmn,fmp,c1,c2   
      common/fracs/oh,ot,ht,qt,ft,deltm    
      ura_1 = (hc**2.d0)/(10.d0*pi**2.d0*fmu)
      ura_2 = rhoi**(ht)*(3.d0*pi2*y_u)**(ft)      
      e_u_ura = ura_1*ura_2
      end

c------------------------------------------------------------------|

      function e_u_anal(rhoi,y_u)
      implicit real*8(a-h,o-z)
      common/const/hc,pi,pi2,fmu,fmu2,fmn,fmp,c1,c2   
      common/fracs/oh,ot,ht,qt,ft,deltm     
      tnorm = 2.d0/(8.d0*pi2)
      hc3=hc*hc*hc
      tfac = tnorm/hc3
      xkfmu = (3.d0*pi2*rhoi*y_u)**(ot)
      xkcon = xkfmu*hc
      xme = (xkcon**2+fmu**2)**oh
      evmu1 = xkcon*xme**3
      evmu2 = (-1.5d0)*fmu2*xme*xkcon
      evmu3 = (-0.5d0)*fmu2*fmu2*dlog((xme+xkcon)/fmu)
      EV_u = tfac*(evmu1+evmu2+evmu3) 
      e_u_anal = EV_u/rhoi
      end

c------------------------------------------------------------------|

      subroutine frac_eval(n,rho,e_s)
      implicit real*8(a-h,o-z)
      common/const/hc,pi,pi2,fmu,fmu2,fmn,fmp,c1,c2   
      common/fracs/oh,ot,ht,qt,ft,deltm     
      common/ys/Y_u(100),Y_e(100),Y_p(100),Y_n(100),ierr(100)  

      dimension rho(100), e_s(100)
      dimension resul(100), intera(100)
      
      x0=0.d0                                                           
      x1=1.d0                                                         
      a=0.d0                                                   
      b=1.d0
      e=1e-6
      do i = 1,n
         ift = 0
         tmp1 = rho(i)
         tmp2 = e_s(i)
         call Bisect(x0,x1,e,x,m,ift,tmp1,tmp2)
         x0 = 0.d0
         x1 = 1.d0
         Y_u(i) = x
         resul(i) = Y(x,tmp1,tmp2)
         intera(i) = m
         ierr(i) = ift       
         Y_e(i)=(((c1*(Y_u(i)*tmp1)**ot)**2+fmu2)**oh/
     1          (c1*tmp1**ot))**3         
         if(ierr(i) .EQ. 1) then
            Y_u(i) = 0.d0
            call Quartile(a,b,e,root,m,ift,tmp1,tmp2)
            a = 0
            b = 1
            resul(i) = Y2(root,tmp1,tmp2)
            intera(i) = m
            ierr(i) = ift
            Y_e(i) = root
         end if
      end do

      do i=1,n
         Y_p(i) = Y_e(i)+Y_u(i)
         Y_n(i) = 1.d0-Y_p(i)
      end do
      return
      end
          
          
c***********************************************
c*        Bisection method subroutine          *
c* ------------------------------------------- *
c* This routine iteratively seeks the zero of  *
c* function Y(x) using the method of interval  *
c* halving until the interval is less than e   *
c* in width. It is assumed that the function   *
c* Y(x) is available from a function routine.  *
c* ------------------------------------------- *
c* Input values: range (x0,x1), and e.         *
c* Output values: root x, Y(x) and number of   *
c* steps m.                                    *
c***********************************************

      subroutine Bisect(x0,x1,e,x,m,ift,tmp1,tmp2)
      implicit real*8(a-h,o-z)
     
      ift = 0
      m=0

      if (Y(x0,tmp1,tmp2) * Y(x1,tmp1,tmp2) .GT. 0.0) then
        ift = 1
        return
      else
  10    y0=Y(x0,tmp1,tmp2)
        x=(x0+x1)/2.d0
        yy=Y(x,tmp1,tmp2)
        m=m+1
        if (yy*y0 .EQ. 0.0) then 
           return
        else if (yy*y0 .LT. 0.0) then
           x1=x
        else if (yy*y0 .GT. 0.0) then
           x0=x
        end if
        if (dabs(x1-x0) .GT. e) then
           goto 10
        end if
      end if
      return
      end

c***********************************************
c         Quartile method subroutine           *
c -------------------------------------------- *
c This routine iteratively seeks the zero of   *
c function Y(x) using the method of interval   *
c quarting until the interval is less than tol *
c in width. It is assumed that the function    *
c Y(x) is available from a function routine.   *
c -------------------------------------------- *
c Input values: range (a, b), and tol.         *
c Output values: found root and number of used *
c steps.                                       *
c***********************************************

      subroutine Quartile(a, b, tol, root, step, ft,tmp1,tmp2)         
      implicit real*8(a-h,o-z)

      ift = 0
      co=0.25
c     co can be 0.3

      istep=0
      if (Y2(a,tmp1,tmp2) * Y2(b,tmp1,tmp2) .GT. 0.0) then
        ift = 1
        return
      else
        do while (dabs(a-b) .GT. tol)
          if (dabs(Y2(a,tmp1,tmp2)) .LT. dabs(Y2(b,tmp1,tmp2))) then
            tm = a + co*(b-a)
          else
            tm = a + (1.0-co)*(b-a)
          end if
          if (Y2(tm,tmp1,tmp2)*Y2(a,tmp1,tmp2) .GT. 0.0) then
            a= tm 
          else 
            b=tm
          end if 
          istep=istep+1
        end do
        root = (a+b)/2.d0
      end if
      return
      end

c---------------------------------------------------------------------------------------------------

c-------Gaussian-Legendre-Quadrature-Numeric-Integration-Routines-BEGIN

c---------------------------------------------------------------------------------------------------


      subroutine lgauss(n)
c
      implicit real*8(a-h,o-z)
      dimension z(200),wz(200)
      common/paspoi/pas(200),poi(200),xfs(200),wfs(200)
      if(n-1) 1,2,3
 1       return

 2       z(1)=0.d0
         wz(1)=2.d0

         return
 3       r=dfloat(n)
      g=-1.d0
      do 147 i=1,n
         test=-2.d0
         ic=n+1-i
         if(ic.lt.i) go to 150
 4       s=g
         t=1.d0
         u=1.d0
         v=0.d0
         do 50 k=2,n
            a=dfloat(k)
            p=((2.d0*a-1.d0)*s*g-(a-1.d0)*t)/a
            dp=((2.d0*a-1.d0)*(s+g*u)-(a-1.d0)*v)/a
            v=u
            u=dp
            t=s
 50      s=p
         if(abs((test-g)/(test+g)).lt.0.5d-09) go to 100 
         sum=0.d0
         if(i.eq.1) go to 52
         do 51 j=2,i
            sum=sum+1.d0/(g-xfs(j-1))
 51      continue
 52      test=g
         g=g-p/(dp-p*sum)
         go to 4
 100     xfs(ic)=-g
         xfs(i)=g
         wfs(i)=2.d0/(r*t*dp)
         wfs(ic)=wfs(i)
 147  g=g-r*t/((r+2.d0)*g*dp+r*v-2.d0*r*t*sum)
 150  do 160 i=1,n
         z(i)=xfs(i)
         wz(i)=wfs(i)
 160  continue
      return
      end

c---------------------------------------------------------------------------------------------------



      subroutine paspoid(xi,xf,ni,nf,xinf)
      implicit real*8(a-h,o-z)
      common/paspoi/pas(200),poi(200),xfs(200),wfs(200)
      coeff=-1.d0
      xs=(xi+xf)/2.d0
      xd=(xf-xi)/2.d0
      do i=ni,nf
         pas(i)=xs+xd*xfs(i)
         poi(i)=xs*wfs(i)
c         write(*,*) pas(i), poi(i)
      end do
      return
      end


c---------------------------------------------------------------------------------------------------

c-------Gaussian-Legendre-Quadrature-Numeric-Integration-Routines-END

c---------------------------------------------------------------------------------------------------

c-------Muon-E/A-Numeric-Integral-Function
      function eps_mu(xkf)
      implicit real*8(a-h,o-z)
      common/const/hc,pi,pi2,fmu,fmu2,fmn,fmp,c1,c2   
      common/fracs/oh,ot,ht,qt,ft,deltm   
      xk=xkf
      xk2=xk*xk      
      eps_mu=dsqrt(xk2+fmu2)*xk2
      end
c-------Function end

c-------Muon-E/A-Numeric-Integral-Function
      function eps_mu0(xkf)
      implicit real*8(a-h,o-z)
      common/const/hc,pi,pi2,fmu,fmu2,fmn,fmp,c1,c2   
      common/fracs/oh,ot,ht,qt,ft,deltm   
      xk=xkf
      xk2=xk*xk
      eps_mu0=dsqrt(1.d0+fmu/xk2)*xk2*xk
      end
c-------Function end


cImported from 'tquart.f90'

      Function Y(x,tmp1,tmp2)
      implicit real*8(a-h,o-z)
      common/const/hc,pi,pi2,fmu,fmu2,fmn,fmp,c1,c2   
      common/fracs/oh,ot,ht,qt,ft,deltm   
      Ye=((((c1*(x*tmp1)**ot)**2+fmu2)**oh)/(c1*tmp1**ot))**3
      Y = -4.d0*tmp2*(1.d0-2.d0*(x+Ye))-deltm+c1*(tmp1*Ye)**(ot)
      return
      end

      Function Y2(x,tmp1,tmp2)
      implicit real*8 (a-h,o-z)
      common/const/hc,pi,pi2,fmu,fmu2,fmn,fmp,c1,c2   
      common/fracs/oh,ot,ht,qt,ft,deltm   
      diff=fmp-fmn
      Y2 = -4.d0*tmp2*(1.d0-2.d0*x)+diff+c1*(x*tmp1)**(ot)
      return
      end
          
