
      PROGRAM server
      implicit real*8 (a-h, o-z)
      common/block1/n,m,prs_match,rho_match,dentemp(100),
     1              den(100),pres(100),eng(100),prssec(100)
       
c SERVER INTIALIZATION
       
      dimension :: rho(100), eps(100)
      dimension :: druke(100), epfrac(100)
      dimension :: break(100), coeff(4,100), depdp(100)
      
      open(100,file='beos.don')
      open(400,file='par.don')
      
      read(400,*) n, m, rho_match, rho_jump, match_type
      n2=n-1
      ro = 0.16d0
      r2o = 2.d0*rho0
     
c     match_type = 0, m = 6 , rho_match = 0.32, rho_jump = 2.0
c     match_type = 1, m = 6 , rho_match = 0.32, rho_jump = 2.0
c     match_type = 2, m = 6 , rho_match = 0.32, rho_jump = 2.0
      write(*,*) n
      do i=1,n
         read(100,*) eng(i), pres(i), den(i)
      end do
      
      call dcsakm(n,den,pres,break,coeff)       
      prs_match  = dcsval(rho_match,n2,break,coeff)
 
      do i=1,m
         if(match_type .EQ. 0) then 
            dentemp(i) = (rho_match)+(rho_match*(rho_jump-1.d0))
     1                   *((i-1)/(1.d0*(m-1)))         
         else if(match_type .EQ. 1) then
            dentemp(i) = (r2o)+(rho_match)
     1                   *((i-1)/(1.d0*(m-1)))
         else 
            dentemp(i) = (r2o)+(0.32)*((i-1)/(1.d0*(m-1)))  
         end if
         prssec(i) = dcsval(dentemp(i),n2,break,coeff)
      end do
       
c SERVER MAIN LOOP

c read the command and parameters
c if command is 0 (or rather "not 1"), evaluate energy and repeat
c if command is 1, stop the server

      icmd = 0
      do 100 while (icmd .ne. 1)         
         read(*,*) icmd, gam           
         if (icmd .ne. 1) then
            dev = sigma(gam)
            write(*, *) dev
         end if       
 100  continue           
           
      stop
      END PROGRAM
      
      
      FUNCTION sigma(gam) 
      implicit real*8 (a-h,o-z)
      common/block1/n,m,prs_match,rho_match,dentemp(100),
     1              den(100),pres(100),eng(100),prssec(100)
      sum=0.d0
      do 20 i=1,m  
         functfit= (prs_match/(rho_match**gam))*(dentemp(i))**gam
         err=dabs((prssec(i)-functfit)/prssec(i)) 
   20    sum=sum+err
      sum=sum/(m*1.d0)   
      sigma=sum               
      return
      end 

