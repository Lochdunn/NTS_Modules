       program matchmaker
       implicit real*8(a-h,o-z)
       common/block1/pr2o,r2o       
       dimension :: den(100), prs(100), eng(100)
       dimension :: rho(100), p(100), e(100)
       dimension :: x(100), y(100)
       dimension :: break(100), coeff(4,100)
       dimension :: dehors(100), facteur(4,100)
       
       open(100,file='gam.srt')
       open(555,file='par.don')
       open(900,file='beos.don')
       open(999,file='eos_beta.don')
       open(444,file='val_out.don')
        
        
       read(555,*) n, m ,rho_match, rho_jump, match_type
              
       nm=n
       k=27
       k2=k+1
       
       m=16
       m2=m+1
       
       r2o=rho_match
       s_i=rho_match
       f_i=0.034584296d0
       f_f=0.3021615185d0
       s_f=1.824d0
       s_step=(s_f-s_i)/k      
       f_step=(f_f-f_i)/m      
       n2=n-1
       
       
       read(100,*) gam
       
       do i = 1,n
          read(900,*) eng(i), prs(i), den(i) 
       end do   
            
       call dcsakm(n,den,prs,break,coeff)       
       pr2o  = dcsval(r2o,n2,break,coeff)       
       alp = (pr2o/(r2o**gam))
        
       call dcsakm(n,den,eng,dehors,facteur)
       er2o = dcsval(r2o,n2,dehors,facteur)
        
       c=(er2o-(((r2o**gam)*alp)/(gam-1.d0)))/r2o
        
       write(444,1554) gam, alp, c

1554   Format(F10.3,2x,F10.3,2x,F8.3)
       
       do i=1,m2
          temp_d = f_i + (i-1)*f_step
          temp_p = dcsval(temp_d,n2,break,coeff)
          temp_e = dcsval(temp_d,n2,dehors,facteur)
          write(999,1550) temp_e, temp_p, temp_d
       end do          
       
       do i=1,k2
          temp_d=s_i+(i-1)*s_step
          temp_p=druke(alp,gam,temp_d)
          temp_e=energie(alp,gam,c,temp_d)
          write(999,1550) temp_e, temp_p, temp_d
       end do
       
1550   Format(F10.3,2x,F10.3,2x,F8.3)
        
       end 
        
       function druke(alp,gam,r)
          implicit real*8(a-h,o-z)
          common/block1/pr2o,r2o
          druke=alp*(r)**gam
       end function
       
       function energie(alp,gam,c,r) 
          implicit real*8(a-h,o-z)
          common/block1/pr2o,r2o
          energie=((alp/(gam-1.d0))*((r)**gam))+c*r
       end function              
