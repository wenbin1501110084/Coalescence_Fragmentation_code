C********************************************************************C
C                                                                    ! 
C     Originally Written by K. C. Han, R. J. Fries, and C.-M. Ko     C 
C     Improved by Wenbin Zhao, 1501110084@pku.edu.cn         !
C     When you use this coalescence code, please cite these          C
C      three papers:                                                 !
C       W.~Zhao, etc. al.[arXiv:2103.14657 [hep-ph]].                C
C       W.~Zhao,etc.al.Phys. Rev. Lett.125, (2020) no.7, 072301.     !
C       K. C. Han, R. J. Fries and C. M. Ko,                         C
C                           Phys. Rev. C93,no. 4, 045207 (2016).     !
C====================================================================C
C     Recombination + Remnant Fragmentatoin Module                   ! 
C     Input:                                                         C 
C     nth: number of thermal partons, idth(i): thermal parton id     ! 
C     xth(i,4): postions and times of thermal partons                C 
C     pth(i,4): 4-momenta of thermal partons                         ! 
C     nsh: number of shower partons, idsh(i): shower parton ids      C 
C     xsh(i,4): postion and time of shower partons                   ! 
C     psh(i,4): 4-momenta of shower partons                          C 
C     Output:                                                        ! 
C     nf:number of final hadrons, idf(i): hadron id                  C 
C     xf(i,4): position and time of the final hadrons                ! 
C     pf(i,4): 4-momenta of the final hadrons                        C 
C     Kst(i): origin of the hadrons(th-th:0, sh-th:1, sh-sh:3)       ! 
C                                                                    C 
C********************************************************************! 
      SUBROUTINE hadronization(nth,idth,xth,pth,nsh,idsh,xsh,psh,
     .     nf,idf,xf,pf,Kst,Npt,Kpx,Ppx,Xpx,
     .    NCTnr, NCTKr,CTPr,CTXr,QQscale,qfscale) 
 
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      common/const/pi,hbarc
      common/parm/sigma,sigmapro,sigmaK,maxf
      
      DIMENSION idth(5000),xth(5000,4),pth(5000,4)
      DIMENSION idsh(5000),xsh(5000,4),psh(5000,4)
      DIMENSION idf(5000),xf(5000,4),pf(5000,4)
      
      DIMENSION Kst(5000),QQscale(5000),rqscale(5000)
      DIMENSION P(5000,5),XV(5000,5),qfscale(5000) 
      DIMENSION Kp(5000,5),qscale(5000)
      DIMENSION KH(5000,2),PH(5000,5),XH(5000,5)
      DIMENSION KfH(5000),PfH(5000,5),XfH(5000,5)
      DIMENSION Kr(5000,5),Pr(5000,5),Xr(5000,5)
      dimension Kpx(5000), Ppx(5000,5),Xpx(5000,5)
      dimension NCTKr(5000), CTPr(5000,4), CTXr(5000,4) 
      numpt = 0
C......Thermal Partons .....................
      DO ith=1, nth
         numpt = numpt + 1
         Kp(numpt,2) = idth(ith) ! parton id
         Kp(numpt,1) = 0         ! thermal
         Kp(numpt,3) = 0         ! origin 
         P(numpt,1) = pth(ith,1)
         P(numpt,2) = pth(ith,2)
         P(numpt,3) = pth(ith,3)
         P(numpt,4) = pth(ith,4)
         XV(numpt,1) = xth(ith,1)
         XV(numpt,2) = xth(ith,2)
         XV(numpt,3) = xth(ith,3)
         XV(numpt,4) = xth(ith,4)
         qscale(numpt) =0.0
      ENDDO

C......Shower Partons .....................                 
      DO ish=1, nsh
         numpt = numpt + 1
         Kp(numpt,2) = idsh(ish) ! parton id                  
         Kp(numpt,1) = 1         ! shower
         Kp(numpt,3) = 0         ! origin                     
         P(numpt,1) = psh(ish,1)
         P(numpt,2) = psh(ish,2)
         P(numpt,3) = psh(ish,3)
         P(numpt,4) = psh(ish,4)
         XV(numpt,1) = xsh(ish,1)
         XV(numpt,2) = xsh(ish,2)
         XV(numpt,3) = xsh(ish,3)
         XV(numpt,4) = xsh(ish,4)
         qscale(numpt)=QQscale(ish)
      ENDDO
      
C.......total number of thermal and shower partons .....
      numtot = nth+nsh
      call recomb(numtot,Kp,P,XV,nH,KH,PH,XH,nr,Kr,Pr,Xr,
     .              NCTnr,NCTKr,CTPr,CTXr,qscale,rqscale)
      nf = 0
C.....the number of recombined hadrons ..............
      do ireco=1, nH
         nf = nf + 1
         
         idf(nf) = KH(ireco,2)
         Kst(nf) = KH(ireco,1)  !...th+sh(1), sh+sh (2), th+th(0)

         pf(nf,1) = PH(ireco,1)
         pf(nf,2) = PH(ireco,2)
         pf(nf,3) = PH(ireco,3)
         pf(nf,4) = PH(ireco,4)
         xf(nf,1) = XH(ireco,1)
         xf(nf,2) = XH(ireco,2)
         xf(nf,3) = XH(ireco,3)
         xf(nf,4) = XH(ireco,4)
      enddo

      call remorg(nr,Kr,Pr,Xr,Npt,Kpx,Ppx,Xpx,rqscale,
     .            qfscale)
       do i=1,10
       enddo
      RETURN
      END

C*****************************************************!
C                                                     C
C   Subroutine for recombination                      !
C   input :                                           C
C   np : number of partons                            !
C   Kp(i,j) :                                         C
C     j=1: thermal(0) or shower(1),sh(fail e cut) -1  !
C           shower(outside kinetic cut) -2            C
C      j=2: particle id, j=3: origin                  !
C   Pp(i,j):  momenta and energies of partons         C
C     j=1(px), 2(py), 3(pz), 4(e)                     !
C   Xp(i,j):  positions and times of partons          C
C     j=1 (x), 2(y), 3(z), 4(t)                       !
C   Output:                                           C
C   nH: number of hadrons from recombination,         !
C   KH(i,2): particle id of the hadrons(1:pion,2:kaon,C
C     3: nucleon, 4: Lambda )                         !
C   KH(i,1): 0(th+th), 1(sh+th), 2(sh+sh)             C
C   PH(i,j): 4-momenta of hadrons (px,py,pz,E)        !
C   XH(i,j): space-time of hadrons (x,y,z,t)          !
C   nr: number of remnant quarks,                     C
C   Kr(i,j): origin and pid of remnant quarks         !
C    j=1 (origin), 2(pid),3:thermal(0),shower(1)      C
C     Kr(i,1)=0 (leading quarks), 1,2,..(numerical    !
C        orders of their mother gluons)               C
C   Pr: 3-momenta of remnant quarks                   !
C   Xr: space-time of remnant quarks                  !
C                                                     C
C*****************************************************!
      subroutine recomb(np,Kp,Pp,Xp,nH,KH,PH,XH,nr,Kr,Pr,Xr,
     .                NCTnr,NCTKr,CTPr,CTXr,qscale,rqscale)
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
      common/parm/sigma,sigmapro,sigmaK,maxf

      dimension Kp(5000,5), Pp(5000,5), Xp(5000,5)
      dimension Kq(5000,5), Pq(5000,5), Xq(5000,5)
      dimension KH(5000,2), PH(5000,5), XH(5000,5)
      dimension Kr(5000,5), Pr(5000,5), Xr(5000,5)
      dimension nchg(5000),qscale(5000),Q0(5000)
      dimension xmdx(4), pMx(3), pdx(4,4)
      dimension xmTx(4), pTx(3), pdBx(4,4)
      dimension NCTKr(5000), CTPr(5000,4), CTXr(5000,4) 
      dimension rqscale(5000)
C..... Constants ................
      pi = 3.1415926535897932384626433832795
      hbarc = 0.197327          ! [ GeV*fm ]  
      hbarc2 = hbarc*hbarc
      gmax = 2.0                ! [GeV] maximum mass of gluon ************* 

!*****widths for recombination ****************** 
      Del = hbarc/sigma
!**** radius derived from Charged ***** 

      Sig2 = sigma*sigma
      SigK2 = sigmaK*sigmaK
      SigN2 = sigmapro*sigmapro
      SigRB = 0.88
      sigmaphi=0.65 
      Sigphi2 = sigmaphi*sigmaphi
      SigRB2 = SigRB*SigRB
      SigLB = sqrt(3.)*SigRB/2.
      SigLB2 = SigLB*SigLB
      SigRO = 1.20!sigma of Omega 
      SigRO2 =SigRO*SigRO        
      SigLBO = sqrt(3.)*SigRO2/2.
      SigLBO2 = SigLBO*SigLBO    

!**** masses of u,d and s quarks **************************
      xmq = 0.25                    !mass_u/d=0.25 at T=T_freezeout
      xmq2 = xmq*xmq
      xms = 0.475                   !mass_s=0.475 at T=T_freezeout
      xms2 = xms*xms
      
!**** masses of hadrons **********************
      xmpi = 0.1395
      xmrho = 0.77
      xmK = 0.49
      xmKs = 0.892
      xmKs1 = 1.41
      xmK1 = 1.27
      xmpro = 0.935
      xmDel = 1.232
      xmN1 = 1.44
      xmDel1 = 1.6
      xmrho1 = 0.64             !*** mass of rho meson in Deltal1->N  rho   
      xmLamb = 1.11568
      xmLamb1 = 1.52            !***
      xmSig = 1.385
C....Initialize number of quarks, number of hadrons, number of remnants
      numq =0
      nH = 0
      nr = 0      
      NCTnr=0
      pttcutsquare=1.5*1.5 ! The pt cut for thermal partons

!**** Let's recombine quarks to get pions ********
!*******************************************************!                       
!     Gluon Decays into q qbar:  gluon -> q qbar      !                       
!*******************************************************!                       
      DO ipt=1, np
 
!******leading quarks **************************** 
         if(abs(Kp(ipt,2)).le.3) then
            numq = numq + 1
            Kq(numq,1) = Kp(ipt,1) ! thermal or shower ******** 
            Kq(numq,2) = Kp(ipt,2) ! id 
            Kq(numq,3) = 0         !*** origin .... 0 for leading quark... 
            
            Pq(numq,1) = Pp(ipt,1) 
            Pq(numq,2) = Pp(ipt,2)
            Pq(numq,3) = Pp(ipt,3)
            Pq(numq,4) = Pp(ipt,4)
            Xq(numq,1) = Xp(ipt,1)
            Xq(numq,2) = Xp(ipt,2)
            Xq(numq,3) = Xp(ipt,3)
            Xq(numq,4) = Xp(ipt,4)
            Q0(numq)   = qscale(ipt)
!******gluons ***********************************                               
         else
            if(Kp(ipt,2).eq.21) then
               xmgx = 2.*xms + (gmax - 2.*xms)*ran() ! virtualities of gluons 
               xgpx = Pp(ipt,1)
               xgpy = Pp(ipt,2)
               xgpz = Pp(ipt,3)
               xge = Pp(ipt,4)
               xgx = Xp(ipt,1)
               xgy = Xp(ipt,2)
               xgz = Xp(ipt,3)
               xgt = Xp(ipt,4)
               
               call gluondec(xmgx,xgx,xgy,xgz,xgt,xgpx,xgpy,xgpz,xgnx,
     .              xgny,xgnz,xgnt,nflq1,qpx1,qpy1,qpz1,qe1,nflq2,
     .              qpx2,qpy2,qpz2,qe2)
                  numq = numq + 1
                  Kq(numq,1) = Kp(ipt,1) !*** thermal or shower *****              
                  Kq(numq,2) = nflq1
                  Kq(numq,3) = ipt
               
                  Pq(numq,1) = qpx1
                  Pq(numq,2) = qpy1
                  Pq(numq,3) = qpz1
                  Pq(numq,4) = qe1
                  Xq(numq,1) = xgnx
                  Xq(numq,2) = xgny
                  Xq(numq,3) = xgnz
                  Xq(numq,4) = xgnt
                  Q0(numq)   = qscale(ipt) 

                  numq = numq + 1
                  Kq(numq,1) = Kp(ipt,1) !*** thermal or shower *****              
                  Kq(numq,2) = nflq2
                  Kq(numq,3) = ipt
               
                  Pq(numq,1) = qpx2
                  Pq(numq,2) = qpy2
                  Pq(numq,3) = qpz2
                  Pq(numq,4) = qe2
                  Xq(numq,1) = xgnx
                  Xq(numq,2) = xgny
                  Xq(numq,3) = xgnz
                  Xq(numq,4) = xgnt
                  Q0(numq)   = qscale(ipt)

            endif
         endif
      enddo
C***************************************************
C     Recombinatoin                              !
C***************************************************
      ktotf = 0
      do ix=1, numq
         nchg(ix) = Kq(ix,2)/abs(Kq(ix,2))
         if(Kq(ix,3).eq.1) ktotf = ktotf + nchg(ix)
         if(abs(Kq(ix,2)).le.2) Pq(ix,5) = xmq
         if(abs(Kq(ix,2)).eq.3) Pq(ix,5) = xms
      enddo
      numqx = numq
     
C_******************************************
      !probbz=ran()
      !pz2=ran()
      do lzz=1,numq*numq*10
         pzzz=ran()
      enddo
C_***************************************** 
!***  first loop **************                                                  
      do i=1, numqx-2
         if(nchg(i).eq.0) goto 30
         px1 = Pq(i,1)
         py1 = Pq(i,2)
         pz1  = Pq(i,3)
         e1 =  Pq(i,4)
         
         x1 = Xq(i,1)
         y1 = Xq(i,2)
         z1 = Xq(i,3)
         t1 = Xq(i,4)
!         print*,"test1"
C**** second loop **************                                                
         do j= i+1, numqx-1
            
C......avoiding double counting of quarks .......                               
            if(nchg(j).eq.0) goto 40
C......Avoiding recombination of q qbar decayed from a same gluon .......       
            !if(Kq(i,3).eq.Kq(j,3)) goto 40
            if((Kq(i,3).eq.Kq(j,3)).and. (Kq(i,3).ge.1)) go to 40 !changed by wenbin 
            px2 = Pq(j,1)
            py2 = Pq(j,2)
            pz2 = Pq(j,3)
            e2 = Pq(j,4)
            
            x2 = Xq(j,1)
            y2 = Xq(j,2)
            z2 = Xq(j,3)
            t2 = Xq(j,4)
            
            if(Kq(i,2)*Kq(j,2).gt.0) then               
!***********third loop *********************                                    
               do j2=j+1, numqx
C.......NOT count th-th-th recombination ...............                       
C.......NOT count th-th-th recombination ...............                       
                if (Kq(i,1).eq.0.and.Kq(j,1).eq.0.and.Kq(j2,1).eq.0)
     .              then! changed by wenbin 2018.11.23             !TTT       
!     .                 goto 50 !               changed by wenbin 2018.11.23                                      !TTT

                  px3 = Pq(j2,1)
                  py3 = Pq(j2,2)
                  pz3 = Pq(j2,3)

                !pt1square=px1*px1+py1*py1
                !pt2square=px2*px2+py2*py2
                !pt3square=px3*px3+py3*py3
                !if(pt1square.lt.pttcutsquare)goto 50
                !if(pt2square.lt.pttcutsquare)goto 50
                !if(pt3square.lt.pttcutsquare)goto 50
              
                else www=1.0
                endif

C.......NOT count ha-ha-ha recombination ...............
!                  if (Kq(i,1)*Kq(j,1)*Kq(j2,1).gt.0)!changed by wenbin 2018.11.23             !TTT
!     .                 goto 50 !               changed by wenbin2018.11.23                                      !TTT

                  if(Kq(i,2)*Kq(j2,2).lt. 0)goto 50  !added by wenbin 2018.11.23

!                  if(((Kq(i,1).eq.0).and.(Kq(j,1).eq.0)
!     .                     .and.(Kq(j2,1).eq.0)))goto 50!don't consider  the th-th-th recombination !Wenbin

                  px3 = Pq(j2,1)
                  py3 = Pq(j2,2)
                  pz3 = Pq(j2,3)
                  e3 = Pq(j2,4)
                  
                  x3 = Xq(j2,1)
                  y3 = Xq(j2,2)
                  z3 = Xq(j2,3)
                  t3 = Xq(j2,4)
                  
                  pxbaryon = px1 + px2 + px3
                  pybaryon = py1 + py2 + py3
                  pzbaryon = pz1 + pz2 + pz3
                  
!...........velocity of center-of-mass .........................                
                  betaBx = pxbaryon/(e1+e2+e3)
                  betaBy = pybaryon/(e1+e2+e3)
                  betaBz = pzbaryon/(e1+e2+e3)
                  
C.........speed of the Lab frame in the center-of-mass frame ......            
                  betaBLx = -betaBx
                  betaBLy = -betaBy
                  betaBLz = -betaBz
!..........Lorentz Boost into center-of-mass frame ..........                   
                  call lorentzboost(betaBx,betaBy,betaBz,e1,px1,py1,
     .                 pz1,e1Bcm,px1Bcm,py1Bcm,pz1Bcm)
                  call lorentzboost(betaBx,betaBy,betaBz,t1,x1,y1,
     .                 z1,t1Bcm,x1Bcm,y1Bcm,z1Bcm)
                  call lorentzboost(betaBx,betaBy,betaBz,e2,px2,py2,
     .                 pz2,e2Bcm,px2Bcm,py2Bcm,pz2Bcm)
                  call lorentzboost(betaBx,betaBy,betaBz,t2,x2,y2,
     .                 z2,t2Bcm,x2Bcm,y2Bcm,z2Bcm)
                  call lorentzboost(betaBx,betaBy,betaBz,e3,px3,
     .                 py3,pz3,e3Bcm,px3Bcm,py3Bcm,pz3Bcm)
                  call lorentzboost(betaBx,betaBy,betaBz,t3,x3,y3,
     .                 z3,t3Bcm,x3Bcm,y3Bcm,z3Bcm)
                  
!******** invariant mass of three quarks  **********!added by wenbin 208.11.23
                  xinvB=e1Bcm+e2Bcm+e3Bcm!added by wenbin 208.11.23
!*******************relative momenta defined in the baryon-rest-frame ************* 
                  xk1x = (Pq(j,5)*px1Bcm - Pq(i,5)*px2Bcm)/
     .                 (Pq(i,5)+Pq(j,5))
                  xk1y = (Pq(j,5)*py1Bcm - Pq(i,5)*py2Bcm)/
     .                 (Pq(i,5)+Pq(j,5))
                  xk1z = (Pq(j,5)*pz1Bcm - Pq(i,5)*pz2Bcm)/
     .                 (Pq(i,5)+Pq(j,5))
                  xk12 = xk1x**2.+xk1y**2.+xk1z**2.
                  
                  xk2x = (Pq(j2,5)*(px1Bcm+px2Bcm)-(Pq(i,5)+
     .                 Pq(j,5))*px3Bcm)/(Pq(i,5)+Pq(j,5)+Pq(j2,5))
                  xk2y = (Pq(j2,5)*(py1Bcm+py2Bcm)-(Pq(i,5)+
     .                 Pq(j,5))*py3Bcm)/(Pq(i,5)+Pq(j,5)+Pq(j2,5))
                  xk2z = (Pq(j2,5)*(pz1Bcm+pz2Bcm)-(Pq(i,5)+
     .                 Pq(j,5))*pz3Bcm)/(Pq(i,5)+Pq(j,5)+Pq(j2,5))
                  xk22 = xk2x**2. + xk2y**2. + xk2z**2.
                  
!***************Speeds of quarks in the baryon-rest frame ***************       

                  vx1Bcm = px1Bcm/e1Bcm
                  vy1Bcm = py1Bcm/e1Bcm
                  vz1Bcm = pz1Bcm/e1Bcm
                  
                  vx2Bcm = px2Bcm/e2Bcm
                  vy2Bcm = py2Bcm/e2Bcm
                  vz2Bcm = pz2Bcm/e2Bcm

                  vx3Bcm = px3Bcm/e3Bcm
                  vy3Bcm = py3Bcm/e3Bcm
                  vz3Bcm = pz3Bcm/e3Bcm
                  
!************Free propgation of earlier created quarks to the time latest produced one *********** 
                  tMaxBcm = max(t1Bcm,t2Bcm,t3Bcm)
                  x1nBcm = x1Bcm + vx1Bcm*(tMaxBcm - t1Bcm)
                  y1nBcm = y1Bcm + vy1Bcm*(tMaxBcm - t1Bcm)
                  z1nBcm = z1Bcm + vz1Bcm*(tMaxBcm - t1Bcm)
                  
                  x2nBcm = x2Bcm + vx2Bcm*(tMaxBcm - t2Bcm)
                  y2nBcm = y2Bcm + vy2Bcm*(tMaxBcm - t2Bcm)
                  z2nBcm = z2Bcm + vz2Bcm*(tMaxBcm - t2Bcm)

                  x3nBcm = x3Bcm + vx3Bcm*(tMaxBcm - t3Bcm)
                  y3nBcm = y3Bcm + vy3Bcm*(tMaxBcm - t3Bcm)
                  z3nBcm = z3Bcm + vz3Bcm*(tMaxBcm - t3Bcm)
                  
C......position of center of mass at the born time of the youngest parton ... 
                  xcmB = (x1nBcm*Pq(i,5)+x2nBcm*Pq(j,5)+x3nBcm*Pq(j2,5))
     .                 /(Pq(i,5)+Pq(j,5)+Pq(j2,5))
                  ycmB = (y1nBcm*Pq(i,5)+y2nBcm*Pq(j,5)+y3nBcm*Pq(j2,5))
     .                 /(Pq(i,5)+Pq(j,5)+Pq(j2,5))
                  zcmB = (z1nBcm*Pq(i,5)+z2nBcm*Pq(j,5)+z3nBcm*Pq(j2,5))
     .                 /(Pq(i,5)+Pq(j,5)+Pq(j2,5))
                  tcmB = tMaxBcm

C............position and time of the formed hadron in the Lab frame ....    
                  call lorentzboost(betaBLx,betaBLy,betaBLz,tcmB,xcmB,
     .                 ycmB,zcmB,tBL,xBL,yBL,zBL)
C..   print*, tBL                                                
                  xr1x = (x1nBcm-x2nBcm)
                  xr1y = (y1nBcm-y2nBcm)
                  xr1z = (z1nBcm-z2nBcm)
                  xr12 = xr1x**2. + xr1y **2. + xr1z**2.
                  
                  xr2x =((Pq(i,5)*x1nBcm+Pq(j,5)*x2nBcm)
     .                 /(Pq(i,5)+Pq(j,5))-x3nBcm)
                  xr2y = ((Pq(i,5)*y1nBcm+Pq(j,5)*y2nBcm)
     .                 /(Pq(i,5)+Pq(j,5))-y3nBcm)
                  xr2z = ((Pq(i,5)*z1nBcm+Pq(j,5)*z2nBcm)
     .                 /(Pq(i,5)+Pq(j,5))-z3nBcm)
                  xr22 = xr2x**2. + xr2y**2. + xr2z**2.
                  
!********* include the Omega's sigma !begin added by wenbin 2018.11.23 *******
                  if(abs(Kq(i,2)*Kq(j,2)
     .                 *Kq(j2,2)).eq. 9)then
                     SigR2=SigRO2
                     SigL2=SigLBO2
                  else 
                     SigR2=SigRB2
                     SigL2=SigLB2
                  endif
                  urhox = 0.5*(xr1x**2./SigR2+(xk1x**2.)*SigR2/hbarc2)
                  urhoy = 0.5*(xr1y**2./SigR2+(xk1y**2.)*SigR2/hbarc2)
                  urhoz = 0.5*(xr1z**2./SigR2+(xk1z**2.)*SigR2/hbarc2)
                  
                  ulambx=0.5*(xr2x**2./SigL2+(xk2x**2.)*SigL2/hbarc2)
                  ulamby=0.5*(xr2y**2./SigL2+(xk2y**2.)*SigL2/hbarc2)
                  ulambz=0.5*(xr2z**2./SigL2+(xk2z**2.)*SigL2/hbarc2)
                  cccut=40.0!CCC the cut for the wigner fucntion !added by wenbin
                  if((urhox+urhoy+urhoz).gt. cccut)goto 50!added by wenbin
                  if((ulambx+ulambx+ulambx).gt. cccut)goto 50!added by wenbin

!********* include the Omega's sigma !end added by wenbin 2018.11.23 *******
                  
                  sumWigB = 0.
                  
!*******1D-ground w. Wig functions .......                                     
                  wig0BRx = exp(-urhox)
                  wig0BRy = exp(-urhoy)
                  wig0BRz = exp(-urhoz)
                  
                  wig0BLx = exp(-ulambx)
                  wig0BLy = exp(-ulamby)
                  wig0BLz = exp(-ulambz)
                  wig0B = wig0BRx*wig0BRy*wig0BRz*wig0BLx*wig0BLy*
     .                 wig0BLz
                  wig0B = wig0B/27.0 ! 1/27 is the color
!.... Summing 3D w.w. functions up to nlev ......... 
                  wigR1x = wig0BRx
                  maxfB=10 !!Proton's highest excited state
                  if(abs(Kq(i,2)*Kq(j,2)*Kq(j2,2)).eq.6)maxfB=5 !Lambda's highest excited state

                  do iRx=0, maxfB
                     wigR1y = wig0BRy
                     do iRy = 0, maxfB-iRx
                        wigR1z = wig0BRz
                        do iRz = 0, maxfB-iRx-iRy
                           wigL1x = wig0BLx
                           do iLx=0, maxfB-iRx-iRy-iRz
                              wigL1y = wig0BLy
                              do iLy=0, maxfB-iRx-iRy-iRz-iLx
                                 wigL1z = wig0BLz
                                 do iLz=0, maxfB-iRx-iRy-iRz-iLx-iLy
                                    sumWigB = sumWigB + wigR1x*
     .                                   wigR1y*wigR1z*wigL1x*wigL1y
     .                                   *wigL1z
                                    wigL1z = wigL1z*ulambz/(iLz+1)
                                 enddo
                                 wigL1y=wigL1y*ulamby/(iLy+1)
                              enddo
                              wigL1x = wigL1x*ulambx/(iLx+1)
                           enddo
                           wigR1z= wigR1z*urhoz/(iRz+1)
                        enddo
                        wigR1y = wigR1y*urhoy/(iRy+1)
                     enddo
                     wigR1x = wigR1x*urhox/(iRx+1)
                  enddo

                  sumWigB=sumWigB/27.0 ! 1/27.0 added by wenbin color freedom 2019.03.10;  *2 for the identifical two quarks
                  Wigpro =  sumWigB!changed by wenbin to account the color degeneracy 2019.02.24
!                  Wigpro =  sumWigB*ZKfactor !changed by wenbin to account the color degeneracy 2019.02.24
!                  Wiglambda = 0.!sumWigB !changed by wenbin 2018.11.13
                  Wiglambda = sumWigB/2.0!Only the Lambda and Sigma0 
!                  print*, Wigpro
                  prob5 = ran()
                  
!****************! begin added by wenbin 2018.11.23 ********************! 
!    Recombine into Omega(3334)in the ground states (n=0) !added by wenbin 2018.11.23
                  if(abs(Kq(i,2)*Kq(j,2)*Kq(j2,2)).eq. 9) then
                     if(prob5.le.wig0B/2.) then 
                        nH = nH + 1
                        KH(nH,2) = 6
!                        if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) KH(nH,1) = 2
!                        if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1
                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25

                        
                        PH(nH,1) = pxbaryon
                        PH(nH,2) = pybaryon
                        PH(nH,3) = pzbaryon
                        PH(nH,4) = sqrt(pxbaryon**2.+pybaryon**2.+
     .                       pzbaryon**2.+xmLamb**2.)
                        
                        XH(nH,1) = xBL
                        XH(nH,2) = yBL
                        XH(nH,3) = zBL
                        XH(nH,4) = tBL
                        
                        nchg(i) = 0
                        nchg(j) = 0
                        nchg(j2)= 0
                        goto 30
                     endif
                  endif
!****************! end added by wenbin 2018.11.23 ********************! 
                  prob4 = ran()
!*********************************************************************! 
!    Recombine into Sigma(1385) and Sigma->Lambda + pion              ! 
!*********************************************************************!
                  
                  if(prob4.le.Wiglambda.and.abs(Kq(i,2)*Kq(j,2)*Kq(j2
     .                 ,2)).eq.6) then
                  if(xinvB.lt.xmLamb)goto 50
!*******************Lambda in the ground states (n=0) *********************  
                     if(prob4.le.wig0B/4.) then
                        nH = nH + 1
                        KH(nH,2) = 4
!                        if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) KH(nH,1) = 2
!                        if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25


                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25
                        
                        PH(nH,1) = pxbaryon
                        PH(nH,2) = pybaryon
                        PH(nH,3) = pzbaryon
                        PH(nH,4) = sqrt(pxbaryon**2.+pybaryon**2.+
     .                       pzbaryon**2.+xmLamb**2.)
                        
                        XH(nH,1) = xBL
                        XH(nH,2) = yBL
                        XH(nH,3) = zBL
                        XH(nH,4) = tBL
                        
                        nchg(i) = 0
                        nchg(j) = 0
                        nchg(j2)= 0
                        goto 30
                     endif
!*******************end of Lambda in the ground states (n=0) *********************  
                     probSigDec = ran()
                     if(probSigDec.le.0.5) then
                        xmSig = xinvB
                        
                        if(xmSig.le.xmLamb+xmpi) xmSig = xmLamb+xmpi
                        call twobody(xmSig,xmLamb,xmpi,pxbaryon,
     .                       pybaryon,pzbaryon,pxSigtoLamb,pySigtoLamb,
     .                       pzSigtoLamb,eSigtoLamb,pxSigtopi,pySigtopi,
     .                       pzSigtopi,eSigtopi)
                        
                        
                        nH = nH + 1
                        KH(nH,2) = 4

                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25
                        
                        PH(nH,1) = pxSigtoLamb
                        PH(nH,2) = pySigtoLamb
                        PH(nH,3) = pzSigtoLamb
                        PH(nH,4) = eSigtoLamb
                        vxSigtoLamb = pxSigtoLamb/eSigtoLamb
                        vySigtoLamb = pySigtoLamb/eSigtoLamb
                        vzSigtoLamb = pzSigtoLamb/eSigtoLamb
                        
                        gamSigtoLamb = hbarc/sqrt(1.-vxSigtoLamb**2.-
     .                       vySigtoLamb**2.-vzSigtoLamb**2.)
                        XH(nH,1) = xBL+(pxSigtoLamb/eSigtoLamb)*
     .                       gamSigtoLamb/xmLamb
                        XH(nH,2) = yBL+(pySigtoLamb/eSigtoLamb)*
     .                       gamSigtoLamb/xmLamb
                        XH(nH,3) = zBL+(pzSigtoLamb/eSigtoLamb)*
     .                       gamSigtoLamb/xmLamb
                        XH(nH,4) = tBL+gamSigtoLamb/xmLamb
                        
                        nH = nH + 1
                        KH(nH,2) = 1

                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25
                        PH(nH,1) = pxSigtopi
                        PH(nH,2) = pySigtopi
                        PH(nH,3) = pzSigtopi
                        PH(nH,4) = eSigtopi
                        vxSigtopi = pxSigtopi/eSigtopi
                        vySigtopi = pySigtopi/eSigtopi
                        vzSigtopi = pzSigtopi/eSigtopi
                        gamSigtopi = hbarc/sqrt(1.-vxSigtopi**2.
     .                       -vySigtopi**2.-vzSigtopi**2.)
                        XH(nH,1) = xBL+(pxSigtopi/eSigtopi)*
     .                       gamSigtopi/xmpi
                        XH(nH,2) = yBL+(pySigtopi/eSigtopi)*
     .                       gamSigtopi/xmpi
                        XH(nH,3) = zBL+(pzSigtopi/eSigtopi)*
     .                       gamSigtopi/xmpi
                        XH(nH,4) = tBL+gamSigtopi/xmpi
                        
                        nchg(i) = 0
                        nchg(j) = 0
                        nchg(j2) = 0
                        goto 30
                     endif

C...... Sigma -> Lambda + pi + pi .......................
                     if(probSigDec.gt.0.5) then
                    if(xinvB.le. 1.192)goto 50 ! Sigma0 -> Lambda +pi
                    if(xinvB.le.xmLamb+xmpi+xmpi)xinvB=xmLamb+xmpi+xmpi
                        call threebody(xinvB,xmLamb,xmpi,xmpi,pxbaryon,
     .                       pybaryon,pzbaryon,pxSigtoLamb1,
     .                       pySigtoLamb1,pzSigtoLamb1,eSigtoLamb1,
     .                       pxSigtopi1,pySigtopi1,pzSigtopi1,
     .                       eSigtopi1,pxSigtopi2,pySigtopi2,
     .                       pzSigtopi2,eSigtopi2)
                        
                        nH = nH + 1
                        KH(nH,2) = 4
 
                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25
                        PH(nH,1) = pxSigtoLamb1
                        PH(nH,2) = pySigtoLamb1
                        PH(nH,3) = pzSigtoLamb1
                        PH(nH,4) = eSigtoLamb1
                        vxSigtoLamb1 = pxSigtoLamb1/eSigtoLamb1
                        vySigtoLamb1 = pySigtoLamb1/eSigtoLamb1
                        vzSigtoLamb1 = pzSigtoLamb1/eSigtoLamb1
                        gamSigtoLamb1 = hbarc/sqrt(1.-vxSigtoLamb1**2.
     .                       -vySigtoLamb1**2.-vzSigtoLamb1**2.)
                        XH(nH,1) = xBL+(pxSigtoLamb1/eSigtoLamb1)*
     .                       gamSigtoLamb1/xmLamb                       
                        XH(nH,2) = yBL+(pySigtoLamb1/eSigtoLamb1)*
     .                       gamSigtoLamb1/xmLamb                       
                        XH(nH,3) = zBL+(pzSigtoLamb1/eSigtoLamb1)*
     .                       gamSigtoLamb1/xmLamb                       
                        XH(nH,4) = tBL+gamSigtoLamb1/xmLamb
                        
                        nH = nH + 1
                        KH(nH,2) = 1

                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25



                        PH(nH,1) = pxSigtopi1
                        PH(nH,2) = pySigtopi1
                        PH(nH,3) = pzSigtopi1
                        PH(nH,4) = eSigtopi1
                        vxSigtopi1 = pxSigtopi1/eSigtopi1
                        vySigtopi1 = pySigtopi1/eSigtopi1
                        vzSigtopi1 = pzSigtopi1/eSigtopi1
                        
                        gamSigtopi1 = hbarc/sqrt(1.-vxSigtopi1**2.
     .                       -vySigtopi1**2.-vzSigtopi1**2.)
                        XH(nH,1) = xBL+(pxSigtopi1/eSigtopi1)*
     .                       gamSigtopi1/xmpi
                        XH(nH,2) = yBL+(pySigtopi1/eSigtopi1)*
     .                       gamSigtopi1/xmpi
                        XH(nH,3) = zBL+(pzSigtopi1/eSigtopi1)*
     .                       gamSigtopi1/xmpi
                        XH(nH,4) = tBL+gamSigtopi1/xmpi
                        
                        nH = nH + 1
                        KH(nH,2) = 1

                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25
                        PH(nH,1) = pxSigtopi2
                        PH(nH,2) = pySigtopi2
                        PH(nH,3) = pzSigtopi2
                        PH(nH,4) = eSigtopi2
                        vxSigtopi2 = pxSigtopi2/eSigtopi2
                        vySigtopi2 = pySigtopi2/eSigtopi2
                        vzSigtopi2 = pzSigtopi2/eSigtopi2
                        gamSigtopi2 = hbarc/sqrt(1.-vxSigtopi2**2.
     .                       -vySigtopi2**2.-vzSigtopi2**2.)
                        XH(nH,1) = xBL+(pxSigtopi2/eSigtopi2)*
     .                       gamSigtopi2/xmpi
                        XH(nH,2) = yBL+(pySigtopi2/eSigtopi2)*
     .                       gamSigtopi2/xmpi
                        XH(nH,3) = zBL+(pzSigtopi2/eSigtopi2)*
     .                       gamSigtopi2/xmpi
                        XH(nH,4) = tBL+gamSigtopi2/xmpi
                        
                        nchg(i) = 0
                        nchg(j) = 0
                        nchg(j2) = 0
                        goto 30
                     endif
                  endif

!*************************************************************************
!          Recombination of nucleons and Delta                           !
!*************************************************************************
                  prob3 = ran()              
!************Recombine proton ***********************************
                  if(prob3.le.Wigpro*2.0.and.(abs(Kq(i,2)).le.2.and.    ! Wigpro*2.0, 2.0 is for the two identifical quarks
     .                 abs(Kq(j,2)).le.2.and.abs(Kq(j2,2)).le.2)) then
                     
!.........Direct N and Delta(N*) -> N + pi in the ground state ...........
                     if(prob3.le.wig0B) then
                        proNg = ran()
                        if(proNg.le.0.25) then
                           
                           nH = nH + 1
                           KH(nH,2) = 3
                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25
                           
                           PH(nH,1) = pxbaryon
                           PH(nH,2) = pybaryon
                           PH(nH,3) = pzbaryon
                           PH(nH,4) = sqrt(pxbaryon**2.+pybaryon**2.+
     .                          pzbaryon**2.+xmpro**2.)
                           XH(nH,1) = xBL
                           XH(nH,2) = yBL
                           XH(nH,3) = zBL
                           XH(nH,4) = tBL
                           
                           nchg(i) = 0
                           nchg(j) = 0
                           nchg(j2) = 0
                           goto 30
                        endif

!................ Delta -> N + pi ...................
                        if(proNg.ge.0.25.and.proNg.lt.0.75) then
                           call twobody(xmDel,xmpro,xmpi,pxbaryon,
     .                          pybaron,pzbaryon,pxDeltoN,pyDeltoN,
     .                          pzDeltoN,eDeltoN,pxDeltopi,pyDeltopi,
     .                          pzDeltopi,eDeltopi)
                           nH = nH + 1
                           KH(nH,2) = 3

                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25

                           PH(nH,1) = pxDeltoN
                           PH(nH,2) = pyDeltoN
                           PH(nH,3) = pzDeltoN
                           PH(nH,4) = eDeltoN
                           vxDeltoN = pxDeltoN/eDeltoN
                           vyDeltoN = pyDeltoN/eDeltoN
                           vzDeltoN = pzDeltoN/eDeltoN
                           gamDeltoN = hbarc/sqrt(1.-vxDeltoN**2.-
     .                          vyDeltoN**2.-vzDeltoN**2.)
                           XH(nH,1) = xBL+(pxDeltoN/eDeltoN)*
     .                          gamDeltoN/xmpro
                           XH(nH,2) = yBL+(pyDeltoN/eDeltoN)*
     .                          gamDeltoN/xmpro
                           XH(nH,3) = zBL+(pzDeltoN/eDeltoN)*
     .                          gamDeltoN/xmpro
                           XH(nH,4) = tBL+gamDeltoN/xmpro
                           
                           nH = nH + 1
                           KH(nH,2) = 1
                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25

                           PH(nH,1) = pxDeltopi
                           PH(nH,2) = pyDeltopi
                           PH(nH,3) = pzDeltopi
                           PH(nH,4) = eDeltopi
                           vxDeltopi = pxDeltopi/eDeltopi
                           vyDeltopi = pyDeltopi/eDeltopi
                           vzDeltopi = pzDeltopi/eDeltopi
                           gamDeltopi = hbarc/sqrt(1.-vxDeltopi**2.
     .                          -vyDeltopi**2.-vzDeltopi**2.)
                           XH(nH,1) = xBL+(pxDeltopi/eDeltopi)*
     .                          gamDeltopi/xmpi
                           XH(nH,2) = yBL+(pyDeltopi/eDeltopi)*
     .                          gamDeltopi/xmpi
                           XH(nH,3) = zBL+(pzDeltopi/eDeltopi)*
     .                          gamDeltopi/xmpi
                           XH(nH,4) = tBL+gamDeltopi/xmpi
                     
                           nchg(i) = 0
                           nchg(j) = 0
                           nchg(j2) = 0
                           kDel = kDel + 1
                           goto 30
                        endif

!............... N*(1440)-> N + pi .....................
                        if(proNg.ge.0.75) then
                           call twobody(xmN1,xmpro,xmpi,pxbaryon,
     .                          pybaryon,pzbaryon,pxN1toN,pyN1toN,
     .                          pzN1toN,eN1toN,pxN1topi,pyN1topi,
     .                          pzN1topi,eN1topi)
                     
                           nH = nH + 1
                           KH(nH,2) = 3

                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25
                           PH(nH,1) = pxN1toN
                           PH(nH,2) = pyN1toN
                           PH(nH,3) = pzN1toN
                           PH(nH,4) = eN1toN
                           vxN1toN = pxN1toN/eN1toN
                           vyN1toN = pyN1toN/eN1toN
                           vzN1toN = pzN1toN/eN1toN
                           gamN1toN = hbarc/sqrt(1.-vxN1toN**2.-
     .                          vyN1toN**2.-vzN1toN**2.)
                           XH(nH,1) = xBL+(pxN1toN/eN1toN)*
     .                          gamN1toN/xmpro
                           XH(nH,2) = yBL+(pyN1toN/eN1toN)*
     .                          gamN1toN/xmpro
                           XH(nH,3) = zBL+(pzN1toN/eN1toN)*
     .                          gamN1toN/xmpro
                           XH(nH,4) = tBL+gamN1toN/xmpro
                           
                           nH = nH + 1
                           KH(nH,2) = 1
                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25
                           PH(nH,1) = pxN1topi
                           PH(nH,2) = pyN1topi
                           PH(nH,3) = pzN1topi
                           PH(nH,4) = eN1topi
                           vxN1topi = pxN1topi/eN1topi
                           vyN1topi = pyN1topi/eN1topi
                           vzN1topi = pzN1topi/eN1topi
                           gamN1topi = hbarc/sqrt(1.-vxN1topi**2.
     .                          -vyN1topi**2.-vzN1topi**2.)
                           XH(nH,1) = xBL+(pxN1topi/eN1topi)*
     .                          gamN1topi/xmpi
                           XH(nH,2) = yBL+(pyN1topi/eN1topi)*
     .                          gamN1topi/xmpi
                           XH(nH,3) = zBL+(pzN1topi/eN1topi)*
     .                          gamN1topi/xmpi
                           XH(nH,4) = tBL+gamN1topi/xmpi                     
                           
                           nchg(i) = 0
                           nchg(j) = 0
                           nchg(j2) = 0
                           goto 30
                        endif
                     endif
                       
                     probBdec = ran()
!************ Recombine Delta and decay Delta -> proton pion *********************
!              if(probBdec.le.0.5) then
                     if(xinvB.ge.xmpro+xmpi.and.xinvB.le.xmpro+xmpi*2.) 
     .                    then
!     if(xinvB.le.xmpro+xmpi) goto 50
                        call twobody(xinvB,xmpro,xmpi,pxbaryon,pybaryon,
     .                       pzbaryon,pxDeltoN,pyDeltoN,pzDeltoN,
     .                       eDeltoN,pxDeltopi,pyDeltopi,pzDeltopi,
     .                       eDeltopi)
                  
                        nH = nH + 1
                        KH(nH,2) = 3
                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25
                        PH(nH,1) = pxDeltoN
                        PH(nH,2) = pyDeltoN
                        PH(nH,3) = pzDeltoN
                        PH(nH,4) = eDeltoN
                        vxDeltoN = pxDeltoN/eDeltoN
                        vyDeltoN = pyDeltoN/eDeltoN
                        vzDeltoN = pzDeltoN/eDeltoN
                        gamDeltoN = hbarc/sqrt(1.-vxDeltoN**2.-
     .                       vyDeltoN**2.-vzDeltoN**2.)
                        XH(nH,1) = xBL+(pxDeltoN/eDeltoN)*
     .                       gamDeltoN/xmpro
                        XH(nH,2) = yBL+(pyDeltoN/eDeltoN)*
     .                       gamDeltoN/xmpro
                        XH(nH,3) = zBL+(pzDeltoN/eDeltoN)*
     .                       gamDeltoN/xmpro
                        XH(nH,4) = tBL+gamDeltoN/xmpro
                        
                        nH = nH + 1
                        KH(nH,2) = 1


                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                   KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
               else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25





                        PH(nH,1) = pxDeltopi
                        PH(nH,2) = pyDeltopi
                        PH(nH,3) = pzDeltopi
                        PH(nH,4) = eDeltopi
                        vxDeltopi = pxDeltopi/eDeltopi
                        vyDeltopi = pyDeltopi/eDeltopi
                        vzDeltopi = pzDeltopi/eDeltopi
                        gamDeltopi = hbarc/sqrt(1.-vxDeltopi**2.
     .                       -vyDeltopi**2.-vzDeltopi**2.)
                        XH(nH,1) = xBL+(pxDeltopi/eDeltopi)*
     .                       gamDeltopi/xmpi
                        XH(nH,2) = yBL+(pyDeltopi/eDeltopi)*
     .                       gamDeltopi/xmpi
                        XH(nH,3) = zBL+(pzDeltopi/eDeltopi)*
     .                       gamDeltopi/xmpi
                        XH(nH,4) = tBL+gamDeltopi/xmpi
                  
                        nchg(i) = 0
                        nchg(j) = 0
                        nchg(j2) = 0
                        goto 30
                     endif

!........    Delta -> proton + pi + pi ........................
!              if(probBdec.gt.0.5) then
                     if(xinvB.gt.xmpro+xmpi*2..and.xinvB.le.1.8) then
!     if(xinvB.le.xmpro+xmpi+xmpi) goto 50
                  call threebody(xinvB,xmpro,xmpi,xmpi,pxbaryon,
     .                 pybaryon,pzbaryon,pxDeltoN,pyDeltoN,pzDeltoN,
     .                 eDeltoN,pxDeltopi1,pyDeltopi1,pzDeltopi1,
     .                 eDeltopi1,pxDeltopi2,pyDeltopi2,pzDeltopi2,
     .                 eDeltopi2)
                  
                  nH = nH + 1
                  KH(nH,2) = 3
                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                   KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
               else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25

                  PH(nH,1) = pxDeltoN
                  PH(nH,2) = pyDeltoN
                  PH(nH,3) = pzDeltoN
                  PH(nH,4) = eDeltoN
                  vxDeltoN = pxDeltoN/eDeltoN
                  vyDeltoN = pyDeltoN/eDeltoN
                  vzDeltoN = pzDeltoN/eDeltoN
                  gamDeltoN = hbarc/sqrt(1.-vxDeltoN**2.-vyDeltoN**2.
     .                 -vzDeltoN**2.)
                  XH(nH,1) = xBL+(pxDeltoN/eDeltoN)*gamDeltoN/xmpro
                  XH(nH,2) = yBL+(pyDeltoN/eDeltoN)*gamDeltoN/xmpro
                  XH(nH,3) = zBL+(pzDeltoN/eDeltoN)*gamDeltoN/xmpro
                  XH(nH,4) = tBL+gamDeltoN/xmpro
                  
                  nH = nH + 1
                  KH(nH,2) = 1

                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25

                  PH(nH,1) = pxDeltopi1
                  PH(nH,2) = pyDeltopi1
                  PH(nH,3) = pzDeltopi1
                  PH(nH,4) = eDeltopi1
                  vxDeltopi1 = pxDeltopi1/eDeltopi1
                  vyDeltopi1 = pyDeltopi1/eDeltopi1
                  vzDeltopi1 = pzDeltopi1/eDeltopi1
                  gamDeltopi1 = hbarc/sqrt(1.-vxDeltopi1**2.
     .                 -vyDeltopi1**2.-vzDeltopi1**2.)
                  XH(nH,1) = xBL+(pxDeltopi1/eDeltopi1)*gamDeltopi1/xmpi
                  XH(nH,2) = yBL+(pyDeltopi1/eDeltopi1)*gamDeltopi1/xmpi
                  XH(nH,3) = zBL+(pzDeltopi1/eDeltopi1)*gamDeltopi1/xmpi
                  XH(nH,4) = tBL+gamDeltopi1/xmpi
                  
                  nH = nH + 1
                  KH(nH,2) = 1

                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25

                  PH(nH,1) = pxDeltopi2
                  PH(nH,2) = pyDeltopi2
                  PH(nH,3) = pzDeltopi2
                  PH(nH,4) = eDeltopi2
                  vxDeltopi2 = pxDeltopi2/eDeltopi2
                  vyDeltopi2 = pyDeltopi2/eDeltopi2
                  vzDeltopi2 = pzDeltopi2/eDeltopi2
                  gamDeltopi2 = hbarc/sqrt(1.-vxDeltopi2**2.
     .                 -vyDeltopi2**2.-vzDeltopi2**2.)
                  XH(nH,1) = xBL+(pxDeltopi2/eDeltopi2)*gamDeltopi2/xmpi
                  XH(nH,2) = yBL+(pyDeltopi2/eDeltopi2)*gamDeltopi2/xmpi
                  XH(nH,3) = zBL+(pzDeltopi2/eDeltopi2)*gamDeltopi2/xmpi
                  XH(nH,4) = tBL+gamDeltopi2/xmpi
                  
                  nchg(i) = 0
                  nchg(j) = 0
                  nchg(j2) = 0
                  goto 30
               endif

!********** Delta -> proton + pi + pi + pi ******************
              if(xinvB.gt.1.8) then
                 xmTx(1) = xmpro
                 xmTx(2) = xmpi
                 xmTx(3) = xmpi
                 xmTx(4) = xmpi
                 pTx(1) = pxbaryon
                 pTx(2) = pybaryon
                 pTx(3) = pzbaryon
!                 print*, kev,xinvB
!                 xinvB = 10.
                 call fourbody(xinvB,xmTx,pTx,pdBx)
                 
                 nH = nH + 1
                 KH(nH,2) = 3
                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                 else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25


                 PH(nH,1) = pdBx(1,1)
                 PH(nH,2) = pdBx(1,2)
                 PH(nH,3) = pdBx(1,3)
                 PH(nH,4) = pdBx(1,4)
                 vxDeltoN = PH(nH,1)/PH(nH,4)
                 vyDeltoN = PH(nH,2)/PH(nH,4)
                 vzDeltoN = PH(nH,3)/PH(nH,4)
                 gamDeltoN = hbarc/sqrt(1.-vxDeltoN**2.-vyDeltoN**2.
     .                -vzDeltoN**2.)
                 XH(nH,1) = xBL+vxDeltoN*gamDeltoN/xmpro
                 XH(nH,2) = yBL+vyDeltoN*gamDeltoN/xmpro
                 XH(nH,3) = zBL+vzDeltoN*gamDeltoN/xmpro
                 XH(nH,4) = tBL+gamDeltoN/xmpro


                 nH = nH + 1
                 KH(nH,2) = 1

                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then 
                  KH(nH,1) = 2
                else 
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25

                 PH(nH,1) = pdBx(2,1)
                 PH(nH,2) = pdBx(2,2)
                 PH(nH,3) = pdBx(2,3)
                 PH(nH,4) = pdBx(2,4)
                 vxDeltopi1 = PH(nH,1)/PH(nH,4)
                 vyDeltopi1 = PH(nH,2)/PH(nH,4)
                 vzDeltopi1 = PH(nH,3)/PH(nH,4)
                 gamDeltopi1 = hbarc/sqrt(1.-vxDeltopi1**2.-
     .                vyDeltopi1**2.-vzDeltopi1**2.)
                 XH(nH,1) = xBL+vxDeltopi1*gamDeltopi1/xmpi
                 XH(nH,2) = yBL+vyDeltopi1*gamDeltopi1/xmpi
                 XH(nH,3) = zBL+vzDeltopi1*gamDeltopi1/xmpi
                 XH(nH,4) = tBL+gamDeltopi1/xmpi

                 nH = nH + 1
                 KH(nH,2) = 1
                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                 else                                              !changed by wenbin 2018.11.25
!                  if(Kq(i,1)*Kq(j,1)*Kq(j2,1).eq.0) KH(nH,1) = 1 !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25


                 PH(nH,1) = pdBx(3,1)
                 PH(nH,2) = pdBx(3,2)
                 PH(nH,3) = pdBx(3,3)
                 PH(nH,4) = pdBx(3,4)
                 vxDeltopi2 = PH(nH,1)/PH(nH,4)
                 vyDeltopi2 = PH(nH,2)/PH(nH,4)
                 vzDeltopi2 = PH(nH,3)/PH(nH,4)
                 gamDeltopi2 = hbarc/sqrt(1.-vxDeltopi2**2.-
     .                vyDeltopi2**2.-vzDeltopi2**2.)
                 XH(nH,1) = xBL+vxDeltopi2*gamDeltopi2/xmpi
                 XH(nH,2) = yBL+vyDeltopi2*gamDeltopi2/xmpi
                 XH(nH,3) = zBL+vzDeltopi2*gamDeltopi2/xmpi
                 XH(nH,4) = tBL+gamDeltopi2/xmpi
                 
                 nH = nH + 1
                 KH(nH,2) = 1
                 if(Kq(i,1)*Kq(j,1)*Kq(j2,1).ne.0) then             !changed by wenbin 2018.11.25
                  KH(nH,1) = 2                                      !changed by wenbin 2018.11.25
                 else                                              !changed by wenbin 2018.11.25
                   if(Kq(i,1).eq. 0 .and.(Kq(j,1).eq. 0 .and.        !changed by wenbin 2018.11.25
     .                Kq(j2,1).eq. 0)) then                          !changed by wenbin 2018.11.25
                      KH(nH,1) = 0                                   !changed by wenbin 2018.11.25
                   else                                            !changed by wenbin 2018.11.25
                      KH(nH,1) = 1                                    !changed by wenbin 2018.11.25
                   endif                                           !changed by wenbin 2018.11.25
                 endif                                           !changed by wenbin 2018.11.25
                 PH(nH,1) = pdBx(4,1)
                 PH(nH,2) = pdBx(4,2)
                 PH(nH,3) = pdBx(4,3)
                 PH(nH,4) = pdBx(4,4)
                 vxDeltopi3 = PH(nH,1)/PH(nH,4)
                 vyDeltopi3 = PH(nH,2)/PH(nH,4)
                 vzDeltopi3 = PH(nH,3)/PH(nH,4)
                 gamDeltopi3 = hbarc/sqrt(1.-vxDeltopi3**2.-
     .                vyDeltopi3**2.-vzDeltopi3**2.)
                 XH(nH,1) = xBL+vxDeltopi3*gamDeltopi3/xmpi
                 XH(nH,2) = yBL+vyDeltopi3*gamDeltopi3/xmpi
                 XH(nH,3) = zBL+vzDeltopi3*gamDeltopi3/xmpi
                 XH(nH,4) = tBL+gamDeltopi3/xmpi

                 nchg(i) = 0
                 nchg(j) = 0
                 nchg(j2) = 0
                 goto 30
              endif
           endif
 50        CONTINUE
        ENDDO
!****** recombination of quark and antiquark *******************
      else
C........NOT consider th-th recombinatoin ......................               
!         if(Kq(i,1).eq.0.and.Kq(j,1).eq.0) 
!     .     then !goto 40 !changed by wenbin 2018.11.23 !TTT
!                pt1square=px1*px1+py1*py1
!                pt2square=px2*px2+py2*py2
!                if(pt1square.lt.pttcutsquare)goto 40
!                if(pt2square.lt.pttcutsquare)goto 40
!         endif
C........NOT consider hard-hard recombinatoin ......................
!         if(Kq(i,1)*Kq(j,1).gt.0) goto 40 !changed by wenbin2018.11.23 !TTT
!*********momentum of the recombined mesons ****************
 
!         if((Kq(i,1).eq.0).and.(Kq(j,1).eq.0))goto 40!don't consider  the th-th-th recombination !Wenbin

         pxpion = px1 + px2
         pypion = py1 + py2           
         pzpion = pz1 + pz2
         epion = e1 + e2
         
         betaMx = pxpion/(e1+e2)
         betaMy = pypion/(e1+e2)
         betaMz = pzpion/(e1+e2)
         betaM2 = betaMx**2.+betaMy**2.+betaMz**2.

C......speeds of the Lab frame in the center-of-mass frame .......             
         betaMLx = -betaMx
         betaMLy = -betaMy
         betaMLz = -betaMz
         
         call lorentzboost(betaMx,betaMy,betaMz,e1,px1,py1,pz1,e1Mcm,
     .        px1Mcm,py1Mcm,pz1Mcm)
         call lorentzboost(betaMx,betaMy,betaMz,e2,px2,py2,pz2,e2Mcm,
     .        px2Mcm,py2Mcm,pz2Mcm)
         call lorentzboost(betaMx,betaMy,betaMz,t1,x1,y1,z1,t1Mcm,
     .        x1Mcm,y1Mcm,z1Mcm)
         call lorentzboost(betaMx,betaMy,betaMz,t2,x2,y2,z2,t2Mcm,
     .        x2Mcm,y2Mcm,z2Mcm)

!********* speeds of quarks and antiquarks in their center-of-mass frame ********************
         vx1Mcm = px1Mcm/e1Mcm
         vy1Mcm = py1Mcm/e1Mcm
         vz1Mcm = pz1Mcm/e1Mcm
         
         vx2Mcm = px2Mcm/e2Mcm
         vy2Mcm = py2Mcm/e2Mcm
         vz2Mcm = pz2Mcm/e2Mcm
         
         dy2 = (t2-t1)**2. - (x2-x1)**2. - (y2-y1)**2. - (z2-z1)**2.
         dyp1 = (t2-t1)*e1 - (x2-x1)*px1 - (y2-y1)*py1 - (z2-z1)*pz1
         dyp2 = (t2-t1)*e2 - (x2-x1)*px2 - (y2-y1)*py2 - (z2-z1)*pz2
         p12 = e1*e2 - px1*px2 - py1*py2 - pz1*pz2
         
         betax1 = px1/e1
         betay1 = py1/e1
         betaz1 = pz1/e1
         beta12 = betax1**2.+betay1**2.+betaz1**2.
         gam1 = 1./sqrt(1.-beta12)
         betax2 = px2/e2
         betay2 = py2/e2
         betaz2 = pz2/e2
         beta22 = betax2**2.+betay2**2.+betaz2**2.
         gam2 = 1./sqrt(1.-beta22)
         
         xm12 = e1**2.-px1**2.-py1**2.-pz1**2.
         xm22 = e2**2.-px2**2.-py2**2.-pz2**2.
         
         dyu1 = gam1*((t2-t1)-(x2-x1)*betax1-(y2-y1)*betay1-
     .        (z2-z1)*betaz1)
         dyu2 =  gam2*((t2-t1)-(x2-x1)*betax2-(y2-y1)*betay2-
     .        (z2-z1)*betaz2)
         u12 = gam1*gam2*(1.-betax1*betax2-betay1*betay2-
     .        betaz1*betaz2)
!********* t1new, t2new which are moments in the rest frame when two 
           
!***********particles are closest in cm frame ************************
         t1new = (-dyu1+dyu2*u12)*gam1/(u12**2.-1.)+t1
         t2new = (dyu2-dyu1*u12)*gam2/(u12**2.-1.)+t2
!*********** 90 % for that ********************************           
         tMaxMcm = max(t1Mcm,t2Mcm)
         x1nMcm = x1Mcm + vx1Mcm*(tMaxMcm-t1Mcm)
         y1nMcm = y1Mcm + vy1Mcm*(tMaxMcm-t1Mcm)
         z1nMcm = z1Mcm + vz1Mcm*(tMaxMcm-t1Mcm)
         
         x2nMcm = x2Mcm + vx2Mcm*(tMaxMcm-t2Mcm)
         y2nMcm = y2Mcm + vy2Mcm*(tMaxMcm-t2Mcm)
         z2nMcm = z2Mcm + vz2Mcm*(tMaxMcm-t2Mcm)
         
         dx2M = (x1nMcm-x2nMcm)**2.
         dy2M = (y1nMcm-y2nMcm)**2.
         dz2M = (z1nMcm-z2nMcm)**2.
         
         dr2M = ((x1nMcm-x2nMcm)**2.+(y1nMcm-y2nMcm)**2.+
     .        (z1nMcm-z2nMcm)**2.)
!********10% for that **************************************
         if(t1new.ge.t1.and.t2new.ge.t2) then
            dr2M = (-dy2-((dyp1**2.)*xm22+(dyp2**2.)*xm12-2.*dyp1*dyp2*
     .           p12)/(p12**2.-xm12*xm22))
         endif
         
         xcmM = (x1nMcm*Pq(i,5)+x2nMcm*Pq(j,5))/(Pq(i,5)+Pq(j,5))
         ycmM = (y1nMcm*Pq(i,5)+y2nMcm*Pq(j,5))/(Pq(i,5)+Pq(j,5))
         zcmM = (z1nMcm*Pq(i,5)+z2nMcm*Pq(j,5))/(Pq(i,5)+Pq(j,5))
         tcmM = tMaxMcm

C........... position and time of the recombined meson in the Lab frame ....    
         call lorentzboost(betaMLx,betaMLy,betaMLz,tcmM,xcmM,ycmM,
     .        zcmM,tML,xML,yML,zML)


!******** invariant mass of quark and antiquark **********
        xinvm = e1Mcm + e2Mcm
        
        dkx2M = (px1Mcm-px2Mcm)**2./4.
        dky2M = (py1Mcm-py2Mcm)**2./4.
        dkz2M = (pz1Mcm-pz2Mcm)**2./4.
        dk2M = dkx2M + dky2M + dkz2M

        drk2M = ((x1nMcm-x2nMcm)*(px1Mcm-px2Mcm)/2.+(y1nMcm-y2nMcm)
     .      *(py1Mcm-py2Mcm)/2.+(z1nMcm-z2nMcm)*(pz1Mcm-pz2Mcm)/2.)**2.
        
        sumWig3D = 0.
!        if(Kq(i,2)*Kq(j,2).eq.-6.or.Kq(i,2)*Kq(j,2).eq.-3) Sig2=SigK2 !changed by wenbin 2018.11.23     
!****** begin added by wenbin 2018.11.23 *********
        if(Kq(i,2)*Kq(j,2).eq.-6.or.Kq(i,2)*Kq(j,2).eq.-3)Sig2=SigK2
        if (Kq(i,2)*Kq(j,2).eq. -9)Sig2=Sigphi2
!****** end added by wenbin 2018.11.23 *********
!****** Sum of 3D w.w functions up to nlev ***********************
        ux = 0.5*(dx2M/Sig2+dkx2M*Sig2/hbarc2)
        uy = 0.5*(dy2M/Sig2+dky2M*Sig2/hbarc2)
        uz = 0.5*(dz2M/Sig2+dkz2M*Sig2/hbarc2)
        uchi = ux+uy+uz
        if(uchi.gt. cccut)goto 40!added by wenbin

!        nchi = anint(uchi)
!        if(nchi.ge.125.) goto 40
!********** maximum and minimum for Poisson Distribution ..........
!        if(nchi.ge.1) then
!           nmax=nchi+anint(sqrt(uchi))*2
!           if(uchi.lt.1.) nmax=2
!        else
!           nmax = 2
!        endif

!        if(nchi.ge.5) then
!           nmin=nchi-1-int(sqrt(uchi))*2
!        else
!           nmin=0
!        endif
        nmin = 0
        nmax = maxf
        if(Kq(i,2)*Kq(j,2).eq.-6.or.Kq(i,2)*Kq(j,2).eq.-3)nmax = 10 !Kaon's highest excited state
        do ichi=nmin, nmax
           pmf=(uchi**ichi)*exp(-uchi)/factor(ichi)
           if(ichi.ge.125) pmf=exp(float(ichi)-uchi)*(uchi/float(ichi))
     .          **ichi/sqrt(2*pi*ichi)                  
           sumWig3D=sumWig3D+pmf
        enddo
!......Ground w. Wig functions .....        
!        wigpig = exp(-dr2M/(2.*Sig2)-dk2M*Sig2/hbarc2) !changed by wenbin 2018.11.23
        wigpig=exp(-uchi)!added by wenbin 2018.11.23
        
        wigpig=wigpig/9.0 !1/9.0 for color degree of freedom! added by wenbin 2019.03.10
        sumWig3D=sumWig3D/9.0 !1/9.0 for color degree of freedom! added by wenbin 2019.03.10
        wpair = sumWig3D
!        wpair = sumWig3D*ZKfactor !modify to account the colour degeneracy !changed by wenbin 2019.02.24
        WigK = wpair
        Wigpion = wpair
        Wigphi = wigpig !added by wenbin 
!        print*,Wigpion
!**********Calculate rates for 2-body, 3-body, and 4-body decays **********************
!*******************ref. to Phys. Rev. C72, 064901 (AMPT paper) *************************
!     Rpi = 2.0
        Rpi = sqrt(dr2M)/2.
        decpro2 = ((xinvm*Rpi/hbarc)**3./(6.*pi**2.))**2.*
     .       factor(4*2-4)*(2*2-1)/(factor(2*2-1)**2.*factor(3*2-4))
        decpro3 = ((xinvm*Rpi/hbarc)**3./(6.*pi**2.))**3.*
     .       factor(4*3-4)*(2*3-1)/(factor(2*3-1)**2.*factor(3*3-4))
        decpro4 = ((xinvm*Rpi/hbarc)**3./(6.*pi**2.))**4.*
     .       factor(4*4-4)*(2.*4.-1.)/(factor(2*4-1)**2.*factor(3*4-4))
        xnorm = decpro2 + decpro3 + decpro4
!     print*, "xnorm=",xnorm
        decpro2 = decpro2/xnorm
        decpro3 = decpro3/xnorm
        decpro4 = decpro4/xnorm
!     print*, "decpro", decpro2,decpro3,decpro4
        
!---------------------------------------------------------------------
!     Recombine K mesons                                 *
!************recombine K* and Decay K*-> K + pi ***************************
        prob2 = ran()
        if(prob2.le.WigK.and.(Kq(i,2)*Kq(j,2).eq.-6.or.
     .       Kq(i,2)*Kq(j,2).eq.-3)) then
           
!********ground state ******************************
           if(prob2.le.wigpig) then
!............direct kaon and K*->K +pion in n=0 ..........................
              if(xinvm.ge.xmK+xmpi)then  !!>>>>>>> Koan decay ground state wenbin 2020.6.18<<<<<<<<<<<
              if(ran().le.0.25) then
                 nH = nH + 1
                 KH(nH,2) = 2
                 PH(nH,1) = pxpion
                 PH(nH,2) = pypion
                 PH(nH,3) = pzpion
                 PH(nH,4) = sqrt(pxpion**2.+pypion**2.+pzpion**2.+
     .                xmK**2.)
                 XH(nH,1) = xML
                 XH(nH,2) = yML
                 XH(nH,3) = zML
                 XH(nH,4) = tML
                 
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
                 
                 nchg(i) = 0
                 nchg(j) = 0
                 goto 30
                 
              else
                 call twobody(xmKs,xmK,xmpi,pxpion,pypion,pzpion,
     .                pxKstoKg,pyKstoKg,pzKstoKg,eKstoKg,pxKstopi,
     .                pyKstopi,pzKstopi,eKstopi)
                 
                 nH = nH + 1
                 KH(nH,2) = 2

!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
                 
                 PH(nH,1) = pxKstoKg
                 PH(nH,2) = pyKstoKg
                 PH(nH,3) = pzKstoKg
                 PH(nH,4) = eKstoKg
                 vxKstoK = pxKstoKg/eKstoK
                 vyKstoK = pyKstoKg/eKstoK
                 vzKstoK = pzKstoKg/eKstoK
                 gamKstoK = hbarc/sqrt(1.-vxKstoK**2.-vyKstoK**2.-
     .                vzKstoK**2.)
                 XH(nH,1) = xML+vxKstoK*gamKstoK/xmK
                 XH(nH,2) = yML+vyKstoK*gamKstoK/xmK
                 XH(nH,3) = zML+vzKstoK*gamKstoK/xmK
                 XH(nH,4) = tML + gamKstoK/xmK
                 
                 nH = nH + 1
                 KH(nH,2) = 1

!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
                 
                 PH(nH,1) = pxKstopi
                 PH(nH,2) = pyKstopi
                 PH(nH,3) = pzKstopi
                 PH(nH,4) = eKstopi
                 vxKstopi = pxKstopi/eKstopi
                 vyKstopi = pyKstopi/eKstopi
                 vzKstopi = pzKstopi/eKstopi
                 gamKstopi = hbarc/sqrt(1.-vxKstopi**2.-vyKstopi**2.
     .                -vzKstopi**2.)
                 XH(nH,1) = xML+vxKstopi*gamKstopi/xmpi
                 XH(nH,2) = yML+vyKstopi*gamKstopi/xmpi
                 XH(nH,3) = zML+vzKstopi*gamKstopi/xmpi
                 XH(nH,4) = tML + gamKstopi/xmpi
                 
                 nchg(i) = 0
                 nchg(j) = 0
                 goto 30
                 
              endif
           else !>>>>>>>>>>>>>>>>>>. ground state Kaon wenbin 2020.6.18 <<<<<<<<<<<<<<<<<<<<<<<<<
  
                 nH = nH + 1
                 KH(nH,2) = 2
                 PH(nH,1) = pxpion
                 PH(nH,2) = pypion
                 PH(nH,3) = pzpion
                 PH(nH,4) = sqrt(pxpion**2.+pypion**2.+pzpion**2.+
     .                xmK**2.)
                 XH(nH,1) = xML
                 XH(nH,2) = yML
                 XH(nH,3) = zML
                 XH(nH,4) = tML

                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25

                 nchg(i) = 0
                 nchg(j) = 0
                 goto 30

           endif  !>>>>>>>>>>>>wenbin K decay groud states Wenbin !20206.18 <<<<<<<<<<<<<<<<<<<<<<<<<<,
           endif     
!========= begin excited states ===============wenbin 
           probKdec = ran()
           if(xinvm.ge.3.0*xmK)then !>>>>>>>>>>>>>>>>>>>>> all chanels <<<<<<<<<<<<<<<<<<< !wenbin 2020.6.18
!************ Recombine Ks  and Decay K* -> K pi *************************************
           if(probKdec.le.0.33) then
              if(xinvm.le.xmK+xmpi) xinvm = xmKs
              call twobody(xinvm,xmK,xmpi,pxpion,pypion,pzpion,
     .             pxKstoK,pyKstoK,pzKstoK,eKstoK,pxKstopi,pyKstopi,
     .             pzKstopi,eKstopi)
              
              nH = nH + 1
              KH(nH,2) = 2


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstoK
              PH(nH,2) = pyKstoK
              PH(nH,3) = pzKstoK
              PH(nH,4) = eKstoK
              vxKstoK = pxKstoK/eKstoK
              vyKstoK = pyKstoK/eKstoK
              vzKstoK = pzKstoK/eKstoK
              gamKstoK = hbarc/sqrt(1.-vxKstoK**2.-vyKstoK**2.-
     .             vzKstoK**2.)
              XH(nH,1) = xML+vxKstoK*gamKstoK/xmK
              XH(nH,2) = yML+vyKstoK*gamKstoK/xmK
              XH(nH,3) = zML+vzKstoK*gamKstoK/xmK
              XH(nH,4) = tML + gamKstoK/xmK
              
              nH = nH + 1
              KH(nH,2) = 1


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstopi
              PH(nH,2) = pyKstopi
              PH(nH,3) = pzKstopi
              PH(nH,4) = eKstopi
              vxKstopi = pxKstopi/eKstopi
              vyKstopi = pyKstopi/eKstopi
              vzKstopi = pzKstopi/eKstopi
              gamKstopi = hbarc/sqrt(1.-vxKstopi**2.-vyKstopi**2.
     .             -vzKstopi**2.)
              XH(nH,1) = xML+vxKstopi*gamKstopi/xmpi
              XH(nH,2) = yML+vyKstopi*gamKstopi/xmpi
              XH(nH,3) = zML+vzKstopi*gamKstopi/xmpi
              XH(nH,4) = tML + gamKstopi/xmpi
              
              nchg(i) = 0
              nchg(j) = 0
              goto 30
           endif


!*********** Decay into K -> K + K + K ***********************************
           if(probKdec.gt.0.33.and.probKdec.le.0.66) then
              if(xinvm.le.xmK*3.) xinvm = 3.*xmK
              call threebody(xinvm,xmK,xmK,xmK,pxpion,pypion,pzpion,
     .             pxKtoK1,pyKtoK1,pzKtoK1,eKtoK1,pxKtoK2,pyKtoK2,
     .             pzKtoK2,eKtoK2,pxKtoK3,pyKtoK3,pzKtoK3,eKtoK3)
              
              nH = nH + 1
              KH(nH,2) = 2
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKtoK1
              PH(nH,2) = pyKtoK1
              PH(nH,3) = pzKtoK1
              PH(nH,4) = eKtoK1
              vxKtoK1 = pxKtoK1/eKtoK1
              vyKtoK1 = pyKtoK1/eKtoK1
              vzKtoK1 = pzKtoK1/eKtoK1
              gamKtoK1 = hbarc/sqrt(1.-vxKtoK1**2.-vyKtoK1**2.-
     .             vzKtoK1**2.)
              XH(nH,1) = xML+vxKtoK1*gamKtoK1/xmK
              XH(nH,2) = yML+vyKtoK1*gamKtoK1/xmK
              XH(nH,3) = zML+vzKtoK1*gamKtoK1/xmK
              XH(nH,4) = tML + gamKtoK1/xmK

              nH = nH + 1
              KH(nH,2) = 2

!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              




              PH(nH,1) = pxKtoK2
              PH(nH,2) = pyKtoK2
              PH(nH,3) = pzKtoK2
              PH(nH,4) = eKtoK2
              vxKtoK2 = pxKtoK2/eKtoK2
              vyKtoK2 = pyKtoK2/eKtoK2
              vzKtoK2 = pzKtoK2/eKtoK2
              gamKtoK2 = hbarc/sqrt(1.-vxKtoK2**2.-vyKtoK2**2.-
     .             vzKtoK2**2.)
              XH(nH,1) = xML+vxKtoK2*gamKtoK2/xmK
              XH(nH,2) = yML+vyKtoK2*gamKtoK2/xmK
              XH(nH,3) = zML+vzKtoK2*gamKtoK2/xmK
              XH(nH,4) = tML + gamKtoK2/xmK

              nH = nH + 1
              KH(nH,2) = 2

!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKtoK3
              PH(nH,2) = pyKtoK3
              PH(nH,3) = pzKtoK3
              PH(nH,4) = eKtoK3
              vxKtoK3 = pxKtoK3/eKtoK3
              vyKtoK3 = pyKtoK3/eKtoK3
              vzKtoK3 = pzKtoK3/eKtoK3
              gamKtoK3 = hbarc/sqrt(1.-vxKtoK3**2.-vyKtoK3**2.-
     .             vzKtoK3**2.)
              XH(nH,1) = xML+vxKtoK3*gamKtoK3/xmK
              XH(nH,2) = yML+vyKtoK3*gamKtoK3/xmK
              XH(nH,3) = zML+vzKtoK3*gamKtoK3/xmK
              XH(nH,4) = tML + gamKtoK3/xmK
              
              nchg(i) = 0
              nchg(j) = 0
              goto 30
           endif

!***** Decay into K + pi + pi ********************************* 
           if(probKdec.gt.0.66) then
              if(xinvm.le.xmK+xmpi*2.) xinvm = xmK + xmpi*2.
              call threebody(xinvm,xmK,xmpi,xmpi,pxpion,pypion,pzpion,
     .             pxKstoK,pyKstoK,pzKstoK,eKstoK,pxKstopi1,pyKstopi1,
     .             pzKstopi1,eKstopi1,pxKstopi2,pyKstopi2,pzKstopi2,
     .             eKstopi2)
                
              nH = nH + 1
              KH(nH,2) = 2


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstoK
              PH(nH,2) = pyKstoK
              PH(nH,3) = pzKstoK
              PH(nH,4) = eKstoK
              vxKstoK = pxKstoK/eKstoK
              vyKstoK = pyKstoK/eKstoK
              vzKstoK = pzKstoK/eKstoK
              gamKstoK = hbarc/sqrt(1.-vxKstoK**2.-vyKstoK**2.-
     .             vzKstoK**2.)
              XH(nH,1) = xML+vxKstoK*gamKstoK/xmK
              XH(nH,2) = yML+vyKstoK*gamKstoK/xmK
              XH(nH,3) = zML+vzKstoK*gamKstoK/xmK
              XH(nH,4) = tML + gamKstoK/xmK
              
              nH = nH + 1
              KH(nH,2) = 1


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstopi1
              PH(nH,2) = pyKstopi1
              PH(nH,3) = pzKstopi1
              PH(nH,4) = eKstopi1
              vxKstopi1 = pxKstopi1/eKstopi1
              vyKstopi1 = pyKstopi1/eKstopi1
              vzKstopi1 = pzKstopi1/eKstopi1
              gamKstopi1 = hbarc/sqrt(1.-vxKstopi1**2.-vyKstopi1**2.-                              
     .             vzKstopi1**2.)
              XH(nH,1) = xML+vxKstopi1*gamKstopi1/xmpi
              XH(nH,2) = yML+vyKstopi1*gamKstopi1/xmpi
              XH(nH,3) = zML+vzKstopi1*gamKstopi1/xmpi
              XH(nH,4) = tML + gamKstopi1/xmpi

              nH = nH + 1
!              KH(nH,2) = 2!changed by wenbin 2018.11.23
              KH(nH,2) = 1!changed by wenbin 2018.11.23

!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstopi2
              PH(nH,2) = pyKstopi2
              PH(nH,3) = pzKstopi2
              PH(nH,4) = eKstopi2
              vxKstopi2 = pxKstopi2/eKstopi2
              vyKstopi2 = pyKstopi2/eKstopi2
              vzKstopi2 = pzKstopi2/eKstopi2
              gamKstopi2 = hbarc/sqrt(1.-vxKstopi2**2.-vyKstopi2**2.-                              
     .             vzKstopi2**2.)
              XH(nH,1) = xML+vxKstopi2*gamKstopi2/xmpi
              XH(nH,2) = yML+vyKstopi2*gamKstopi2/xmpi
              XH(nH,3) = zML+vzKstopi2*gamKstopi2/xmpi
              XH(nH,4) = tML + gamKstopi2/xmpi

              nchg(i) = 0
              nchg(j) = 0
              goto 30
           endif
           else !!>>>>>>>>>>>>>>>>>>>>> all chanels <<<<<<<<<<<<<<<<<<< 
           if(xinvm.lt.3.0*xmK.and.xinvm.ge.xmk+2.0*xmpi)then !>>>>>>>>>>>>>>>>>>>>> two chanels <<<<<<<<<<<<<<<<<<< !wenbin 2020.6.18
!************* K*->K+pi *************** !wenbin 2020.06.18
             if(probKdec.le.0.50) then
              if(xinvm.le.xmK+xmpi) xinvm = xmKs
              call twobody(xinvm,xmK,xmpi,pxpion,pypion,pzpion,
     .             pxKstoK,pyKstoK,pzKstoK,eKstoK,pxKstopi,pyKstopi,
     .             pzKstopi,eKstopi)
              
              nH = nH + 1
              KH(nH,2) = 2


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstoK
              PH(nH,2) = pyKstoK
              PH(nH,3) = pzKstoK
              PH(nH,4) = eKstoK
              vxKstoK = pxKstoK/eKstoK
              vyKstoK = pyKstoK/eKstoK
              vzKstoK = pzKstoK/eKstoK
              gamKstoK = hbarc/sqrt(1.-vxKstoK**2.-vyKstoK**2.-
     .             vzKstoK**2.)
              XH(nH,1) = xML+vxKstoK*gamKstoK/xmK
              XH(nH,2) = yML+vyKstoK*gamKstoK/xmK
              XH(nH,3) = zML+vzKstoK*gamKstoK/xmK
              XH(nH,4) = tML + gamKstoK/xmK
              
              nH = nH + 1
              KH(nH,2) = 1


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstopi
              PH(nH,2) = pyKstopi
              PH(nH,3) = pzKstopi
              PH(nH,4) = eKstopi
              vxKstopi = pxKstopi/eKstopi
              vyKstopi = pyKstopi/eKstopi
              vzKstopi = pzKstopi/eKstopi
              gamKstopi = hbarc/sqrt(1.-vxKstopi**2.-vyKstopi**2.
     .             -vzKstopi**2.)
              XH(nH,1) = xML+vxKstopi*gamKstopi/xmpi
              XH(nH,2) = yML+vyKstopi*gamKstopi/xmpi
              XH(nH,3) = zML+vzKstopi*gamKstopi/xmpi
              XH(nH,4) = tML + gamKstopi/xmpi
              
              nchg(i) = 0
              nchg(j) = 0
              goto 30
           endif
!***** Decay into K + pi + pi ********************************* 
           if(probKdec.gt.0.5) then
              !if(xinvm.le.xmK+xmpi*2.) xinvm = xmK + xmpi*2.
              call threebody(xinvm,xmK,xmpi,xmpi,pxpion,pypion,pzpion,
     .             pxKstoK,pyKstoK,pzKstoK,eKstoK,pxKstopi1,pyKstopi1,
     .             pzKstopi1,eKstopi1,pxKstopi2,pyKstopi2,pzKstopi2,
     .             eKstopi2)
                
              nH = nH + 1
              KH(nH,2) = 2


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstoK
              PH(nH,2) = pyKstoK
              PH(nH,3) = pzKstoK
              PH(nH,4) = eKstoK
              vxKstoK = pxKstoK/eKstoK
              vyKstoK = pyKstoK/eKstoK
              vzKstoK = pzKstoK/eKstoK
              gamKstoK = hbarc/sqrt(1.-vxKstoK**2.-vyKstoK**2.-
     .             vzKstoK**2.)
              XH(nH,1) = xML+vxKstoK*gamKstoK/xmK
              XH(nH,2) = yML+vyKstoK*gamKstoK/xmK
              XH(nH,3) = zML+vzKstoK*gamKstoK/xmK
              XH(nH,4) = tML + gamKstoK/xmK
              
              nH = nH + 1
              KH(nH,2) = 1


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstopi1
              PH(nH,2) = pyKstopi1
              PH(nH,3) = pzKstopi1
              PH(nH,4) = eKstopi1
              vxKstopi1 = pxKstopi1/eKstopi1
              vyKstopi1 = pyKstopi1/eKstopi1
              vzKstopi1 = pzKstopi1/eKstopi1
              gamKstopi1 = hbarc/sqrt(1.-vxKstopi1**2.-vyKstopi1**2.-                              
     .             vzKstopi1**2.)
              XH(nH,1) = xML+vxKstopi1*gamKstopi1/xmpi
              XH(nH,2) = yML+vyKstopi1*gamKstopi1/xmpi
              XH(nH,3) = zML+vzKstopi1*gamKstopi1/xmpi
              XH(nH,4) = tML + gamKstopi1/xmpi

              nH = nH + 1
!              KH(nH,2) = 2!changed by wenbin 2018.11.23
              KH(nH,2) = 1!changed by wenbin 2018.11.23

!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstopi2
              PH(nH,2) = pyKstopi2
              PH(nH,3) = pzKstopi2
              PH(nH,4) = eKstopi2
              vxKstopi2 = pxKstopi2/eKstopi2
              vyKstopi2 = pyKstopi2/eKstopi2
              vzKstopi2 = pzKstopi2/eKstopi2
              gamKstopi2 = hbarc/sqrt(1.-vxKstopi2**2.-vyKstopi2**2.-                              
     .             vzKstopi2**2.)
              XH(nH,1) = xML+vxKstopi2*gamKstopi2/xmpi
              XH(nH,2) = yML+vyKstopi2*gamKstopi2/xmpi
              XH(nH,3) = zML+vzKstopi2*gamKstopi2/xmpi
              XH(nH,4) = tML + gamKstopi2/xmpi

              nchg(i) = 0
              nchg(j) = 0
              goto 30
           endif
 
           else !>>>>>>>>>>>>>>>>>>>>> two chanels <<<<<<<<<<<<<<<<<<< !wenbin 2020.6.18
           if(xinvm.lt.xmk+2.0*xmpi)then !>>>>>>>>>>>>>>>>>>>>> one chanels <<<<<<<<<<<<<<<<<<< !wenbin 2020.6.18
           if(probKdec.le.1.0) then
              if(xinvm.le.xmK+xmpi) xinvm = xmKs
              call twobody(xinvm,xmK,xmpi,pxpion,pypion,pzpion,
     .             pxKstoK,pyKstoK,pzKstoK,eKstoK,pxKstopi,pyKstopi,
     .             pzKstopi,eKstopi)
              
              nH = nH + 1
              KH(nH,2) = 2


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstoK
              PH(nH,2) = pyKstoK
              PH(nH,3) = pzKstoK
              PH(nH,4) = eKstoK
              vxKstoK = pxKstoK/eKstoK
              vyKstoK = pyKstoK/eKstoK
              vzKstoK = pzKstoK/eKstoK
              gamKstoK = hbarc/sqrt(1.-vxKstoK**2.-vyKstoK**2.-
     .             vzKstoK**2.)
              XH(nH,1) = xML+vxKstoK*gamKstoK/xmK
              XH(nH,2) = yML+vyKstoK*gamKstoK/xmK
              XH(nH,3) = zML+vzKstoK*gamKstoK/xmK
              XH(nH,4) = tML + gamKstoK/xmK
              
              nH = nH + 1
              KH(nH,2) = 1


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxKstopi
              PH(nH,2) = pyKstopi
              PH(nH,3) = pzKstopi
              PH(nH,4) = eKstopi
              vxKstopi = pxKstopi/eKstopi
              vyKstopi = pyKstopi/eKstopi
              vzKstopi = pzKstopi/eKstopi
              gamKstopi = hbarc/sqrt(1.-vxKstopi**2.-vyKstopi**2.
     .             -vzKstopi**2.)
              XH(nH,1) = xML+vxKstopi*gamKstopi/xmpi
              XH(nH,2) = yML+vyKstopi*gamKstopi/xmpi
              XH(nH,3) = zML+vzKstopi*gamKstopi/xmpi
              XH(nH,4) = tML + gamKstopi/xmpi
              
              nchg(i) = 0
              nchg(j) = 0
              goto 30
           endif


           endif!>>>>>>>>>>>>>>>>>>>>> one chanels <<<<<<<<<<<<<<<<<<< !wenbin 2020.6.18
           endif!>>>>>>>>>>>>>>>>>>>>> two chanels <<<<<<<<<<<<<<<<<<< !wenbin 2020.6.18  
           endif!>>>>>>>>>>>>>>>>>>>>> all chanels <<<<<<<<<<<<<<<<<<< !wenbin 2020.6.18
        endif
       
!-----------------------------------------------------------------------------
!     Recombine mesons formed from light quarks                          *
!-----------------------------------------------------------------------------
        prob1 = ran()
!********** recombine pion ******************************************************
        if((Kq(i,2)*Kq(j,2).eq.-2.or.Kq(i,2)*Kq(j,2).eq.-4.or.
     .       Kq(i,2)*Kq(j,2).eq.-1).and.(prob1.le.Wigpion)) then

           probdec = ran()
!*************** direction pion and rho in the ground statte ******************
           if(prob1.le.wigpig) then
!*****************Direct pion **************
              if(ran().le.0.25) then
                 nH = nH + 1
                 KH(nH,2) = 1


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
                 
                 PH(nH,1) = pxpion
                 PH(nH,2) = pypion
                 PH(nH,3) = pzpion
                 PH(nH,4) = sqrt(pxpion**2.+pypion**2.+pzpion**2.
     .                +xmpi**2.)
                 XH(nH,1) = xML
                 XH(nH,2) = yML
                 XH(nH,3) = zML
                 XH(nH,4) = tML

                 nchg(i) = 0
                 nchg(j) = 0
                 goto 30

!*************** rho -> pi + pi in the ground state ********************
              else
                 call twobody(xmrho,xmpi,xmpi,pxpion,pypion,pzpion,
     .                pxrhotopi1,pyrhotopi1,pzrhotopi1,erhotopi1,
     .                pxrhotopi2,pyrhotopi2,pzrhotopi2,erhotopi2)

                 nH = nH + 1
                 KH(nH,2) = 1


!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
                 
                 PH(nH,1) = pxrhotopi1
                 PH(nH,2) = pyrhotopi1
                 PH(nH,3) = pzrhotopi1
                 PH(nH,4) = erhotopi1
                 vxrhotopi1 = pxrhotopi1/erhotopi1
                 vyrhotopi1 = pyrhotopi1/erhotopi1
                 vzrhotopi1 = pzrhotopi1/erhotopi1
                 
                 gamrhotopi1 = hbarc/sqrt(1.-vxrhotopi1**2.-
     .                vyrhotopi1**2.-vzrhotopi1**2.)
                 XH(nH,1) = xML+vxrhotopi1*gamrhotopi1/xmpi
                 XH(nH,2) = yML+vyrhotopi1*gamrhotopi1/xmpi
                 XH(nH,3) = zML+vzrhotopi1*gamrhotopi1/xmpi
                 XH(nH,4) = tML+gamrhotopi1/xmpi
                 
                 nH = nH + 1
                 KH(nH,2) = 1

!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
                 
                 PH(nH,1) = pxrhotopi2
                 PH(nH,2) = pyrhotopi2
                 PH(nH,3) = pzrhotopi2
                 PH(nH,4) = erhotopi2
                 vxrhotopi2 = pxrhotopi2/erhotopi2
                 vyrhotopi2 = pyrhotopi2/erhotopi2
                 vzrhotopi2 = pzrhotopi2/erhotopi2
                 gamrhotopi2 = hbarc/sqrt(1.-vxrhotopi2**2.-
     .                vyrhotopi2**2.-vzrhotopi2**2.)
                 XH(nH,1) = xML+vxrhotopi2*gamrhotopi2/xmpi
                 XH(nH,2) = yML+vyrhotopi2*gamrhotopi2/xmpi
                 XH(nH,3) = zML+vzrhotopi2*gamrhotopi2/xmpi
                 XH(nH,4) = tML+gamrhotopi2/xmpi
              endif

              nchg(i) = 0
              nchg(j) = 0
              goto 30

           endif
!************* four, three, two  pion decay *************!WB 
           if(xinvm.gt.xmpi*4.) then !wenbin
!**************** two-pion decay **************************************************************
             if(probdec.le.decpro2) then
!             if(xinvm.le.xmpi*3.) then
              call twobody(xinvm,xmpi,xmpi,pxpion,pypion,pzpion,
     .             pxpitopi1,pypitopi1,pzpitopi1, epitopi1,pxpitopi2,
     .             pypitopi2,pzpitopi2,epitopi2)

              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              

              PH(nH,1) = pxpitopi1
              PH(nH,2) = pypitopi1
              PH(nH,3) = pzpitopi1
              PH(nH,4) = epitopi1
              vxpitopi1 = pxpitopi1/epitopi1
              vypitopi1 = pypitopi1/epitopi1
              vzpitopi1 = pzpitopi1/epitopi1
              gampitopi1 = hbarc/sqrt(1.-vxpitopi1**2.-
     .             vypitopi1**2.-vzpitopi1**2.)
              XH(nH,1) = xML+vxpitopi1*gampitopi1/xmpi
              XH(nH,2) = yML+vypitopi1*gampitopi1/xmpi
              XH(nH,3) = zML+vzpitopi1*gampitopi1/xmpi
              XH(nH,4) = tML+gampitopi1/xmpi

              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxpitopi2
              PH(nH,2) = pypitopi2
              PH(nH,3) = pzpitopi2
              PH(nH,4) = epitopi2
              vxpitopi2 = pxpitopi2/epitopi2
              vypitopi2 = pypitopi2/epitopi2
              vzpitopi2 = pzpitopi2/epitopi2
              gampitopi2 = hbarc/sqrt(1.-vxpitopi2**2.-
     .             vypitopi2**2.-vzpitopi2**2.)
              XH(nH,1) = xML+vxpitopi2*gampitopi2/xmpi
              XH(nH,2) = yML+vypitopi2*gampitopi2/xmpi
              XH(nH,3) = zML+vzpitopi2*gampitopi2/xmpi
              XH(nH,4) = tML+gampitopi2/xmpi
                
              nchg(i) = 0
              nchg(j) = 0
              goto 30
!           endif
           endif
!***************** three-pion decay ******************************************************
            if(probdec.gt.decpro2.and.probdec.le.decpro2+decpro3) then
!            if(xinvm.gt.xmpi*3.and.xinvm.le.xmpi*4.) then
              call threebody(xinvm,xmpi,xmpi,xmpi,pxpion,pypion,
     .             pzpion,pxpitopi1,pypitopi1,pzpitopi1,epitopi1,
     .             pxpitopi2,pypitopi2,pzpitopi2,epitopi2,pxpitopi3,
     .             pypitopi3,pzpitopi3,epitopi3)
                
              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxpitopi1
              PH(nH,2) = pypitopi1
              PH(nH,3) = pzpitopi1
              PH(nH,4) = epitopi1
              vxpitopi1 = pxpitopi1/epitopi1
              vypitopi1 = pypitopi1/epitopi1
              vzpitopi1 = pzpitopi1/epitopi1
              gampitopi1 = hbarc/sqrt(1.-vxpitopi1**2.-
     .                       vypitopi1**2.-vzpitopi1**2.)
              XH(nH,1) = xML+vxpitopi1*gampitopi1/xmpi
              XH(nH,2) = yML+vypitopi1*gampitopi1/xmpi
              XH(nH,3) = zML+vzpitopi1*gampitopi1/xmpi
              XH(nH,4) = tML+gampitopi1/xmpi

              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxpitopi2
              PH(nH,2) = pypitopi2
              PH(nH,3) = pzpitopi2
              PH(nH,4) = epitopi2
              vxpitopi2 = pxpitopi2/epitopi2
              vypitopi2 = pypitopi2/epitopi2
              vzpitopi2 = pzpitopi2/epitopi2
              gampitopi2 = hbarc/sqrt(1.-vxpitopi2**2.-
     .             vypitopi2**2.-vzpitopi2**2.)
              XH(nH,1) = xML+vxpitopi2*gampitopi2/xmpi
              XH(nH,2) = yML+vypitopi2*gampitopi2/xmpi
              XH(nH,3) = zML+vzpitopi2*gampitopi2/xmpi
              XH(nH,4) = tML+gampitopi2/xmpi
              
              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxpitopi3
              PH(nH,2) = pypitopi3
              PH(nH,3) = pzpitopi3
              PH(nH,4) = epitopi3
              vxpitopi3 = pxpitopi3/epitopi3
              vypitopi3 = pypitopi3/epitopi3
              vzpitopi3 = pzpitopi3/epitopi3
              gampitopi3 = hbarc/sqrt(1.-vxpitopi3**2.-
     .             vypitopi3**2.-vzpitopi3**2.)
              XH(nH,1) = xML+vxpitopi3*gampitopi3/xmpi
              XH(nH,2) = yML+vypitopi3*gampitopi3/xmpi
              XH(nH,3) = zML+vzpitopi3*gampitopi3/xmpi
              XH(nH,4) = tML+gampitopi3/xmpi
              
              nchg(i) = 0
              nchg(j) = 0
              goto 30
!           endif
           endif
!************* four-pion decay ************************************************************
             if(probdec.gt.decpro2+decpro3) then
!             if(xinvm.gt.xmpi*4.) then
                
              xmdx(1) = xmpi
              xmdx(2) = xmpi
              xmdx(3) = xmpi
              xmdx(4) = xmpi
              pMx(1) = pxpion
              pMx(2) = pypion
              pMx(3) = pzpion

              call fourbody(xinvm,xmdx,pMx,pdx)
              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              PH(nH,1) = pdx(1,1)
              PH(nH,2) = pdx(1,2)
              PH(nH,3) = pdx(1,3)
              PH(nH,4) = pdx(1,4)
              vxpitopi1 = PH(nH,1)/PH(nH,4)
              vypitopi1 = PH(nH,2)/PH(nH,4)
              vzpitopi1 = PH(nH,3)/PH(nH,4)
              gampitopi1 = hbarc/sqrt(1.-vxpitopi1**2.-vypitopi1**2.                               
     .             -vzpitopi1**2.)
              XH(nH,1) = xML+vxpitopi1*gampitopi1/xmpi
              XH(nH,2) = yML+vypitopi1*gampitopi1/xmpi
              XH(nH,3) = zML+vzpitopi1*gampitopi1/xmpi
              XH(nH,4) = tML+gampitopi1/xmpi

              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              PH(nH,1) = pdx(2,1)
              PH(nH,2) = pdx(2,2)
              PH(nH,3) = pdx(2,3)
              PH(nH,4) = pdx(2,4)
              vxpitopi2 = PH(nH,1)/PH(nH,4)
              vypitopi2 = PH(nH,2)/PH(nH,4)
              vzpitopi2 = PH(nH,3)/PH(nH,4)
              gampitopi2 = hbarc/sqrt(1.-vxpitopi2**2.-vypitopi2**2.                               
     .             -vzpitopi2**2.)
              XH(nH,1) = xML+vxpitopi2*gampitopi2/xmpi
              XH(nH,2) = yML+vypitopi2*gampitopi2/xmpi
              XH(nH,3) = zML+vzpitopi2*gampitopi2/xmpi
              XH(nH,4) = tML+gampitopi2/xmpi
              
              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              PH(nH,1) = pdx(3,1)
              PH(nH,2) = pdx(3,2)
              PH(nH,3) = pdx(3,3)
              PH(nH,4) = pdx(3,4)
              vxpitopi3 = PH(nH,1)/PH(nH,4)
              vypitopi3 = PH(nH,2)/PH(nH,4)
              vzpitopi3 = PH(nH,3)/PH(nH,4)
              gampitopi3 = hbarc/sqrt(1.-vxpitopi3**2.-vypitopi3**2.                               
     .                    -vzpitopi3**2.)
              XH(nH,1) = xML+vxpitopi3*gampitopi3/xmpi
              XH(nH,2) = yML+vypitopi3*gampitopi3/xmpi
              XH(nH,3) = zML+vzpitopi3*gampitopi3/xmpi
              XH(nH,4) = tML+gampitopi3/xmpi

              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              PH(nH,1) = pdx(4,1)
              PH(nH,2) = pdx(4,2)
              PH(nH,3) = pdx(4,3)
              PH(nH,4) = pdx(4,4)
              vxpitopi4 = PH(nH,1)/PH(nH,4)
              vypitopi4 = PH(nH,2)/PH(nH,4)
              vzpitopi4 = PH(nH,3)/PH(nH,4)
              gampitopi4 = hbarc/sqrt(1.-vxpitopi4**2.-vypitopi4**2.                               
     .             -vzpitopi4**2.)
              XH(nH,1) = xML+vxpitopi4*gampitopi4/xmpi
              XH(nH,2) = yML+vypitopi4*gampitopi4/xmpi
              XH(nH,3) = zML+vzpitopi4*gampitopi4/xmpi
              XH(nH,4) = tML+gampitopi4/xmpi
              
              nchg(i) = 0
              nchg(j) = 0
              goto 30
!           endif                ! close four-body decay
           endif                ! close four-body decay
           endif         ! close the if(xinvm.gt.xmpi*4.) then
!********** end of four,three,two decay ************** !WB

!*************  three, two  pion decay *************!WB 
          if(xinvm.gt.xmpi*3.and.xinvm.le.xmpi*4.) 
     .       then!WB
           decpro22=decpro2*1.0/(decpro2+decpro3)
           decpro33=decpro3*1.0/(decpro2+decpro3)
!!**************** two-pion decay **************************************************************
             if(probdec.le.decpro22) then
!!             if(xinvm.le.xmpi*3.) then
              call twobody(xinvm,xmpi,xmpi,pxpion,pypion,pzpion,
     .             pxpitopi1,pypitopi1,pzpitopi1, epitopi1,pxpitopi2,
     .             pypitopi2,pzpitopi2,epitopi2)

              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              

              PH(nH,1) = pxpitopi1
              PH(nH,2) = pypitopi1
              PH(nH,3) = pzpitopi1
              PH(nH,4) = epitopi1
              vxpitopi1 = pxpitopi1/epitopi1
              vypitopi1 = pypitopi1/epitopi1
              vzpitopi1 = pzpitopi1/epitopi1
              gampitopi1 = hbarc/sqrt(1.-vxpitopi1**2.-
     .             vypitopi1**2.-vzpitopi1**2.)
              XH(nH,1) = xML+vxpitopi1*gampitopi1/xmpi
              XH(nH,2) = yML+vypitopi1*gampitopi1/xmpi
              XH(nH,3) = zML+vzpitopi1*gampitopi1/xmpi
              XH(nH,4) = tML+gampitopi1/xmpi

              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxpitopi2
              PH(nH,2) = pypitopi2
              PH(nH,3) = pzpitopi2
              PH(nH,4) = epitopi2
              vxpitopi2 = pxpitopi2/epitopi2
              vypitopi2 = pypitopi2/epitopi2
              vzpitopi2 = pzpitopi2/epitopi2
              gampitopi2 = hbarc/sqrt(1.-vxpitopi2**2.-
     .             vypitopi2**2.-vzpitopi2**2.)
              XH(nH,1) = xML+vxpitopi2*gampitopi2/xmpi
              XH(nH,2) = yML+vypitopi2*gampitopi2/xmpi
              XH(nH,3) = zML+vzpitopi2*gampitopi2/xmpi
              XH(nH,4) = tML+gampitopi2/xmpi
                
              nchg(i) = 0
              nchg(j) = 0
              goto 30
!           endif
           endif
!***************** three-pion decay ******************************************************
            if(probdec.gt.decpro22.and.probdec.le.decpro22+decpro33)
     .          then                              
!            if(xinvm.gt.xmpi*3.and.xinvm.le.xmpi*4.) then
              call threebody(xinvm,xmpi,xmpi,xmpi,pxpion,pypion,
     .             pzpion,pxpitopi1,pypitopi1,pzpitopi1,epitopi1,
     .             pxpitopi2,pypitopi2,pzpitopi2,epitopi2,pxpitopi3,
     .             pypitopi3,pzpitopi3,epitopi3)
                
              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxpitopi1
              PH(nH,2) = pypitopi1
              PH(nH,3) = pzpitopi1
              PH(nH,4) = epitopi1
              vxpitopi1 = pxpitopi1/epitopi1
              vypitopi1 = pypitopi1/epitopi1
              vzpitopi1 = pzpitopi1/epitopi1
              gampitopi1 = hbarc/sqrt(1.-vxpitopi1**2.-
     .                       vypitopi1**2.-vzpitopi1**2.)
              XH(nH,1) = xML+vxpitopi1*gampitopi1/xmpi
              XH(nH,2) = yML+vypitopi1*gampitopi1/xmpi
              XH(nH,3) = zML+vzpitopi1*gampitopi1/xmpi
              XH(nH,4) = tML+gampitopi1/xmpi

              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxpitopi2
              PH(nH,2) = pypitopi2
              PH(nH,3) = pzpitopi2
              PH(nH,4) = epitopi2
              vxpitopi2 = pxpitopi2/epitopi2
              vypitopi2 = pypitopi2/epitopi2
              vzpitopi2 = pzpitopi2/epitopi2
              gampitopi2 = hbarc/sqrt(1.-vxpitopi2**2.-
     .             vypitopi2**2.-vzpitopi2**2.)
              XH(nH,1) = xML+vxpitopi2*gampitopi2/xmpi
              XH(nH,2) = yML+vypitopi2*gampitopi2/xmpi
              XH(nH,3) = zML+vzpitopi2*gampitopi2/xmpi
              XH(nH,4) = tML+gampitopi2/xmpi
              
              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxpitopi3
              PH(nH,2) = pypitopi3
              PH(nH,3) = pzpitopi3
              PH(nH,4) = epitopi3
              vxpitopi3 = pxpitopi3/epitopi3
              vypitopi3 = pypitopi3/epitopi3
              vzpitopi3 = pzpitopi3/epitopi3
              gampitopi3 = hbarc/sqrt(1.-vxpitopi3**2.-
     .             vypitopi3**2.-vzpitopi3**2.)
              XH(nH,1) = xML+vxpitopi3*gampitopi3/xmpi
              XH(nH,2) = yML+vypitopi3*gampitopi3/xmpi
              XH(nH,3) = zML+vzpitopi3*gampitopi3/xmpi
              XH(nH,4) = tML+gampitopi3/xmpi
              
              nchg(i) = 0
              nchg(j) = 0
              goto 30
!           endif
           endif
!************* four-pion decay ************************************************************
!none
           endif         ! close the if(xinvm.gt.xmpi*3.and.xinvm.le.xmpi*4.) 

!********** end of three,two decay ********** WB 
!
!
!*************  two  pion decay *************!WB 
              if(xinvm.ge.xmpi*2. .and.xinvm.le.xmpi*3.) then !WB
!**************** two-pion decay **************************************************************
!             if(probdec.le.decpro22) then
!             if(xinvm.le.xmpi*3.) then
              call twobody(xinvm,xmpi,xmpi,pxpion,pypion,pzpion,
     .             pxpitopi1,pypitopi1,pzpitopi1, 
     .             epitopi1,pxpitopi2,
     .             pypitopi2,pzpitopi2,epitopi2)

              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              

              PH(nH,1) = pxpitopi1
              PH(nH,2) = pypitopi1
              PH(nH,3) = pzpitopi1
              PH(nH,4) = epitopi1
              vxpitopi1 = pxpitopi1/epitopi1
              vypitopi1 = pypitopi1/epitopi1
              vzpitopi1 = pzpitopi1/epitopi1
              gampitopi1 = hbarc/sqrt(1.-vxpitopi1**2.-
     .             vypitopi1**2.-vzpitopi1**2.)
              XH(nH,1) = xML+vxpitopi1*gampitopi1/xmpi
              XH(nH,2) = yML+vypitopi1*gampitopi1/xmpi
              XH(nH,3) = zML+vzpitopi1*gampitopi1/xmpi
              XH(nH,4) = tML+gampitopi1/xmpi

              nH = nH + 1
              KH(nH,2) = 1
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
                   KH(nH,1) = 2
                 else
                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
                      KH(nH,1) = 0
                   else
                      KH(nH,1) = 1
                   endif
                 endif                                      !changed by wenbin 2018.11.25
              
              PH(nH,1) = pxpitopi2
              PH(nH,2) = pypitopi2
              PH(nH,3) = pzpitopi2
              PH(nH,4) = epitopi2
              vxpitopi2 = pxpitopi2/epitopi2
              vypitopi2 = pypitopi2/epitopi2
              vzpitopi2 = pzpitopi2/epitopi2
              gampitopi2 = hbarc/sqrt(1.-vxpitopi2**2.-
     .             vypitopi2**2.-vzpitopi2**2.)
              XH(nH,1) = xML+vxpitopi2*gampitopi2/xmpi
              XH(nH,2) = yML+vypitopi2*gampitopi2/xmpi
              XH(nH,3) = zML+vzpitopi2*gampitopi2/xmpi
              XH(nH,4) = tML+gampitopi2/xmpi
                
              nchg(i) = 0
              nchg(j) = 0
              goto 30
!           endif
!           endif
!***************** three-pion decay ******************************************************
!none
!************* four-pion decay ************************************************************
!none
           endif         ! close the if(xinvm.le.xmpi*3.) then

!********** end of two decay ********** WB 




        endif                   ! close light quarks
       endif                     ! quark and anti-quark
!************ begin of added by wenbin 2018.11.24 *************
!************recombine phi ***************************
!        prob6 = ran()
!        if(prob6.le.Wigphi.and.(Kq(i,2)*Kq(j,2).eq. -9)) then
           
!********ground state of phi ******************************
!              if(ran().le.0.25) then
!                 nH = nH + 1
!                 KH(nH,2) = 5
!                 PH(nH,1) = pxpion
!                 PH(nH,2) = pypion
!                 PH(nH,3) = pzpion
!                 PH(nH,4) = sqrt(pxpion**2.+pypion**2.+pzpion**2.+
!     .                xmphi**2.)
!                 XH(nH,1) = xML
!                 XH(nH,2) = yML
!                 XH(nH,3) = zML
!                 XH(nH,4) = tML
                 
!                 if(Kq(i,1)*Kq(j,1).ne.0) KH(nH,1) = 2        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).eq.0) KH(nH,1) = 1        !chaned by wenbin 2018.11.25
!                 if(Kq(i,1)*Kq(j,1).ne.0) then               !changed by wenbin 2018.11.25
!                   KH(nH,1) = 2
!                 else
!                   if(Kq(i,1).eq.0 .and. Kq(j,1).eq.0)then 
!                      KH(nH,1) = 0
!                   else
!                      KH(nH,1) = 1
!                   endif
!                 endif                                      !changed by wenbin 2018.11.25
                               
!                 nchg(i) = 0
!                 nchg(j) = 0
!                 goto 30
!              endif
!        endif
!************ end of added by wenbin 2018.11.24 *************

 40   continue
      ENDDO                     ! j=i+1, numq !!second loop Wenbin
 30   CONTINUE   
      ENDDO

!*************************************!                                         
!     Remnant  Quarks                  !                                        
!*************************************!        
      kremtot = 0
      do irem=1, numqx
         if(nchg(irem).ne.0.and.(Kq(irem,1).eq.1)) then!changed by wenbin 2018.11.23!Kq(irem,1)=1 for shower parton Kq(irem,1)=0 for thermal parton !noted by wenbin 2018.11.23 
!         if(nchg(irem).ne.0.) then!Both the remnant thermal and shower partons fragment  !changed by wenbin 2018.11.23
            nr = nr + 1
            Kr(nr,1) = Kq(irem,3)
            Kr(nr,2) = Kq(irem,2)
            Pr(nr,1) = Pq(irem,1)
            Pr(nr,2) = Pq(irem,2)
            Pr(nr,3) = Pq(irem,3)
            Pr(nr,4) = Pq(irem,4)
            Xr(nr,1) = Xq(irem,1)
            Xr(nr,2) = Xq(irem,2)
            Xr(nr,3) = Xq(irem,3)
            Xr(nr,4) = Xq(irem,4)
            rqscale(nr)=Q0(irem)
            !write(*,*)"111 ",Pr(nr,2),nr
            kremtot = kremtot + Kr(nr,2)/abs(Kr(nr,2))
         endif
!*************************************!
!     coalesced thermal quarks  Quarks                  !
!*************************************!
         if(nchg(irem).eq.0.and.(Kq(irem,1).eq.0)) then!changed by wenbin 2018.11.23!Kq(irem,1)=1 for shower parton Kq(irem,1)=0 for thermal parton !noted by wenbin 2018.11.23
!         if(nchg(irem).ne.0.) then!Both the remnant thermal and shower
!         partons fragment  !changed by wenbin 2018.11.23
            NCTnr = NCTnr + 1
           NCTKr(NCTnr) = Kq(irem,2)
           CTPr(NCTnr,1) = Pq(irem,1)
           CTPr(NCTnr,2) = Pq(irem,2)
           CTPr(NCTnr,3) = Pq(irem,3)
           CTPr(NCTnr,4) = Pq(irem,4)
           CTXr(NCTnr,1) = Xq(irem,1)
           CTXr(NCTnr,2) = Xq(irem,2)
           CTXr(NCTnr,3) = Xq(irem,3)
           CTXr(NCTnr,4) = Xq(irem,4)
            !kremtot = kremtot + Kr(nr,2)/abs(Kr(nr,2))
         endif


      enddo


!*************************************!
!     coalesced thermal quarks  Quarks                  !
!*************************************!
      !kremtot = 0
      !write(*,*)"numqx= ",numqx
!      do irem=1, numqx
!         if(nchg(irem).eq.0.and.(Kq(irem,1).eq.0)) then!changed by wenbin 2018.11.23!Kq(irem,1)=1 for shower parton Kq(irem,1)=0 for thermal parton !noted by wenbin 2018.11.23
!         if(nchg(irem).ne.0.) then!Both the remnant thermal and shower
!         partons fragment  !changed by wenbin 2018.11.23
!            NCTnr = NCTnr + 1
!           NCTKr(NCTnr) = Kq(irem,2)
!           CTPr(NCTnr,1) = Pq(irem,1)
!           CTPr(NCTnr,2) = Pq(irem,2)
!           CTPr(NCTnr,3) = Pq(irem,3)
!           CTPr(NCTnr,4) = Pq(irem,4)
!           CTXr(NCTnr,1) = Xq(irem,1)
!           CTXr(NCTnr,2) = Xq(irem,2)
!           CTXr(NCTnr,3) = Xq(irem,3)
!           CTXr(NCTnr,4) = Xq(irem,4)
            !kremtot = kremtot + Kr(nr,2)/abs(Kr(nr,2))
!         endif
!      enddo
      !write(*,*)"a3= ",NCTnr,nr
C.... Adding fake parton .....                                                  
C...  print*, kremtot                                                       
!      if(mod(kremtot,3).ne.0) then
!         nr = nr + 1
!         Kr(nr,1) = 0
!         if(abs(mod(kremtot,3)).eq.1) Kr(nr,2) = -mod(kremtot,3)
!         if(mod(kremtot,3).eq.2) Kr(nr,2) = 1
!         if(mod(Kremtot,3).eq.-2) Kr(nr,2) = -1
!         Pr(nr,1) = 0.
!         Pr(nr,2) = 0.
!         Pr(nr,3) = 0.
!         Pr(nr,4) = 0.
!         Xr(nr,1) = 0.
!         Xr(nr,2) = 0.
!         Xr(nr,3) = 0.
!         Xr(nr,4) = 0.
!      endif
      return
      end

!********************************
!     factorial function        ! 
!********************************
      function factor(Nf)
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      factor = 1.
      do i=2, Nf
         factor=factor*i
      enddo
      end function factor


!=============================================                                  
!    Subroutine for Gluon Decay 1:1:1             !  added by wenbin 2018.12.11                                
!=============================================                                  
      subroutine gluondec111(xmg,xg,yg,zg,tg,px,py,pz,xn,yn,
     .  zn,tn,nfl1,px1,py1,pz1,e1,nfl2,px2,py2,pz2,e2,xmu,xms)
!     xmg : mass of gluon, (px,py,pz)  : momentum of gluon                          
!     (xg, yg, zg, tg) : space and time of gluon                                    
!     (xn, yn, zn, tn) : space and time of quark and antiquark                      
!     nfl1, nfl2 : flavors of created quarks from the decay                         
!     px1, py1, pz1; px2, py2, pz2 : momenta of the quarks 

      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
      hbarc = 0.19732
!      xmu = 0.3
!      xms = 0.45
!*****ratio = Gamma(g->ssbar)/Gamma(g->uubar, ddbar)******* 
         prob=ran()
         if(prob.gt. 0.5) then
            xmq = xms
            nfl1 = 3
            nfl2 = -3
         endif

         if(prob.ge. 0.25.and.prob.lt. 0.5) 
     .        then
            xmq = xmu
            nfl1 = 1
            nfl2 = -1
         endif

         if(prob.lt. 0.25) then
            xmq = xmu
            nfl1 = 2
            nfl2 = -2
         endif


      eg = sqrt(xmg**2.+px**2.+py**2.+pz**2.)

!**** velocities of gluons in a rest frame **********                           
      betagx = -px/eg
      betagy = -py/eg
      betagz = -pz/eg
      gamma = 1./sqrt(betagx**2.+betagy**2.+betagz**2.)
      tau = hbarc*gamma/xmg
      
!***  momenta of quark and antiquark in the gluon-rest-frame.                    
      pq = sqrt((xmg**2.)/4.- xmq**2.)
      tn = tg + tau
      xn = xg - betagx*tau
      yn = yg - betagy*tau
      zn = zg - betagz*tau
      
      theta = acos(1.-2.*ran())
      phi = 2.*pi*ran()
      pqx = pq*sin(theta)*cos(phi)
      apqx = -pqx
      pqy = pq*sin(theta)*sin(phi)
      apqy = -pqy
      pqz = pq*cos(theta)
      apqz = -pqz
      eq = sqrt(xmq**2.+ pq**2.)
      
!**** Matrix Elements of Lorentz Boost **********                               
      call lorentzboost(betagx,betagy,betagz,eq,pqx,pqy,pqz,e1,px1,
     .     py1,pz1)
      call lorentzboost(betagx,betagy,betagz,eq,apqx,apqy,apqz,e2,
     .     px2,py2,pz2)
      
      return
      end
!=============================================                                  
!    Subroutine for Gluon Decay              !                                  
!=============================================                                  
      subroutine gluondec(xmg,xg,yg,zg,tg,px,py,pz,xn,yn,zn,tn,nfl1,
     .     px1,py1,pz1,e1,nfl2,px2,py2,pz2,e2)
!     xmg : mass of gluon, (px,py,pz)  : momentum of gluon                          
!     (xg, yg, zg, tg) : space and time of gluon                                    
!     (xn, yn, zn, tn) : space and time of quark and antiquark                      
!     nfl1, nfl2 : flavors of created quarks from the decay                         
!     px1, py1, pz1; px2, py2, pz2 : momenta of the quarks 

      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
      hbarc = 0.19732
      xmu = 0.33
      xms = 0.43
!*****ratio = Gamma(g->ssbar)/Gamma(g->uubar, ddbar)******* 
      ratio = 0.5*sqrt((xmg**2.-4.*xms**2.)/(xmg**2.-4.*xmu**2.))*
     .     ((xmg**2.+2.*xms**2.)/(xmg**2.+2.*xmu**2.))

      if(xmg.ge.0.9) then
         prob=ran()
         if(prob.lt.ratio/(1.+ratio)) then
            xmq = xms
            nfl1 = 3
            nfl2 = -3
         endif
            
         if(prob.ge.ratio/(1.+ratio).and.prob.lt.(ratio+0.5)/(1.+ratio)) 
     .        then
            xmq = xmu
            nfl1 = 1
            nfl2 = -1
         endif
     
         if(prob.ge.(ratio+0.5)/(1.+ratio)) then
            xmq = xmu
            nfl1 = 2
            nfl2 = -2
         endif
     
      else
         xmq = xmu
         if(ran().ge.0.5) then
            nfl1 = 2
            nfl2 = -2
         else
            nfl1 = 1
            nfl2 = -1
         endif
      endif

      eg = sqrt(xmg**2.+px**2.+py**2.+pz**2.)

!**** velocities of gluons in a rest frame **********                           
      betagx = -px/eg
      betagy = -py/eg
      betagz = -pz/eg
      gamma = 1./sqrt(betagx**2.+betagy**2.+betagz**2.)
      tau = hbarc*gamma/xmg
      
!***  momenta of quark and antiquark in the gluon-rest-frame.                    
      pq = sqrt((xmg**2.)/4.- xmq**2.)
      tn = tg + tau
      xn = xg - betagx*tau
      yn = yg - betagy*tau
      zn = zg - betagz*tau
      
      theta = acos(1.-2.*ran())
      phi = 2.*pi*ran()
      pqx = pq*sin(theta)*cos(phi)
      apqx = -pqx
      pqy = pq*sin(theta)*sin(phi)
      apqy = -pqy
      pqz = pq*cos(theta)
      apqz = -pqz
      eq = sqrt(xmq**2.+ pq**2.)
      
!**** Matrix Elements of Lorentz Boost **********                               
      call lorentzboost(betagx,betagy,betagz,eq,pqx,pqy,pqz,e1,px1,
     .     py1,pz1)
      call lorentzboost(betagx,betagy,betagz,eq,apqx,apqy,apqz,e2,
     .     px2,py2,pz2)
      
      return
      end



!=============================
!     Lorentz Boost           !                                                  
!=============================                                                  
      subroutine lorentzboost(betax,betay,betaz,tlab,xlab,ylab,
     .     zlab,tn,xn,yn,zn)
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
      
      beta2 = betax**2.+betay**2.+betaz**2.
      gamma = 1./sqrt(1.-beta2)
      xlam00= gamma
      xlam01 = -gamma*betax
      xlam02 = -gamma*betay
      xlam03 = -gamma*betaz
      xlam10 = xlam01
      xlam11 = 1.+(gamma-1.)*(betax**2.)/beta2
      xlam12 = (gamma-1.)*betax*betay/beta2
      xlam13 = (gamma-1.)*betax*betaz/beta2
      xlam20 = xlam02
      xlam21 = xlam12
      xlam22 = 1.+(gamma-1.)*(betay**2.)/beta2
      xlam23 = (gamma-1.)*betay*betaz/beta2
      xlam30 = xlam03
      xlam31 = xlam13
      xlam32 = xlam23
      xlam33 = 1.+(gamma-1.)*(betaz**2.)/beta2
      tn = tlab*xlam00+xlab*xlam01+ylab*xlam02+zlab*xlam03
      xn = tlab*xlam10+xlab*xlam11+ylab*xlam12+zlab*xlam13
      yn = tlab*xlam20+xlab*xlam21+ylab*xlam22+zlab*xlam23
      zn = tlab*xlam30+xlab*xlam31+ylab*xlam32+zlab*xlam33
      return
      end

      
!=======================================================
!     subroutine for calculation of two body decay     |
!=======================================================
      subroutine twobody(xM, xm1, xm2,px,py,pz,px1,py1,pz1,e1,
     .     px2,py2,pz2,e2)
!***  input xM, xm1,xm2,px,py,pz **************************
!***  output px1,py1,pz1,e1,px2,py2,pz2,e2 ****************
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
!     nseed = 6
      e = sqrt(xM**2.+px**2.+py**2.+pz**2.)
!******back speed of mother particle
      betax = -px/e
      betay = -py/e
      betaz = -pz/e
      beta2 = betax**2.+betay**2.+betaz**2.
      gam = 1./sqrt(1.- beta2)
!**** momeuntum of daughter particle in mother-rest-frame***
      pd = sqrt((xM**2.-(xm1-xm2)**2.)*(xM**2.-(xm1+xm2)**2.))/(2.*xM)
!**** randomly isotropic ***********************
      theta = acos(1.-2.*ran())
      phi = 2.*pi*ran()
      
      pdx = pd*sin(theta)*cos(phi)
      apdx = -pdx
      pdy = pd*sin(theta)*sin(phi)
      apdy = -pdy
      pdz = pd*cos(theta)
      apdz = -pdz
      ed1 = sqrt(xm1**2. + pd**2.)
      ed2 = sqrt(xm2**2. + pd**2.)
      
!**** Matrix Elements of Lorentz Boost **********
      xlam00 = gam
      xlam01 = -gam*betax
      xlam02 = -gam*betay
      xlam03 = -gam*betaz
      xlam11 = 1.+(gam-1.)*(betax**2.)/beta2
      xlam12 = (gam-1.)*(betax*betay)/beta2
      xlam13 = (gam-1.)*(betax*betaz)/beta2
      xlam22 = 1.+(gam-1.)*(betay**2.)/beta2
      xlam23 = (gam-1.)*(betay*betaz)/beta2
      xlam10 = xlam01
      xlam20 = xlam02
      xlam21 = xlam12
      xlam30 = xlam03
      xlam31 = xlam13
      xlam32 = xlam23
      xlam33 = 1.+(gam-1.)*(betaz**2.)/beta2
      
!******momenta of daughter particles in the lab frame ************
      px1 = xlam10*ed1 + xlam11*pdx + xlam12*pdy + xlam13*pdz
      py1 = xlam20*ed1 + xlam21*pdx + xlam22*pdy + xlam23*pdz
      pz1 = xlam30*ed1 + xlam31*pdx + xlam32*pdy + xlam33*pdz
      e1 =  xlam00*ed1 + xlam01*pdx + xlam02*pdy + xlam03*pdz
      
      px2 = xlam10*ed2 + xlam11*apdx + xlam12*apdy + xlam13*apdz
      py2 = xlam20*ed2 + xlam21*apdx + xlam22*apdy + xlam23*apdz
      pz2 = xlam30*ed2 + xlam31*apdx + xlam32*apdy + xlam33*apdz
      e2 =  xlam00*ed2 + xlam01*apdx + xlam02*apdy + xlam03*apdz
      return
      end


 
!=========================================================!
!     subroutine for calculation of three body decay      !
!=========================================================!
      subroutine threebody(xM, xm1, xm2, xm3, px,py,pz,px1,py1,pz1,e1,
     .     px2,py2,pz2,e2,px3,py3,pz3,e3)
!***  input       xM, xm1,xm2,xm3,px,py,pz      **************************
!***  output px1,py1,pz1,e1,px2,py2,pz2,e2,px3,py3,pz3,e3 ****************
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
      xs = xM
      call  twocoupledmass(xs, xm1, xm2, xm3, xM12)
      call twobody(xM,xM12,xm3,px,py,pz,px12,py12,pz12,e12,xpx3,
     .     xpy3,xpz3,xe3)
      px3 = xpx3
      py3 = xpy3
      pz3 = xpz3
      e3 = xe3
      call twobody(xM12,xm1,xm2,px12,py12,pz12,xpx1,xpy1,xpz1,xe1,
     .     xpx2,xpy2,xpz2,xe2)
      px1 = xpx1
      py1 = xpy1
      pz1 = xpz1
      e1 = xe1
      px2 = xpx2
      py2 = xpy2
      pz2 = xpz2
      e2 = xe2
      return
      end
      
      
!*****************************************************
!     coupled mass in three body decay                !
!     ref. to PDG booklet p.273                      !
!*****************************************************
      subroutine twocoupledmass(xs, xm1,xm2,xm3,xM12)
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
 12   continue
!.....m12**2.      ................................
      x = (xm1+xm2)**2.+ran()*((xs-xm3)**2.-(xm1+xm2)**2.)
      
!.....m23**2. .......................................
      y = (xm2+xm3)**2.+ran()*((xs-xm1)**2.-(xm2+xm3)**2.)
      
!.... energies of 2 and 3 in the xM12-rest frame .......
      e2s = (x-xm1**2.+xm2**2.)/(2.*sqrt(x))
      e3s = (xs**2.-x-xm3**2.)/(2.*sqrt(x))

!..... maximum and minimum of y for a given value of xm12**2 ..............
      ymax = (e2s+e3s)**2.-(sqrt(e2s**2.-xm2**2.)-sqrt(e3s**2.- 
     .     xm3**2.))**2.
      ymin = (e2s+e3s)**2.-(sqrt(e2s**2.-xm2**2.)+sqrt(e3s**2.- 
     .     xm3**2.))**2.
!      write(*,*)"e2s=, xm2= ",e2s,xm2,ymax,ymin
      if(y.le.ymax.and.y.ge.ymin) goto 11
      goto 12
 11   xM12 = sqrt(x)
      return
      end


!=====================================================
!     Subroutine for four-body decay                   !
!=====================================================
!***********************************************!
!     input : xM, pM, xmd                           !
!     output : pd momenta-energies of daughters     !
!***********************************************!
      subroutine fourbody(xM,xmd,pM,pd)
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
      dimension xmd(4),pM(3),pd(4,4)
      pMx = pM(1)
      pMy = pM(2)
      pMz = pM(3)
      eM = sqrt(pMx**2.+pMy**2.+pMz**2.+xM**2.)
      xm1 = xmd(1)
      xm2 = xmd(2)
      xm3 = xmd(3)
      xm4 = xmd(4)
      call threecoupledmass(xM,xm1,xm2,xm3,xm4,xM123)
      call twobody(xM,xM123,xm4,pMx,pMy,pMz,px123,py123,
     .     pz123,e123,px4,py4,pz4,e4)
      pd(4,1) = px4
      pd(4,2) = py4
      pd(4,3) = pz4
      pd(4,4) = e4
      
      call threebody(xM123,xm1,xm2,xm3,px123,py123,pz123,px1,
     .     py1,pz1,e1,px2,py2,pz2,e2,px3,py3,pz3,e3)
      pd(1,1) = px1
      pd(1,2) = py1
      pd(1,3) = pz1
      pd(1,4) = e1
      
      pd(2,1) = px2
      pd(2,2) = py2
      pd(2,3) = pz2
      pd(2,4) = e2
      
      pd(3,1) = px3
      pd(3,2) = py3
      pd(3,3) = pz3
      pd(3,4) = e3
      
      return
      end
      
!*************************************************
!     Three coupled mass for 4-body decay        !
!*************************************************
      subroutine threecoupledmass(xs,xm1,xm2,xm3,xm4,xM123)
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)      
      common/const/pi,hbarc
      
 21   continue
      x = (xm1+xm2+xm3) + ran()*(xs-xm1-xm2-xm3-xm4)
!***  trial (m1,m2) coupled mass *************************
      xMc = (xm1+xm2) + ran()*(x-xm1-xm2-xm3)
      
      y = alambda(xs,x**2.,xm4**2.)*alambda(x**2.,xMc**2.,xm3**2.)* 
     .     alambda(xMc**2.,xm1**2.,xm2**2.)/(x*xMc)
      call fourmaxmin(xs,xm1,xm2,xm3,xm4,xMax,xMin)
      g = xMin + ran()*(xMax-xMin)
      if(g.le.y) goto 22
      goto 21
 22   xM123 = x
      return
      end
      
      
!*****************************************
!     Define Lambda function                !
!*****************************************
      function alambda(x,y,z)
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      alambda = x*x+y*y+z*z-2.*x*y-2.*y*z-2.*x*z
      return
      end function alambda
      
!*************************************************
!     Maximum and Minimum of weights in 4-body     !
!*************************************************
!     input: xs, xm1, xm2, xm3, xm4                !
!     output : xMax, xMin                          !
!*************************************************
      subroutine fourmaxmin(xs,xm1,xm2,xm3,xm4,xMax,xMin)
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
      n = 1000
      xMt = xm1+xm2+xm3+xm4+ran()*(xs-xm1-xm2-xm3-xm4)
      xMd = xm1+xm2+ran()*(xMt-xm1-xm2-xm3)
      xMax = alambda(xs,xMt**2.,xm4**2.)*alambda(xMt**2.,xMd**2.,
     .     xm3**2.)*alambda(xMd**2.,xm1**2.,xm2**2.)/(xMt*xMd)
      xMin = alambda(xs,xMt**2.,xm4**2.)*alambda(xMt**2.,xMd**2.,
     .     xm3**2.)*alambda(xMd**2.,xm1**2.,xm2**2.)/(xMt*xMd)
      
      do i=1, n
!**** trial mass of (m1,m2,m3) system *****************
         xMt = xm1+xm2+xm3+xm4+ran()*(xs-xm1-xm2-xm3-xm4)
!**** trial mass of (m1,m2,m3) system ****************
         xMd = xm1+xm2+ran()*(xMt-xm1-xm2-xm3)
         
         y =  alambda(xs,xMt**2.,xm4**2.)*alambda(xMt**2.,xMd**2.,
     .        xm3**2.)*alambda(xMd**2.,xm1**2.,xm2**2.)/(xMt*xMd)
         if(y.ge.xMax) xMax = y
         if(y.le.xMin) xMin = y
      enddo
      return
      end
  
   
C*******************************************************************!
C                                                                   C
C   Subroutine for forming gluons from q qbar decayed from same     !
C     gluons and diquark (or anti-diquarks) in the remnant partons  C
C    Input:                                                         !
C    Nq: number of remnant quarks, Kq(i,j): origin and pid of       C
C    remnant quarks. j=1 (origin), 2(pid)                           ! 
C     Kq(i,1)=0 (leading quarks), 1,2,..(numerical                  C
C        orders of their mother gluons)                             !
C    Pq(i,j) : 3-momentum of remnant partons (px:j=1,py:j=2,pz:j=3) C
C    Xq(i,j) : position and time of the remnant partons (x,y,z,t)   !
C   Output:                                                         C
C   Npt: number of partons (quark,gluon,diquarks)....               !
C   Kpt: id of partons                                              C
C   Ppt: 3-momenta of the remnant partons                           !
C   Xpt: position and time of the ordered remnant partons           C
C                                                                   C
C*******************************************************************!
      subroutine remorg(Nq,Kq,Pq,Xq,Npt,Kpt,Ppt,Xpt,rqscale,
     .                  qscale)
      
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      dimension Kq(5000,5),Pq(5000,5),Xq(5000,5)
      dimension Ppg(5000,5),Xpg(5000,5)
      dimension Kx(5000),rqscale(5000),qscale(5000)
      dimension Kpt(5000), Ppt(5000,5),Xpt(5000,5)
      dimension nchg(5000),ndhg(5000)


      do ir=1, Nq
         nchg(ir) = Kq(ir,2)/abs(Kq(ir,2))
      enddo
      
      kg = 0
      sumsharge=0.0
!***********************************************
!     Restroing gluons from q and qbar          !
!***********************************************
      do i=1, Nq
         if(nchg(i).eq.0) goto 41
C........ leading quarks ........................
         if(Kq(i,1).eq.0) then
            kg = kg + 1
            Ppt(kg,1) = Pq(i,1)
            Ppt(kg,2) = Pq(i,2)
            Ppt(kg,3) = Pq(i,3)
            Ppt(kg,4) = Pq(i,4)
            Xpt(kg,1) = Xq(i,1)
            Xpt(kg,2) = Xq(i,2)
            Xpt(kg,3) = Xq(i,3)
            Xpt(kg,4) = Xq(i,4)
            !write(*,*)"222 ",Ppt(kg,2)
            Kpt(kg) = Kq(i,2)
            qscale(kg)=rqscale(i)
            ndhg(kg) = 1
            nchg(i) = 0
            goto 41
         endif

!..... form gluons ..................
         do j=i+1, Nq
            if(nchg(j).eq.0) goto 42
            if(Kq(i,1).eq.Kq(j,1).and.Kq(j,1).ne.0) then
               kg = kg + 1
               Ppt(kg,1) = Pq(i,1) + Pq(j,1)
               Ppt(kg,2) = Pq(i,2) + Pq(j,2)
               Ppt(kg,3) = Pq(i,3) + Pq(j,3)
               Ppt(kg,4) = Pq(i,4) + Pq(j,4)
               Xpt(kg,1) = (Xq(i,1)+Xq(j,1))/2.
               Xpt(kg,2) = (Xq(i,2)+Xq(j,2))/2.
               Xpt(kg,3) = (Xq(i,3)+Xq(j,3))/2.
               Xpt(kg,4) = (Xq(i,4)+Xq(j,4))/2.
               qscale(kg)= (rqscale(i)+rqscale(j))/2.
               Kpt(kg) = 21
               ndhg(kg) = 1
               nchg(i) = 0
               nchg(j) = 0
               goto 41
            endif
 42         continue
         enddo
!..... form half gluons ..................!added by wenbin 
!         if(nchg(i) .ne. 0 .and. Kq(i,1).ne.0 ) then
!              probgg = ran()!changed by wenbin
!              if(probgg.le. 0.50) then
!                 kg = kg + 1
!                 Ppt(kg,1) = Pq(i,1)*2.
!                 Ppt(kg,2) = Pq(i,2)*2.
!                 Ppt(kg,3) = Pq(i,3)*2.
!                 Ppt(kg,4) = Pq(i,4)*2.
!                 Xpt(kg,1) = (Xq(i,1)+Xq(i,1))/2.
!                 Xpt(kg,2) = (Xq(i,2)+Xq(i,2))/2.
!                 Xpt(kg,3) = (Xq(i,3)+Xq(i,3))/2.
!                 Xpt(kg,4) = (Xq(i,4)+Xq(i,4))/2.
!                 Kpt(kg) = 21 !Kq(i,2)
!                 qscale(kg)=rqscale(i)
!                 ndhg(kg) = 1
!                 nchg(i) = 0
!                 goto 41
!              else
!                 nchg(i) = 0
!                 goto 41
!              endif
!         endif
!...............................................

 41      continue
      enddo
      Npt=kg
 71   continue
      return
      end

