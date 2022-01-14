!****************************************************
!     Running recombination module                  !
!****************************************************
      program main
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
      common/parm/sigma,sigmapro,sigmaK,maxf      
       CHARACTER*200 NAME1
      DIMENSION P(5000,4),XV(5000,4)
      DIMENSION Kp(5000,4),idsh(5000),Kpx(5000)
      DIMENSION PH(5000,4),Kst(5000),Qscale(5000)
      DIMENSION KfH(5000), Pth(5000,4),XH(5000,4)
      DIMENSION Ppx(5000,4), Xpx(5000,4)
      dimension Xth(5000,4),Kth(5000),qfscale(5000)
      dimension xnbin(100,2), xnbinK(100,2)
      dimension xnbinpro(100,2)
      DIMENSION NCTKr(5000), CTPr(5000,4),CTXr(5000,4)

C.... input files .........................
      call GETARG(1,NAME1)
      open (unit=30, file='shower_parton', status='unknown')
      open (unit=45, file='thermal_parton.dat', status='unknown')

C.... output files ........................
      open (unit=51, file='coaleced_hadron.dat', status='unknown')
      open (unit=52, file='remnant_jet_parton.dat', status='unknown')
      open (unit=53, file='coalesced_thermal_partons.dat'
     .             ,status='unknown')

      open (unit=88, file='input.txt', status='unknown')
!.....Read input parameters ....... 
      READ(88,*) NEV
      READ(88,*) sigma
      READ(88,*) sigmaK
      READ(88,*) sigmapro
      READ(88,*) maxf

!.........................      
      step = 30./40.
      pi = 2.*asin(1.)
      
      ntherm = 1340
      xnevt = float(NEV)
      nremq = 0
      kpi = 0
      kK = 0
      kpro = 0
C.... initialization of binning ..........
      do ini=1, 40
         do jni=1, 2
            xnbin(ini,jni) = 0.
            xnbinK(ini,jni) = 0.
            xnbinpro(ini,jni) = 0.
         enddo
      enddo

      DO IEV=1, NEV             
!..... read thermal parton from hydro ...
         READ(45,*)ntherm
         numq = 0
      do itherm=1, ntherm
         READ(45,*)Kth(itherm),Pth(itherm,1),
     .        Pth(itherm,2),Pth(itherm,3),Pth(itherm,4),
     .  Xth(itherm,1),Xth(itherm,2),Xth(itherm,3),Xth(itherm,4)       !KTh thermal parton ID
      enddo
!..... read LBT parton from LBT result 
!.....number of partons .........................

         nchth = 0
!.....  Shower partons ..........................
         READ(30,*)numparton
         DO ish=1, numparton
            READ(30,*)idsh(ish),P(ish,1),P(ish,2),
     .     P(ish,3),parton_mass,XV(ish,1),XV(ish,2),XV(ish,3),
     .           XV(ish,4)
          Qscale(ish)=0.0
              P(ish,4)=sqrt(P(ish,1)*P(ish,1)+P(ish,2)*P(ish,2)
     .         + P(ish,3)*P(ish,3)+parton_mass*parton_mass)
         enddo
         call hadronization(ntherm,Kth,Xth,Pth,numparton,idsh,XV,P,nH,
     .        KfH,XH,PH,Kst,Npt,Kpx,Ppx,Xpx, !
     .        NCTnr, NCTKr,CTPr,CTXr,Qscale,qfscale)
            write(51,17)IEV,nH
C.... Output the produced hadrons ........................
         do ihad=1, nH
            if(ihad.lt.nH) then
               if ( PH(ihad,1).eq.PH(ihad+1,1) .and.
     .              PH(ihad,2).eq.PH(ihad+1,2) .and.
     .              PH(ihad,3).eq.PH(ihad+1,3) .and.
     .              KfH(ihad).eq.KfH(ihad+1)) then
                    ptmag = sqrt(PH(ihad,1)*PH(ihad,1) +
     .                           PH(ihad,2)*PH(ihad,2) +
     .                           PH(ihad,3)*PH(ihad,3))
                   theta1 = atan2(PH(ihad,2), PH(ihad,1))
                   thetap1 = theta1+ran()*0.314 - 0.157 
                   theta2 = asin(PH(ihad,3)/ptmag)
                   thetap2 = theta2+ran()*0.314 - 0.157 
                   PH(ihad,1) = ptmag*cos(thetap)*cos(thetap2)
                   PH(ihad,2) = ptmag*sin(thetap)*cos(thetap2)
                   PH(ihad,3) = ptmag*sin(thetap2)
                endif
            endif
            write(51,16) IEV,KfH(ihad),PH(ihad,1),PH(ihad,2),
     .           PH(ihad,3),PH(ihad,4),XH(ihad,1),XH(ihad,2),XH(ihad,3),
     .           XH(ihad,4),Kst(ihad)!Kst is origin of the hadrons(th-th:0, sh-th:1, sh-sh:2,frag:3)

         enddo

!********* write the Remnant partons ******** 
            write(52,17)IEV,Npt
C.... Binning the writing ........................
         do ihad=1, Npt
            if (isnan(Xpx(ihad,1)))Xpx(ihad,1) = 40.
            if (isnan(Xpx(ihad,2)))Xpx(ihad,2) = 40.
            if (isnan(Xpx(ihad,3)))Xpx(ihad,3) = 40.
            if (isnan(Xpx(ihad,4)))Xpx(ihad,4) = 40.
            
            write(52,18) IEV,Kpx(ihad),Ppx(ihad,1),Ppx(ihad,2),
     .           Ppx(ihad,3),Ppx(ihad,4),Xpx(ihad,1),Xpx(ihad,2),
     .           Xpx(ihad,3),Xpx(ihad,4)
         enddo

!********* write the coaleced thermal  partons ********
            write(53,*)"# ",IEV,NCTnr
C.... Binning the produced hadrons ........................
         do ihad=1, NCTnr
           ttau=sqrt(CTXr(ihad,4)*CTXr(ihad,4)-CTXr(ihad,3)*CTXr(ihad,3)
     .              )
           eetas=0.5*log((CTXr(ihad,4)+CTXr(ihad,3))/(CTXr(ihad,4)
     .                   -CTXr(ihad,3)))
            write(53,18) IEV,NCTKr(ihad),CTPr(ihad,1),CTPr(ihad,2),
     .           CTPr(ihad,3),CTPr(ihad,4),CTXr(ihad,1),CTXr(ihad,2),
     .           eetas,ttau
         enddo


       ENDDO

 15   format(1(f6.2),4(e21.8))
 16   format(1(I9),1(I7),8(f21.8),1(I9))
 17   format(1(I9),1(I9))
 18   format(1(I9),1(I7),8(f21.8))
      end

