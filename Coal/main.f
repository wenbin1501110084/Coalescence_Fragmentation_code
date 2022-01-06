!****************************************************
!     Running recombination module                  !
!****************************************************
      program main
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
      common/parm/sigma,sigmapro,sigmaK,maxf      

      DIMENSION P(5000,4),XV(5000,4)
      DIMENSION Kp(5000,4),idsh(5000),Kpx(5000)
      DIMENSION PH(5000,4),Kst(5000)
      DIMENSION KfH(5000), Pth(5000,4),XH(5000,4)
      DIMENSION Ppx(5000,4), Xpx(5000,4)
      dimension Xth(5000,4),Kth(5000)
      dimension xnbin(100,2), xnbinK(100,2)
      dimension xnbinpro(100,2)
      DIMENSION NCTKr(5000), CTPr(5000,4),CTXr(5000,4)

C.... input files .........................
!      open (unit=10, file='shower.dat', status='unknown')
!      open (unit=9, file='numsh.dat', status='unknown')
!      open (unit=15, file='thermal.dat', status='unknown')

      open (unit=30, file='parton/LBT_result', status='unknown')
      open (unit=45, file='parton/oscar.dat', status='unknown')

C.... output files ........................
!      open (unit=31, file='ana_pion.dat', status='unknown')
!      open (unit=32, file='ana_kaon.dat', status='unknown')
!      open (unit=33, file='ana_proton.dat', status='unknown')
      open (unit=51, file='thermal_jet.dat', status='unknown')!output file !noted by wenbin 2018.11.30
      open (unit=52, file='remnant_jet_parton.dat', status='unknown')!output file !noted by wenbin 2018.11.30
      open (unit=53, file='coalesced_thermal_partons.dat'
     .             ,status='unknown')!output file !noted by wenbin 2018.11.30

!      open (unit=217, file='numofreco.dat', status='unknown')
      

      open (unit=88, file='input.txt', status='unknown')
!.....Read input parameters ....... 
      READ(88,*) NEV
      READ(88,*) sigma
      READ(88,*) sigmaK
      READ(88,*) sigmapro
      READ(88,*) maxf

!      NEV = 50000
!.........................      
      step = 30./40.
      pi = 2.*asin(1.)
      
      ntherm = 1340
      xnevt = float(NEV)
      nremq = 0
      numq = 0
      kpi = 0
      kK = 0
      kpro = 0
      massud=0.33
      masss=0.450
C.... initialization of binning ..........
      do ini=1, 40
         do jni=1, 2
            xnbin(ini,jni) = 0.
            xnbinK(ini,jni) = 0.
            xnbinpro(ini,jni) = 0.
         enddo
      enddo

!      do itherm=1, ntherm
!         READ(15,*) Kth(itherm),Pth(itherm,1),Pth(itherm,2),
!     .        Pth(itherm,3),pth(itherm,4),Xth(itherm,1),
!     .        Xth(itherm,2),Xth(itherm,3),Xth(itherm,4)       !KTh thermal parton ID  
!      enddo
      
      DO IEV=1, NEV             !...event loop ..............
!..... read thermal parton from hydro !added by wenbin 2018.11.26...
         READ(45,*) mid,ntherm,mid,mid
      do itherm=1, ntherm
         READ(45,*)mid,Kth(itherm),Pth(itherm,1),
     .        Pth(itherm,2),Pth(itherm,3),Pth(itherm,4),amid,
     .  Xth(itherm,1),Xth(itherm,2),Xth(itherm,3),Xth(itherm,4)       !KTh thermal parton ID  
      enddo
!..... read LBT parton from LBT result !added by wenbin 2018.11.26
!.....number of partons .........................
         READ(30,*)kev,  numparton !kev the event id, numparton: number of partons wenbin 2018.11.23
         nchth = 0
!.....  Shower partons ..........................
         DO ish=1, numparton
            READ(30,*) idsh(ish),P(ish,1),P(ish,2),
     .           P(ish,3),amid,mid,XV(ish,1),XV(ish,2),XV(ish,3),
     .           XV(ish,4)!,amid,amid,amid,amid,amid,amid
          !write(*,*)idsh(ish),P(ish,1),P(ish,2)
         if((abs(idsh(ish)).eq.1).or.(abs(idsh(ish)).eq.2))then
              P(ish,4)=sqrt(P(ish,1)*P(ish,1)+P(ish,2)*P(ish,2)
     .         + P(ish,3)*P(ish,3)+massud*massud)
         !endif
         else
              if(abs(idsh(ish)).eq.3)then
                P(ish,4)=sqrt(P(ish,1)*P(ish,1)+P(ish,2)*P(ish,2)
     .           + P(ish,3)*P(ish,3)+masss*masss)
              else
                P(ish,4)=amid
              endif
         endif
        !write(*,*)ntherm,numparton
!!.....number of partons .........................
!         READ(9,*) kev, numparton !kev the event id, numparton: number of partons wenbin 2018.11.23
!         nchth = 0
!!.....  Shower partons ..........................
!         DO ish=1, numparton
!            READ(10,*) kev,idsh(ish),Kp(ish,3),P(ish,1),P(ish,2),
!     .           P(ish,3),P(ish,4),XV(ish,1),XV(ish,2),XV(ish,3),
!     .           XV(ish,4)
!idsh:id, 

!..... counting number of quarks ..................
            if(idsh(ish).eq.21) then
               numq = numq + 2
            else
               numq = numq + 1
            endif
         enddo
        !write(*,*)ntherm,numparton

          !write(*,*)"zzzzzzzz"
!         call hadronization(ntherm,Kth,Xth,Pth,numparton,idsh,XV,P,nH,
!     .        KfH,XH,PH,Kst)
         call hadronization(ntherm,Kth,Xth,Pth,numparton,idsh,XV,P,nH,
     .        KfH,XH,PH,Kst,Npt,Kpx,Ppx,Xpx, !added by wenbin for outputing remnant jet partons
     .        NCTnr, NCTKr,CTPr,CTXr)!CTnr, NCTKr,CTPr,CTXr the information of the coalesced hadrons Wenbin
!********* write the hadrons ********
           !write(*,*)"NCTnr,CTPr=1 ",NCTnr,CTPr(NCTnr,1) 
            write(51,17)IEV,nH
C.... Binning the produced hadrons ........................
         do ihad=1, nH
!            pT = sqrt(PH(ihad,1)**2.+PH(ihad,2)**2.)
            write(51,16) IEV,KfH(ihad),PH(ihad,1),PH(ihad,2),
     .           PH(ihad,3),PH(ihad,4),XH(ihad,1),XH(ihad,2),XH(ihad,3),
     .           XH(ihad,4),Kst(ihad)!Kst is origin of the hadrons(th-th:0, sh-th:1, sh-sh:2,frag:3)

         enddo

!********* write the Remnant partons ******** 
            write(52,17)IEV,Npt
C.... Binning the writing ........................
         do ihad=1, Npt
!            pT = sqrt(PH(ihad,1)**2.+PH(ihad,2)**2.)
            write(52,18) IEV,Kpx(ihad),Ppx(ihad,1),Ppx(ihad,2),
     .           Ppx(ihad,3),Ppx(ihad,4),Xpx(ihad,1),Xpx(ihad,2),
     .           Xpx(ihad,3),Xpx(ihad,4)
            !write(*,*) Ppx(ihad,3),Ppx(ihad,4)
         enddo

!********* write the coaleced thermal  partons ********
            write(53,*)"# ",IEV,NCTnr
C.... Binning the produced hadrons ........................
        if(NCTnr.eq.0)then
             write(53,*)1,0,0,0,0,0,0,0,0,0
             write(53,*)1,0,0,0,0,0,0,0,0,0

        endif   

         do ihad=1, NCTnr
!            pT = sqrt(PH(ihad,1)**2.+PH(ihad,2)**2.)
         !....... transfer the t,z to tau, etas...wenbin
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

