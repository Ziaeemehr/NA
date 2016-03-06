C======================================================================
C PROGRAM: ising.f 
C TYPE   : main
C PURPOSE: simulate 2D Ising model with Glauber or Metropolis dynamic
C VERSION: 05 February 2001
!echo 'plot "magnet" u 1:2 w l title " " ' | gnuplot --persist

C======================================================================
        parameter(Ntot=10,nMC=1000,beta=0.3,Ndyn=1)
        dimension nspin(0:Ntot-1,0:Ntot-1)
        nran(i) = mod(int(i*ran5(idum)),i)
        open(1,file='magnet')
C   Initialisation of the spins
        nmagnet = 0
        do i=0,Ntot-1
           do j=0,Ntot-1
              Nspin(i,j) = 1 
              nmagnet = nmagnet + 1
           end do
        end do
        idum = -4856
C   Monte Carlo steps
        do iter=1,nMC
           do iter2=1,Ntot**3
C   Random choice of the spin
              imod = nran(Ntot)
              jmod = nran(Ntot)
C   Energy difference calculation
              imodp = mod(imod+1,Ntot)
              imodm = mod(imod-1+Ntot,Ntot)
              jmodp = mod(jmod+1,Ntot)
              jmodm = mod(jmod-1+Ntot,Ntot)
              nrjdif = 2*Nspin(imod,jmod)*(Nspin(imodp,jmod)+
     %                  Nspin(imodm,jmod)+Nspin(imod,jmodp)+
     %                  Nspin(imod,jmodm)) 
	      if (Ndyn.eq.0) then
C   Glauber dynamic
                 xprob = exp(-beta*nrjdif)
                 xprob = xprob/(1.+xprob)
              else if (Ndyn.eq.1) then
C   Metropolis dynamic
                 if (exp(-beta*nrjdif).lt.1.) then
                    xprob = exp(-beta*nrjdif)
                 else
                    xprob = 1.
                 end if
              end if
              if (ran5(idum).lt.xprob) then
                 Nspin(imod,jmod) = -Nspin(imod,jmod)
                 nmagnet = nmagnet + 2*Nspin(imod,jmod)
              end if
           end do
           xmagnet = real(nmagnet)/real(Ntot**2)
           write(1,*)iter*Ntot,xmagnet
        end do
	end              
C======================================================================
C PROGRAM: ran5.f
C TYPE   : function
C PURPOSE: generate random numbers
C COMMENT: Initialize idum with negative integer
C======================================================================
      REAL FUNCTION RAN5(IDUM)
      INTEGER IDUM
      INTEGER IA, IM, IQ, IR, NTAB
      REAL    AM, ATAB
      PARAMETER (IA=16807, IM=2147483647, AM=1.0/IM)
      PARAMETER (IQ=127773, IR=2836, NTAB=32, ATAB=NTAB-1)
      INTEGER J, K
      REAL V(NTAB), Y
      SAVE V, Y
      DATA V/NTAB*0.0/, Y/0.5/
      IF (IDUM.LE.0) THEN
          IDUM = MAX(-IDUM,1)
          DO 12 J=NTAB,1,-1
              K = IDUM/IQ
              IDUM = IA*(IDUM-K*IQ) - IR*K
              IF (IDUM.LT.0) IDUM = IDUM+IM
              V(J) = AM*IDUM
   12     CONTINUE
          Y = V(1)
      END IF
    1 CONTINUE
          K = IDUM/IQ
          IDUM = IA*(IDUM-K*IQ) - IR*K
          IF (IDUM.LT.0) IDUM = IDUM+IM
          J = 1 + INT(ATAB*Y)
          Y = V(J)
          RAN5 = Y
          V(J) = AM*IDUM
      IF (RAN5.EQ.0.0 .OR. RAN5.EQ.1.0) GOTO 1
      RETURN
      END
