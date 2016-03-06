c TWO-DIMENSIONAL ISING MODEL: "Maybe it works, maybe it doesn't."
       program ising2

c INITIALIZE: nn*nn is size, Temp is temperature in units of J/k_B, warm
c is # of time steps for warm-up, mcs is total # time steps, ss(nn,nn)
c are the spins, ibr keeps track of spin neighbors.  Then open
c file(s) to write stuff on.  Then initial conditions for spins. Then
c initialize junk variables for averages.  Then neighbour table to do
c periodic boundary conditions.

       real ran3,aran,Temp,prob,energy,mag,fooe,foom
       integer iseed,itime,nn,nnn,ii,ix,iy,warm,mcs
       parameter (nn=20)
       parameter (Temp=1.0)
       parameter (warm=1)
       parameter (mcs=100)
       parameter (iseed=-12888333)
       integer ss(nn,nn),ibr(nn,2)

       open (77,*)

       do 6 ix=1,nn 
       do 5 iy=1,nn 
       ss(ix,iy) = 1 
5      continue
6      continue

       fooe=0.0
       foom=0.0

       do 5000 ix=1,nn
       ibr(ix,1)=ix-1
       ibr(ix,2)=ix+1
       if(ix.eq.1)ibr(ix,1)=nn
       if(ix.eq.nn)ibr(ix,2)=1
5000      continue

c LOOP: Total of "mcs" monte carlo steps with "warm" to warm up, 
c then calculate average energy and magnetization

       do 111 itime=1,warm
       call mcmove(ss,nn,Temp,iseed,ibr)
111    continue


       do 112 itime=1,mcs-warm
       call mcmove(ss,nn,Temp,iseed,ibr)
          do 90 ix = 1,nn
          do 89 iy = 1,nn
              fooe = fooe-ss(ix,iy)*ss(ix,ibr(iy,1))
     &                 -ss(ix,iy)*ss(ibr(ix,2),iy)
          foom = foom + ss(ix,iy)
89        continue
90        continue
112    continue
      
       energy = fooe/(float(nn*nn*(mcs-warm)))
       mag = foom/(float(nn*nn*(mcs-warm)))

       if (Temp.ge.2.269) then
       foom=0.0
       else
       foom = (1.0 - (sinh (2.0/Temp) )**(-4.0) )**(1.0/8.0)
       endif

       write (6,*) 'Energy/spin, Mag/spin',energy,mag
       write (6,*) 'Analytically: -2 to 0, and', foom
       write (6,*) 'Last configuration on file 77'
          do 909 ix = 1,nn
          do 899 iy = 1,nn
       write (77,*) ix,iy,ss(ix,iy)
899       continue
909       continue
      
       stop
       end

c ONE MONTE CARLO STEP by Metropolis: Flip probability 1 if Enew < Eold, 
c else prob is exp -(Enew-Eold)/T.  Simplified here since only there 
c are five cases in d=1 for external field = 0.
c FLIP WITH prob1 prob2  1.0   1.0   1.0   (Below spins called)
c             +     -     -     -     -           ss2
c            +++   +++   ++-   ++-   -+-      ss1 ss0 ss3
c             +     +     +     -     -           ss4

      subroutine mcmove(ss,nn,Temp,iseed,ibr)
      integer nn,iseed,ix,iy,iix,iiy,ss0,ss1,ss2,ss3,ss4,de
      integer ix1,ix2,iy1,iy2
      integer ss(nn,nn),ibr(nn,2)
      real ran3,prob1,prob2,Temp

      prob1 = exp(-8.0/Temp)
      prob2 = exp(-4.0/Temp)

      do 100 iix=1,nn
      do 17  iiy=1,nn
      ix = 1 + nn * ran3(iseed)
      iy = 1 + nn * ran3(iseed)
      ss0=ss(ix,iy)
      ss1=ss(ibr(ix,1),iy)
      ss2=ss(ix,ibr(iy,1))
      ss3=ss(ibr(ix,2),iy)
      ss4=ss(ix,ibr(iy,2))
      de=2*ss0*(ss1+ss2+ss3+ss4)
         if (de.eq.8.and.ran3(iseed).gt.prob1) then 
78       continue
         else if (de.eq.4.and.ran3(iseed).gt.prob2) then 
87       continue
         else
         ss(ix,iy) = -1*ss(ix,iy)
         endif
17    continue
100   continue

      return
      end

c  RANDOM # GENERATOR
      function ran3(idum)
c        implicit real*4(m)
c        parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      dimension ma(55)
      save ma
      save iff,inext,inextp

      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      return
      end
