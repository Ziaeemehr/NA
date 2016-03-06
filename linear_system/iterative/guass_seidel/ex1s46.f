c
c     Numerical Analysis:
c     Mathematics of Scientific Computing
c     Third Edition
c     D.R. Kincaid & E.W. Cheney
c     Brooks/Cole Publ., 2002
c     Copyright (c) 1996
c
c     Section 4.6
c
c     Example of Jacobi method and Gauss-Seidel method
c
c    
c     file: ex1s46.f
c
      data M/50/
c
      print *
      print *,' Jacobi method example'
      print *,' Section 4.6, Kincaid-Cheney'
      print *
c
      x1 = 0.
      x2 = 0.
      print *,'  k        x1             x2'
      print 4,0,x1,x2
      do 2 k=1,M
         y1 = (6./7.)*x2 + (3./7.)
         y2 = (8./9.)*x1 - (4./9.)
         x1 = y1
         x2 = y2
         if (0 .eq. mod(k,10)) print 4,k,x1,x2
 2    continue
c
c     Gauss-Seidel method
c
      print *
      print *
      print *,' Gauss-Seidel method example'
      print *,' Section 4.6, Kincaid-Cheney'
      print *
c
      x1 = 0.
      x2 = 0.
      print *,'  k        x1             x2'
      print 4,0,x1,x2
      do 3 k=1,M
         x1 = (6./7.)*x2 + (3./7.)
         x2 = (8./9.)*x1 - (4./9.)
         if (0 .eq. mod(k,10)) print 4,k,x1,x2
 3    continue
c
 4    format (3x,i2,2x,2(e13.6,2x))
      stop
      end
