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
c     Example of Gauss-Seidel method
c
c
c     file: ex3s46.f 
c
      dimension x(3)
      data x/0.,0.,0./
c
      print *
      print *,' Gauss-Seidel example'
      print *,' Section 4.6, Kincaid-Cheney'
      print *
c
      print *,'  k       x1             x2             x3'
      print 3,0,x
      do 2 k=1,20
         x(1) = 0.5*x(2) + 1.0
         x(2) = -(1.0/6.0)*x(1) + (1.0/3.0)*x(3) - (2.0/3.0)
         x(3) = -(0.5)*x(1) + (3.0/8.0)*x(2) + (5.0/8.0)
         print 3,k,x
 2    continue
c
 3    format (3x,i2,2x,3(e13.6,2x))
      stop
      end
