!*******************************************************************************
!   Copyright(C) 1999-2013 Intel Corporation. All Rights Reserved.
!   
!   The source code, information  and  material ("Material") contained herein is
!   owned  by Intel Corporation or its suppliers or licensors, and title to such
!   Material remains  with Intel Corporation  or its suppliers or licensors. The
!   Material  contains proprietary information  of  Intel or  its  suppliers and
!   licensors. The  Material is protected by worldwide copyright laws and treaty
!   provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way  without Intel's  prior  express written  permission. No  license
!   under  any patent, copyright  or  other intellectual property rights  in the
!   Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
!   implication, inducement,  estoppel or  otherwise.  Any  license  under  such
!   intellectual  property  rights must  be express  and  approved  by  Intel in
!   writing.
!   
!   *Third Party trademarks are the property of their respective owners.
!   
!   Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
!   this  notice or  any other notice embedded  in Materials by Intel or Intel's
!   suppliers or licensors in any way.
!
!*******************************************************************************
!  Content:
!      D D O T  Example Program Text
!*******************************************************************************

      program  DDOT_MAIN
*
      integer          n, incx, incy
      integer          xmax, ymax
      parameter       (xmax=20, ymax=20)
      double precision sx(xmax), sy(ymax), res
      integer          i
*       Intrinsic Functions
      intrinsic        abs
*       External Subroutines
      external         DDOT, PrintVectorD
      double precision DDOT
*
*       Executable Statements
*
      print*
      print*,'   D D O T  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n
      read*, incx, incy
      if ( ((1+(n-1)*abs(incx)).gt.xmax).or.
     $     ((1+(n-1)*abs(incy)).gt.ymax) ) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, (sx(i),i=1,1+(n-1)*abs(incx))
      read*, (sy(i),i=1,1+(n-1)*abs(incy))
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      call PrintVectorD(0,n,sx,incx,'X ')
      call PrintVectorD(0,n,sy,incy,'Y ')
*
*       Call DDOT subroutine
      res = DDOT(n,sx,incx,sy,incy)
*
      print*
      print*, '     OUTPUT DATA'
      print 101, res

      stop
 100  format(7x,'N=',i2)
 101  format(7x,'DDOT = ',f8.3)
 999  stop 1
      end
