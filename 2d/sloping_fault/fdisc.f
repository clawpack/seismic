c
c
c
c     =================================================
      function fdisc(x,y)
c     =================================================

      implicit double precision (a-h,o-z)

c
c     # for computing cell averages for initial data that has a
c     # discontinuity along some curve.  fdisc should be negative to the 
c     # left of the curve and positive to the right

c     # three layers with one linear and one quadratic interface
c     # fdisc is positive in the middle layer, negative in the other two

      if (y .gt. 0.5d0) then
           fdisc = 0.6d0 + 0.1d0*x - y
        else
           fdisc = y - 0.4d0 - 0.3*(x-0.5d0) + 0.2d0*(x-0.5d0)**2 
        endif
c
      return
      end

