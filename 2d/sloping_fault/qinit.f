
c
c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit double precision (a-h,o-z)
       dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
       dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

       do 20 i=1,mx
          xi = xlower + (i-0.5d0)*dx
          do 20 j=1,my
             yj = ylower + (j-0.5d0)*dy
             q(1,i,j) = 0.d0
             q(2,i,j) = 0.d0
             q(3,i,j) = 0.d0
             q(4,i,j) = 0.d0
             q(5,i,j) = 0.d0
 20          continue

       return
       end
