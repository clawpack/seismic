!
! ------------------------------------------------------------------
!
      subroutine bc3amr(val,aux,nrow,ncol,nfil,meqn,naux, &
                       hx, hy, hz,level, time, &
                       xleft,  xright,  yfront, yrear, &
                       zbot, ztop, &
                       xlower,ylower,zlower, &
                       xupper,yupper,zupper, &
                       xperiodic, yperiodic,zperiodic)


!
!
! :::::::::: BC3AMR ::::::::::::::::::::::::::::::::::::::::::::::;
!
!     Take a grid patch with mesh widths hx,hy,hz, of dimensions nrow by
!     ncol by nfil,  and set the values of any piece of
!     of the patch which extends outside the physical domain
!     using the boundary conditions.
!     ------------------------------------------------
!     # Standard boundary condition choices for amr3ez in clawpack
!
!     # At each boundary  k = 1 (left),  2 (right),  3 (front), 4 (rear),
!                         5 (bottom) 6 (top)
!     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
!     #            =  1  for zero-order extrapolation
!     #            =  2  for periodic boundary coniditions
!     #            =  3  for solid walls, assuming this can be implemented
!     #                  by reflecting the data about the boundary and then
!     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
!     #                  or 4'th (for k = 5,6) component of q.
!     ------------------------------------------------
!
!     The corners of the grid patch are at
!        (xleft,yfront,zbot)  --  lower front left corner
!        (xright,yrear,ztop) --  upper rear right corner

!     The physical domain itself is a rectangular parallelopiped bounded by
!        (xlower,ylower,zlower)  -- lower front left corner
!        (xupper,yupper,zupper)  -- upper rear right corner
!
!     the picture is the following:
!
!                            __________________________(xupper,yupper,zupper)
!                           /                         /|
!                          /                         / |
!                         /                         /  |
!                        /_________________________/   |
!                        |                         |   |
!                        |                         |   |
!                     ___|_____(xright,yrear,ztop) |   |
!                    /___|____/|                   |   |
!                    |   |    ||                   |   |
!                    |   |    ||                   |   |
!                    |   |    ||                   |   |
!                    |___|____|/                   |   |
!  (xleft,yfront,zbot)   |                         |  /
!                        |                         | /
!                        |_________________________|/
!  (xlower,ylower,zlower)
!
!     Any cells that lie outside the physical domain are ghost cells whose
!     values should be set in this routine.  This is tested for by comparing
!     xleft with xlower to see if values need to be set at the left, as in
!     the figure above, and similarly at the other boundaries.
!
!     Patches are guaranteed to have at least 1 row of cells filled
!     with interior values so it is possible to  extrapolate.
!     Fix trimbd if you want more than 1 row pre-set.
!
!     Make sure the order the boundaries are specified is correct
!     so that diagonal corner cells are also properly taken care of.
!
!     Periodic boundaries are set before calling this routine, so if the
!     domain is periodic in one direction only you
!     can safely extrapolate in the other direction.
!
!     Don't overwrite ghost cells in periodic directions!
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      use amr_module, only:  mthbc
      implicit none

      integer nrow,ncol,nfil,meqn,naux,level
      integer i,j,k,m,nxl,nxr,ibeg,nyf,nyr,jbeg,nzb,nzt,kbeg

      real (kind=8)  val(meqn,nrow,ncol,nfil)
      real (kind=8)  aux(naux,nrow,ncol,nfil)
      logical xperiodic, yperiodic, zperiodic
      real (kind=8)  hx,hy,hz,hxmarg,hymarg,hzmarg,time
      real (kind=8)  xleft,xright,yfront,yrear,zbot,ztop
      real (kind=8)  xlower,ylower,zlower,xupper,yupper,zupper

      REAL (kind=8) s, xcell, ycell, zcell

      REAL (kind=8) :: rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer
      common /cparam/  rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer

      REAL (kind=8) :: t0wall, amplitude, pulse_span, trans_halfdepth, trans_halfwidth
      common /combc/ t0wall, amplitude, pulse_span, trans_halfdepth, trans_halfwidth

      hxmarg = hx*.01d0
      hymarg = hy*.01d0
      hzmarg = hz*.01d0

      if (xperiodic .and. yperiodic .and. zperiodic) go to 699
!
!
!-------------------------------------------------------
!     # xlower boundary:
!-------------------------------------------------------
      if (xleft .ge. xlower-hxmarg) then
!        # not a physical boundary -- ghost cells lie within another
!        # grid and values are set elsewhere in amr code.
         go to 199
         endif
!
!     # number of ghost cells lying outside physical domain:
      nxl = (xlower+hxmarg-xleft)/hx
!
      go to (100,110,120,130) mthbc(1)+1
!
  100 continue
  !     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc3amr'
      stop
      go to 199

  110 continue
!     # zero-order extrapolation:
         do 115 k = 1,nfil
            do 115 j = 1,ncol
               do 115 i=1,nxl
                  do m=1,meqn
                     val(m,i,j,k) = val(m,nxl+1,j,k)
                  end do
           115 continue
      go to 199

  120 continue
!     # periodic:   handled elsewhere in amr
      go to 199

  130 continue
!     # solid wall:
      do 135 k = 1,nfil
       do 135 j = 1,ncol
        do 135 i=1,nxl
         do 135 m=1,meqn
               val(m,i,j,k) = val(m,2*nxl+1-i,j,k)
  135    continue
!     # negate the normal velocity:
      do 136 i=1,nxl
         do 136 j = 1,ncol
            do 136 k = 1,nfil
               val(7,i,j,k) = -val(7,i,j,k)
  136       continue
      go to 199

  199 continue
!
!-------------------------------------------------------
!     # xupper boundary:
!-------------------------------------------------------
      if (xright .le. xupper+hxmarg) then
!        # not a physical boundary -- ghost cells lie within another
!        # grid and values are set elsewhere in amr code.
         go to 299
         endif
!
!     # number of ghost cells lying outside physical domain:
      nxr = (xright - xupper + hxmarg)/hx
      ibeg = max0(nrow-nxr+1, 1)
!
      go to (200,210,220,230) mthbc(2)+1
!
  200 continue
!     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc3amr'
      stop
      go to 299

  210 continue
!     # zero-order extrapolation:
       do 215 k = 1,nfil
        do 215 j = 1,ncol
         do 215 i=ibeg,nrow
          do 215 m=1,meqn
             val(m,i,j,k) = val(m,ibeg-1,j,k)
  215     continue
      go to 299

  220 continue
!     # periodic:   handled elsewhere in amr
      go to 299

  230 continue
!     # solid wall:
      do 235 k = 1,nfil
       do 235 j = 1,ncol
        do 235 i=ibeg,nrow
         do 235 m=1,meqn
            val(m,i,j,k) = val(m,2*ibeg-1-i,j,k)
  235    continue
!     # negate the normal velocity:
       do 236 k = 1,nfil
       do 236 j = 1,ncol
       do 236 i=ibeg,nrow
          val(7,i,j,k) = -val(7,i,j,k)
  236  continue
      go to 299

  299 continue
!
!-------------------------------------------------------
!     # ylower boundary:
!-------------------------------------------------------
      if (yfront .ge. ylower-hymarg) then
!        # not a physical boundary -- ghost cells lie within another
!        # grid and values are set elsewhere in amr code.
         go to 399
         endif
!
!     # number of ghost cells lying outside physical domain:
      nyf = (ylower+hymarg-yfront)/hy
!
      go to (300,310,320,330) mthbc(3)+1
!
  300 continue
!     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc3amr'
      stop
      go to 399
!
  310 continue
!     # zero-order extrapolation:
      do 315 k = 1,nfil
       do 315 j=1,nyf
        do 315 i=1,nrow
         do 315 m=1,meqn
            val(m,i,j,k) = val(m,i,nyf+1,k)
  315 continue
      go to 399

  320 continue
!     # periodic:   handled elsewhere in amr
      go to 399

  330 continue
!     # solid wall:
      do 335 k = 1,nfil
       do 335 j=1,nyf
        do 335 i=1,nrow
         do 335 m=1,meqn
            val(m,i,j,k) =  val(m,i,2*nyf+1-j,k)
  335    continue
!     # negate the normal velocity:
      do 336 k = 1,nfil
       do 336 j=1,nyf
        do 336 i=1,nrow
           val(8,i,j,k) = -val(8,i,j,k)
  336   continue
      go to 399

  399 continue
!
!-------------------------------------------------------
!     # yupper boundary:
!-------------------------------------------------------
      if (yrear .le. yupper+hymarg) then
!        # not a physical boundary -- ghost cells lie within another
!        # grid and values are set elsewhere in amr code.
         go to 499
         endif
!
!     # number of ghost cells lying outside physical domain:
      nyr = (yrear - yupper + hymarg)/hy
      jbeg = max0(ncol-nyr+1, 1)
!
      go to (400,410,420,430) mthbc(4)+1
!
  400 continue
!     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc3amr'
      stop
      go to 499

  410 continue
!     # zero-order extrapolation:
       do 415 k = 1,nfil
        do 415 j=jbeg,ncol
         do 415 i=1,nrow
          do 415 m=1,meqn
             val(m,i,j,k) =  val(m,i,jbeg-1,k)
  415     continue
      go to 499

  420 continue
!     # periodic:   handled elsewhere in amr
      go to 499

  430 continue
!     # solid wall:
      do 435 m=1,meqn
         do 435 j=jbeg,ncol
            do 435 i=1,nrow
               do 435 k = 1,nfil
                  val(m,i,j,k) =  val(m,i,2*jbeg-1-j,k)
  435          continue
!     # negate the normal velocity:
      do 436 j=jbeg,ncol
         do 436 i=1,nrow
            do 436 k = 1,nfil
               val(8,i,j,k) = -val(8,i,j,k)
  436       continue
      go to 499

  499 continue

!
!-------------------------------------------------------
!     # zlower boundary:
!-------------------------------------------------------
      if (zbot .ge. zlower-hzmarg) then
!        # not a physical boundary -- ghost cells lie within another
!        # grid and values are set elsewhere in amr code.
         go to 599
         endif
!
!     # number of ghost cells lying outside physical domain:
      nzb = (zlower+hzmarg-zbot)/hz
!
      go to (500,510,520,530) mthbc(5)+1
!
  500 continue
!     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(5)=0 and no BCs specified in bc3amr'
      stop
      go to 599
!
  510 continue
!     # zero-order extrapolation:
      do 515 k=1,nzb
       do 515 j = 1,ncol
        do 515 i=1,nrow
         do 515 m=1,meqn
            val(m,i,j,k) = val(m,i,j,nzb+1)
  515    continue
      go to 599

  520 continue
!     # periodic:   handled elsewhere in amr
      go to 599

  530 continue
!     # solid wall:
      do 535 k=1,nzb
       do 535 j = 1,ncol
        do 535 i=1,nrow
         do 535 m=1,meqn
            val(m,i,j,k) =  val(m,i,j,2*nzb+1-k)
  535    continue
!     # negate the normal velocity:
      do 536 k = 1,nzb
       do 536 j = 1,ncol
        do 536 i = 1,nrow
           val(9,i,j,k) = -val(9,i,j,k)
  536   continue
      go to 599

  599 continue
!
!-------------------------------------------------------
!     # zupper boundary:
!-------------------------------------------------------
      if (ztop .le. zupper+hzmarg) then
!        # not a physical boundary -- ghost cells lie within another
!        # grid and values are set elsewhere in amr code.
         go to 699
         endif
!
!     # number of ghost cells lying outside physical domain:
      nzt = (ztop - zupper + hzmarg)/hz
      kbeg = max0(nfil-nzt+1, 1)
!
      go to (600,610,620,630) mthbc(6)+1
!
  600 continue

!     # first-order extrapolation:
      do 605 k = kbeg,nfil
        do 605 j = 1,ncol
          do 605 i = 1,nrow
            do 605 m = 1,meqn
              val(m,i,j,k) =  val(m,i,j,2*kbeg-1-k)
  605       continue

!    # negate the normal stress:
     do 606 k = kbeg,nfil
       do 606 j = 1,ncol
         do 606 i = 1,nrow
           val(3,i,j,k) = -val(3,i,j,k)
           val(5,i,j,k) = -val(5,i,j,k)
           val(6,i,j,k) = -val(6,i,j,k)
 606     continue

      go to 699

  610 continue
!     # zero-order extrapolation:
      do 615 k = kbeg,nfil
         do 615 j = 1,ncol
            do 615 i = 1,nrow
               do 615 m = 1,meqn
                  val(m,i,j,k) =  val(m,i,j,kbeg-1)
  615          continue
      go to 699

  620 continue
!     # periodic:   handled elsewhere in amr
      go to 699

  630 continue
!     # solid wall:
       do 635 k = kbeg,nfil
        do 635 j = 1,ncol
         do 635 i = 1,nrow
          do 635 m=1,meqn
             val(m,i,j,k) =  val(m,i,j,2*kbeg-1-k)
  635     continue
!     # negate the normal velocity:
      do 636 k=kbeg,nfil
       do 636 j = 1,ncol
        do 636 i=1,nrow
           val(9,i,j,k) = -val(9,i,j,k)
  636   continue
      go to 699

  699 continue

      return
      end
