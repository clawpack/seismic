subroutine setaux(mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,maux,aux)

    ! Called at start of computation before calling qinit, and
    ! when AMR is used, also called every time a new grid patch is created.
    ! Use to set auxiliary arrays
    !   aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc).
    ! Note that ghost cell values may need to be set if the aux arrays
    ! are used by the Riemann solver(s).

    implicit none
    integer, intent(in) :: mbc,mx,my,mz,maux
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz
    real(kind=8), intent(out) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

    real(kind=8) :: xcell, ycell, zcell, cp, cs, rho_cell, lambda_cell, mu_cell
    real(kind=8) :: xpcorn(4), ypcorn(4), zpcorn(4), mag, zmin, zmax
    integer :: i,j,k

    real (kind=8) :: center(3), theta, xcb(2), ycb(2), mindepth
    common /fault/  center, theta, xcb, ycb, mindepth

    !# NEED AUX DEF HERE

    lambda_cell = 60.d9  ! Pa
    mu_cell = 30.d9      ! Pa
    rho_cell = 2500.d0   ! kg/m**3
    cp = dsqrt((lambda_cell + 2*mu_cell)/rho_cell)
    cs = dsqrt(mu_cell/rho_cell)


    ! Loop over all cells
    do k=1-mbc,mz + mbc
      zcell = zlower + (k-0.5d0)*dz
      do j=1-mbc,my + mbc
        ycell = ylower + (j-0.5d0)*dy
        do i=1-mbc,mx + mbc
          xcell = xlower + (i-0.5d0)*dx

          aux(1,i,j,k) = rho_cell
          aux(2,i,j,k) = lambda_cell
          aux(3,i,j,k) = mu_cell

          ! Calculate pressure and shear wave speeds
          aux(4,i,j,k) = cp
          aux(5,i,j,k) = cs

          ! compute mapping info for lower face in x direction
          call mapc2p(xcell - 0.5d0*dx, ycell - 0.5d0*dy, zcell - 0.5d0*dz, xpcorn(1), ypcorn(1), zpcorn(1))
          call mapc2p(xcell - 0.5d0*dx, ycell - 0.5d0*dy, zcell + 0.5d0*dz, xpcorn(2), ypcorn(2), zpcorn(2))
          call mapc2p(xcell - 0.5d0*dx, ycell + 0.5d0*dy, zcell + 0.5d0*dz, xpcorn(3), ypcorn(3), zpcorn(3))
          call mapc2p(xcell - 0.5d0*dx, ycell + 0.5d0*dy, zcell - 0.5d0*dz, xpcorn(4), ypcorn(4), zpcorn(4))
          ! normal information and area ratio are computed
          ! note normal will be in xz plane
          mag = ((ypcorn(3) - ypcorn(1))*(zpcorn(2) - zpcorn(4)) - (ypcorn(2) - ypcorn(4))*(zpcorn(3) - ypcorn(1)))**2
          mag = mag + ((xpcorn(3) - xpcorn(1))*(ypcorn(2) - ypcorn(4)) - (xpcorn(2) - xpcorn(4))*(ypcorn(3) - ypcorn(1)))**2
          mag = dsqrt(mag)
          aux(6,i,j,k) = ((ypcorn(3) - ypcorn(1))*(zpcorn(2) - zpcorn(4)) - (ypcorn(2) - ypcorn(4))*(zpcorn(3) - ypcorn(1)))/mag
          aux(7,i,j,k) = ((xpcorn(3) - xpcorn(1))*(ypcorn(2) - ypcorn(4)) - (xpcorn(2) - xpcorn(4))*(ypcorn(3) - ypcorn(1)))/mag
          aux(8,i,j,k) = 0.5d0*mag/(dy*dz)

          ! compute mapping info for lower face in y direction
          call mapc2p(xcell - 0.5d0*dx, ycell - 0.5d0*dy, zcell - 0.5d0*dz, xpcorn(1), ypcorn(1), zpcorn(1))
          call mapc2p(xcell - 0.5d0*dx, ycell - 0.5d0*dy, zcell + 0.5d0*dz, xpcorn(2), ypcorn(2), zpcorn(2))
          call mapc2p(xcell + 0.5d0*dx, ycell - 0.5d0*dy, zcell + 0.5d0*dz, xpcorn(3), ypcorn(3), zpcorn(3))
          call mapc2p(xcell + 0.5d0*dx, ycell - 0.5d0*dy, zcell - 0.5d0*dz, xpcorn(4), ypcorn(4), zpcorn(4))
          ! only need area ratio from cross-product of diagonals, which will point in the y direction
          mag = (xpcorn(2) - xpcorn(4))*(zpcorn(3) - zpcorn(1)) - (xpcorn(3) - xpcorn(1))*(zpcorn(2) - zpcorn(4))
          aux(9,i,j,k) = 0.5d0*mag/(dx*dz)

          ! compute mapping info for lower face in z direction
          call mapc2p(xcell - 0.5d0*dx, ycell - 0.5d0*dy, zcell - 0.5d0*dz, xpcorn(1), ypcorn(1), zpcorn(1))
          call mapc2p(xcell - 0.5d0*dx, ycell + 0.5d0*dy, zcell - 0.5d0*dz, xpcorn(2), ypcorn(2), zpcorn(2))
          call mapc2p(xcell + 0.5d0*dx, ycell + 0.5d0*dy, zcell - 0.5d0*dz, xpcorn(3), ypcorn(3), zpcorn(3))
          call mapc2p(xcell + 0.5d0*dx, ycell - 0.5d0*dy, zcell - 0.5d0*dz, xpcorn(4), ypcorn(4), zpcorn(4))
          ! for this face, the normal information is needed in addition to the area ratio
          ! the normal will reside in the xz plane
          mag = ((ypcorn(3) - ypcorn(1))*(zpcorn(2) - zpcorn(4)) - (ypcorn(2) - ypcorn(4))*(zpcorn(3) - ypcorn(1)))**2
          mag = mag + ((xpcorn(3) - xpcorn(1))*(ypcorn(2) - ypcorn(4)) - (xpcorn(2) - xpcorn(4))*(ypcorn(3) - ypcorn(1)))**2
          mag = dsqrt(mag)
          aux(10,i,j,k) = ((ypcorn(3) - ypcorn(1))*(zpcorn(2) - zpcorn(4)) - (ypcorn(2) - ypcorn(4))*(zpcorn(3) - ypcorn(1)))/mag
          aux(11,i,j,k) = ((xpcorn(3) - xpcorn(1))*(ypcorn(2) - ypcorn(4)) - (xpcorn(2) - xpcorn(4))*(ypcorn(3) - ypcorn(1)))/mag
          aux(12,i,j,k) = 0.5d0*mag/(dx*dy)

          ! compute capacity function value
          zmin = dmin1(zpcorn(1),zpcorn(2),zpcorn(3),zpcorn(4))
          mag = dmax1(zpcorn(1),zpcorn(2),zpcorn(3),zpcorn(4))
          call mapc2p(xcell - 0.5d0*dx, ycell - 0.5d0*dy, zcell + 0.5d0*dz, xpcorn(1), ypcorn(1), zpcorn(1))
          call mapc2p(xcell - 0.5d0*dx, ycell + 0.5d0*dy, zcell + 0.5d0*dz, xpcorn(2), ypcorn(2), zpcorn(2))
          call mapc2p(xcell + 0.5d0*dx, ycell + 0.5d0*dy, zcell + 0.5d0*dz, xpcorn(3), ypcorn(3), zpcorn(3))
          call mapc2p(xcell + 0.5d0*dx, ycell - 0.5d0*dy, zcell + 0.5d0*dz, xpcorn(4), ypcorn(4), zpcorn(4))
          zmax = dmax1(zpcorn(1),zpcorn(2),zpcorn(3),zpcorn(4))
          mag = 0.5d0*(zmax - zmin + dmin1(zpcorn(1),zpcorn(2),zpcorn(3),zpcorn(4)) - mag)
          aux(13,i,j,k) = mag/dz

          ! set fault slip:
          if ((abs(zcell+0.5d0*dz - center(3)) < 0.5d0*dz) .and. &
              (xcb(1) <= xcell) .and. (xcell <= xcb(2)) .and. &
              (ycb(1) <= ycell) .and. (ycell <= ycb(2))) then
                aux(14,i,j,k) = exp(-( &
                              ((xcell-center(1))/(center(1)-xcb(1)))**2 &
                              + ((ycell-center(2))/(center(2)-ycb(1)))**2 &
                              ))
          else
              aux(14,i,j,k) = 0.d0
          endif

        end do
      end do
    end do

end subroutine setaux
