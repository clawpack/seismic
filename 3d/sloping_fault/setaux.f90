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

    real(kind=8) :: xcell, ycell, zcell
    real(kind=8) :: xpcorn(4), ypcorn(4), zpcorn(4), mag
    integer :: i,j,k

    ! Auxiliary variables:
    !       1 slip
    !       2 capacity function value

    ! Loop over all cells
    do k=1-mbc,mz + mbc
      zcell = zlower + (k-0.5d0)*dz
      do j=1-mbc,my + mbc
        ycell = ylower + (j-0.5d0)*dy
        do i=1-mbc,mx + mbc
          xcell = xlower + (i-0.5d0)*dx

          ! initialize slip to zero
          aux(1,i,j,k) = 0.d0

          ! compute capacity function value ( volume = (vol of lower face in y direction)*dy )
          call mapc2p(xcell - 0.5d0*dx, ycell - 0.5d0*dy, zcell - 0.5d0*dz, xpcorn(1), ypcorn(1), zpcorn(1))
          call mapc2p(xcell - 0.5d0*dx, ycell - 0.5d0*dy, zcell + 0.5d0*dz, xpcorn(2), ypcorn(2), zpcorn(2))
          call mapc2p(xcell + 0.5d0*dx, ycell - 0.5d0*dy, zcell + 0.5d0*dz, xpcorn(3), ypcorn(3), zpcorn(3))
          call mapc2p(xcell + 0.5d0*dx, ycell - 0.5d0*dy, zcell - 0.5d0*dz, xpcorn(4), ypcorn(4), zpcorn(4))
          ! only need area ratio from cross-product of diagonals, which will point in the y direction
          mag = dabs((xpcorn(3) - xpcorn(1))*(zpcorn(2) - zpcorn(4)) - (xpcorn(2) - xpcorn(4))*(zpcorn(3) - zpcorn(1)))
          aux(2,i,j,k) = 0.5d0*mag/(dx*dz)


        end do
      end do
    end do

end subroutine setaux
