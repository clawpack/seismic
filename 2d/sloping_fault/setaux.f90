subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)

    implicit none
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(out) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8) :: xcell, ycell, cp, cs
    integer :: i,j

    real(kind=8) :: rho_cell, lambda_cell, mu_cell

    ! Arrays to temporarily store computational and physical corners of grid cells
    real(kind=8) :: xccorn(4),yccorn(4),xpcorn(4),ypcorn(4)
    real(kind=8) :: norm, xn, yn, areap, b2c

    real (kind=8) :: center(2), theta, xcb(2), mindepth
    common /fault/  center, theta, xcb, mindepth


! c     #   (lambda = nu*E/((1+nu)(1-2nu))), E=young modulus, nu=poisson ration
! c     #   aux(1,i,j) is the density of the elastic material
! c     #   aux(2,i,j) is the lambda Lame parameter for elasticity
! c     #   aux(3,i,j) is the mu Lame parameter for elasticity (shear)
!
! c     #   aux(4,i,j) is cp is the P-wave speed (pressure wave)
! c     #   aux(5,i,j) is cs is the S-wave speed (shear wave)
!
! c     #   aux(6,i,j) is x component of normal at "left" boundary of grid point (i,j)
! c     #   aux(7,i,j) is y component of normal at "left" boundary of grid point (i,j)
! c     #   aux(8,i,j) =  is ratio of left edge to dy
!
! c     #   aux(9,i,j) is x component of normal at "bottom" boundary of grid point (i,j)
! c     #   aux(10,i,j) is y component of normal at "bottom" boundary of grid point (i,j)
! c     #   aux(11,i,j) =  is ratio of bottom edge to dx
!
! c     #   aux(12,i,j) = kappa  is ratio of cell area to dx*dy
!
! c     #   aux(13,i,j) = slip:
    !

    lambda_cell = 60.d9  ! Pa
    mu_cell = 30.d9      ! Pa
    rho_cell = 2500.d0   ! kg/m**3
    cp = dsqrt((lambda_cell + 2*mu_cell)/rho_cell)
    cs = dsqrt(mu_cell/rho_cell)

    ! Loop over all cells
      do j=1-mbc,my + mbc
        ycell = ylower + (j-0.5d0)*dy
        do i=1-mbc,mx + mbc

          aux(1,i,j) = rho_cell
          aux(2,i,j) = lambda_cell
          aux(3,i,j) = mu_cell

          ! Calculate pressure and shear wave speeds
          aux(4,i,j) = cp
          aux(5,i,j) = cs

          ! ============================================
          !BEGINS CALCULATING AND SAVING NORMALS, LENGTH RATIOS AND AREA RATIO FOR MAPPED GRID

          ! Computes corners of grid cell
! c           # lower left corner:
          xccorn(1) = xlower + float(i-1)*dx
          yccorn(1) = ylower + float(j-1)*dy
          call mapc2p(xccorn(1),yccorn(1),xpcorn(1),ypcorn(1))

! c           # upper left corner:
          xccorn(2) = xccorn(1)
          yccorn(2) = yccorn(1) + dy
          call mapc2p(xccorn(2),yccorn(2),xpcorn(2),ypcorn(2))

! c           # upper right corner:
          xccorn(3) = xccorn(1) + dx
          yccorn(3) = yccorn(1) + dy
          call mapc2p(xccorn(3),yccorn(3),xpcorn(3),ypcorn(3))

! c           # lower right corner:
          xccorn(4) = xccorn(1) + dx
          yccorn(4) = yccorn(1)
          call mapc2p(xccorn(4),yccorn(4),xpcorn(4),ypcorn(4))

          ! Compute inner normals
          !left
          xn = ypcorn(2) - ypcorn(1)
          yn = -(xpcorn(2) - xpcorn(1))
          if (dabs(yn) < 1e-10) then
            yn = 0.0
          end if
          norm = dsqrt(xn*xn + yn*yn)
          aux(6,i,j) = xn/norm
          aux(7,i,j) = yn/norm
          aux(8,i,j) = norm/dy

          !bottom
          xn = -(ypcorn(4) - ypcorn(1))
          yn = xpcorn(4) - xpcorn(1)
          if (dabs(xn) < 1e-10) then
            xn = 0.0
          end if
          norm = dsqrt(xn*xn + yn*yn)
          aux(9,i,j) = xn/norm
          aux(10,i,j) = yn/norm
          aux(11,i,j) = norm/dx

          ! Computes area of grid cell using cross product of diagonal
          areap = (xpcorn(3) - xpcorn(1))*(ypcorn(2) - ypcorn(4))
          areap = areap - (ypcorn(3)-ypcorn(1))*(xpcorn(2)-xpcorn(4))
          areap = 0.5*dabs(areap)
          aux(12,i,j) = areap/(dx*dy)

          ! set slip:
          xccorn(1) = xlower + float(i-1)*dx
          yccorn(1) = ylower + float(j-1)*dy
          if ((abs(yccorn(1)-center(2)) < 0.5d0*dy) .and. &
              (xccorn(1) >= xcb(1)) .and. (xccorn(1) <= xcb(2))) then
                aux(13,i,j) = exp(-((xccorn(1)-center(1))/(center(1)-xcb(1))**2)
            else
                aux(13,i,j) = 0.d0
            endif

        end do
      end do

end subroutine setaux
