
subroutine b4step3(mbc,mx,my,mz,meqn,q,xlower,ylower,zlower, &
    dx,dy,dz,t,dt,maux,aux)

    ! Called before each call to step3.
    ! Use to set time-dependent aux arrays or perform other tasks.
    !
    ! This default version does nothing.

    implicit none
    integer, intent(in) :: mbc,mx,my,mz,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

    real (kind=8) :: center(3), theta, xcb(2), ycb(2), mindepth
    common /fault/  center, theta, xcb, ycb, mindepth

    integer :: i, j, k
    real(kind=8) :: xcell, ycell, zcell

    do k = 1-mbc,mz+mbc
      zcell = zlower + (k-0.5d0)*dz
      if (abs(zcell - 0.50*dz - center(3)) < 0.5d0*dz) then
        do j = 1-mbc,my+mbc
          ycell = ylower + (j-0.5d0)*dy
          if (ycb(1) <= ycell - 0.5d0*dy .and. ycell + 0.5d0*dy <= ycb(2)) then
            do i = 1-mbc,mx+mbc
              xcell = xlower + (i-0.5d0)*dx
              aux(1,i,j,k) = 0.d0
              if (t <= 1.d0 .and. xcb(1) <= xcell - 0.5d0*dx .and. xcell + 0.5d0*dx <= xcb(2)) then
                aux(1,i,j,k) = 1.d0
              end if
            end do
          end if
        end do
      end if
    end do

end subroutine b4step3
