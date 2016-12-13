
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)

    ! Called before each call to step2.
    ! Use to set time-dependent aux arrays or perform other tasks.
    !
    ! This default version does nothing.

    use fault_module, only: center, xcb, nsubfaults, subfaults

    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    real(kind=8) :: xcell, ycell, subfaultx

    integer :: i,j,subfaulti

    do j=1-mbc,my+mbc
      ycell = ylower + (j-0.5d0)*dy
      if (abs(ycell - 0.5d0*dy - center(2)) < 0.5d0*dy) then
        subfaulti = 1
        subfaultx = xcb(1)
        do i=1-mbc,mx+mbc
          xcell = xlower + (i-0.5d0)*dx
          aux(13,i,j) = 0.d0
          if (xcb(1) <= xcell - 0.5d0*dx .and. xcell + 0.5d0*dx <= xcb(2)) then
            ! check if current subfault needs to be updated
            do while(xcell > subfaultx + subfaults(subfaulti)%width)
              subfaultx = subfaultx + subfaults(subfaulti)%width
              subfaulti = subfaulti + 1
            end do

            if (subfaults(subfaulti)%rupture_time <= t .and. &
                  t <= subfaults(subfaulti)%rupture_time + subfaults(subfaulti)%rise_time) then
              aux(13,i,j) = subfaults(subfaulti)%slip/subfaults(subfaulti)%rise_time
            end if

          end if
        end do
      end if
    end do

end subroutine b4step2
