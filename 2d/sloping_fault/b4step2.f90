
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)

!   set slip before call to step2

    use fault_module, only: center, xcb, nsubfaults, subfaults, final_rupture_rise_time, LAT2METER

    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer :: i, j, k 
    real(kind=8) :: xcell, ycell

    aux(13,:,:) = 0.d0
    if (t <= final_rupture_rise_time) then

      do j=1-mbc,my+mbc
        ycell = ylower + (j-0.5d0)*dy
        if (abs(ycell - 0.5d0*dy - center(2)) < 0.5d0*dy) then

          do i=1-mbc,mx+mbc
            xcell = xlower + (i-0.5d0)*dx
            if (xcb(1) <= xcell - 0.5d0*dx .and. xcell + 0.5d0*dx <= xcb(2)) then

              do k=1,nsubfaults
                if (subfaults(k)%longitude*LAT2METER <= xcell .and. &
                    xcell <= subfaults(k)%longitude*LAT2METER + subfaults(k)%width .and. &
                    subfaults(k)%rupture_time <= t .and. &
                    t <= subfaults(k)%rupture_time + subfaults(k)%rise_time) then
                  aux(13,i,j) = subfaults(k)%slip/subfaults(k)%rise_time
                  exit
                end if
              end do

            end if
          end do
 
        end if
      end do

    end if

end subroutine b4step2
