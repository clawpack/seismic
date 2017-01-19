
subroutine b4step3(mbc,mx,my,mz,meqn,q,xlower,ylower,zlower, &
    dx,dy,dz,t,dt,maux,aux)

!   set slip before call to step2

    use fault_module, only: center, xcb, ycb, nsubfaults, subfaults, final_rupture_rise_time, LAT2METER

    implicit none
    integer, intent(in) :: mbc,mx,my,mz,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

    integer :: i, j, k, l
    real(kind=8) :: xcell, ycell, zcell

    aux(1,:,:,:) = 0.d0
    if (t <= final_rupture_rise_time) then

      do k=1-mbc,mz+mbc
        zcell = zlower + (k-0.5d0)*dz
        if (abs(zcell - 0.5d0*dz - center(3)) < 0.5d0*dz) then

          do j=1-mbc,my+mbc
            ycell = ylower + (j-0.5d0)*dy
            if (ycb(1) <= ycell - 0.5d0*dy .and. ycell + 0.5d0*dy <= ycb(2)) then

              do i=1-mbc,mx+mbc
                xcell = xlower + (i-0.5d0)*dx
                if (xcb(1) <= xcell - 0.5d0*dx .and. xcell + 0.5d0*dx <= xcb(2)) then
              
                  do l=1,nsubfaults
                    if (subfaults(l)%longitude*LAT2METER <= xcell .and. &
                        xcell <= subfaults(l)%longitude*LAT2METER + subfaults(l)%width .and. &
                        subfaults(l)%latitude - 0.5d0*subfaults(l)%length <= ycell .and. &
                        ycell <= subfaults(l)%latitude + 0.5d0*subfaults(l)%length .and. &
                        subfaults(l)%rupture_time <= t .and. &
                        t <= subfaults(l)%rupture_time + subfaults(l)%rise_time) then
                      aux(1,i,j,k) = subfaults(l)%slip/subfaults(l)%rise_time
                      exit
                    end if
                  end do

                end if
              end do

            end if
          end do
 
        end if
      end do

    end if

end subroutine b4step3
