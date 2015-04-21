subroutine src3(meqn,mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation 
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: meqn,mbc,mx,my,mz,maux
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz,t,dt
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc, 1-mbc:mz+mbc)
    integer :: i, j, k
    real(kind=8) :: xcell, ycell, zcell, vol_cell, two_pi, current_amp

    real (kind=8) :: src_x, src_y, src_z, amplitude, t_span
    common /csrc/    src_x, src_y, src_z, amplitude, t_span

    ! Assuming non-mapped cartesian grid for this example
    vol_cell = dx*dy*dz
    if (t .le. t_span) then
        two_pi = 8.d0*datan(1.d0)
        current_amp = amplitude*(1.d0 - dcos(two_pi*t/t_span))/2.d0
    else
        current_amp = 0.d0
    end if

    do k=1-mbc,mz+mbc
        zcell = zlower + (k-0.5d0)*dz
        if (zcell - 0.5d0*dz .le. src_z .and. src_z .le. zcell + 0.5d0*dx) then
            do j=1-mbc,my+mbc
                ycell = ylower + (j-0.5d0)*dy
                if (ycell - 0.5d0*dy .le. src_y .and. src_y .le. ycell + 0.5d0*dy) then
                    do i=1-mbc,mx+mbc
                        xcell = xlower + (i-0.5d0)*dx
                        if (xcell - 0.5d0*dx .le. src_x .and. src_x .le. xcell + 0.5d0*dx) then

                            q(7,i,j,k) = q(7,i,j,k) + dt/vol_cell*current_amp
                            q(8,i,j,k) = q(8,i,j,k) + dt/vol_cell*current_amp
                            q(9,i,j,k) = q(9,i,j,k) + dt/vol_cell*current_amp

                        end if
                    end do
                end if
            end do
        end if
    end do

end subroutine src3
