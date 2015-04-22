subroutine src3(meqn,mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation 
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version does nothing. 
 
    use amr_module, only: lfine, hxposs, hyposs, hzposs

    implicit none
    integer, intent(in) :: meqn,mbc,mx,my,mz,maux
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz,t,dt
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc, 1-mbc:mz+mbc)
    integer :: i, j, k, m, counter
    real(kind=8) :: relpos(3)
    real(kind=8) :: xcell, ycell, zcell, vol_cell, two_pi, current_amp, adj_amp

    real (kind=8) :: src_x, src_y, src_z, amplitude, t_span
    common /csrc/    src_x, src_y, src_z, amplitude, t_span

    real (kind=8) :: tol = 1.0d-10

    ! Only add source to finest grid
    if (t .le. t_span .and. dabs(hxposs(lfine)-dx) < tol .and. dabs(hyposs(lfine)-dy) < tol .and. &
        dabs(hzposs(lfine)-dz) < tol) then
        vol_cell = dx*dy*dz
        two_pi = 8.d0*datan(1.d0)
        current_amp = amplitude*(1.d0 - dcos(two_pi*t/t_span))/2.d0

        do k=1-mbc,mz+mbc
            zcell = zlower + (k-0.5d0)*dz
            relpos(3) = (src_z - (zcell - 0.5d0*dz))/dz
            if (-tol < relpos(3) .and. relpos(3) < 1.d0+tol) then
                do j=1-mbc,my+mbc
                    ycell = ylower + (j-0.5d0)*dy
                    relpos(2) = (src_y - (ycell - 0.5d0*dy))/dy
                    if (-tol < relpos(2) .and. relpos(2) < 1.d0+tol) then
                        do i=1-mbc,mx+mbc
                            xcell = xlower + (i-0.5d0)*dx
                            relpos(1) = (src_x - (xcell - 0.5d0*dx))/dx

                            if (-tol < relpos(1) .and. relpos(1) < 1.d0+tol) then
                                counter = 0
                                do m=1,3
                                    if(dabs(relpos(m)) < tol .or. dabs(relpos(m)-1.d0) < tol) then
                                        counter = counter + 1
                                    end if
                                end do

                                adj_amp = current_amp/2.d0**counter
                                q(1,i,j,k) = q(1,i,j,k) + dt/vol_cell*adj_amp
                                q(2,i,j,k) = q(2,i,j,k) + dt/vol_cell*adj_amp
                                q(3,i,j,k) = q(3,i,j,k) + dt/vol_cell*adj_amp

                            end if
                        end do
                    end if
                end do
            end if
        end do

        write(6,*) dt/vol_cell*current_amp

    end if

end subroutine src3
