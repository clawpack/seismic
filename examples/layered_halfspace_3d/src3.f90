subroutine src3(meqn,mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation 
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: meqn,mbc,mx,my,mz,maux
    real (kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz,t,dt
    real (kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real (kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc, 1-mbc:mz+mbc)
    integer :: i, j, k, m
    real (kind=8) :: xcell, ycell, zcell, pi, current_amp, delta, src_radius
    real (kind=8) :: weights1d(3),  locmult1d(3), force(3)
    real (kind=8) :: weights(27), xlocmult(27), ylocmult(27), zlocmult(27)

    real (kind=8) :: src_x, src_y, src_z, src_cell_radius, amplitude, t_span
    common /csrc/    src_x, src_y, src_z, src_cell_radius, amplitude, t_span



    if (t .le. t_span) then
        pi = 4.d0*datan(1.d0)
        src_radius = src_cell_radius*max(dx,dy,dz)

        weights1d(1) = 1.d0/6.d0
        weights1d(2) = 4.d0/6.d0
        weights1d(3) = 1.d0/6.d0
        locmult1d(1) = -1.d0
        locmult1d(2) = 0.d0
        locmult1d(3) = 1.d0
        do i=1,3
            do j=1,3
                do k=1,3
                    weights(k + 3*(j-1) + 9*(i-1)) = weights1d(i)*weights1d(j)*weights1d(k)
                    xlocmult(k + 3*(j-1) + 9*(i-1)) = locmult1d(i)
                    ylocmult(k + 3*(j-1) + 9*(i-1)) = locmult1d(j)
                    zlocmult(k + 3*(j-1) + 9*(i-1)) = locmult1d(k)
                end do
            end do
        end do

        current_amp = amplitude*(1.d0 - dcos(2.d0*pi*t/t_span))/2.d0

        do k=1-mbc,mz+mbc
            zcell = zlower + (k-0.5d0)*dz
            if (dabs(src_z - (zcell-0.5d0*dz)) < src_radius .or. &
                dabs(src_z - (zcell+0.5d0*dz)) < src_radius) then

                do j=1-mbc,my+mbc
                    ycell = ylower + (j-0.5d0)*dy
                    if (dabs(src_y - (ycell-0.5d0*dy)) < src_radius .or. &
                        dabs(src_y - (ycell+0.5d0*dy)) < src_radius) then
                        do i=1-mbc,mx+mbc
                            xcell = xlower + (i-0.5d0)*dx
                            if (dabs(src_x - (xcell-0.5d0*dx)) < src_radius .or. &
                                dabs(src_x - (xcell+0.5d0*dx)) < src_radius) then

                                ! compute cell average of source term using simpson rule
                                force(1) = 0.d0
                                force(2) = 0.d0
                                force(3) = 0.d0
                                do m = 1,27
                                    force(1) = force(1) + weights(m)*current_amp* &
                                                (-pi)*dsin(pi*(xcell + 0.5d0*xlocmult(m)*dx - src_x)/src_radius)/&
                                                    (2.d0*src_radius**2.d0)* &
                                                (1.0 + dcos(pi*(ycell + 0.5d0*ylocmult(m)*dy - src_y)/src_radius))/&
                                                    (2.d0*src_radius)* &
                                                (1.0 + dcos(pi*(zcell + 0.5d0*zlocmult(m)*dz - src_z)/src_radius))/&
                                                    (2.d0*src_radius)
                                    force(2) = force(2) + weights(m)*current_amp* &
                                                (1.0 + dcos(pi*(xcell + 0.5d0*xlocmult(m)*dx - src_x)/src_radius))/&
                                                    (2.d0*src_radius)* &
                                                (-pi)*dsin(pi*(ycell + 0.5d0*ylocmult(m)*dy - src_y)/src_radius)/&
                                                    (2.d0*src_radius**2.d0)* &
                                                (1.0 + dcos(pi*(zcell + 0.5d0*zlocmult(m)*dz - src_z)/src_radius))/&
                                                    (2.d0*src_radius)
                                    force(3) = force(3) + weights(m)*current_amp* &
                                                (1.0 + dcos(pi*(xcell + 0.5d0*xlocmult(m)*dx - src_x)/src_radius))/&
                                                    (2.d0*src_radius)* &
                                                (1.0 + dcos(pi*(ycell + 0.5d0*ylocmult(m)*dy - src_y)/src_radius))/&
                                                    (2.d0*src_radius)* &
                                                (-pi)*dsin(pi*(zcell + 0.5d0*zlocmult(m)*dz - src_z)/src_radius)/&
                                                    (2.d0*src_radius**2.d0)
                                end do

                                q(7,i,j,k) = q(7,i,j,k) - dt*force(1)/aux(1,i,j,k)
                                q(8,i,j,k) = q(8,i,j,k) - dt*force(2)/aux(1,i,j,k)
                                q(9,i,j,k) = q(9,i,j,k) - dt*force(3)/aux(1,i,j,k)

                            end if
                        end do
                    end if
                end do
            end if
        end do
    end if

end subroutine src3
