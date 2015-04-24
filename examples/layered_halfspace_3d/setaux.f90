subroutine setaux(mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,maux,aux)

! Called at start of computation before calling qinit, and
! when AMR is used, also called every time a new grid patch is created.
! Use to set auxiliary arrays:
!
!  aux(1,:,:,:) - density
!  aux(2,:,:,:) - lambda (lame parameter)
!  aux(3,:,:,:) - mu (lame parameter)
!  aux(4,:,:,:) - p-wave speed
!  aux(5,:,:,:) - s-wave speed
!
! Note that ghost cell values may need to be set if the aux arrays
! are used by the Riemann solver(s).

    implicit none
    integer, intent(in) :: mbc,mx,my,mz,maux
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz
    real(kind=8), intent(out) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    integer :: i,j,k
    real(kind=8) :: zcell, cp, cs, rho, lambda, mu
    
    REAL (kind=8) :: rho1, lambda1, mu1, rho2, lambda2, mu2, layer_boundary
    common /cparam/  rho1, lambda1, mu1, rho2, lambda2, mu2, layer_boundary

    ! Loop over all cells
    do k=1-mbc,mz + mbc
        zcell = zlower + (k-0.5d0)*dz
        do j=1-mbc,my + mbc
            do i=1-mbc,mx + mbc
          
                if (zcell .ge. layer_boundary) then
                    rho = rho1
                    lambda = lambda1
                    mu = mu1
                else
                    rho = rho2
                    lambda = lambda2
                    mu = mu2
                end if
            
                aux(1,i,j,k) = rho
                aux(2,i,j,k) = lambda
                aux(3,i,j,k) = mu

                ! Calculate p and s wave speeds
                cp = dsqrt((lambda + 2.d0*mu)/rho)
                cs = dsqrt(mu/rho)
                aux(4,i,j,k) = cp
                aux(5,i,j,k) = cs

            end do
        end do
    end do

end subroutine setaux
