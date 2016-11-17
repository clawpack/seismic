subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version does nothing.

    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer :: i, j
    real(kind=8) :: lambda_plate, mu_plate, rho_plate, lambda_water, mu_water, rho_water, g
    common /material/ lambda_plate, mu_plate, rho_plate, lambda_water, mu_water, rho_water, g

    do j=1-mbc,my+mbc
      do i=1-mbc,mx+mbc
        q(6,i,j) = q(6,i,j) + dt*q(5,i,j)
      end do
    end do

end subroutine src2
