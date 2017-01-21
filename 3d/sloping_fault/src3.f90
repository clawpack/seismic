subroutine src3(meqn,mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation 
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! Advances vertical displacement term
 
    implicit none
    integer, intent(in) :: meqn,mbc,mx,my,mz,maux
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz,t,dt
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

    integer :: i,j,k

    do k=1-mbc,mz+mbc
      do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
        
          q(10,i,j,k) = q(10,i,j,k) + dt*q(9,i,j,k)

        end do
      end do
    end do

end subroutine src3
