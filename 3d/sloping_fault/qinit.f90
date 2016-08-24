subroutine qinit(meqn,mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,q,maux,aux)

    ! Set initial conditions for the q array.
    ! Interpolate from reloaded solution if available
    ! Otherwise, start with trivial ICs

    implicit none

    integer, intent(in) :: mbc,mx,my,mz,maux,meqn
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    integer :: i,j,k,m, ind

    q(:,:,:,:) = 0.d0

    return
end
