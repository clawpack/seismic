
subroutine b4step3(mbc,mx,my,mz,meqn,q,xlower,ylower,zlower, &
    dx,dy,dz,t,dt,maux,aux)

    ! Set slip to zero after t = 1.0
    implicit none
    integer, intent(in) :: mbc,mx,my,mz,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

    if (t > 1.d0) then
        aux(14,:,:,:) = 0.d0
    end if


end subroutine b4step3
