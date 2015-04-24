subroutine qinit(meqn,mbc,mx,my,mz,xlow,ylow,zlow,dx,dy,dz,q,maux,aux)

! Set (trivial) initial conditions for the q array:
!
!  q(1) - sigma_xx
!  q(2) - sigma_yy
!  q(3) - sigma_zz
!  q(4) - sigma_xy
!  q(5) - sigma_xz
!  q(6) - sigma_yz
!  q(7) - u
!  q(8) - v
!  q(9) - w


    implicit none

    integer, intent(in) :: mbc,mx,my,mz,maux,meqn
    real(kind=8), intent(in) :: xlow,ylow,zlow,dx,dy,dz
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    integer :: i,j,k,m

    do k=1-mbc,mz+mbc
        do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
                do m=1,meqn
                   q(m,i,j,k) = 0.d0
                end do
            end do
        end do
    end do

    return

end
