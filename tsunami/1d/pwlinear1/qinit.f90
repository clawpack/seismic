subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.
    ! This default version simply sets eta = max(h + b,0)

    ! For more specific initial conditions
    !  copy this to an application directory and
    !  loop over all grid cells to set values of q(1:meqn, 1:mx).

    !use geoclaw_module, only: dry_tolerance !uncomment if needed
    use geoclaw_module, only: grav  !uncomment if needed
    use setprob_module, only: xgrid,zgrid,mx_grid

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i
    real(kind=8) :: xcell,eta


    do i=1,mx
      xcell = 0.5*(xgrid(i) + xgrid(i+1))
      eta = 0.d0  ! initialize surface to sea level

      ! square pulse for testing:
      if ((xcell > -120e3) .and. (xcell < -80e3)) eta = 2.d0

      q(1,i) = max(0.0, eta - aux(1,i))
      q(2,i) = 0.d0

   enddo


end subroutine qinit
