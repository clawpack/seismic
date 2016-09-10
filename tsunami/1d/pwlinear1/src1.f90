
subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version integrates manning friction or other friction terms if present


    use geoclaw_module, only: dry_tolerance, grav, DEG2RAD
    use geoclaw_module, only: ifrictiontype => friction_forcing
    use geoclaw_module, only: frictioncoeff => friction_coefficient
    use setprob_module, only: xgrid,zgrid,mx_grid


    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)

    !Locals
    real(kind=8) :: gamma
    integer :: i
    real(kind=8) :: u, rcell
    integer :: ndim


    if (frictioncoeff.eq.0.d0 .or. ifrictiontype.eq.0) go to 800
      ! integrate source term based on Manning formula
    if (ifrictiontype.eq.1) then
      do i=1,mx
         if (q(1,i)<=dry_tolerance) then
            q(2,i) = 0.0
         else
            gamma= dsqrt(q(2,i)**2)*(grav*frictioncoeff**2)/(q(1,i)**(7.0/3.0))
            q(2,i)= q(2,i)/(1.d0 + dt*gamma)
        endif
      enddo
    elseif (ifrictiontype.eq.2) then
      do i=1,mx
         if (q(1,i)<=dry_tolerance) then
            q(2,i) = 0.0
         else
            gamma= q(1,i)*grav*dtan(frictioncoeff*DEG2RAD)
            gamma = max(0.d0, abs(q(2,i)) - dt*abs(gamma))
            q(2,i) = sign(gamma, q(2,i))
        endif
      enddo
    endif

 800 continue

    ! add radial source term if ndim==2:

    ndim = 1
    do i=1,mx
         rcell = 0.5*(xgrid(i) + xgrid(i+1)) - xgrid(1)  ! r = dist from center
         q(1,i) = q(1,i) - dt*(ndim-1)/rcell * q(2,i)
         if (q(1,i) .gt. dry_tolerance) then
             u = q(2,i)/q(1,i)
           else
             u = 0.d0
           endif
         q(2,i) = q(2,i) - dt*(ndim-1)/rcell * q(1,i)*u**2
         enddo

end subroutine src1


