subroutine b4step1(mbc,mx,meqn,q,xlower,dx,t,dt,maux,aux)

    ! Called before each call to step1.
    ! Use to set time-dependent aux arrays or perform other tasks.

    ! this version checks for negative depths and outputs gauge information

    use gauges_module
    use geoclaw_module, only: dry_tolerance
    use setprob_module, only: hmax, xgrid, mx_grid
    use dtopo_module, only: dtopo, x_dtopo, t_dtopo, mx_dtopo, mt_dtopo, &
                            t0_dtopo, tf_dtopo, dx_dtopo, dt_dtopo

    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc)

    !local variables
    integer :: i,ig,j,m,mvars,k
    real(kind=8) :: dz,alpha,beta,xcell


      do i=1-mbc,mx+mbc
         if (q(1,i)<=dry_tolerance) then
            q(1,i) = max(q(1,i),0.0)
            do m=2,meqn
               q(m,i)=0.d0
            enddo
         endif
      enddo


      ! update hmax, max depth seen at each point:
      do i=1,mx
          hmax(i) = dmax1(hmax(i), q(1,i))
          enddo


      if (allocated(igauge)) then
         mvars = meqn+maux
         if (.not.allocated(sol)) then
            allocate(sol(1:mvars))
         endif

         do ig = 1,mgauges
            if (t.ge.t0gauge(ig).and.t.le.tFgauge(ig)) then
               call return_gauge(meqn,maux,mvars,mx,dx,xlower, &
     &               xgauge(ig),sol(1:mvars),q(1:meqn,1:mx),aux(1:maux,1:mx))

               write(OUTGAUGEUNIT,*) igauge(ig),1,t,(sol(j),j=1,mvars)
            endif
         enddo
      endif

      ! adjust topography if dtopo is active:

      if (allocated(dtopo) .and. (t .lt. tf_dtopo)) then
         k = floor((t-t0_dtopo)/dt_dtopo)
         alpha = (t - t_dtopo(k+1))/dt_dtopo
         do i=1,mx
             xcell = 0.5*(xgrid(i) + xgrid(i+1))
             if ((xcell.le.x_dtopo(1)) .or. (xcell.ge.x_dtopo(mx_dtopo))) then
                dz = 0.d0
               else
                 j = floor((xcell-x_dtopo(1))/mx_dtopo)
                 beta = (xcell - x_dtopo(j+1))/dx_dtopo
                 dz = alpha*(beta*dtopo(k+1,j+1) +(1.d0-beta)*dtopo(k+1,j+2)) &
                      + (1.d0-alpha)*(beta*dtopo(k+2,j+1) +(1.d0-beta)*dtopo(k+2,j+2))
               endif
             aux(1,i) = aux(3,i) + dz  ! adjust initial topo from aux(3,:)
             enddo
         endif

end subroutine b4step1

