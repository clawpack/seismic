subroutine setaux(mbc,mx,xlower,dx,maux,aux)

    ! Called at start of computation before calling qinit, and
    ! when AMR is used, also called every time a new grid patch is created.
    ! Use to set auxiliary arrays aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc).
    ! Note that ghost cell values may need to be set if the aux arrays
    ! are used by the Riemann solver(s).
    !
    ! This version sets aux(1,:) to b(x) for shallow flow
    ! and aux(2,:) to the ratio of cell width to dxc (capacity function).

    !use geoclaw_module, only: dry_tolerance !uncomment if needed
    !use geoclaw_module, only: grav  !uncomment if needed
    use setprob_module, only: xgrid,zgrid,mx_grid

    implicit none
    integer, intent(in) :: mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(out) ::  aux(maux,1-mbc:mx+mbc)

    !locals
    integer :: i,i0,i1,j
    real(kind=8) :: xcell,zcell,a

    if (mx .ne. mx_grid) then
        write(6,*) 'mx_grid from grid.data must agree with mx'
        stop
        endif

    do i=1,mx
        aux(1,i) = 0.5d0*(zgrid(i) + zgrid(i+1))
        aux(2,i) = (xgrid(i+1) - xgrid(i))/dx
        if (aux(2,i) <= 0.d0) then
            write(6,*) '+++ i,xgrid(i),xgrid(i+1): ',i,xgrid(i),xgrid(i+1)
            endif
    enddo

    aux(:,0) = aux(:,1)
    aux(:,-1) = aux(:,1)
    aux(:,mx+1) = aux(:,mx)
    aux(:,mx+2) = aux(:,mx)

end subroutine setaux
