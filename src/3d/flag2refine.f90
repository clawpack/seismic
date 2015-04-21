! -------------------------------------------------------------------
subroutine flag2refine(mx,my,mz,mbc,meqn,maux,xlower,ylower, &
                        zlower,dx,dy,dz,t, &
                        level,tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)
! -------------------------------------------------------------------

! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! The logical function allowflag(x,y,t,level) is called to check whether
! further refinement at this level is allowed at this particular location
! and time.  The default library version of this routine returns .true.
! for all arguments.  Copy that routine to the application directory and
! modify it if needed to restrict the region where refinement is allowed.
!
! The trace of the stress is computed at points where refinement is
! allowed.  If this value is over tolsp, the point is flagged for
! refinement.  This function assumes the following stress tensor (sigma)
! structure for q:
!
!  q(1) - sigma_xx
!  q(2) - sigma_yy
!  q(3) - sigma_zz
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
    use regions_module, only: regions, num_regions

    implicit none


    real (kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real (kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real (kind=8), intent(inout) :: amrflags(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    logical     allowflag
    external    allowflag
    real (kind=8), intent(in) :: DOFLAG, DONTFLAG, xlower, ylower, zlower, dx, dy, dz, tolsp, t, level

    real (kind=8) :: xcell, ycell, zcell
    integer :: mx, my, mz, mbc, meqn, maux, i, j, k, m, mreg, min_level, max_level

!   # loop over interior points on this grid:
    do k = 1,mz
        zcell = zlower + (k-0.5d0)*dz
        do j = 1,my
            ycell = ylower + (j-0.5d0)*dy
            do i = 1,mx
                xcell = xlower + (i-0.5d0)*dx
                amrflags(i,j,k) = DONTFLAG

!               # check which regions, if any, the point is in
                mreg = -1
                m = 1
                min_level = 0
                max_level = 1000
                do while (mreg .eq. -1 .and. m .le. num_regions)
                    if (regions(m)%t_low .le. t .and. t .le. regions(m)%t_hi .and. &
                        regions(m)%x_low .le. xcell .and. xcell .le. regions(m)%x_hi .and. &
                        regions(m)%y_low .le. ycell .and. ycell .le. regions(m)%y_hi .and. &
                        regions(m)%z_low .le. zcell .and. zcell .le. regions(m)%z_hi) then
                        mreg = m
                        min_level = max(min_level, regions(m)%min_level)
                        max_level = min(max_level, regions(m)%max_level)
                    else
                        m = m + 1
                    end if
                end do

!               # if a region is found, check if the current level is valid
!               # if no region is found, use allowflag and check the specified tolerance
                if (mreg < num_regions .and. min_level .le. level .and. level .le. max_level) then
                    amrflags(i,j,k) = DOFLAG
                else if (mreg .eq. num_regions .and. allowflag(xcell,ycell,zcell,t,level) &
                         .and. dabs(q(1,i,j,k) + q(2,i,j,k) + q(3,i,j,k))/3.d0 .ge. tolsp) then
                    amrflags(i,j,k) = DOFLAG
                end if
            
            end do
        end do
    end do

    return
end
