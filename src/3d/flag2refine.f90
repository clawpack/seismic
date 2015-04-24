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
! First, each point is checked against the min_level and max_level
! requirements of any regions present.  If no changes need to be made,
! the infinity norm of the stress tensor is checked against the user
! specified tolsp value.  This function assumes the first 6 components of
! q are the 6 stress tensor components.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
    use regions_module, only: regions, num_regions

    implicit none

    integer, intent(in) :: mx, my, mz, mbc, meqn, maux
    real (kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real (kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real (kind=8), intent(inout) :: amrflags(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    logical     allowflag
    external    allowflag
    real (kind=8), intent(in) :: DOFLAG, DONTFLAG, xlower, ylower, zlower, dx, dy, dz, tolsp, t, level

    real (kind=8) :: xcell, ycell, zcell, max_stress
    integer :: i, j, k, m, min_level, max_level
    integer :: infinity = 1e3

!   # loop over interior points on this grid:
    do k = 1,mz
        zcell = zlower + (k-0.5d0)*dz
        do j = 1,my
            ycell = ylower + (j-0.5d0)*dy
            do i = 1,mx
                xcell = xlower + (i-0.5d0)*dx

!               # obtain the overall min and max levels from any regions containing the point
                min_level = 0
                max_level = infinity
                do m =1,num_regions
                    if (regions(m)%t_low .le. t .and. t .le. regions(m)%t_hi .and. &
                        regions(m)%x_low .le. xcell .and. xcell .le. regions(m)%x_hi .and. &
                        regions(m)%y_low .le. ycell .and. ycell .le. regions(m)%y_hi .and. &
                        regions(m)%z_low .le. zcell .and. zcell .le. regions(m)%z_hi) then
                        min_level = max(min_level, regions(m)%min_level)
                        max_level = min(max_level, regions(m)%max_level)
                    end if
                end do

!               # if point is in region, make sure that region is refined as specified
!               # if nothing needs to be changed, use specified tolerance and stress
                if (min_level > 0 .and. level < min_level) then
                    amrflags(i,j,k) = DOFLAG
                else if (min_level > 0 .and. max_level < level) then
                    amrflags(i,j,k) = DONTFLAG
                else if (allowflag(xcell,ycell,zcell,t,level)) then
                    max_stress = 0.d0
                    do m = 1,6
                        max_stress = max(max_stress, dabs(q(m,i,j,k)))  
                    end do
                    if (max_stress .ge. tolsp) then
                        amrflags(i,j,k) = DOFLAG
                    else
                        amrflags(i,j,k) = DONTFLAG
                    end if
                else
                    amrflags(i,j,k) = DONTFLAG
                end if
            
            end do
        end do
    end do

    return
end
