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
    implicit none


    real (kind=8) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real (kind=8) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real (kind=8) :: amrflags(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    logical     allowflag
    external    allowflag
    real (kind=8) :: DOFLAG, DONTFLAG, xlower, ylower, zlower, dx, dy, dz, tolsp, t
    integer :: mx, my, mz, mbc, meqn, maux, i, j, k

!   # loop over interior points on this grid:
    do k = 1,mz
        do j = 1,my
            do i = 1,mx

                amrflags(i,j,k) = DONTFLAG
                if (allowflag(x,y,z,t,level) .and. dabs(q(1,i,j,k) + q(2,i,j,k) &
                    + q(3,i,j,k))/3.d0 .ge. tolsp) then
                    amrflags(i,j,k) = DOFLAG
                end if
            
            end do
        end do
    end do

    return
end
