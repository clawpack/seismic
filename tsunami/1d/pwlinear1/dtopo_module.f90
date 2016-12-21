
module dtopo_module
    implicit none
    integer :: mx_dtopo,mt_dtopo
    real(kind=8) :: xlow_dtopo,t0_dtopo,dx_dtopo,dt_dtopo,xhi_dtopo,tf_dtopo
    real(kind=8), allocatable :: dtopo(:,:), t_dtopo(:), x_dtopo(:)
    real(kind=8) :: dtopo_xshift
    logical :: dtopo_finalized
    save

contains

subroutine read_dtopo()
    use setprob_module, only: xgrid, mx_grid
    implicit none
    integer :: iunit, i, j, k, dtopo_type, num_dtopofiles
    character*150 :: dtopofname

   iunit = 72
   dtopofname = 'dtopo.data'
!  # open the unit with new routine from Clawpack 4.4 to skip over
!  # comment lines starting with #:
   call opendatafile(iunit, dtopofname)

   read(iunit,*) num_dtopofiles
   if (num_dtopofiles .eq. 0) then
       dtopo_finalized = .true.
       return
     else if (num_dtopofiles .eq. 1) then
       read(iunit,*) dtopofname
     else
       write(6,*) '*** Can support at most 1 dtopo file'
       stop
     endif

    close(iunit)

    dtopo_finalized = .false.

    dtopo_type = 3
    open(unit=iunit, file=dtopofname, status='unknown',form='formatted')

    ! Read in header directly
    read(iunit,*) mx_dtopo
    read(iunit,*) mt_dtopo
    read(iunit,*) xlow_dtopo
    read(iunit,*) t0_dtopo
    read(iunit,*) dx_dtopo
    read(iunit,*) dt_dtopo

    xhi_dtopo = xlow_dtopo + dx_dtopo*(mx_dtopo-1)
    tf_dtopo = t0_dtopo + dt_dtopo*(mt_dtopo-1)

    allocate(dtopo(mt_dtopo,mx_dtopo))
    allocate(x_dtopo(mx_dtopo))
    allocate(t_dtopo(mt_dtopo))

    xlow_dtopo = xlow_dtopo + dtopo_xshift

    do j=1,mx_dtopo
        x_dtopo(j) = xlow_dtopo + (j-1)*dx_dtopo
        enddo

    do k=1,mt_dtopo
        t_dtopo(k) = t0_dtopo + (k-1)*dt_dtopo
        enddo

    
    select case(abs(dtopo_type))
         case(2)
            ! read the data
            do k = 1,mt_dtopo
               do i = 1,mx_dtopo
                  read(iunit,*) dtopo(k,i)
               enddo
            enddo
         case(3)
            do k = 1,mt_dtopo
              read(iunit,*) (dtopo(k,i), i=1,mx_dtopo)
            enddo
      end select

     close(iunit)

end subroutine read_dtopo

end module dtopo_module
