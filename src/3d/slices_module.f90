module slices_module
!     # Output the results for a general system of conservation laws
!     # in 3 dimensions as a 2d amrclaw output file, along one
!     # coordinate-aligned slice.


    implicit none
    save

    ! Number of slices and their info (x,y,z,vx,vy,vz)
    integer :: num_slices
    real (kind=8), allocatable, dimension(:) :: slice_x,slice_y,slice_z
    real (kind=8), allocatable, dimension(:) :: slice_vx,slice_vy,slice_vz

contains

    subroutine set_slices(fname)

        implicit none

        ! Input
        character(len=*), intent(in), optional :: fname

        ! Local variables
        integer :: i
        integer :: iunit = 7

        ! Open file
        if (present(fname)) then
            call opendatafile(iunit, fname)
        else
            call opendatafile(iunit, 'slices.data')
        endif

        read(iunit,*) num_slices

        allocate(slice_x(num_slices))
        allocate(slice_y(num_slices))
        allocate(slice_z(num_slices))
        allocate(slice_vx(num_slices))
        allocate(slice_vy(num_slices))
        allocate(slice_vz(num_slices))

        ! Read in slices
        do i=1,num_slices
            read(iunit,*) slice_x(i),slice_y(i),slice_z(i)
            read(iunit,*) slice_vx(i),slice_vy(i),slice_vz(i)
        end do

        close(iunit)

    end subroutine set_slices


! #######################################################################


    subroutine print_slices(lst, lend, time, nvar, naux)

        use amr_module

        integer :: lst, lend, naux, nvar, nnz
        real (kind=8) :: time

        character*15 :: fname1, fname2, fname3, fname4
        real (kind=8) :: val(nvar)
        real (kind=8) :: alpha, xlow, xup, ylow, yup, zlow, zup
        real (kind=8) :: xim, xip, yjm, yjp, zkm, zkp

        integer :: i, j, k, l, ipos, ivar, iaux, idigit, indprint, level, loc, locaux
        integer :: matunit1, matunit2, mitot, mjtot, mktot, mptr, ndim
        integer :: ngrids(num_slices)
        integer :: nx, ny, nz, nstp
        integer iadd
        integer iaddaux

        iadd(ivar,i,j,k)   = loc     +    (ivar-1) &
                                  +    (i-1)*nvar &
                                  +    (j-1)*nvar*mitot &
                                  +    (k-1)*nvar*mitot*mjtot
        iaddaux(iaux,i,j,k) = locaux +    (iaux-1) &
                                  +    (i-1)*naux &
                                  +    (j-1)*naux*mitot &
                                  +    (k-1)*naux*mitot*mjtot


!     ###  make the file names and open output files
        fname1 = 'slice_x.qxxxx'
        fname2 = 'slice_x.txxxx'
!        fname3 = 'fort.axxxx'
!        fname4 = 'fort.bxxxx'
        matunit1 = 50
        matunit2 = 80
 !       matunit3 = 70
 !       matunit4 = 71

        nstp     = matlabu
        do ipos = 13, 10, -1
            idigit = mod(nstp,10)
            fname1(ipos:ipos) = char(ichar('0') + idigit)
            fname2(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
        end do

        do i = 1,num_slices
            fname1(7:7) = char(ichar('0') + i)
            open(unit=matunit1+i,file=fname1,status='unknown',form='formatted')

            ngrids(i) = 0
        end do

        if (output_format == 3) then
!         # binary output
            write(6,*) 'Binary slice output not yet supported'
            stop
!          open(unit=matunit4,file=fname4,status='unknown',
!     &            access='stream')
          endif

        level = lst

        do while (level .le. lend)
            mptr = lstart(level)

            do while (mptr .ne. 0)
                nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
                ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
                nz      = node(ndkhi,mptr) - node(ndklo,mptr) + 1
                loc     = node(store1, mptr)
                locaux  = node(storeaux,mptr)
                mitot   = nx + 2*nghost
                mjtot   = ny + 2*nghost
                mktot   = nz + 2*nghost
                xlow = rnode(cornxlo,mptr)
                ylow = rnode(cornylo,mptr)
                zlow = rnode(cornzlo,mptr)
                xup = xlow + nx*hxposs(level)
                yup = ylow + ny*hyposs(level)
                zup = zlow + nz*hzposs(level)

                do l = 1,num_slices

! ############################# YZ SLICES #######################################
                  if ((dabs(slice_x(l)) > 0.d0) .and. &
                  (slice_x(l) .ge. xlow - 0.01d0*hxposs(level)) .and. &
                  (slice_x(l) .le. xup + 0.01d0*hxposs(level))) then

                          ngrids(l)  = ngrids(l) + 1
                          write(matunit1+l,1301) mptr, level, ny, nz
                      1301         format(i5,'                 grid_number',/, &
                                 i5,'                 AMR_level',/, &
                                 i5,'                 my',/, &
                                 i5,'                 mz')

                          write(matunit1+l,1302) ylow, zlow, hyposs(level), hzposs(level)
                      1302        format(e18.8,'    ylow', /, &
                                 e18.8,'    zlow', /, &
                                 e18.8,'    dy', /, &
                                 e18.8,'    dz')

                          if (output_format == 1) then
                              write(matunit1+l,*) ' '
                              indprint = -1
                              do i = nghost, mitot-nghost
                                  xim = xlow + (i-nghost - 0.5d0)*hxposs(level)
                                  xip = xlow + (i-nghost + 0.5d0)*hxposs(level)
                                  if ((slice_x(l) .ge. xim) .and. (slice_x(l) .lt. xip)) then
                                      if (i==nghost) then
                                          alpha = 1.d0
                                      else if (i==mitot-nghost) then
                                          alpha = 0.d0
                                      else
                                          alpha = (slice_x(l) - xim)/hxposs(level)
                                      end if

                                      indprint = i
                                      do k = nghost+1, mktot-nghost
                                          do j = nghost+1, mjtot-nghost
                                              do ivar=1,nvar
                                                  val(ivar) = alpha*alloc(iadd(ivar,i+1,j,k)) &
                                                       + (1.d0-alpha)*alloc(iadd(ivar,i,j,k))
                                                  if (dabs(val(ivar)) < 1.d-90) then
                                                      val(ivar) = 0.d0
                                                  end if
                                              end do
                                              write(matunit1+l,109) (val(ivar), ivar=1,nvar)
                          109                 format(50e26.16)
                                          end do ! j loop
                                          write(matunit1+l,*) ' '
                                      end do ! k loop
                                  endif
                                  write(matunit1+l,*) ' '
                              end do ! i loop

                              if (indprint < 0) then
                                  write(6,*) '*** ERROR in slices_module'
                                  write(6,*) 'xlow, xup, dx: ',xlow,xup,hxposs(level)
                                  stop
                              end if
                              write(matunit1+l,*) ' '
                          end if
                      end if
! ############################# YZ SLICES #######################################

! ############################# XZ SLICES #######################################
                    if ((dabs(slice_y(l)) > 0.d0) .and. &
                    (slice_y(l) .ge. ylow - 0.01d0*hyposs(level)) .and. &
                    (slice_y(l) .le. yup + 0.01d0*hyposs(level))) then

                        ngrids(l)  = ngrids(l) + 1
                        write(matunit1+l,1303) mptr, level, nx, nz
                1303         format(i5,'                 grid_number',/, &
                               i5,'                 AMR_level',/, &
                               i5,'                 mx',/, &
                               i5,'                 mz')

                        write(matunit1+l,1304) xlow, zlow, hxposs(level), hzposs(level)
                1304        format(e18.8,'    xlow', /, &
                               e18.8,'    zlow', /, &
                               e18.8,'    dx', /, &
                               e18.8,'    dz')

                        if (output_format == 1) then
                            write(matunit1+l,*) ' '
                            indprint = -1
                            do j = nghost, mjtot-nghost
                                yjm = ylow + (j - nghost - 0.5d0)*hyposs(level)
                                yjp = ylow + (j - nghost + 0.5d0)*hyposs(level)
                                if ((slice_y(l) .ge. yjm) .and. (slice_y(l) .lt. yjp)) then
                                    if (j==nghost) then
                                        alpha = 1.d0
                                    else if (j==mjtot-nghost) then
                                        alpha = 0.d0
                                    else
                                        alpha = (slice_y(l) - yjm)/hyposs(level)
                                    end if

                                    indprint = j
                                    do k = nghost+1, mktot-nghost
                                        do i = nghost+1, mitot-nghost
                                            do ivar=1,nvar
                                                val(ivar) = alpha*alloc(iadd(ivar,i,j+1,k)) &
                                                     + (1.d0-alpha)*alloc(iadd(ivar,i,j,k))
                                                if (dabs(val(ivar)) < 1.d-90) then
                                                    val(ivar) = 0.d0
                                                end if
                                            end do
                                            write(matunit1+l,110) (val(ivar), ivar=1,nvar)
                        110                 format(50e26.16)
                                        end do ! i loop
                                        write(matunit1+l,*) ' '
                                    end do ! k loop
                                endif
                                write(matunit1+l,*) ' '
                            end do ! j loop

                            if (indprint < 0) then
                                write(6,*) '*** ERROR in slices_module'
                                write(6,*) 'ylow, yup, dy: ',ylow,yup,hyposs(level)
                                stop
                            end if
                            write(matunit1+l,*) ' '
                        end if
                    end if
! ############################# XZ SLICES #######################################

! ############################# XY SLICES #######################################
                  if ((dabs(slice_vz(l)) > 0.d0) .and. &
                  (slice_z(l) .ge. zlow - 0.01d0*hzposs(level)) .and. &
                  (slice_z(l) .le. zup + 0.01d0*hzposs(level))) then
                        ngrids(l)  = ngrids(l) + 1
                          write(matunit1+l,1305) mptr, level, nx, ny
             1305         format(i5,'                 grid_number',/, &
                                 i5,'                 AMR_level',/, &
                                 i5,'                 mx',/, &
                                 i5,'                 my')

                          write(matunit1+l,1306) xlow, ylow, hxposs(level), hyposs(level)
              1306        format(e18.8,'    xlow', /, &
                                 e18.8,'    ylow', /, &
                                 e18.8,'    dx', /, &
                                 e18.8,'    dy')
                          if (output_format == 1) then
                              write(matunit1+l,*) ' '
                              indprint = -1
                              do k = nghost, mktot-nghost
                                  zkm = zlow + (k-nghost - 0.5d0)*hzposs(level)
                                  zkp = zlow + (k-nghost + 0.5d0)*hzposs(level)
                                  alpha = -1.d0  ! dummy value to indicate row is not in slice
                                  if ((slice_z(l) .ge. zkm) .and. (slice_z(l) .lt. zkp)) then
                                      if (k==nghost) then
                                          alpha = 1.d0
                                      else if (k == mktot-nghost) then
                                          alpha = 0.d0
                                      else
                                          alpha = (slice_z(l) - zkm)/hzposs(level)
                                      end if

                                      indprint = k
                                      do j = nghost+1, mjtot-nghost
                                          do i = nghost+1, mitot-nghost
                                              do ivar=1,nvar
                                                  val(ivar) = alpha*alloc(iadd(ivar,i,j,k+1)) &
                                                       + (1.d0-alpha)*alloc(iadd(ivar,i,j,k))
                                                  if (dabs(val(ivar)) < 1d-90) then
                                                      val(ivar) = 0.d0
                                                  end if
                                              end do
                                              write(matunit1+l,111) (val(ivar), ivar=1,nvar)
                          111                 format(50e26.16)
                                          end do ! i loop
                                          write(matunit1+l,*) ' '
                                      end do ! j loop
                                  endif
                                  write(matunit1+l,*) ' '
                              end do ! k loop

                              if (indprint < 0) then
                                  write(6,*) '*** ERROR in slices_module'
                                  stop
                              end if
                              write(matunit1+l,*) ' '
                          end if
                        end if
! ###############################################################################
                  end do
                  mptr = node(levelptr, mptr)
              end do
              level = level + 1

          end do

      !     --------------
      !     # fort.t file:
      ! --------------



          do i = 1,num_slices
              fname2(7:7) = char(ichar('0') + i)
              open(unit=matunit2+i,file=fname2,status='unknown',form='formatted')
              ndim = 3
          !     # NOTE: we need to print out nghost too in order to strip
          !     #       ghost cells from q when reading in pyclaw.io.binary
              write(matunit2+i,1000) time,nvar,ngrids(i),naux,2,nghost
              1000 format(e18.8,'    time', /, &
                  i5,'                 meqn'/, &
                  i5,'                 ngrids'/, &
                  i5,'                 naux'/, &
                  i5,'                 ndim'/, &
                  i5,'                 nghost'/)

              close(unit=matunit1+i)
              close(unit=matunit2+i)
          end do

        write(6,601) matlabu,time
    601 format('AMRCLAW: Frame ',i4, &
            ' slice files done at time t = ', d12.6,/)

        matlabu = matlabu + 1

    end subroutine print_slices







end module slices_module
