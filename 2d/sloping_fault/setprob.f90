!   ==================
    subroutine setprob
!   ==================

    implicit none

    character*12 fname
    integer iunit

    REAL (kind=8) :: ylower_p, yupper_p, yf1, yf2, ycf, xf1, xf2
    common /mapped/  ylower_p, yupper_p, yf1, yf2, ycf, xf1, xf2
!
!
      iunit = 7
      fname = 'setprob.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
      call opendatafile(iunit, fname)
                

!
    read(7,*) ylower_p
    read(7,*) yupper_p
    read(7,*) yf1
    read(7,*) yf2
    read(7,*) xf1
    read(7,*) xf2

    ycf = 1 + (yf1+yf2)/(2*(yupper_p-ylower_p))
    ycf = 0.9d0
    write(6,*) '+++ ycf = ',ycf

    return
    end

