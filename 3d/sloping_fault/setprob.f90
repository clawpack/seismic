!==================
subroutine setprob
!==================

    implicit none

    character*12 fname
    integer iunit

    real (kind=8) :: center(3), theta, xcb(2), ycb(2), mindepth
    common /fault/  center, theta, xcb, ycb, mindepth

    real (kind=8) :: width, length
!
!
    iunit = 7
    fname = 'setprob.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
    call opendatafile(iunit, fname)


!
    read(7,*) width ! domain_depth
    read(7,*) width ! domain width
    read(7,*) width ! domain_length
    read(7,*) center(1)
    read(7,*) center(2)
    read(7,*) width
    read(7,*) length
    read(7,*) theta
    read(7,*) center(3)

    center(3) = -center(3)
    xcb(1) = center(1) - 0.5d0*width
    xcb(2) = center(1) + 0.5d0*width
    ycb(1) = center(2) - 0.5d0*length
    ycb(2) = center(2) + 0.5d0*length

    mindepth = dmin1(dabs(center(3) - 0.5d0*width*dsin(theta)), &
                      dabs(center(3) + 0.5d0*width*dsin(theta)))

    return
end
