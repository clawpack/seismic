module fault_module

    implicit none
    real(kind=8), parameter :: LAT2METER = 111133.84012073893 ! from clawpack.geoclaw.data
    integer :: nsubfaults
    type subfault
      real(kind=8) :: width, length, depth, slip, longitude, latitude, rupture_time, rise_time
    end type subfault
    type(subfault), allocatable :: subfaults(:)
    real(kind=8) :: center(3), theta, xcb(2), ycb(2), final_rupture_rise_time

contains

    subroutine load_fault(fname)

        implicit none

        character*12, intent(in) :: fname

        integer :: i
        real(kind=8) :: input_line(12), xp1, yp1, zp1, xp2, yp2, zp2, total_width, total_length

        call opendatafile(7, fname)

        read(7,*)
        read(7,*) nsubfaults
        read(7,*)

        allocate(subfaults(nsubfaults))

        ! Read in subfaults
        final_rupture_rise_time = 0.d0
        do i=1,nsubfaults
          read(7,*) input_line
          theta = input_line(2)/180.0*4.d0*datan(1.d0)
          subfaults(i)%width = input_line(3)
          subfaults(i)%depth = input_line(4)
          subfaults(i)%slip = input_line(5)
          subfaults(i)%length = input_line(8)
          subfaults(i)%longitude = input_line(9)
          subfaults(i)%latitude = input_line(10)
          subfaults(i)%rupture_time = input_line(11)
          subfaults(i)%rise_time = input_line(12)
 
          final_rupture_rise_time = max(final_rupture_rise_time, subfaults(i)%rupture_time + subfaults(i)%rise_time)
        end do
         
        xp1 = subfaults(1)%longitude*LAT2METER
        yp1 = subfaults(1)%latitude - 0.5d0*subfaults(1)%length
        zp1 = -subfaults(1)%depth

        xp2 = subfaults(nsubfaults)%longitude*LAT2METER + dcos(theta)*subfaults(nsubfaults)%width
        yp2 = subfaults(nsubfaults)%latitude + 0.5d0*subfaults(nsubfaults)%length
        zp2 = -subfaults(nsubfaults)%depth - dsin(theta)*subfaults(nsubfaults)%width

        center(1) = 0.5d0*(xp1 + xp2)
        center(2) = 0.5d0*(yp1 + yp2)
        center(3) = 0.5d0*(zp1 + zp2)
        total_width = dsqrt((xp2-xp1)**2 + (zp2-zp1)**2)
        total_length = yp2-yp1

        xcb(1) = center(1) - 0.5d0*total_width
        xcb(2) = center(1) + 0.5d0*total_width
        ycb(1) = center(2) - 0.5d0*total_length
        ycb(2) = center(2) + 0.5d0*total_length

        close(7)

    end subroutine load_fault

end module fault_module
