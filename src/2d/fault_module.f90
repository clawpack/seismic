module fault_module

    implicit none
    integer :: nsubfaults
    type subfault
      real(kind=8) :: width, slip, rupture_time, rise_time
    end type subfault
    type(subfault), allocatable :: subfaults(:)
    real(kind=8) :: center(2), theta, xcb(2)

contains

    subroutine load_fault(fname)

        implicit none

        character*12, intent(in) :: fname

        integer :: i
        real(kind=8) :: total_width
        real(kind=8) :: input_line(12)

        call opendatafile(7, fname)

        read(7,*)
        read(7,*) nsubfaults
        read(7,*)

        allocate(subfaults(nsubfaults))

        ! Read in subfaults
        read(7,*) input_line
        theta = input_line(2)/180.0*4.d0*datan(1.d0)
        subfaults(1)%width = input_line(3)
        subfaults(1)%slip = input_line(5)
        subfaults(1)%rupture_time = input_line(11)
        subfaults(1)%rise_time = input_line(12)
        total_width = subfaults(1)%width
        center(1) = 0.5d0*input_line(9)*111.d3
        center(2) = -0.5d0*input_line(4)

        do i=2,nsubfaults-1
          read(7,*) input_line
          subfaults(i)%width = input_line(3)
          subfaults(i)%slip = input_line(5)
          subfaults(i)%rupture_time = input_line(11)
          subfaults(i)%rise_time = input_line(12)
          total_width = total_width + subfaults(i)%width
        end do

        read(7,*) input_line
        subfaults(nsubfaults)%width = input_line(3)
        subfaults(nsubfaults)%slip = input_line(5)
        subfaults(nsubfaults)%rupture_time = input_line(11)
        subfaults(nsubfaults)%rise_time = input_line(12)
        total_width = total_width + subfaults(nsubfaults)%width
        center(1) = center(1) + 0.5*(input_line(9)*111.d3  &
                                      + cos(theta)*subfaults(nsubfaults)%width)
        center(2) = center(2) - 0.5*(input_line(4) &
                                      + sin(theta)*subfaults(nsubfaults)%width)

        xcb(1) = center(1) - 0.5d0*total_width
        xcb(2) = center(1) + 0.5d0*total_width

        close(7)

    end subroutine load_fault

end module fault_module
