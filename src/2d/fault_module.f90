module fault_module

    implicit none
    integer :: nsubfaults
    real(kind=8), allocatable :: widths(:), slips(:), rupture_times(:), rise_times(:)
    real(kind=8) :: center(2), theta, xcb(2), mindepth

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

        allocate(widths(nsubfaults),slips(nsubfaults))
        allocate(rupture_times(nsubfaults),rise_times(nsubfaults))

        ! Read in subfaults
        read(7,*) input_line
        theta = input_line(2)/180.0*4.d0*datan(1.d0)
        widths(1) = input_line(3)
        slips(1) = input_line(5)
        rupture_times(1) = input_line(11)
        rise_times(1) = input_line(12)
        total_width = widths(1)
        center(1) = 0.5d0*input_line(9)*111.d3
        center(2) = -0.5d0*input_line(4)

        do i=2,nsubfaults-1
          read(7,*) input_line
          widths(i) = input_line(3)
          slips(i) = input_line(5)
          rupture_times(i) = input_line(11)
          rise_times(i) = input_line(12)
          total_width = total_width + widths(i)
        end do

        read(7,*) input_line
        widths(nsubfaults) = input_line(3)
        slips(nsubfaults) = input_line(5)
        rupture_times(nsubfaults) = input_line(11)
        rise_times(nsubfaults) = input_line(12)
        total_width = total_width + widths(nsubfaults)
        center(1) = center(1) + 0.5*(input_line(9)*111.d3  &
                                      + cos(theta)*widths(nsubfaults))
        center(2) = center(2) - 0.5*(input_line(4) &
                                      + sin(theta)*widths(nsubfaults))

        xcb(1) = center(1) - 0.5d0*total_width
        xcb(2) = center(1) + 0.5d0*total_width

        mindepth = dmin1(dabs(center(2) - 0.5d0*total_width*dsin(theta)), &
                        dabs(center(2) + 0.5d0*total_width*dsin(theta)))

        close(7)

    end subroutine load_fault

end module fault_module
