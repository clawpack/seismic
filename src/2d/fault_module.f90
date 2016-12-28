module fault_module

    implicit none
    real(kind=8), parameter :: LAT2METER = 111133.84012073893 ! from clawpack.geoclaw.data
    integer :: nsubfaults
    type subfault
      real(kind=8) :: width, depth, slip, longitude, rupture_time, rise_time
    end type subfault
    type(subfault), allocatable :: subfaults(:)
    real(kind=8) :: center(2), theta, xcb(2)

contains

    subroutine load_fault(fname)

        implicit none

        character*12, intent(in) :: fname

        integer :: i
        real(kind=8) :: input_line(12)

        call opendatafile(7, fname)

        read(7,*)
        read(7,*) nsubfaults
        read(7,*)

        allocate(subfaults(nsubfaults))

        ! Read in subfaults
        do i=1,nsubfaults
          read(7,*) input_line
          theta = input_line(2)/180.0*4.d0*datan(1.d0)
          subfaults(i)%width = input_line(3)
          subfaults(i)%depth = input_line(4)
          subfaults(i)%slip = input_line(5)
          subfaults(i)%longitude = input_line(9)
          subfaults(i)%rupture_time = input_line(11)
          subfaults(i)%rise_time = input_line(12)
        end do
  
        xcb(1) = subfaults(1)%longitude*LAT2METER
        xcb(2) = subfaults(nsubfaults)%longitude*LAT2METER + dcos(theta)*subfaults(nsubfaults)%width
        center(1) = 0.5d0*(xcb(1) + xcb(2))
        center(2) = -0.5d0*(subfaults(1)%depth + subfaults(nsubfaults)%depth + &
                    dsin(theta)*subfaults(nsubfaults)%width)

        close(7)

    end subroutine load_fault

end module fault_module
