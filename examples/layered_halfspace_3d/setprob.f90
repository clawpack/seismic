subroutine setprob

    implicit none
      
    real (kind=8) :: rho1, lambda1, mu1, rho2, lambda2, mu2, layer_boundary
    common /cparam/  rho1, lambda1, mu1, rho2, lambda2, mu2, layer_boundary

    real (kind=8) :: src_x, src_y, src_z, src_cell_radius, amplitude, t_span
    common /csrc/    src_x, src_y, src_z, src_cell_radius, amplitude, t_span

    ! local variables:
    character*12 :: fname
    integer :: iunit, i, nrst

    iunit = 7
    fname = 'setprob.data'
!   # open the unit with new routine from Clawpack 4.4 to skip over
!   # comment lines starting with #:
      call opendatafile(iunit, fname)

!
    read(7,*) rho1
    read(7,*) lambda1
    read(7,*) mu1
    read(7,*) rho2
    read(7,*) lambda2
    read(7,*) mu2
    read(7,*) layer_boundary
    read(7,*) src_x
    read(7,*) src_y
    read(7,*) src_z
    read(7,*) src_cell_radius
    read(7,*) amplitude
    read(7,*) t_span

    return
end
