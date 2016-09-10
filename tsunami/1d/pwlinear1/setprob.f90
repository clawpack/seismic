
subroutine setprob

   !set gauges and other parameters for 1D geoclaw
   ! for other problem set-up features, copy to your directory and modify

   use gauges_module
   use geoclaw_module
   use setprob_module

   implicit none
   integer :: ndim, iunit
   character*25 fname
   common /comsrc/ ndim

   integer :: i

   !common /comgrid/ xgrid, zgrid, mx_grid

   call set_gauges()
   call set_geo()

   open(unit=58, file='grid.data', status='old',form='formatted')
   read(58,*) mx_grid
   if (mx_grid+1 > mx_grid_max) then
      write(6,*) '*** too many topo values'
      stop
      endif 

   do i=1,mx_grid+1
      read(58,*) xgrid(i),zgrid(i)
      hmax(i) = 0.d0
      enddo


   iunit = 7
   fname = 'setprob.data'
!  # open the unit with new routine from Clawpack 4.4 to skip over
!  # comment lines starting with #:
   call opendatafile(iunit, fname)

end subroutine setprob
