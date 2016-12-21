module setprob_module

    implicit none
    save

    integer, parameter :: mx_grid_max=50000
    integer :: mx_grid
    real(kind=8), dimension(mx_grid_max) ::  xgrid, zgrid

    real(kind=8), dimension(mx_grid_max) ::  hmax

end module setprob_module
