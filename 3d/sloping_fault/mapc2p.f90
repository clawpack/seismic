subroutine mapc2p(xc, yc, zc, xp, yp, zp)

    use fault_module, only: center, theta, xcb, ycb

    implicit none

    real(kind=8), intent(in) :: xc, yc, zc
    real(kind=8), intent(out) :: xp, yp, zp

    ! Local variables
    real (kind=8) :: ls, alpha, xrot, zrot

    if (xc < xcb(1)) then
      ls = dsqrt((xc-xcb(1))**2 + (zc-center(3))**2)
    elseif (xc > xcb(2)) then
      ls = dsqrt((xc-xcb(2))**2 + (zc-center(3))**2)
    else
      ls = dabs(zc - center(3))
    end if

    alpha = ls/(-center(3))
    xrot = center(1) + dcos(theta)*(xc-center(1)) + dsin(theta)*(zc-center(3))
    zrot = center(3) - dsin(theta)*(xc-center(1)) + dcos(theta)*(zc-center(3))

    yp = yc
    if (alpha < 1.d0) then
      xp = (1.d0-alpha)*xrot + alpha*xc
      zp = (1.d0-alpha)*zrot + alpha*zc
    else
      xp = xc
      zp = zc
    end if

    return
end
