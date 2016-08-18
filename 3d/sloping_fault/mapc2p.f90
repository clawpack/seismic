subroutine mapc2p(xc, yc, zc, xp, yp, zp)
    implicit none

    real(kind=8), intent(in) :: xc, yc, zc
    real(kind=8), intent(out) :: xp, yp, zp

    ! Variables from setprob:
    real (kind=8) :: center(3), theta, xcb(2), ycb(2), mindepth
    common /fault/  center, theta, xcb, ycb, mindepth

    ! Local variables
    real (kind=8) :: ls, alpha, zrot

    if (xc < xcb(1)) then
      ls = dsqrt((xc-xcb(1))**2 + (zc-center(3))**2)
    elseif (xc > xcb(2)) then
      ls = dsqrt((xc-xcb(2))**2 + (zc-center(3))**2)
    elseif (yc < ycb(1)) then
        ls = dsqrt((yc-ycb(1))**2 + (zc-center(3))**2)
    elseif (yc > ycb(2)) then
      ls = dsqrt((yc-ycb(2))**2 + (zc-center(3))**2)
    elseif (xc < xcb(1) .and. yc < ycb(1)) then
      ls = dsqrt((xc - xcb(1))**2 + (yc-ycb(1))**2 + (zc-center(3))**2)
    elseif (xc > xcb(2) .and. yc < ycb(1)) then
      ls = dsqrt((xc - xcb(2))**2 + (yc-ycb(1))**2 + (zc-center(3))**2)
    elseif (xc < xcb(1) .and. yc > ycb(2)) then
      ls = dsqrt((xc - xcb(1))**2 + (yc-ycb(2))**2 + (zc-center(3))**2)
    elseif (xc > xcb(2) .and. yc > ycb(2)) then
      ls = dsqrt((xc - xcb(2))**2 + (yc-ycb(2))**2 + (zc-center(3))**2)
    else
      ls = dabs(zc - center(3))
    end if

    alpha = ls/mindepth
    zrot = center(3) - (xc-center(1))*dsin(theta)

    xp = xc
    yp = yc
    if (alpha < 1.d0) then
      zp = (1.d0-alpha)*zrot + alpha*zc
    else
      zp = zc
    end if

    return
end
