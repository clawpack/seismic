  !=====================================================
  subroutine mapc2p(xc,yc,xp,yp)
  !=====================================================
    ! Maps for sloping fault
    ! on input,  (xc,yc) is a computational grid point
    ! on output, (xp,yp) is corresponding point in physical space

    use fault_module, only: center, theta, xcb

    implicit none
    real (kind=8), intent(in) :: xc,yc
    real (kind=8), intent(out) :: xp,yp

    ! Local variables
    real (kind=8) :: ls, alpha, xrot, yrot

    if (xc < xcb(1)) then
      ls = dsqrt((xc-xcb(1))**2 + (yc-center(2))**2)
    elseif (xc > xcb(2)) then
      ls = dsqrt((xc-xcb(2))**2 + (yc-center(2))**2)
    else
      ls = dabs(yc - center(2))
    end if

    alpha = ls/(-center(2))
    xrot = center(1) + dcos(theta)*(xc-center(1)) + dsin(theta)*(yc-center(2))
    yrot = center(2) - dsin(theta)*(xc-center(1)) + dcos(theta)*(yc-center(2))

    if (alpha < 1.d0) then
      xp = (1.d0-alpha)*xrot + alpha*xc
      yp = (1.d0-alpha)*yrot + alpha*yc
    else
      xp = xc
      yp = yc
    end if

    return
    end
