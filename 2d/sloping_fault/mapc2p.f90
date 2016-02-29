  !=====================================================
  subroutine mapc2p(xc,yc,xp,yp)
  !=====================================================
     ! Maps for sloping fault
     ! on input,  (xc,yc) is a computational grid point
     ! on output, (xp,yp) is corresponding point in physical space

     implicit none
     REAL (kind=8) :: xc,yc,xp,yp

     ! local variables
     REAL (kind=8)   :: yf
     
     ! Variables from setprob:
     REAL (kind=8) :: ylower_p, yupper_p, yf1, yf2, ycf, xf1, xf2
        
     common /mapped/  ylower_p, yupper_p, yf1, yf2, ycf, xf1, xf2

          
    xp = xc

    if (xc <= xf1) then
        yf = yf1
      else if (xc <= xf2) then
        yf = yf1+(xc-xf1)*(yf2-yf1)/(xf2-xf1)
      else
        yf = yf2
      endif

    if (yc <= ycf) then
        yp = ylower_p + yc*(yf-ylower_p)/ycf
      else
        yp = yf + (yc-ycf)*(yupper_p-yf)/(1.-ycf)
      endif
     
    return
    end
