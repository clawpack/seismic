function [xp, yp, zp] = mapc2p(xc, yc, zc)

    global fault_width fault_length theta xcenter ycenter zcenter
    
    xcl = xcenter - 0.5*fault_width;
    xcu = xcenter + 0.5*fault_width;
    ycl = ycenter - 0.5*fault_length;
    ycu = ycenter + 0.5*fault_length;
    
    zp1 = zcenter - (xcl-xcenter)*sin(theta);
    zp2 = zcenter - (xcu-xcenter)*sin(theta);
    tol = norm([zp1 zp2],inf);
    
    ls = abs(zc - zcenter);

    ind = xc < xcl;
    ls(ind) = sqrt((xc(ind)-xcl).^2 + (zc(ind)-zcenter).^2);
    ind = xc > xcu;
    ls(ind) = sqrt((xc(ind)-xcu).^2 + (zc(ind)-zcenter).^2);    
    ind = yc < ycl;
    ls(ind) = sqrt((yc(ind)-ycl).^2 + (zc(ind)-zcenter).^2);
    ind = yc > ycu;
    ls(ind) = sqrt((yc(ind)-ycu).^2 + (zc(ind)-zcenter).^2);
    
    ind = xc < xcl & yc < ycl;
    ls(ind) = sqrt((xc(ind) - xcl).^2 + (yc(ind) - ycl).^2 + (zc(ind) - zcenter).^2);
    ind = xc > xcu & yc < ycl;
    ls(ind) = sqrt((xc(ind) - xcu).^2 + (yc(ind) - ycl).^2 + (zc(ind) - zcenter).^2);
    ind = xc < xcl & yc > ycu;
    ls(ind) = sqrt((xc(ind) - xcl).^2 + (yc(ind) - ycu).^2 + (zc(ind) - zcenter).^2);
    ind = xc > xcu & yc > ycu;
    ls(ind) = sqrt((xc(ind) - xcu).^2 + (yc(ind) - ycu).^2 + (zc(ind) - zcenter).^2);

    zrot = zcenter - (xc - xcenter)*sin(theta);
    
    xp = xc;
    yp = yc;
    zp = zc;
    
    ind = ls < tol;
    zp(ind) = (tol-ls(ind))/tol.*zrot(ind) + ls(ind)/tol.*zc(ind);
    
end
