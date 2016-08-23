function [xp, yp, zp] = mapc2p(xc, yc, zc)

    % Parameters from run
    fault_width = 50735;
    fault_length = 50735;
    theta = 0.0;
    xcenter = 25e3;
    ycenter = 0.0;
    zcenter = -19.3e3;    
    
    xcl = xcenter - 0.5*fault_width;
    xcu = xcenter + 0.5*fault_width;
    
    zp1 = zcenter - (xcl-xcenter)*sin(theta);
    zp2 = zcenter - (xcu-xcenter)*sin(theta);
    tol = min(abs([zp1 zp2]));
    
    ls = abs(zc - zcenter);

    ind = xc < xcl;
    ls(ind) = sqrt((xc(ind)-xcl).^2 + (zc(ind)-zcenter).^2);
    ind = xc > xcu;
    ls(ind) = sqrt((xc(ind)-xcu).^2 + (zc(ind)-zcenter).^2);    

    zrot = zc - (xc - xcenter)*sin(theta);
    
    xp = xc;
    yp = yc;
    zp = zc;
    
    ind = ls < tol;
    zp(ind) = (tol-ls(ind))/tol.*zrot(ind) + ls(ind)/tol.*zc(ind);
    
end
