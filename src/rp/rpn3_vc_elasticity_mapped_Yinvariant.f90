! ==================================================================
subroutine rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! ==================================================================

! Riemann solver for the elasticity equations in 3d, with varying
! material properties.
!
! waves: 6
! equations: 9
! aux fields: 6
!
! Conserved quantities:
!       1 sigma_xx
!       2 sigma_yy
!       3 sigma_zz
!       4 sigma_xy
!       5 sigma_xz
!       6 sigma_yz
!       7 u
!       8 v
!       9 w
!
! Auxiliary variables:
!       1 rho
!       2 lambda
!       3 mu
!       4 cp
!       5 cs
!       6 --- not used ---
!       7 nx at lower wall in e1 direction (see setaux1d)
!       8 ny at lower wall in e1 direction
!       9 nz at lower wall in e1 direction
!       10 area ratio of lower wall in e1 direction
!       11 nx at lower wall in e2 direction
!       12 ny at lower wall in e2 direction
!       13 nz at lower wall in e2 direction
!       14 area ratio of lower wall in e2 direction
!       15 nx at lower wall in e3 direction
!       16 ny at lower wall in e3 direction
!       17 nz at lower wall in e3 direction
!       18 area ratio of lower wall in e3 direction

! Note that although there are 9 eigenvectors, 3 eigenvalues are
! always zero and so we only need to compute 6 waves.

! Solve Riemann problems along one slice of data.
! This data is along a slice in the x-direction if ixyz=1
!                               the y-direction if ixyz=2.
!                               the z-direction if ixyz=3.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr

! Note waves 1-2 are the P-waves, waves 3-6 are the S-waves

    implicit none
    integer, intent(in) :: ixyz, maxm,meqn,mwaves,mbc,mx, maux
    double precision, intent(in) ::   ql(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) ::   qr(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: auxl(maux,1-mbc:maxm+mbc)
    double precision, intent(in) :: auxr(maux,1-mbc:maxm+mbc)
    double precision, intent(out) :: wave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision, intent(out) ::    s(mwaves,1-mbc:maxm+mbc)
    double precision, intent(out) :: amdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(out) :: apdq(meqn,1-mbc:maxm+mbc)
    integer :: i, j, sig_xx, sig_yy, sig_zz, sig_xy, sig_xz, sig_yz, u, v, w
    double precision :: dsig_xx, dsig_yy, dsig_zz, dsig_xy, dsig_xz, dsig_yz, du, dv, dw
    double precision :: laml, mul, bulkl, cpl, csl, lamr, mur, bulkr, cpr, csr
    double precision :: det, a1, a2, a3, a4, a5, a6

    ! Variables for the mapping in the xy plane
    double precision :: nx, ny, nz, nxz, arearatio
    real(kind=8) :: tx, ty, tz, ttx, tty, ttz
    real(kind=8) :: dsig_n, dsig_t, dsig_tt, du_n, du_t, du_tt



!       These are just for readability
        sig_xx = 1
        sig_yy = 2
        sig_zz = 3
        sig_xy = 4
        sig_xz = 5
        sig_yz = 6
        u = 7
        v = 8
        w = 9

!     # split the jump in q at each interface into waves
!     # The jump is split into 1 leftgoing wave traveling at speed -cp
!     # and 2 leftgoing waves traveling at speed -cs
!     # relative to the material properties to the left of the interface,
!     # and 1 rightgoing wave traveling at speed cp
!     # and 2 rightgoing waves traveling at speed cs
!     # relative to the material properties to the right of the interface,



!     # find a1-a6, the coefficients of the 6 eigenvectors:
    do i = 2-mbc, mx+mbc

        ! Compute Delta Q values
        dsig_xx = ql(sig_xx,i) - qr(sig_xx,i-1)
        dsig_yy = ql(sig_yy,i) - qr(sig_yy,i-1)
        dsig_zz = ql(sig_zz,i) - qr(sig_zz,i-1)
        dsig_xy = ql(sig_xy,i) - qr(sig_xy,i-1)
        dsig_xz = ql(sig_xz,i) - qr(sig_xz,i-1)
        dsig_yz = ql(sig_yz,i) - qr(sig_yz,i-1)
        du = ql(u,i) - qr(u,i-1)
        dv = ql(v,i) - qr(v,i-1)
        dw = ql(w,i) - qr(w,i-1)

        nx = auxl(7,i)
        ny = auxl(8,i)
        nz = auxl(9,i)
        arearatio = auxl(10,i)
        
        ! ***the remainder only works for y-invariant mapping ***
        if (ixyz .eq. 2) then
            tx = 0.d0 
            ty = 0.d0
            tz = 1.d0
            ttx = 1.d0
            tty = 0.d0
            ttz = 0.d0
        else
            tx = 0.d0
            ty = 1.d0
            tz = 0.d0
            ttx = -nz
            tty = 0.d0
            ttz = nx
        end if
            
        ! Compute normal/tangent jumps in stress/velocity
        dsig_n = (dsig_xx*nx + dsig_xy*ny + dsig_xz*nz)*nx &
                +(dsig_xy*nx + dsig_yy*ny + dsig_yz*nz)*ny &
                +(dsig_xz*nx + dsig_yz*ny + dsig_zz*nz)*nz 
        du_n = du*nx + dv*ny + dw*nz
        
        dsig_t = (dsig_xx*nx + dsig_xy*ny + dsig_xz*nz)*tx &
                +(dsig_xy*nx + dsig_yy*ny + dsig_yz*nz)*ty &
                +(dsig_xz*nx + dsig_yz*ny + dsig_zz*nz)*tz 
        du_t = du*tx + dv*ty + dw*tz


        dsig_tt = (dsig_xx*nx + dsig_xy*ny + dsig_xz*nz)*ttx & 
                 +(dsig_xy*nx + dsig_yy*ny + dsig_yz*nz)*tty &
                 +(dsig_xz*nx + dsig_yz*ny + dsig_zz*nz)*ttz 
        du_tt = du*ttx + dv*tty + dw*ttz

        ! material properties in cells i (on right) and i-1 (on left):
        lamr = auxl(2,i)
        mur = auxl(3,i)
        bulkr = lamr + 2.d0*mur
        cpr = auxl(4,i)
        csr = auxl(5,i)

        laml = auxr(2,i-1)
        mul = auxr(3,i-1)
        bulkl = laml + 2.d0*mul
        cpl = auxr(4,i-1)
        csl = auxr(5,i-1)

        ! Compute the P-wave stengths
        det = bulkl*cpr + bulkr*cpl
        if (det < 1.e-10) then
            write(6,*) 'det=0 in rpn3'
            stop
        end if
        a1 = (cpr*dsig_n + bulkr*du_n) / det
        a2 = (cpl*dsig_n - bulkl*du_n) / det
        
        ! Compute the S-wave strengths
        det = mul*csr + mur*csl
        if (det .eq. 0.d0) then
            a3 = 0.d0
            a4 = 0.d0
            a5 = 0.d0
            a6 = 0.d0
        else
            a3 = (csr*dsig_t + mur*du_t) / det 
            a4 = (csr*dsig_tt + mur*du_tt) / det
            a5 = (csl*dsig_t - mul*du_t) / det 
            a6 = (csl*dsig_tt - mul*du_tt) / det
        end if         

        wave(:,1,i) = 0.d0
        wave(sig_xx,1,i) = a1 * (laml + 2.d0*mul*nx*nx)
        wave(sig_yy,1,i) = a1 * (laml + 2.d0*mul*ny*ny)
        wave(sig_zz,1,i) = a1 * (laml + 2.d0*mul*nz*nz)
        wave(sig_xz,1,i) = a1 * 2.d0*mul*nx*nz
        wave(sig_yz,1,i) = a1 * 2.d0*mul*ny*nz
        wave(u,1,i) = a1 * cpl*nx
        wave(v,1,i) = a1 * cpl*ny
        wave(w,1,i) = a1 * cpl*nz
        s(1,i) = -cpl

        wave(:,2,i) = 0.d0
        wave(sig_xx,2,i) = a2 * (lamr + 2.d0*mur*nx*nx)
        wave(sig_yy,2,i) = a2 * (lamr + 2.d0*mur*ny*ny)
        wave(sig_zz,2,i) = a2 * (lamr + 2.d0*mur*nz*nz)
        wave(sig_xz,2,i) = a2 * 2.d0*mur*nx*nz
        wave(u,2,i) = -a2 * cpr*nx
        wave(v,2,i) = -a2 * cpr*ny
        wave(w,2,i) = -a2 * cpr*nz
        s(2,i) = cpr

        wave(:,3,i) = 0.d0
        wave(sig_xy,3,i) = a3 * mul*(nx*ty + ny*tx)
        wave(sig_yz,3,i) = a3 * mul*(nz*ty + ny*tz)
        wave(sig_xz,3,i) = a3 * mul*(nx*tz + nz*tx)
        wave(u,3,i) = a3 * csl*tx
        wave(v,3,i) = a3 * csl*ty
        wave(w,3,i) = a3 * csl*tz
        s(3,i) = -csl

        wave(:,4,i) = 0.d0
        wave(sig_xx,4,i) = a4 * 2.d0*mul*nx*ttx
        wave(sig_zz,4,i) = a4 * 2.d0*mul*nz*ttz
        wave(sig_xy,4,i) = a4 * mul*(nx*tty + ny*ttx) 
        wave(sig_yz,4,i) = a4 * mul*(ny*ttz + nz*tty)
        wave(sig_xz,4,i) = a4 * mul*(nx*ttz + nz*ttx)
        wave(u,4,i) = a4 * csl*ttx
        wave(v,4,i) = a4 * csl*tty
        wave(w,4,i) = a4 * csl*ttz
        s(4,i) = -csl

        wave(:,5,i) = 0.d0
        wave(sig_xy,5,i) = a5 * mur*(nx*ty + ny*tx)
        wave(sig_yz,5,i) = a5 * mur*(nz*ty + ny*tz)
        wave(sig_xz,5,i) = a5 * mur*(nx*tz + nz*tx)
        wave(u,5,i) = -a5 * csr*tx
        wave(v,5,i) = -a5 * csr*ty
        wave(w,5,i) = -a5 * csr*tz
        s(5,i) = csr

        wave(:,6,i) = 0.d0
        wave(sig_xx,6,i) = a6 * 2.d0*mur*nx*ttx
        wave(sig_zz,6,i) = a6 * 2.d0*mur*nz*ttz
        wave(sig_xy,6,i) = a6 * mur*(nx*tty + ny*ttx) 
        wave(sig_yz,6,i) = a6 * mur*(ny*ttz + nz*tty)
        wave(sig_xz,6,i) = a6 * mur*(nx*ttz + nz*ttx)
        wave(u,6,i) = -a6 * csr*ttx
        wave(v,6,i) = -a6 * csr*tty
        wave(w,6,i) = -a6 * csr*ttz
        s(6,i) = csr

        ! Scale speeds by area ratio
        do j=1,mwaves
            s(j,i) = s(j,i)*arearatio
        end do

        ! compute the leftgoing and rightgoing flux differences:
        ! Note s1,s3,s4 < 0   and   s2,s5,s6 > 0.
        do j=1,meqn
            amdq(j,i) = s(1,i)*wave(j,1,i) + s(3,i)*wave(j,3,i) + s(4,i)*wave(j,4,i)
            apdq(j,i) = s(2,i)*wave(j,2,i) + s(5,i)*wave(j,5,i) + s(6,i)*wave(j,6,i)
        end do
    end do

    return
end subroutine rpn3
