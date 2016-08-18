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
!       1 lambda
!       2 mu
!       3 cp
!       4 cs
!       5 nx at left-wall in x-direction
!       6 ny at left-wall in x-direction
!       7 area ratio of left-wall in x-direction
!       8 nx at left-wall in y-direction
!       9 ny at left-wall in y-direction
!      10 area ratio of left-wall in y-direction
!      11 area ratio of left-wall in z-direction

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
    double precision :: nx, ny, nxy, nx2, ny2, arearatio



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


        ! determine x-y plane normal info
        if (ixyz == 1) then
            nx = auxr(5,i-1)
            ny = auxr(6,i-1)
            arearatio = auxr(7,i-1)
        else if (ixyz == 2) then
            nx = auxr(8,i-1)
            ny = auxr(9,i-1)
            arearatio = auxr(10,i-1)
        else if (ixyz == 3) then
            nx = 0.d0
            ny = 0.d0
            arearatio = auxr(11,i-1)
        end if

        nx2 = nx*nx
        ny2 = ny*ny
        nxy = nx*ny

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

        ! material properties in cells i (on right) and i-1 (on left):
        lamr = auxl(1,i)
        mur = auxl(2,i)
        bulkr = lamr + 2.d0*mur
        cpr = auxl(3,i)
        csr = auxl(4,i)

        laml = auxr(1,i-1)
        mul = auxr(2,i-1)
        bulkl = laml + 2.d0*mul
        cpl = auxr(3,i-1)
        csl = auxr(4,i-1)

        ! Compute the P-waves
        do j = 1, meqn
            wave(j,1,i) = 0.d0
            wave(j,2,i) = 0.d0
        end do
        s(1,i) = -cpl
        s(2,i) = cpr

        det = bulkl*cpr + bulkr*cpl
        if (det < 1.e-10) then
            write(6,*) 'det=0 in rpn3'
            write (6,*) 'cpr', cpr, 'bulkl', bulkl, 'cpl', cpl, 'bulkr', bulkr, 'i', i
            write (6,*) 'test', auxl(1,i)

            stop
        else
            if (ixyz == 3) then
                a1 = (cpr*dsig_zz + bulkr*dw) / det
                a2 = (cpl*dsig_zz - bulkl*dw) / det

                wave(sig_zz,1,i) = a1 * bulkl
                wave(sig_xx,1,i) = a1 * laml
                wave(sig_yy,1,i) = a1 * laml
                wave(w,1,i) = a1 * cpl

                wave(sig_zz,2,i) = a2 * bulkr
                wave(sig_xx,2,i) = a2 * lamr
                wave(sig_yy,2,i) = a2 * lamr
                wave(w,2,i) = -a2 * cpr
            else
                a1 = (cpr*(dsig_xx*nx2 + dsig_yy*ny2 + 2.d0*nxy*dsig_xy) + bulkr*(nx*du + ny*dv)) / det
                a2 = (cpl*(dsig_xx*nx2 + dsig_yy*ny2 + 2.d0*nxy*dsig_xy) - bulkl*(nx*du + ny*dv)) / det

                wave(sig_xx,1,i) = a1 * (laml + 2.d0*mul*nx2)
                wave(sig_yy,1,i) = a1 * (laml + 2.d0*mul*ny2)
                wave(sig_zz,1,i) = a1 * laml
                wave(sig_xy,1,i) = a1 * 2.d0*mul*nxy
                wave(u,1,i) = a1 * cpl * nx
                wave(v,1,i) = a1 * cpl * ny

                wave(sig_xx,2,i) = a2 * (lamr + 2.d0*mur*nx2)
                wave(sig_yy,2,i) = a2 * (lamr + 2.d0*mur*ny2)
                wave(sig_zz,2,i) = a2 * lamr
                wave(sig_xy,2,i) = a2 * 2.d0*mur*nxy
                wave(u,2,i) = -a2 * cpr * nx
                wave(v,2,i) = -a2 * cpr * ny
            end if
        end if


        ! Compute the S-waves
        do j = 1, meqn
            wave(j,3,i) = 0.d0
            wave(j,4,i) = 0.d0
            wave(j,5,i) = 0.d0
            wave(j,6,i) = 0.d0
        end do
        s(3,i) = -csl
        s(4,i) = -csl
        s(5,i) = csr
        s(6,i) = csr

        det = mul*csr + mur*csl
        if (det > 1.e-10) then
            if (ixyz == 3) then
                a3 = (csr*dsig_xz + mur*du) / det
                a4 = (csr*dsig_yz + mur*dv) / det
                a5 = (csl*dsig_xz - mul*du) / det
                a6 = (csl*dsig_yz - mul*dv) / det

                wave(sig_xz,3,i) = a3 * mul
                wave(u,3,i) = a3 * csl

                wave(sig_yz,4,i) = a4 * mul
                wave(v,4,i) = a4 * csl

                wave(sig_xz,5,i) = a5 * mur
                wave(u,5,i) = -a5 * csr

                wave(sig_yz,6,i) = a6 * mur
                wave(v,6,i) = -a6 * csr
            else
                a3 = (csr*(dsig_xy*(nx2 - ny2) + nxy*(dsig_yy - dsig_xx)) + mur*(nx*dv - ny*du)) / det
                a4 = (csr*(dsig_xz*nx + dsig_yz*ny) + mur*dw) / det
                a5 = (csl*(dsig_xy*(nx2 - ny2) + nxy*(dsig_yy - dsig_xx)) - mul*(nx*dv - ny*du)) / det
                a6 = (csl*(dsig_xz*nx + dsig_yz*ny) - mul*dw) / det

                wave(sig_xx,3,i) = -a3 * 2.d0*mul*nxy
                wave(sig_yy,3,i) = a3 * 2.d0*mul*nxy
                wave(sig_xy,3,i) = a3 * mul*(nx2 - ny2)
                wave(u,3,i) = -a3 * csl*ny
                wave(v,3,i) = a3 * csl*nx

                wave(sig_xz,4,i) = a4 * mul*nx
                wave(sig_yz,4,i) = a4 * mul*ny
                wave(w,4,i) = a4 * csl

                wave(sig_xx,5,i) = -a5 * 2.d0*mur*nxy
                wave(sig_yy,5,i) = a5 * 2.d0*mur*nxy
                wave(sig_xy,5,i) = a5 * mur*(nx2 - ny2)
                wave(u,5,i) = a5 * csr*ny
                wave(v,5,i) = -a5 * csr*nx

                wave(sig_xz,6,i) = a6 * mur*nx
                wave(sig_yz,6,i) = a6 * mur*ny
                wave(w,6,i) = -a6 * csr
            end if
        end if

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

