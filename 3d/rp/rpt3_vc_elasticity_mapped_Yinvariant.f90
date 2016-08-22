! ==================================================================
subroutine rpt3(ixyz,icoor,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! ==================================================================

!     # Riemann solver in the transverse direction for the elasticity equations
!     # with varying material properties

!     #
!     # On input,

!     #    ql,qr is the data along some one-dimensional slice, as in rpn3
!     #         This slice is
!     #             in the x-direction if ixyz=1,
!     #             in the y-direction if ixyz=2, or
!     #             in the z-direction if ixyz=3.
!     #    asdq is an array of flux differences (A^*\Dq).
!     #         asdq(i,:) is the flux difference propagating away from
!     #         the interface between cells i-1 and i.
!     #    Note that asdq represents B^*\Dq if ixyz=2 or C^*\Dq if ixyz=3.

!     #    ixyz indicates the direction of the original Riemann solve,
!     #         called the x-like direction in the table below:

!     #               x-like direction   y-like direction   z-like direction
!     #      ixyz=1:        x                  y                  z
!     #      ixyz=2:        y                  z                  x
!     #      ixyz=3:        z                  x                  y

!     #    icoor indicates direction in which the transverse solve should
!     #         be performed.
!     #      icoor=2: split in the y-like direction.
!     #      icoor=3: split in the z-like direction.

!     #    For example,
!     #      ixyz=1, icoor=2 means asdq=A^*\Dq, and should be split in y into
!     #                        bmasdq = B^-A^*\Dq,
!     #                        bpasdq = B^+A^*\Dq.
!     #
!     #      ixyz=2, icoor=2 means asdq=B^*\Dq, and should be split in z into
!     #                        bmasdq = C^-B^*\Dq,
!     #                        bpasdq = C^+B^*\Dq.

!     #    The parameter imp is generally needed only if aux
!     #    arrays are being used, in order to access the appropriate
!     #    variable coefficients:

!     #    imp = 1 if asdq = A^- \Dq,  the left-going flux difference
!     #          2 if asdq = A^+ \Dq, the right-going flux difference

!     #    aux2(:,:,2) is a 1d slice of the aux array along the row
!     #                 where the data ql, qr lie.
!     #    aux1(:,:,2) and aux3(:,:,2) are neighboring rows in the
!     #                 y-like direction
!     #    aux2(:,:,1) and aux2(:,:,3) are neighboring rows in the
!     #                z-like direction


    implicit none
    integer, intent(in) :: ixyz, icoor, imp, maxm,meqn,mwaves,mbc,mx, maux
    double precision, intent(in) :: ql(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: qr(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: asdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: aux1(maux,1-mbc:maxm+mbc,3)
    double precision, intent(in) :: aux2(maux,1-mbc:maxm+mbc,3)
    double precision, intent(in) :: aux3(maux,1-mbc:maxm+mbc,3)
    double precision, intent(out) :: bmasdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(out) :: bpasdq(meqn,1-mbc:maxm+mbc)
    integer :: i, iadj, j, sig_xx, sig_yy, sig_zz, sig_xy, sig_xz, sig_yz, u, v, w
    double precision :: wave(meqn,mwaves)
    double precision :: s(mwaves)
    double precision :: dsig_xx, dsig_yy, dsig_zz, dsig_xy, dsig_xz, dsig_yz, du, dv, dw
    double precision :: lama, mua, bulka, cpa, csa, lamb, mub, bulkb, cpb, csb
    double precision :: lam, mu, bulk, cp, cs
    double precision :: det, a1, a2, a3, a4, a5, a6

    ! Variables for the mapping in the xz plane
    double precision :: nxa, nza, nxb, nzb, arearatioa, arearatiob

!   These are just for readability
    sig_xx = 1
    sig_yy = 2
    sig_zz = 3
    sig_xy = 4
    sig_xz = 5
    sig_yz = 6
    u = 7
    v = 8
    w = 9

!     # split the flux difference asdq into 3 downward parts,
!     # one traveling at speed -cp and 2 traveling at speed -cs
!     # relative to the material properties to below the interface,
!     # and 3 upward parts, one traveling at speed cp
!     # and two traveling at speed cs
!     # relative to the material properties above the interface,


    do i=2-mbc,mx+mbc

!        # imp is used to flag whether the original wave is going to left or right.
        iadj = i-2+imp    !#  =  i-1 for amdq,  i for apdq

        dsig_xx = asdq(sig_xx,i)
        dsig_yy = asdq(sig_yy,i)
        dsig_zz = asdq(sig_zz,i)
        dsig_xy = asdq(sig_xy,i)
        dsig_xz = asdq(sig_xz,i)
        dsig_yz = asdq(sig_yz,i)
        du = asdq(u,i)
        dv = asdq(v,i)
        dw = asdq(w,i)

        if (icoor == 2) then
!           # transverse direction is y-like direction so
!           # auxN(:,:,2) holds data in appropriate plane and N=(1,2,3)
!           # for row (below,at,above) the slice of q data

            ! determine x-z plane normal info
            if (ixyz + icoor == 5) then
                ! transverse direction is x
                nxb = aux2(6,iadj,2)
                nzb = aux2(7,iadj,2)
                arearatiob = aux2(8,iadj,2)

                nxa = aux3(6,iadj,2)
                nza = aux3(7,iadj,2)
                arearatioa = aux3(8,iadj,2)
            else if (ixyz + icoor == 3 .or. ixyz + icoor == 6) then
                ! transverse direction is y
                nxb = 0.d0
                nzb = 0.d0
                arearatiob = aux2(9,iadj,2)

                nxa = 0.d0
                nza = 0.d0
                arearatioa = aux3(9,iadj,2)
            else if (ixyz + icoor == 4) then
                ! transverse direction is z
                nxb = aux2(10,iadj,2)
                nzb = aux2(11,iadj,2)
                arearatiob = aux2(12,iadj,2)

                nxa = aux3(10,iadj,2)
                nza = aux3(11,iadj,2)
                arearatioa = aux3(12,iadj,2)
            end if

            ! Assign material parameters
            lamb = aux1(2,iadj,2)
            mub = aux1(3,iadj,2)
            bulkb = lamb + 2.d0*mub
            cpb = aux1(4,iadj,2)
            csb = aux1(5,iadj,2)

            lam = aux2(2,iadj,2)
            mu = aux2(3,iadj,2)
            bulk = lam + 2.d0*mu
            cp = aux2(4,iadj,2)
            cs = aux2(5,iadj,2)

            lama = aux3(2,iadj,2)
            mua = aux3(3,iadj,2)
            bulka = lama + 2.d0*mua
            cpa = aux3(4,iadj,2)
            csa = aux3(5,iadj,2)
        else !! (icoor .eq. 3)
!           # transverse direction is z-like direction so
!           # aux2(:,:,N) holds data in appropriate plane and N=(1,2,3)
!           # for row (below,at,above) the slice of q data

            ! determine x-y plane normal info
            if (ixyz + icoor == 5) then
                ! transverse direction is x
                nxb = aux2(6,iadj,2)
                nzb = aux2(7,iadj,2)
                arearatiob = aux2(8,iadj,2)

                nxa = aux2(6,iadj,3)
                nza = aux2(7,iadj,3)
                arearatioa = aux2(8,iadj,3)
            else if (ixyz + icoor == 3 .or. ixyz + icoor == 6) then
                ! transverse direction is y
                nxb = 0.d0
                nzb = 0.d0
                arearatiob = aux2(9,iadj,2)

                nxa = 0.d0
                nza = 0.d0
                arearatioa = aux2(9,iadj,3)
            else if (ixyz + icoor == 4) then
                ! transverse direction is z
                nxb = aux2(10,iadj,2)
                nzb = aux2(11,iadj,2)
                arearatiob = aux2(12,iadj,2)

                nxa = aux2(10,iadj,3)
                nza = aux2(11,iadj,3)
                arearatioa = aux2(12,iadj,3)
            end if

            ! Assign material parameters
            lamb = aux2(2,iadj,1)
            mub = aux2(3,iadj,1)
            bulkb = lamb + 2.d0*mub
            cpb = aux2(4,iadj,1)
            csb = aux2(5,iadj,1)

            lam = aux2(2,iadj,2)
            mu = aux2(3,iadj,2)
            bulk = lam + 2.d0*mu
            cp = aux2(4,iadj,2)
            cs = aux2(5,iadj,2)

            lama = aux2(2,iadj,3)
            mua = aux2(3,iadj,3)
            bulka = lama + 2.d0*mua
            cpa = aux2(4,iadj,3)
            csa = aux2(5,iadj,3)
        endif

        ! Compute the P-wave parts (a1 downward, a2 upward)
        do j = 1, meqn
            wave(j,1) = 0.d0
            wave(j,2) = 0.d0
        end do
        s(1) = -cpb
        s(2) = cpa

        if (ixyz + icoor == 3 .or. ixyz + icoor == 6) then
            ! transverse direction is y
            a1 = (cp*dsig_yy + bulk*dv) / (bulk*cpb + bulkb*cp)
            a2 = (cp*dsig_yy - bulk*dv) / (bulk*cpa + bulka*cp)

            wave(sig_yy,1) = a1 * bulkb
            wave(sig_xx,1) = a1 * lamb
            wave(sig_zz,1) = a1 * lamb
            wave(v,1) = a1 * cpb

            wave(sig_yy,2) = a2 * bulka
            wave(sig_xx,2) = a2 * lama
            wave(sig_zz,2) = a2 * lama
            wave(v,2) = -a2 * cpa
        else
            ! transverse direction is x or z
            a1 = (cp*(dsig_xx*nxb*nxb + dsig_zz*nzb*nzb + 2.d0*nxb*nzb*dsig_xz) + bulk*(nxb*du + nzb*dw)) / (bulk*cpb + bulkb*cp)
            a2 = (cp*(dsig_xx*nxa*nxa + dsig_zz*nza*nza + 2.d0*nxa*nza*dsig_xz) - bulk*(nxa*du + nza*dw)) / (bulk*cpa + bulka*cp)

            wave(sig_xx,1) = a1 * (lamb + 2.d0*mub*nxb*nxb)
            wave(sig_zz,1) = a1 * (lamb + 2.d0*mub*nzb*nzb)
            wave(sig_yy,1) = a1 * lamb
            wave(sig_xz,1) = a1 * 2.d0*mub*nxb*nzb
            wave(u,1) = a1 * cpb * nxb
            wave(w,1) = a1 * cpb * nzb
            s(1) = -cpb

            wave(sig_xx,2) = a2 * (lama + 2.d0*mua*nxa*nxa)
            wave(sig_zz,2) = a2 * (lama + 2.d0*mua*nza*nza)
            wave(sig_yy,2) = a2 * lama
            wave(sig_xz,2) = a2 * 2.d0*mua*nxa*nza
            wave(u,2) = -a2 * cpa * nxa
            wave(w,2) = -a2 * cpa * nza
        end if

        ! Compute the S-wave parts (a3,a4 downward, a5,a6 upward)
        do j = 1, meqn
            wave(j,3) = 0.d0
            wave(j,4) = 0.d0
            wave(j,5) = 0.d0
            wave(j,6) = 0.d0
        end do
        s(3) = -csb
        s(4) = -csb
        s(5) = csa
        s(6) = csa

        if (ixyz + icoor == 3 .or. ixyz + icoor == 6) then
            ! transverse direction is y

            det = mub*cs + mu*csb
            if (det > 1.e-10) then
                a3 = (cs*dsig_xy + mu*du) / det
                a4 = (cs*dsig_yz + mu*dw) / det

                wave(sig_xy,3) = a3 * mub
                wave(u,3) = a3 * csb

                wave(sig_yz,4) = a4 * mub
                wave(w,4) = a4 * csb
            end if

            det = mua*cs + mu*csa
            if (det > 1.e-10) then
                a5 = (cs*dsig_xy - mu*du) / det
                a6 = (cs*dsig_yz - mu*dw) / det

                wave(sig_xy,5) = a5 * mua
                wave(u,5) = -a5 * csa

                wave(sig_yz,6) = a6 * mua
                wave(w,6) = -a6 * csa
            end if
        else
            ! transverse direction is x or z

            det = mub*cs + mu*csb
            if (det > 1.e-10) then
                a3 = (cs*(dsig_xz*(nxb*nxb - nzb*nzb) + nxb*nzb*(dsig_zz - dsig_xx)) + mu*(nxb*dw - nzb*du)) / det
                a4 = (cs*(dsig_xy*nxb + dsig_yz*nzb) + mu*dv) / det

                wave(sig_xx,3) = -a3 * 2.d0*mub*nxb*nzb
                wave(sig_zz,3) = a3 * 2.d0*mub*nxb*nzb
                wave(sig_xz,3) = a3 * mub*(nxb*nxb - nzb*nzb)
                wave(u,3) = -a3 * csb*nzb
                wave(w,3) = a3 * csb*nxb

                wave(sig_xy,4) = a4 * mub*nxb
                wave(sig_yz,4) = a4 * mub*nzb
                wave(v,4) = a4 * csb
            end if

            det = mua*cs + mu*csa
            if (det > 1.e-10) then
                a5 = (cs*(dsig_xz*(nxa*nxa - nza*nza) + nxa*nza*(dsig_zz - dsig_xx)) - mu*(nxa*dw - nza*du)) / det
                a6 = (cs*(dsig_xy*nxa + dsig_yz*nza) - mua*dv) / det

                wave(sig_xx,5) = -a5 * 2.d0*mua*nxa*nza
                wave(sig_zz,5) = a5 * 2.d0*mua*nxa*nza
                wave(sig_xz,5) = a5 * mua*(nxa*nxa - nza*nza)
                wave(u,5) = a5 * csa*nza
                wave(w,5) = -a5 * csa*nxa

                wave(sig_xy,6) = a6 * mua*nxa
                wave(sig_yz,6) = a6 * mua*nza
                wave(v,6) = -a6 * csa
             end if
         end if

        ! Scale speeds by appropriate area ratios
        s(1) = arearatiob*s(1)
        s(3) = arearatiob*s(3)
        s(4) = arearatiob*s(4)

        s(2) = arearatioa*s(2)
        s(5) = arearatioa*s(5)
        s(6) = arearatioa*s(6)

!        # Compute downward and upward flux difference:
!       Remember that s1,s3,s4 < 0 and s2,s5,s6 > 0
        do j=1,meqn
            bmasdq(j,i) = s(1)*wave(j,1) + s(3)*wave(j,3) + s(4)*wave(j,4)
            bpasdq(j,i) = s(2)*wave(j,2) + s(5)*wave(j,5) + s(6)*wave(j,6)
        end do

    end do

    return
end subroutine rpt3
