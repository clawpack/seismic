! ==================================================================
subroutine rptt3(ixyz,icoor,imp,impt,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,bsasdq,cmbsasdq,cpbsasdq)
! ==================================================================

!     # Double transverse Riemann solver for elasticity equations
!     # with varying material properties.

!     #
!     # On input,

!     #    ql,qr is the data along some one-dimensional slice, as in rpn3
!     #         This slice is
!     #             in the x-direction if ixyz=1,
!     #             in the y-direction if ixyz=2, or
!     #             in the z-direction if ixyz=3.

!     #    bsasdq is an array of flux differences that result from a
!     #         transverse splitting (a previous call to rpt3).
!     #         This stands for B^* A^* \Dq but could represent any of
!     #         6 possibilities, e.g.  C^* B^* \Dq, as specified by ixyz
!     #         and icoor (see below).
!     #         Moreover, each * represents either + or -, as specified by
!     #         imp and impt.

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
!     #        ixyz=1, icoor=3 means bsasdq=B^*A^*\Dq, and should be
!     #                        split in z into
!     #                           cmbsasdq = C^-B^*A^*\Dq,
!     #                           cpbsasdq = C^+B^*A^*\Dq.
!     #
!     #        ixyz=2, icoor=3 means bsasdq=C^*B^*\Dq, and should be
!     #                        split in x into
!     #                           cmbsasdq = A^-C^*B^*\Dq,
!     #                           cpbsasdq = A^+C^*B^*\Dq.

!     #    The parameters imp and impt are generally needed only if aux
!     #    arrays are being used, in order to access the appropriate
!     #    variable coefficients:

!     #    imp =  1 if bsasdq = B^*A^- \Dq, a left-going flux difference
!     #           2 if bsasdq = B^*A^+ \Dq, a right-going flux difference
!     #    impt = 1 if bsasdq = B^-A^* \Dq, a down-going flux difference
!     #           2 if bsasdq = B^+A^* \Dq, an up-going flux difference

!     #    aux2(:,:,2) is a 1d slice of the aux array along the row
!     #                 where the data ql, qr lie.
!     #    aux1(:,:,2) and aux3(:,:,2) are neighboring rows in the
!     #                 y-like direction
!     #    aux2(:,:,1) and aux2(:,:,3) are neighboring rows in the
!     #                z-like direction


    implicit none
    integer, intent(in) :: ixyz, icoor, imp, impt, maxm,meqn,mwaves,mbc,mx, maux
    double precision, intent(in) :: ql(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: qr(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: bsasdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(in) :: aux1(maux,1-mbc:maxm+mbc,3)
    double precision, intent(in) :: aux2(maux,1-mbc:maxm+mbc,3)
    double precision, intent(in) :: aux3(maux,1-mbc:maxm+mbc,3)
    double precision, intent(out) :: cmbsasdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(out) :: cpbsasdq(meqn,1-mbc:maxm+mbc)
    integer :: i, iadj, j, sig_xx, sig_yy, sig_zz, sig_xy, sig_xz, sig_yz, u, v, w
    double precision :: wave(meqn,mwaves)
    double precision :: s(mwaves)
    double precision :: dsig_xx, dsig_yy, dsig_zz, dsig_xy, dsig_xz, dsig_yz, du, dv, dw
    double precision :: lama, mua, bulka, cpa, csa, lamb, mub, bulkb, cpb, csb
    double precision :: lam, mu, bulk, cp, cs
    double precision :: det, a1, a2, a3, a4, a5, a6

    ! Variables for the mapping in the xy plane
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


!     # split the flux difference bsasdq into 3 downward parts,
!     # one traveling at speed -cp and 2 traveling at speed -cs
!     # relative to the material properties to below the interface,
!     # and 3 upward parts, one traveling at speed cp
!     # and two traveling at speed cs
!     # relative to the material properties above the interface,

    do i=2-mbc,mx+mbc

!        # imp is used to flag whether the original wave is going to left or right.
        iadj = i-2+imp    !#  =  i-1 for bsamdq,  i for bsapdq

        dsig_xx = bsasdq(sig_xx,i)
        dsig_yy = bsasdq(sig_yy,i)
        dsig_zz = bsasdq(sig_zz,i)
        dsig_xy = bsasdq(sig_xy,i)
        dsig_xz = bsasdq(sig_xz,i)
        dsig_yz = bsasdq(sig_yz,i)
        du = bsasdq(u,i)
        dv = bsasdq(v,i)
        dw = bsasdq(w,i)


        if (impt == 1) then
!           # bsasdq propagates to "left" in transverse direction used in rpt3
!           # so we either use auxN(:,:,1) or aux1(:,:,N) for N=1,2,3
!           # depending on icoor:

            if (icoor == 2) then
!           # new double-transverse direction is e2 direction

                ! obtain mapped-grid parameters
                nxb = aux2(11,iadj,1)
                nzb = aux2(13,iadj,1)
                arearatiob = aux2(14,iadj,1)

                nxa = aux3(11,iadj,1)
                nza = aux3(13,iadj,1)
                arearatioa = aux3(14,iadj,1)

                ! Assign material parameters
                lamb = aux1(2,iadj,1)
                mub = aux1(3,iadj,1)
                bulkb = lamb + 2.d0*mub
                cpb = aux1(4,iadj,1)
                csb = aux1(5,iadj,1)

                lam = aux2(2,iadj,1)
                mu = aux2(3,iadj,1)
                bulk = lam + 2.d0*mu
                cp = aux2(4,iadj,1)
                cs = aux2(5,iadj,1)

                lama = aux3(2,iadj,1)
                mua = aux3(3,iadj,1)
                bulka = lama + 2.d0*mua
                cpa = aux3(4,iadj,1)
                csa = aux3(5,iadj,1)
            else !! (icoor .eq. 3)
!           # new double-transverse direction e3 direction

                ! obtain mapped-grid parameters
                nxb = aux1(15,iadj,2)
                nzb = aux1(17,iadj,2)
                arearatiob = aux1(18,iadj,2)

                nxa = aux1(15,iadj,3)
                nza = aux1(17,iadj,3)
                arearatioa = aux1(18,iadj,3)

                ! Assign material parameters
                lamb = aux1(2,iadj,1)
                mub = aux1(3,iadj,1)
                bulkb = lamb + 2.d0*mub
                cpb = aux1(4,iadj,1)
                csb = aux1(5,iadj,1)

                lam = aux1(2,iadj,2)
                mu = aux1(3,iadj,2)
                bulk = lam + 2.d0*mu
                cp = aux1(4,iadj,2)
                cs = aux1(5,iadj,2)

                lama = aux1(2,iadj,3)
                mua = aux1(3,iadj,3)
                bulka = lama + 2.d0*mua
                cpa = aux1(4,iadj,3)
                csa = aux1(5,iadj,3)
            endif
        else
!           # bsasdq propagates to "right" in transverse direction used in rpt3
!           # so we either use auxN(:,:,3) or aux3(:,:,N) for N=1,2,3
!           # depending on icoor:
            if (icoor == 2) then
!           # new double-transverse direction is e2 direction

                ! obtain mapped-grid parameters
                nxb = aux2(11,iadj,3)
                nzb = aux2(13,iadj,3)
                arearatiob = aux2(14,iadj,3)

                nxa = aux3(11,iadj,3)
                nza = aux3(13,iadj,3)
                arearatioa = aux3(14,iadj,3)

                ! Assign material parameters
                lamb = aux1(2,iadj,3)
                mub = aux1(3,iadj,3)
                bulkb = lamb + 2.d0*mub
                cpb = aux1(4,iadj,3)
                csb = aux1(5,iadj,3)

                lam = aux2(2,iadj,3)
                mu = aux2(3,iadj,3)
                bulk = lam + 2.d0*mu
                cp = aux2(4,iadj,3)
                cs = aux2(5,iadj,3)

                lama = aux3(2,iadj,3)
                mua = aux3(3,iadj,3)
                bulka = lama + 2.d0*mua
                cpa = aux3(4,iadj,3)
                csa = aux3(5,iadj,3)
            else !! (icoor .eq. 3)
!           # new double-transverse direction is e3 direction

                ! obtain mapped-grid parameters
                nxb = aux3(15,iadj,2)
                nzb = aux3(17,iadj,2)
                arearatiob = aux3(18,iadj,2)

                nxa = aux3(15,iadj,3)
                nza = aux3(17,iadj,3)
                arearatioa = aux3(18,iadj,3)

                ! Assign material parameters
                lamb = aux3(2,iadj,1)
                mub = aux3(3,iadj,1)
                bulkb = lamb + 2.d0*mub
                cpb = aux3(4,iadj,1)
                csb = aux3(5,iadj,1)

                lam = aux3(2,iadj,2)
                mu = aux3(3,iadj,2)
                bulk = lam + 2.d0*mu
                cp = aux3(4,iadj,2)
                cs = aux3(5,iadj,2)

                lama = aux3(2,iadj,3)
                mua = aux3(3,iadj,3)
                bulka = lama + 2.d0*mua
                cpa = aux3(4,iadj,3)
                csa = aux3(5,iadj,3)
            endif
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
            cmbsasdq(j,i) = s(1)*wave(j,1) + s(3)*wave(j,3) + s(4)*wave(j,4)
            cpbsasdq(j,i) = s(2)*wave(j,2) + s(5)*wave(j,5) + s(6)*wave(j,6)
        end do

    end do


    return
end subroutine rptt3