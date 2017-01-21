! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
  implicit none
!
!     # Riemann solver in the transverse direction for the elastic equations
!     # with varying material properties in a mapped grid
!
!
!     # Contents of ql and qr:
!     #
!     # q(1,:) = sigma^{11} if ixy=1   or   sigma^{22} if ixy=2
!     # q(2,:) = sigma^{22} if ixy=1   or   sigma^{11} if ixy=2
!     # q(3,:) = sigma^{12} = sigma^{21}
!     # q(4,:) = u          if ixy=1   or   v          if ixy=2
!     # q(5,:) = v          if ixy=1   or   u          if ixy=2
!     #
!     # auxN holds corresponding slices of the aux array:
!     #  N = 1 for row below
!     #      2 for this row
!     #      3 for row above
!     #
!     #  auxN(1,i) = rho
!     #  auxN(2,i) = lambda
!     #  auxN(3,i) = mu
!     #  auxN(4,i) = cp
!     #  auxN(5,i) = cs
!
!
!
!     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.
!
!     # imp=1  means  asdq=amdq,    imp=2 means asdq=apdq
!

  integer, intent(in) :: ixy, imp, maxm, meqn, mwaves, mbc, mx, maux
  double precision, intent(in) :: ql, qr, aux1, aux2, aux3, asdq
  double precision, intent(out) :: bmasdq, bpasdq
  dimension     ql(meqn, 1-mbc:maxm+mbc)
  dimension     qr(meqn, 1-mbc:maxm+mbc)
  dimension   asdq(meqn, 1-mbc:maxm+mbc)
  dimension bmasdq(meqn, 1-mbc:maxm+mbc)
  dimension bpasdq(meqn, 1-mbc:maxm+mbc)
  dimension   aux1(maux, 1-mbc:maxm+mbc)
  dimension   aux2(maux, 1-mbc:maxm+mbc)
  dimension   aux3(maux, 1-mbc:maxm+mbc)

  integer :: i, i1
  double precision :: dsigxx, dsigyy, dsigxy, dux, duy
  double precision :: alamm, amum, bulkm, cpm, csm
  double precision :: alam, amu, bulk, cp, cs
  double precision :: alamp, amup, bulkp, cpp, csp
  double precision :: det, a1, a2, a3, a4

  ! Variables required for mapped grid version
  integer :: map, m
  double precision :: nxm, nym, nx2m, ny2m, nxym, nxp, nyp, nx2p, ny2p, nxyp
  double precision :: dsignm, dsigtm, dunm, dutm
  double precision :: dsignp, dsigtp, dunp, dutp
  double precision :: wave(meqn,mwaves), s(mwaves)
!
!
!
!     # set ku to point to  the component of the system that corresponds
!     # to velocity in the direction of this slice, kv to the orthogonal
!     # velocity.  Similarly ksig11 and ksig22 point to normal stresses.
!     # 3rd component is always shear stress sig12.
!
  map = 9 - 3*(ixy-1)

  do i = 2-mbc, mx+mbc
!
!        # imp is used to flag whether wave is going to left or right,
!        # since material properties are different on the two sides
!
    if (imp.eq.1) then
    ! # asdq = amdq, moving to left
      i1 = i-1
    else
    ! # asdq = apdq, moving to right
      i1 = i
    endif

  ! # jumps in asdq:
    dsigxx = asdq(1,i)
    dsigyy = asdq(2,i)
    dsigxy = asdq(3,i)
    dux = asdq(4,i)
    duy = asdq(5,i)

  ! Define direction of normal to grid edge normals for downgoing fluctuation
    nxm = aux2(map,i1)
    nym = aux2(map+1,i1)
    nx2m = nxm*nxm
  	ny2m = nym*nym
  	nxym = nxm*nym

  !Define direction of normal to grid edge normals for upgoing fluctuation
    nxp = aux3(map,i1)
    nyp = aux3(map+1,i1)
    nx2p = nxp*nxp
    ny2p = nyp*nyp
    nxyp = nxp*nyp

    ! Compute jumps in normal/tangential traction/velocity
    dsignm = dsigxx*nx2m + dsigyy*ny2m + 2.d0*dsigxy*nxym
    dsigtm = (dsigyy - dsigxx)*nxym + dsigxy*(nx2m-ny2m)
    dunm = dux*nxm + duy*nym
    dutm = dux*nym - duy*nxm

    dsignp = dsigxx*nx2p + dsigyy*ny2p + 2.d0*dsigxy*nxyp
    dsigtp = (dsigyy - dsigxx)*nxyp + dsigxy*(nx2p-ny2p)
    dunp = dux*nxp + duy*nyp
    dutp = dux*nyp - duy*nxp

  !
  ! # The flux difference asdq is split into downward moving parts
  ! # traveling at speeds -cp and -cs relative to the medium below and
  ! # upward moving parts traveling
  ! # at speeds +cp and +cs relative to the medium above.
  !
  ! # Note that the sum of these parts does not give all of asdq
  ! # since there is also reflection at the interfaces which decreases
  ! # the flux.
  !
  !
  !
  ! # Material parameters in each row of cells:
    alamm = aux1(2,i1)
    alam  = aux2(2,i1)
    alamp = aux3(2,i1)
    amum  = aux1(3,i1)
    amu   = aux2(3,i1)
    amup  = aux3(3,i1)
    bulkm = alamm + 2.d0*amum
    bulk  = alam  + 2.d0*amu
    bulkp = alamp + 2.d0*amup
    cpm = aux1(4,i1)
    cp = aux2(4,i1)
    cpp = aux3(4,i1)
    csm = aux1(5,i1)
    cs = aux2(5,i1)
    csp = aux3(5,i1)

    ! P-wave strengths:
    a1 = (cp*dsignm + bulk*dunm) / (bulkm*cp + bulk*cpm)
    a2 = (cp*dsignp - bulk*dunp) / (bulkp*cp + bulk*cpp)

    ! S-wave strengths depending on if slip is imposed:

    det = amum*cs + amu*csm
    if (det.eq.0.d0) then
      a3 = 0.d0
    else
      a3 = (cs*dsigtm - amu*dutm) / det
    end if

    det = amup*cs + amu*csp
    if (det.eq.0.d0) then
      a4 = 0.d0
    else
      a4 = (cs*dsigtp + amu*dutp) / det
    end if

    ! Compute the waves.
    wave(:,1) = 0.d0
    wave(1,1) = a1 * (alamm + 2.d0*amum*nx2m)
    wave(2,1) = a1 * (alamm + 2.d0*amum*ny2m)
    wave(3,1) = a1 * (2.d0*amum*nxym)
    wave(4,1) = a1 * cpm * nxm
    wave(5,1) = a1 * cpm * nym
    s(1) = -cpm

    wave(:,2) = 0.d0
    wave(1,2) = a2 * (alamp + 2.d0*amup*nx2p)
    wave(2,2) = a2 * (alamp + 2.d0*amup*ny2p)
    wave(3,2) = a2 * (2.d0*amup*nxyp)
    wave(4,2) = - a2 * cpp * nxp
    wave(5,2) = - a2 * cpp * nyp
    s(2) = cpp

    wave(:,3) = 0.d0
    wave(1,3) = - a3 * (2.d0*nxym*amum)
    wave(2,3) = a3 * (2.d0*nxym*amum)
    wave(3,3) = a3 * amum*(nx2m - ny2m)
    wave(4,3) = - a3 * csm * nym
    wave(5,3) = a3 * csm * nxm
    s(3) = -csm

    wave(:,4) = 0.d0
    wave(1,4) = - a4 * (2.d0*nxyp*amup)
    wave(2,4) = a4 * (2.d0*nxyp*amup)
    wave(3,4) = a4 * amup*(nx2p - ny2p)
    wave(4,4) =  a4 * csp * nyp
    wave(5,4) = -a4 * csp * nxp
    s(4) = csp

    ! # Scale spped by relative length of edge of mapped grid
    s(1) = s(1)*aux2(map+2,i1)
    s(2) = s(2)*aux3(map+2,i1)
    s(3) = s(3)*aux2(map+2,i1)
    s(4) = s(4)*aux3(map+2,i1)
!
!        # The down-going flux difference bmasdq is the product  -c * wave
!        # summed over down-going P-wave and S-wave:
!
	! compute the leftgoing and rightgoing flux differences:
        ! Note cpm_s,csm_s < 0   and   cpp_s, csp_s > 0.
    do m=1,meqn
      bmasdq(m,i) = s(1)*wave(m,1) + s(3)*wave(m,3)
      bpasdq(m,i) = s(2)*wave(m,2) + s(4)*wave(m,4)
    enddo

  end do
!
!
  return
end
