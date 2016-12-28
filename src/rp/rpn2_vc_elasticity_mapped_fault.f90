subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)

  ! Riemann solver for the elasticity equations in 2d, with varying
  ! material properties rho, lambda, and mu, in a mapped grid
  !
  ! This Riemann solver is for mapped grids. It implements the Riemann solver
  ! in the normal direction (it doesn't do a rotation of coordinates)
  !
  ! Note that although there are 5 eigenvectors, one eigenvalue
  ! is always zero and so we only need to compute 4 waves.
  !
  ! solve Riemann problems along one slice of data.
  !
  ! On input, ql contains the state vector at the left edge of each cell
  !           qr contains the state vector at the right edge of each cell
  !
  ! Note that the i'th Riemann problem has left state qr(:,i-1)
  !                                    and right state ql(:,i)
  ! From the basic clawpack routines, this routine is called with ql = qr
  !
  ! This data is along a slice in the x-direction if ixy=1
  !                            or the y-direction if ixy=2.
  !
  ! Contents of ql and qr:
  !
  ! q(1,:) = sigma^{11} if ixy=1   or   sigma^{22} if ixy=2
  ! q(2,:) = sigma^{22} if ixy=1   or   sigma^{11} if ixy=2
  ! q(3,:) = sigma^{12} = sigma^{21}
  ! q(4,:) = u          if ixy=1   or   v          if ixy=2
  ! q(5,:) = v          if ixy=1   or   u          if ixy=2
  !
  ! auxl and auxr hold corresponding slice of the aux array:
  ! Here it is assumed that auxl=auxr gives the cell values
  ! for this slice.
  !
  !  auxl(1,i) = rho, density
  !  auxl(2,i) = lambda
  !  auxl(3,i) = mu
  !  auxl(4,i) = cp, P-wave speed
  !  auxl(5,i) = cs, S-wave speed
  !
  !
  ! On output, wave contains the waves,
  !            s the speeds,
  !            amdq the  left-going flux difference  A^- \Delta q
  !            apdq the right-going flux difference  A^+ \Delta q
  !
  ! Note that the waves are *not* in order of increasing lambda.
  ! Instead the 1- and 2-waves are the P-waves and the 3- and 4-waves
  ! are the S-waves.   (The 5th wave has speed zero and is not used.)

  implicit none

  integer, intent(in) :: ixy, maxm, meqn, mwaves, mbc, mx, maux
  double precision, intent(in) :: ql, qr, auxl, auxr
  double precision, intent(out) :: wave, s, amdq, apdq

  dimension wave( meqn, mwaves, 1-mbc:maxm+mbc)
  dimension    s(mwaves, 1-mbc:maxm+mbc)
  dimension   ql(meqn, 1-mbc:maxm+mbc)
  dimension   qr(meqn, 1-mbc:maxm+mbc)
  dimension apdq(meqn, 1-mbc:maxm+mbc)
  dimension amdq(meqn, 1-mbc:maxm+mbc)
  dimension auxl(maux, 1-mbc:maxm+mbc)
  dimension auxr(maux, 1-mbc:maxm+mbc)

  integer :: i, m
  double precision :: dsigxx, dsigyy, dsigxy, dux, duy
  double precision :: alamr, amur, bulkr, cpr, csr
  double precision :: alaml, amul, bulkl, cpl, csl
  double precision :: det, a1, a2, a3, a4

  ! Variables required for mapped grid version
  integer :: map, mw
  double precision :: nx, ny, nx2, ny2, nxy, dsign, dsigt, dun, dut, slip

  ! set index for normal information
  map = 6 + 3*(ixy-1)

  ! split the jump in q at each interface into waves
  ! The jump is split into leftgoing waves traveling at speeds -cp, -cs
  ! relative to the material properties to the left of the interface,
  ! and rightgoing waves traveling at speeds +cp, +cs
  ! relative to the material properties to the right of the interface,

  do i = 2-mbc, mx+mbc

    ! Define direction of normal to grid edge normals
    nx = auxl(map,i)
    ny = auxl(map+1,i)
    nx2 = nx*nx
    ny2 = ny*ny
    nxy = nx*ny

    ! Compute jumps across cells
    dsigxx = ql(1,i) - qr(1,i-1)
    dsigyy = ql(2,i) - qr(2,i-1)
    dsigxy = ql(3,i) - qr(3,i-1)
    dux = ql(4,i) - qr(4,i-1)
    duy = ql(5,i) - qr(5,i-1)

    ! Compute jumps in normal/tangential traction/velocity
    dsign = dsigxx*nx2 + dsigyy*ny2 + 2.d0*dsigxy*nxy
    dsigt = (dsigyy - dsigxx)*nxy + dsigxy*(nx2-ny2)
    dun = dux*nx + duy*ny
    dut = dux*ny - duy*nx

    ! material properties in cells i (on right) and i-1 (on left):
    alamr = auxl(2,i)
    amur = auxl(3,i)
    bulkr = alamr + 2.d0*amur
    cpr = auxl(4,i)
    csr = auxl(5,i)

    alaml = auxr(2,i-1)
    amul = auxr(3,i-1)
    bulkl = alaml + 2.d0*amul
    cpl = auxr(4,i-1)
    csl = auxr(5,i-1)

    ! P-wave strengths:
    det = bulkl*cpr + bulkr*cpl
    if (det.eq.0.d0) then
      write(6,*) 'det=0 in rpn2'
      stop
    endif
    a1 = (cpr*dsign + bulkr*dun) / det
    a2 = (cpl*dsign - bulkl*dun) / det

    ! S-wave strengths depending on if slip is imposed:
    slip = auxl(13,i)
    if (ixy .eq. 2 .and. dabs(slip) > 1.d-10) then
      a3 = (qr(4,i-1)*ny - qr(5,i-1)*nx - 0.5d0*slip)/csl
      a4 = (ql(4,i)  *ny - ql(5,i)  *nx + 0.5d0*slip)/csr
    else
      det = amul*csr + amur*csl
      if (det.eq.0.d0) then
        ! no s-waves
        a3 = 0.d0
        a4 = 0.d0
      else
        a3 = (csr*dsigt - amur*dut) / det
        a4 = (csl*dsigt + amul*dut) / det
      end if
    end if

    ! 5th wave has velocity 0 so is not computed or propagated.


    ! Compute the waves.
    wave(:,1,i) = 0.d0
    wave(1,1,i) = a1 * (alaml + 2.d0*amul*nx2)
    wave(2,1,i) = a1 * (alaml + 2.d0*amul*ny2)
    wave(3,1,i) = a1 * (2.d0*amul*nxy)
    wave(4,1,i) = a1 * cpl * nx
    wave(5,1,i) = a1 * cpl * ny
    s(1,i) = -cpl

    wave(:,2,i) = 0.d0
    wave(1,2,i) = a2 * (alamr + 2.d0*amur*nx2)
    wave(2,2,i) = a2 * (alamr + 2.d0*amur*ny2)
    wave(3,2,i) = a2 * (2.d0*amur*nxy)
    wave(4,2,i) = - a2 * cpr * nx
    wave(5,2,i) = - a2 * cpr * ny
    s(2,i) = cpr

    wave(:,3,i) = 0.d0
    wave(1,3,i) = - a3 * (2.d0*nxy*amul)
    wave(2,3,i) = a3 * (2.d0*nxy*amul)
    wave(3,3,i) = a3 * amul*(nx2 - ny2)
    wave(4,3,i) = - a3 * csl * ny
    wave(5,3,i) = a3 * csl * nx
    s(3,i) = -csl

    wave(:,4,i) = 0.d0
    wave(1,4,i) = - a4 * (2.d0*nxy*amur)
    wave(2,4,i) = a4 * (2.d0*nxy*amur)
    wave(3,4,i) = a4 * amur*(nx2 - ny2)
    wave(4,4,i) =  a4 * csr * ny
    wave(5,4,i) = -a4 * csr * nx
    s(4,i) = csr

    ! Scales speed by relative length of edge of mapped grid
    do mw=1,mwaves
      s(mw,i) = s(mw,i)*auxl(map+2,i)
    end do


    ! compute the leftgoing and rightgoing flux differences:
    ! Note s(i,1),s(i,3) < 0   and   s(i,2),s(i,4) > 0.
    do m=1,meqn
        amdq(m,i) = s(1,i)*wave(m,1,i) + s(3,i)*wave(m,3,i)
        apdq(m,i) = s(2,i)*wave(m,2,i) + s(4,i)*wave(m,4,i)
    end do
  end do

  return
end subroutine rpn2
