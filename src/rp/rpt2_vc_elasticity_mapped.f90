! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
      implicit double precision (a-h,o-z)
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
      dimension     ql(meqn, 1-mbc:maxm+mbc)
      dimension     qr(meqn, 1-mbc:maxm+mbc)
      dimension   asdq(meqn, 1-mbc:maxm+mbc)
      dimension bmasdq(meqn, 1-mbc:maxm+mbc)
      dimension bpasdq(meqn, 1-mbc:maxm+mbc)
      dimension   aux1(maux, 1-mbc:maxm+mbc)
      dimension   aux2(maux, 1-mbc:maxm+mbc)
      dimension   aux3(maux, 1-mbc:maxm+mbc)

      ! Variables required for mapped grid version
      integer :: map, mw
      double precision :: nxm, nym, nx2m, ny2m, nxym
      double precision :: nxp, nyp, nx2p, ny2p, nxyp
      double precision :: cpm_s, cpp_s, csm_s, csp_s
      double precision :: wave(meqn,mwaves)
!
!
!
!     # set ku to point to  the component of the system that corresponds
!     # to velocity in the direction of this slice, kv to the orthogonal
!     # velocity.  Similarly ksig11 and ksig22 point to normal stresses.
!     # 3rd component is always shear stress sig12.
!
!
      if (ixy.eq.1) then
         ksig11 = 1
         ksig22 = 2
         ku = 4
         kv = 5
	 map = 9
      else
         ksig11 = 2
         ksig22 = 1
         ku = 5
         kv = 4
	 map = 6
      endif
!
!
      do i = 2-mbc, mx+mbc
!
!        # imp is used to flag whether wave is going to left or right,
!        # since material properties are different on the two sides
!
         if (imp.eq.1) then
!            # asdq = amdq, moving to left
             i1 = i-1
         else
!            # asdq = apdq, moving to right
             i1 = i
         endif

	 !Define direction of normal to grid edge normals for downgoing fluctuation
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

!
!        # The flux difference asdq is split into downward moving parts
!        # traveling at speeds -cp and -cs relative to the medium below and
!        # upward moving parts traveling
!        # at speeds +cp and +cs relative to the medium above.
!
!        # Note that the sum of these parts does not give all of asdq
!        # since there is also reflection at the interfaces which decreases
!        # the flux.
!
!        # jumps in asdq:
         dsig11 = asdq(1,i)
         dsig22 = asdq(2,i)
         dsig12 = asdq(3,i)
         du     = asdq(4,i)
         dv     = asdq(5,i)
!
!
!        # Material parameters in each row of cells:
         alamm = aux1(2,i1)
         alam  = aux2(2,i1)
         alamp = aux3(2,i1)
         amum  = aux1(3,i1)
         amu   = aux2(3,i1)
         amup  = aux3(3,i1)
         bulkm = alamm + 2.d0*amum
         bulk  = alam  + 2.d0*amu
         bulkp = alamp + 2.d0*amup

!        # P-wave and S-wave speeds in each row of cells:
         cpm = aux1(4,i1)
         cp  = aux2(4,i1)
         cpp = aux3(4,i1)
         csm = aux1(5,i1)
         cs  = aux2(5,i1)
         csp = aux3(5,i1)

!        # transmitted part of down-going P-wave:
         det = bulkm*cp + bulk*cpm
         if (det .eq. 0.d0) then
            write(6,*) 'det1 = 0 in rpt2'
            stop
            endif
         a1 = (cp*(dsig11*nx2m + dsig22*ny2m + 2*nxym*dsig12) + bulk*(nxm*du + nym*dv)) / det

!        # transmitted part of up-going P-wave:
         det = bulk*cpp + bulkp*cp
         if (det .eq. 0.d0) then
            write(6,*) 'det2 = 0 in rpt2'
            stop
            endif
         a2 = (cp*(dsig11*nx2p + dsig22*ny2p + 2*nxyp*dsig12) - bulk*(nxp*du + nyp*dv)) / det
!
!        # transmitted part of down-going S-wave:
         det = amum*cs + amu*csm
         if (det .eq. 0.d0) then
             a3 = 0.d0
         else
             a3 = (cs*(dsig12*(nx2m - ny2m) + nxym*(dsig22 - dsig11)) + amu*(nxm*dv - nym*du)) / det
         endif

!        # transmitted part of up-going S-wave:
         det = amu*csp + amup*cs
         if (det .eq. 0.d0) then
             a4 = 0.d0
         else
	     a4 = (cs*(dsig12*(nx2p - ny2p) + nxyp*(dsig22 - dsig11)) + amu*(nyp*du - nxp*dv)) / det
         endif

         ! Calculate waves: eigenvector*alphas
         ! Calculate downgoing P-wave transverse fluctuation (cpm speed)
    wave(:,1) = 0.d0
	  wave(1,1) = a1 * (alamm + 2*amum*nx2m)
	  wave(2,1) = a1 * (alamm + 2*amum*ny2m)
	  wave(3,1) = a1 * (2*amum*nxym)
	  wave(4,1) = a1 * cpm * nxm
	  wave(5,1) = a1 * cpm *nym

	  ! Calculate upgoing P-wave transverse fluctuation (cpp speed)
    wave(:,2) = 0.d0
	  wave(1,2) = a2 * (alamp + 2*amup*nx2p)
	  wave(2,2) = a2 * (alamp + 2*amup*ny2p)
	  wave(3,2) = a2 * (2*amup*nxyp)
	  wave(4,2) = - a2 * cpp * nxp
	  wave(5,2) = - a2 * cpp *nyp

	  ! Calculate downgoing S-wave transverse fluctuation (csm speed)
    wave(:,3) = 0.d0
	  wave(1,3) = - a3 * (2*nxym*amum)
	  wave(2,3) = a3 * (2*nxym*amum)
	  wave(3,3) = a3 * amum*(nx2m - ny2m)
	  wave(4,3) = - a3 * csm *nym
	  wave(5,3) = a3 * csm * nxm

	  ! Calculate upgoing S-wave transverse fluctuation (csp speed)
    wave(:,4) = 0.d0
	  wave(1,4) = - a4 * (2*nxyp*amup)
	  wave(2,4) = a4 * (2*nxyp*amup)
	  wave(3,4) = a4 * amup*(nx2p - ny2p)
	  wave(4,4) =  a4 * csp * nyp
	  wave(5,4) = -a4 * csp * nxp

!        # Scale P-wave and S-wave speeds for mapped grid with corresponding edge scaling and sign
         cpm_s = -cpm*aux2(map+2,i1)
         cpp_s = cpp*aux3(map+2,i1)
         csm_s = -csm*aux2(map+2,i1)
         csp_s = csp*aux3(map+2,i1)
!
!        # The down-going flux difference bmasdq is the product  -c * wave
!        # summed over down-going P-wave and S-wave:
!
	! compute the leftgoing and rightgoing flux differences:
        ! Note cpm_s,csm_s < 0   and   cpp_s, csp_s > 0.
        do m=1,meqn
            bmasdq(m,i) = cpm_s*wave(m,1) + csm_s*wave(m,3)
            bpasdq(m,i) = cpp_s*wave(m,2) + csp_s*wave(m,4)
        enddo
!
      enddo
!
      return
      end
