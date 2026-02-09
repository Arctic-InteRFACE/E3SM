! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_sediments
!
!> \brief submodules for ocean sediments
!> \author Nicole Jeffery
!> \date   06/09/25
!> \details
!>  This module is based on the subroutine sed_settling_tile the CROCO model
!>  The code and tools are open source. The parts of the code inherited from 
!>  the ROMS model are made available under the conditions of the MIT/X license, 
!>  those inherited from the AGRIF library under the terms of the CeCILL-C license
!>
!=======================================================================
!

module ocn_sediment_settling

   ! module includes
   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_log
   use mpas_constants

   use ocn_config
   use ocn_mesh
   use ocn_constants

   implicit none
   private
   save

   ! private module variables
   real (kind=RKIND), parameter ::  &
      small = 1.e-12_RKIND

   ! bottom erosion parameters
   real (kind=RKIND), parameter ::  &
      phi  = 0.7_RKIND,       & ! porosity of the upper seafloor sediments
      fraction_i = 1.0_RKIND, & ! Fraction of sediment class i in seafloor
      particle_density = 2650.0_RKIND, & ! density of particles in seafloor
      gamma = 1.0e-8_RKIND,   & ! 1.0e-6_RKIND,   & ! ranges from 10^-3 * (10^-3 to 10^-5)
      !this is way too small --->
      !tau_c = 6.4e-7_RKIND,   & ! critical bottom stress N/m2 or 6.4e-6. Should be function of bed type
      tau_c = 0.2_RKIND,      & ! critical bottom stress N/m2
      ERATE = 4.0e-5_RKIND      ! CROCO, with the erosion constant independent of 
                                ! gamma, particle density and sinking rate 
      
   ! public method interfaces
   public :: ocn_sediment_settling_flx, &
             ocn_sediment_init, &
             ocn_compute_bottom_shear_stress

!***********************************************************************

   contains

!***********************************************************************
!
!  routine ocn_sediment_settling
!
!> \brief Computes the settling tendency for sediment tracers
!> \author Nicole Jeffery
!> \date   06/09/2025
!> \details
!
!>  The routine computes the vertical settling (sinking) of suspended
!>  sediment via a semi-Lagrangian advective flux algorithm. It uses a
!>  parabolic,  vertical reconstructuion of the suspended  sediment in
!>  the water column with PPT/WENO constraints to avoid oscillations.
!>
!>  References:
!>
!>  Colella, P. and P. Woodward, 1984: The piecewise parabolic method
!>    (PPM) for gas-dynamical simulations, J. Comp. Phys., 54, 174-201.
!>
!>  Liu, X.D., S. Osher, and T. Chan, 1994: Weighted essentially
!>    nonoscillatory shemes, J. Comp. Phys., 115, 200-212.
!>
!>  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.
!>    Arango, 2008:  Development of a three-dimensional,  regional,
!>    coupled wave, current, and sediment-transport model, Computers
!>    & Geosciences, 34, 1284-1306.
!-----------------------------------------------------------------------

      subroutine ocn_sediment_settling_flx (sedimentConc, sinkingRate, &
                                            workTend, layerThickness,  &
                                            dt, kmin, kmax,     &
                                            sedimentBottomFlux, &
                                            netSedimentStorage, &
                                            bottomStressCell,   &
                                            netErosionFlux)

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------
                                            
      real (kind=RKIND), dimension(:), intent(inout) :: &
          workTend              ! settling tendency

      real (kind=RKIND), intent(inout) :: &
          sedimentBottomFlux, & ! kg/m2/s
          netErosionFlux,     & ! kg/m2/s
          netSedimentStorage    ! kg/m2
      
      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------
 
      real (kind=RKIND), dimension(:), intent(in) :: &
          sedimentConc,      &
          layerThickness
          
      real (kind=RKIND), intent(in) :: &
          sinkingRate, &   !m/s
          dt,    &
          bottomStressCell ! kg/m/s2

      integer, intent(in):: &
          kmin, kmax
      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer ::          &
         k, m, ks, err
      
      integer, dimension(:), allocatable :: &
         ksource

      real (kind=RKIND) ::  & ! for conservation check
         sedCol_init,       &
         sedCol_final,      &
         sedSeafloor_init,  &
         totSediment,       &
         roundOff,          &
         sedError

      real (kind=RKIND), dimension(:), allocatable :: &
         layerDepth_inv , &
         layerDepth_inv2, &
         layerDepth_inv3, &
         sedConcQ, &
         sedConcTmp,&
         fluxConc, &
         WL, WR,   &
         concR,    &
         concL,    &
         vertGridLayer, &
         erosion

      real (kind=RKIND) :: &
         dltR, dltL, cff, cffR, cffL, cu, erosion_o
   !-----------------------------------------------------------------------
   !  Compute bottom skin-friction stress due to combined maximum wave
   !  and current interaction
   !-----------------------------------------------------------------------

   err = 0

   allocate(layerDepth_inv   (kmax ), &
            layerDepth_inv2  (kmax ), &
            layerDepth_inv3  (kmax ), &
            sedConcQ         (kmax ), &
            sedConcTmp       (kmax ), &
            fluxConc         (kmax+1), &
            WL               (kmax), &
            WR               (kmax), &
            concR            (kmax), &
            concL            (kmax), &
            vertGridLayer    (kmax+1), &
            ksource          (kmax), &
            erosion          (kmax))

   vertGridLayer(:) = 0.0_RKIND
   do m = kmax,kmin,-1
      k = kmax+kmin-m
      vertGridLayer(k+1) = vertGridLayer(k) + layerThickness(k)
      layerDepth_inv(k) = 1/layerThickness(m)
      if (m .gt. kmin) then
        layerDepth_inv2(k) = 1/(layerThickness(m) + layerThickness(m-1))
        if (m .lt. kmax) &
           layerDepth_inv3(k) = 1/(layerThickness(m+1)+layerThickness(m) + layerThickness(m-1))
      endif !m > kmin
   enddo

   ! Vertical sinking of suspended particles:   Copy concentration of
   !========= ======= == ========= ==========   suspended sediment into
   ! scratch array "sedConcQ" (q-central, restrict it to be positive) which
   ! are hereafter interpreted as a set of grid-box averaged values for
   ! concentration. then reconstruct vertical profile of "sedConcQQ" in terms
   ! of a set of parabolic segments within each grid box; and, finally,
   ! compute semi-Lagrangian flux due to sinking.

   sedConcQ(:) = 0.0_RKIND
   fluxConc(:) = 0.0_RKIND
   workTend(:) = 0.0_RKIND
   erosion(:)  = 0.0_RKIND
   erosion_o = ERATE !gamma * abs(sinkingRate) * particle_density
   if (bottomStressCell > tau_c) then
      erosion(kmin) = erosion_o * (1.0_RKIND - phi) * fraction_i &
                 * (bottomStressCell/tau_c - 1.0_RKIND)
   endif
   do k = kmin,kmax
      m = kmax+kmin-k
      sedConcQ(k) = sedimentConc(m)
   enddo

   call ocn_sediment_column_total(sedimentConc, layerThickness, kmin, kmax, sedCol_init)
   sedSeafloor_init = netSedimentStorage
   totSediment = sedCol_init + netSedimentStorage 

   do k = kmin,kmax-1
      fluxConc(k+1) = (sedConcQ(k+1)-sedConcQ(k)) * layerDepth_inv2(k)
   enddo
   do k = kmin+1,kmax-1
      m = kmax+kmin-k
      dltR = layerThickness(m) * fluxConc(k+1)
      dltL = layerThickness(m) * fluxConc(k)
      cff = layerThickness(m+1) + 2.0_RKIND * layerThickness(m) + layerThickness(m-1);
      cffR = cff* fluxConc(k+1)
      cffL = cff* fluxConc(k)
      
   ! Apply PPM monotonicity constraint
      if (dltR*dltL .le. 0.0_RKIND) then
         dltR = 0.0_RKIND
         dltL=0.0_RKIND
      elseif (abs(dltR) .gt. abs(cffL)) then
         dltR=cffL
      elseif (abs(dltL) .gt. abs(cffR)) then
         dltL=cffR
      endif
      ! Compute right and left
      ! side values qR,qL of
      cff=(dltR-dltL)*layerDepth_inv3(k) ! parabolic segments
      dltR=dltR-cff*layerThickness(m-1)
      dltL=dltL+cff*layerThickness(m+1)
      concR(k)=sedConcQ(k)+dltR
      concL(k)=sedConcQ(k)-dltL
      WR(k)=( 2.0_RKIND*dltR-dltL )**2
      WL(k)=( dltR-2.0_RKIND*dltL )**2
   enddo ! k
   do k = kmin+1,kmax-2
      dltL=max(WL(k),  small)
      dltR = max(WR(k+1), small)
      concR(k) = (dltR * concR(k) + dltL * concL(k+1))/(dltR + dltL)
      concL(k+1) = concR(k)
   enddo !k

   ! Bottom boundary conditions
   ! strictly monotone
   concR(kmax) = sedConcQ(kmax)
   concL(kmax) = sedConcQ(kmax)
   concR(kmax-1) = sedConcQ(kmax)

   concL(kmin+1)=sedConcQ(kmin)
   concR(kmin)=sedConcQ(kmin)
   concL(kmin)=sedConcQ(kmin)

   do k=kmin,kmax            ! Since the reconciled interfacial values may
      dltR=concR(k)-sedConcQ(k) ! cause non-monotonic
      dltL=sedConcQ(k)-concL(k)    ! behavior of parabolic
      cffR=2.0_RKIND*dltR                  ! segments inside grid
      cffL=2.0_RKIND*dltL                  ! box apply monotonicity
                                          ! constraint again.
      if (dltR*dltL .lt. 0.0_RKIND) then
         dltR=0.0_RKIND
         dltL=0.0_RKIND
      elseif(abs(dltR) .gt. abs(cffL)) then
         dltR=cffL
      elseif (abs(dltL) .gt. abs(cffR)) then
	 dltL=cffR
      endif
      concR(k)=sedConcQ(k)+dltR
      concL(k)=sedConcQ(k)-dltL
   enddo     ! k --> discard everything, except qR,qL
   !
   ! After this moment reconstruction is considered complete. The next
   ! stage is to compute vertical advective fluxes FC. It is expected
   ! that sinking may occur relatively fast, the algorithm is designed
   ! to be free of CFL criterion, which is achieved by allowing
   ! integration bounds for semi-Lagrangian advective flux to use as
   ! many grid boxes in upstream direction as necessary.

   cff=dt*abs(sinkingRate)
   do k=kmin,kmax
      m = kmax+1+kmin-k
      fluxConc(k)=0.0_RKIND
      WL(k)=-vertGridLayer(m)+cff
      WR(k)=layerThickness(m-1) * sedConcQ(k)  ! FC is finite volume flux;
      ksource(k)=k                             ! ksource(:,k) is index of
   enddo                     ! vertical grid box which contains the departure
   do k=kmin,kmax            ! point (restricted by N);
      do ks=k,kmax-1
	 m = kmax+kmin-ks
         if (WL(k) > -vertGridLayer(m)) then
            ksource(k)=ks+1
            fluxConc(k)=fluxConc(k)+WR(ks)
         endif
      enddo                     ! add in content of whole grid boxes participating
   enddo                     ! k --> discard WR  in fluxConc

   do k=kmin,kmax            ! Finalize computation of flux: add fractional part
      ks=ksource(k)
      m=kmax+kmin+1-ks
      cu=min(1.0_RKIND,(WL(k)+vertGridLayer(m))*layerDepth_inv(ks))
      fluxConc(k)=fluxConc(k) + layerThickness(m-1)*cu*( concL(ks) &
            +cu*( 0.5_RKIND*(concR(ks)-concL(ks)) &
            -(1.5_RKIND-cu)*(concR(ks)+concL(ks)-2.0_RKIND*sedConcQ(ks))  ))
   enddo
   do k=kmin,kmax
      m = kmax+kmin-k
! MULTIPLY workTend by layerThickness
!      workTend(m) = (fluxConc(k+1)-fluxConc(k)) * layerDepth_inv(k)/dt * layerThickness(m)
      workTend(m) = ((fluxConc(k+1)-fluxConc(k))/dt + erosion(k))
!      sedimentConc(m) = sedConcQ(k) ...
!	   + workTend(m) * dt /layerThickness(m)
!
      sedConcTmp(m) = sedConcQ(k) + workTend(m) * dt / layerThickness(m)
   enddo
   fluxConc(kmin) = fluxConc(kmin) - erosion(kmin)*dt
   sedimentBottomFlux=fluxConc(kmin)/dt       ! save this
   netSedimentStorage = netSedimentStorage + fluxConc(kmin) ! save this
   netErosionFlux = erosion(kmin)

   ! column conservation check
   call ocn_sediment_column_total(sedConcTmp, layerThickness, kmin, kmax, sedCol_final)
   sedError = totSediment - sedCol_final  - netSedimentStorage 
   roundOff = MAX(totSediment*small, sedCol_final*small, abs(netSedimentStorage)*small, small)

   if (abs(sedError) > roundOff) then
      call mpas_log_write( &
         'Conservation failure in ocn_sediment_settling_flx', MPAS_LOG_WARN)
         call mpas_log_write('totSediment $r.', realArgs=(/totSediment/))
         call mpas_log_write('sedCol_final  $r.', realArgs=(/sedCol_final/))
         call mpas_log_write('netSedimentStorage $r.', realArgs=(/netSedimentStorage/))
         call mpas_log_write('sedError  $r.', realArgs=(/sedError/))
         call mpas_log_write('roundOff $r.', realArgs=(/roundOff/))
         call mpas_log_write('netErosionFlux $r.', realArgs=(/netErosionFlux/))
         call mpas_log_write('sedimentBottomFlux $r.', realArgs=(/sedimentBottomFlux/))
         err = 1

   endif
   if (err == 1) &
      call mpas_log_write( &
         'ERROR: Conservation failure in ocn_sediment_settling_flx', MPAS_LOG_CRIT)

   deallocate(layerDepth_inv, &
      layerDepth_inv2, &
      layerDepth_inv3, &
      sedConcQ, &
      sedConcTmp, &
      fluxConc, &
      WL, &
      WR, &
      concR, &
      concL, &
      vertGridLayer, &
      ksource, &
      erosion)

   end subroutine ocn_sediment_settling_flx !}}

!***********************************************************************
!
!  routine ocn_sediment_init
!
!> \brief   Initializes ocean sediments
!> \author  Nicole Jeffery
!> \date    06/12/2025
!> \details
!>  This routine initializes fields required for sinking tracers
!
!-----------------------------------------------------------------------

   subroutine ocn_sediment_init(domain,err)!{{{

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: ecosysAuxiliary

! input flux components in ecosysAuxiliary
      real (kind=RKIND), dimension(:), pointer :: &
         sedimentBottomFlux,    &
         netSedimentStorage,    &
         netErosionFlux

      err = .false.

      call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

      call mpas_pool_get_array(ecosysAuxiliary, 'sedimentBottomFlux', & 
                                     sedimentBottomFlux)  
      call mpas_pool_get_array(ecosysAuxiliary, 'netErosionFlux', & 
                                     netErosionFlux)  
      call mpas_pool_get_array(ecosysAuxiliary,  'netSedimentStorage', & 
                                     netSedimentStorage)  

      netErosionFlux(:) = 0.0_RKIND
      sedimentBottomFlux(:) = 0.0_RKIND
      netSedimentStorage(:) = 0.0_RKIND  ! needs to go in restart


   end subroutine ocn_sediment_init!}}}

!***********************************************************************
!
!  routine ocn_sediment_column_total
!
!> \brief   sums the total column sediments
!> \author  Nicole Jeffery
!> \date    09/15/2025
!> \details
!>  This routine computes the total column sediment mass
!
!-----------------------------------------------------------------------

   subroutine ocn_sediment_column_total(sedimentConc,layerThickness, &
                                        kmin, kmax, sedColTot)!{{{


      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------
                                            
      real (kind=RKIND), dimension(:), intent(in) :: &
          sedimentConc, & ! kg/m3
          layerThickness  ! m
          
      integer, intent(in):: &
          kmin, kmax

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------
                                                                             
      real (kind=RKIND), intent(out) :: &
         sedColTot             ! kg/m2

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer ::          &
         k

      sedColTot = 0.0_RKIND
     
      do k = kmin,kmax
         sedColTot = sedColTot + sedimentConc(k) * layerThickness(k)
      enddo

   end subroutine ocn_sediment_column_total!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_explicit_bottom_drag_tend
!
!> \brief   Computes bottom shear stress
!> \author  Nicole Jeffery
!> \date    17 August 2025
!> \details
!>  This routine uses the explicit bottom drag formulation
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_bottom_shear_stress(normVelocity, &
                                 KECell, layerThickEdgeDrag, bottomStressCell) !{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normVelocity,      &!< [in] normal velocity
         KECell,            &!< [in] kinetic energy at cell
         layerThickEdgeDrag  !< [in] mean layer thickness at edge

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
         bottomStressCell       !< [inout] bottom shear stress magnitude (N/m2)

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell,      &! loop index for cells
         nCells
         
      call mpas_timer_start('explicit bottom stress for erosion')

      ! Explicit bottom drag term:
      ! du/dt = ... - c |u| u / h
      ! appied to bottom layer only.
      ! This term comes from the bottom boundary condition in the vertical
      ! momentum mixing, and is explicit if both |u| and u are chosen to be at
      ! time level n.  Magnitude of the bottom stress is c * rho * 2 * KE
      !--------------------------------------------------------------------

      ! Need 0,1,2 halo cells
      nCells = nCellsHalo(2)

      ! compute bottom drag magnitude at cell centers
#ifdef MPAS_OPENACC
      !$acc parallel loop present(KECell, maxLevelCell, bottomStressCell)
#else
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCells
         ! the magnitude of the bottom drag is rho*C*u**2 = rho*C*(2*KE)
         bottomStressCell(iCell) = config_explicit_bottom_drag_coeff * rho_sw &
                           * 2.0_RKIND * KECell(maxLevelCell(iCell),iCell)
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      call mpas_timer_stop('explicit bottom stress for erosion')

   !--------------------------------------------------------------------

   end subroutine ocn_compute_bottom_shear_stress!}}}

!***********************************************************************
end module ocn_sediment_settling

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
