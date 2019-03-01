!======================================================================================================================!
!                                                                                                                      !
!     Program:      ACCESS                                                                                             !
!                   Atmospheric Chemistry and Canopy Exchange Simulation System                                        !
!                   Full BVOC chemistry version                                                                        !
!                                                                                                                      !
!     Version:      3.1.0                                                                                              !
!                                                                                                                      !
!======================================================================================================================!
!                                                                                                                      !
!     Last Update:  Dec 2017                                                                                           !
!                                                                                                                      !
!     Contact:      Rick D. Saylor, PhD                                                                                !
!                   Physical Scientist                                                                                 !
!                   U. S. Department of Commerce                                                                       !
!                   National Oceanic and Atmospheric Administration                                                    !
!                   Air Resources Laboratory                                                                           !
!                   Atmospheric Turbulence and Diffusion Division                                                      !
!                   456 S. Illinois Ave                                                                                !
!                   Oak Ridge, TN 37830                                                                                !
!                   email: Rick.Saylor@noaa.gov                                                                        !
!                                                                                                                      !
!**********************************************************************************************************************!
!                   NOTE: See Legal Notice in Main.f90                                                                 !
!**********************************************************************************************************************!
!                                                                                                                      !
!     Module:       GlobalData                                                                                         !
!                                                                                                                      !
!     Description:  Contains all global variables used by the model                                                    !
!                                                                                                                      !
!======================================================================================================================!
module GlobalData
  implicit none
  public
  save

  ! double precision kind parameter
  integer, parameter :: dp=kind(0.d0)
  ! 4 bytes integer
  integer, parameter :: i4=kind(1)

  ! tstart - time of simulation start
  real(kind=dp) :: tstart
  ! tend = time of simulation end
  real(kind=dp) :: tend
  ! dt - integration time step (seconds) 
  real(kind=dp) :: dt
  ! dthalf - half of integration time step (seconds)
  real(kind=dp) :: dthalf
  ! dtenv - time step for input environmental data (seconds)
  real(kind=dp) :: dtenv
  ! ntenv - total number of environmental time steps for simulation
  integer(kind=i4) :: ntenv
  ! dtout - time step for output (seconds)
  !       MUST be even multiple of dt
  real(kind=dp) :: dtout
  ! ntout - total number of output time steps
  integer(kind=i4) :: ntout
  ! zadeg - zenith angle at current meteorological data time step (degrees)
  real(kind=dp) :: zadeg
  ! zarad - zenith angle at current meteorological data time step (radians)
  real(kind=dp) :: zarad

  ! t - current simulation time
  real(kind=dp) :: t
  ! nt - current output time step number
  integer(kind=i4) :: nt
  ! nte - current met data time step number
  integer(kind=i4) :: nte
  ! zcurrent - current z value (only valid within chemistry integration!) (m)
  real(kind=dp) :: zcurrent
  ! zpt - current vertical grid point (only valid within chemistry integration!)
  integer(kind=i4) :: zpt

  ! number of grid nodes
  integer(kind=i4), parameter :: npts = 87
  ! number of grid nodes in canopy
  integer(kind=i4)            :: ncnpy

  ! vertical grid definition (cm)
  real(kind=dp), dimension(npts) :: z
  ! within canopy grid resolution (cm)
  real(kind=dp)                  :: dzhc
  ! stretched grid parameter ()
  real(kind=dp)                  :: alfa

  ! column bottom and top (cm)
  real(kind=dp)                  :: z0
  real(kind=dp)                  :: zi

  ! met measurement reference height (cm)
  real(kind=dp)  :: zref

  ! species measurement reference height (cm)
  real(kind=dp)  :: zcnsref

  ! canopy height ()
  real(kind=dp)  :: hc
  ! canopy height (cm)
  real(kind=dp)  :: hccm

  ! roughness length for momentum (cm)
  real(kind=dp)  :: z0m
  ! roughness length for heat (cm)
  real(kind=dp)  :: z0h
  ! zero-plane displacement height (cm)
  real(kind=dp)  :: d

  ! layer leaf area density of canopy (cm2 leaf/cm3)
  real(kind=dp), dimension(npts) :: lad
  ! layer leaf area index of canopy (cm2 leaf/cm2)
  real(kind=dp), dimension(npts) :: lai
  ! top-down cumulative leaf area index of canopy (cm2 leaf/cm2)
  real(kind=dp), dimension(npts) :: clai
  ! integrated leaf area index of canopy (cm2 leaf/cm2)
  real(kind=dp)                  :: laitot

!$INSERT MECHANISM-SPECIFIC SPECIES PARAMETERS

  ! Gas-phase photolytic rate constant data
  ! nphoto - number of J values calculated in CalcJPhoto
  integer(kind=4), parameter :: nphoto = 51
  ! kphoto - gas-phase photolytic rate constants
  real(kind=dp), dimension(nphoto) :: kphoto

  ! Photolysis reaction indices
  !
  ! ozone
  !   O3 + hv -> O(1D) + O2
  integer(kind=i4), parameter :: jo3a      = 1
  !   O3 + hv -> O(3P) + O2
  integer(kind=i4), parameter :: jo3b      = 2
  ! hydrogen peroxide
  !   H2O2 + hv -> 2OH
  integer(kind=i4), parameter :: jh2o2     = 3
  ! nitrogen dioxide
  !   NO2 + hv -> NO + O(3P)
  integer(kind=i4), parameter :: jno2      = 4
  ! nitrate radical
  !   NO3 + hv -> NO + O2
  integer(kind=i4), parameter :: jno3a     = 5
  !   NO3 + hv -> NO2 + O(3P)
  integer(kind=i4), parameter :: jno3b     = 6
  ! nitrous acid
  !   HONO + hv -> NO + OH
  integer(kind=i4), parameter :: jhono     = 7
  ! nitric acid
  !   HNO3 + hv -> NO2 + OH
  integer(kind=i4), parameter :: jhno3     = 8
  ! peroxynitric acid
  !   HO2NO2 + hv -> NO2 + HO2
  integer(kind=i4), parameter :: jho2no2   = 9
  ! formaldehyde
  !   HCHO + hv -> HCO + H
  integer(kind=i4), parameter :: jhchoa    = 10
  !   HCHO + hv -> CO + H2
  integer(kind=i4), parameter :: jhchob    = 11
  ! acetaldehyde
  !   CH3CHO + hv -> CH3 + HCO
  integer(kind=i4), parameter :: jch3cho   = 12
  ! propionaldehyde
  !   C2H5CHO + hv -> C2H5 + HCO
  integer(kind=i4), parameter :: jc2h5cho  = 13
  ! n-butyraldehyde
  !   n-C3H7CHO + hv -> n-C3H7 + HCO
  integer(kind=i4), parameter :: jc3h7choa = 14
  !   n-C3H7CHO + hv -> CH3CHO + CH2=CH2
  integer(kind=i4), parameter :: jc3h7chob = 15
  ! isobutyraldehyde
  !   i-C3H7CHO + hv -> i-C3H7 + HCO
  integer(kind=i4), parameter :: jic3h7cho = 16
  ! methacrolein
  !   MACR + hv -> CH2=C(CH3) + HCO
  integer(kind=i4), parameter :: jmacra    = 17
  !   MACR + hv -> CH2+C(CH3)CO + H
  integer(kind=i4), parameter :: jmacrb    = 18
  ! acetone
  !   CH3COCH3 + hv -> CH3CO + CH3
  integer(kind=i4), parameter :: jacet     = 19
  ! methyl ethyl ketone
  !   MEK + hv -> CH3CO + CH2CH3
  integer(kind=i4), parameter :: jmek      = 20
  ! methyl vinyl ketone
  !   MVK + hv -> CH3CH=CH2 + CO
  integer(kind=i4), parameter :: jmvka     = 21
  !   MVK + hv -> CH3CO + CH2=CH2
  integer(kind=i4), parameter :: jmvkb     = 22
  ! glyoxal
  !   CHOCHO + hv -> 2CO + H2
  integer(kind=i4), parameter :: jglyoxa   = 23
  !   CHOCHO + hv -> HCHO + CO
  integer(kind=i4), parameter :: jglyoxb   = 24
  !   CHOCHO + hv -> 2HCO
  integer(kind=i4), parameter :: jglyoxc   = 25
  ! methyl gloxal
  !   MGLYOX + hv -> CH3CO + HCO
  integer(kind=i4), parameter :: jmglyox   = 26
  ! biacetyl
  !   BIACET + hv -> 2CH3CO
  integer(kind=i4), parameter :: jbiacet   = 27
  ! methyl hydroperoxide
  !   CH3OOH + hv -> CH3O + OH
  integer(kind=i4), parameter :: jch3ooh   = 28
  ! methyl nitrate
  !   CH3NO3 + hv -> CH3O + NO2
  integer(kind=i4), parameter :: jch3no3   = 29
  ! ethyl nitrate
  !   C2H5NO3 + hv -> C2H5O + NO2
  integer(kind=i4), parameter :: jc2h5no3  = 30
  ! n-propyl nitrate
  !   n-C3H7NO3 + hv -> n-C3H7O + NO2
  integer(kind=i4), parameter :: jnc3h7no3 = 31
  ! isopropyl nitrate
  !   i-C3H7NO3 + hv -> i-C3H7O + NO2
  integer(kind=i4), parameter :: jic3h7no3 = 32
  ! t-butyl nitrate
  !   t-C4H9NO3 + hv -> t-C4H9O + NO2
  integer(kind=i4), parameter :: jtc4h9no3 = 33
  ! 1-hydroxy-2-propanone nitrate
  !   CH3C(O)CH2NO3 + hv -> CH3C(O)CH2O + NO2
  integer(kind=i4), parameter :: jnoaa     = 34
  !   CH3C(O)CH2NO3 + hv -> CH3CO + HCHO + NO2
  integer(kind=i4), parameter :: jnoab     = 35
  ! glycolaldehyde
  !   HOCH2CHO + hv -> HCO + HCHO + HO2
  integer(kind=i4), parameter :: jhoch2cho = 36
  ! hydroxyacetone (acetol)
  !   ACETOL + hv -> CH3CO + HCHO + HO2
  integer(kind=i4), parameter :: jacetol   = 37
  ! peroxyacetic acid
  !   CH3C(O)OOH + hv -> CH3O2 + OH
  integer(kind=i4), parameter :: jpaa      = 38
  ! pyruvic acid
  !   CH3COCOOH + hv -> CH3CO3 + HO2
  integer(kind=i4), parameter :: jpyra     = 39
  ! hydroperoxyaldehyde
  !   HPALD + hv -> products
  integer(kind=i4), parameter :: jhpald    = 40
  ! acrolein
  !   CH2=CHCHO + hv -> products
  integer(kind=i4), parameter :: jacro     = 41
  ! hydroxy methyl hydroperoxide
  !   HOCH2OOH + hv -> HOCH2O + OH
  integer(kind=i4), parameter :: jhmhp     = 42
  ! nitroxy ethanol
  !   CH2(OH)CH2(ONO2) + hv -> CH2(OH)CH2(O) + NO2
  integer(kind=i4), parameter :: jneth     = 43
  ! nitroxy acetone
  !   CH3COCH2(ONO2) + hv -> CH3COCH2(O) + NO2
  integer(kind=i4), parameter :: jnacet    = 44
  ! peroxy acetyl nitrate
  !   PAN + hv -> CH3CO(OO) + NO2
  integer(kind=i4), parameter :: jpana     = 45
  ! peroxy acetyl nitrate
  !   PAN + hv -> CH3CO(O) + NO3
  integer(kind=i4), parameter :: jpanb     = 46
  ! peroxy propionyl nitrate
  !   PPN + hv -> CH3CH2CO(OO) + NO2
  integer(kind=i4), parameter :: jppna     = 47
  ! peroxy propionyl nitrate
  !   PPN + hv -> CH3CH2CO(O) + NO3
  integer(kind=i4), parameter :: jppnb     = 48
  ! dinitrogen pentoxide
  !   N2O5 + hv -> NO3 + NO + O(3P)
  integer(kind=i4), parameter :: jn2o5a    = 49
  ! dinitrogen pentoxide
  !   N2O5 + hv -> NO3 + NO2
  integer(kind=i4), parameter :: jn2o5b    = 50
  ! benzaldehyde
  !   BenzALD + hv -> HCO + HO2 + CO
  integer(kind=i4), parameter :: jbald     = 51

  ! Soil exchange data
  ! vs - soil exchange coefficients (cm/s)
  real(kind=dp), dimension(ninteg) :: vs
  ! vsh2o - soil exchange coefficient for water vapor (cm/s)
  real(kind=dp)                    :: vsh2o
  ! qsoil - effective soil molar humidity (mol/cm3)
  real(kind=dp)                    :: qsoil
  ! ksoil - thermal conductivity of soil (W/m-K)
  real(kind=dp)                    :: ksoil
  ! tsoilk - soil/litter temperature (K)
  real(kind=dp)                    :: tsoilk
  ! tk0 - lower boundary condition for air temperature
  !       derived from the surface energy balance (K)
  real(kind=dp)                    :: tk0
  ! dtdzsoil - temperature gradient in soil surface (K/m)
  real(kind=dp)                    :: dtdzsoil
  ! gbg - ground boundary layer conductance (mol/m2-s)
  real(kind=dp)                    :: gbg
  ! effrhsoil - effective soil fractional relative humidity
  real(kind=dp)                    :: effrhsoil

  ! nsoil - number of species with non-zero soil compensation point concentrations
  integer(kind=i4) :: nsoil
  ! csoil - soil compensation point concentration (molecules/cm3)
  real(kind=dp), dimension(ninteg) :: csoil

  ! rbg - ground boundary layer resistance (s/cm)
  real(kind=dp)                    :: rbg
  ! rsoil - resistance to diffusion of water vapor thru soil pore space (s/cm)
  real(kind=dp)                    :: rsoil
  ! rsoil - resistance to diffusion thru soil pore space for chemical species (s/cm)
  real(kind=dp), dimension(ninteg) :: rsoill

  ! isoiltype - soil type
  integer(kind=i4) :: isoiltype
  ! stheta - volumetric soil water content (m3/m3)
  real(kind=dp)                    :: stheta
  ! sattheta - saturation volumetric soil water content (m3/m3)
  real(kind=dp)                    :: sattheta
  ! rtheta - residual volumetric soil water content (m3/m3)
  real(kind=dp)                    :: rtheta
  ! sbcoef - Clapp and Hornberger exponent
  real(kind=dp)                    :: sbcoef
  ! satphi - saturation matric potential (suction) (m)
  real(kind=dp)                    :: satphi
  ! dsoil - depth of topsoil (cm)
  real(kind=dp)                    :: dsoil

  ! saved soil physics data
  ! vsh2oout - soil exchange coefficient for water vapor (cm/s)
  real(kind=dp), allocatable     :: vsh2oout(:)
  ! qsoilout - effective soil molar humidity (mol/cm3)
  real(kind=dp), allocatable     :: qsoilout(:)
  ! effrhsoilout - effective soil fractional relative humidity ()
  real(kind=dp), allocatable     :: effrhsoilout(:)
  ! rbgout - ground boundary layer resistance (s/cm)
  real(kind=dp), allocatable     :: rbgout(:)
  ! gbgout - ground aerodynamic conductance (mol/m2-s)
  real(kind=dp), allocatable     :: gbgout(:)
  ! rsoilout - resistance to diffusion of water vapor through soil pore space (s/cm)
  real(kind=dp), allocatable     :: rsoilout(:)
  ! tsoilkout - soil temperature over simulation (K)
  real(kind=dp), allocatable     :: tsoilkout(:)
  ! tk0out - lower boundary condition for air temperature over simulation (K)
  real(kind=dp), allocatable     :: tk0out(:)

  ! Soil type data
  ! Soil type strings
   character(len=16), parameter, dimension(11) ::  &
                     ssoiltype     = (/         &
                                   'sand            ', &
                                   'loamy sand      ', &
                                   'sandy loam      ', &
                                   'loam            ', &
                                   'silt loam       ', &
                                   'sandy clay loam ', &
                                   'clay loam       ', &
                                   'silty clay loam ', &
                                   'sandy clay      ', &
                                   'silty clay      ', &
                                   'clay            ' &
                                             /)

  ! xsattheta - volumetric soil water content at saturation (m3/m3)
  real(kind=dp), dimension(11)     :: xsattheta
  data xsattheta /0.437,0.437,0.453,0.463,0.501,0.398,0.464,0.471,0.430,0.479,0.475/
  ! xrtheta - residual volumetric soil water content (m3/m3)
  real(kind=dp), dimension(11)     :: xrtheta
  data xrtheta /0.020,0.035,0.041,0.027,0.015,0.068,0.075,0.040,0.109,0.056,0.090/
  ! xsbcoef - Clapp-Hornberger exponents
  real(kind=dp), dimension(11)     :: xsbcoef
  data xsbcoef /4.05,4.38,4.90,5.39,5.30,7.12,8.52,7.75,10.4,10.4,11.4/
  ! xsatphi - soil matric potential (suction) at saturation (cm)
  real(kind=dp), dimension(11)     :: xsatphi
  data xsatphi /12.1,9.0,21.8,78.6,47.8,29.9,35.6,63.0,15.3,49.0,40.5/

  ! Soil thermal conductivity data
  ! Estimated from:
  ! Kersten, M. S. (1949) Thermal Properties of Soils. University of Minnesota, Engineering
  !  Experiment Station Bulletin No. 28 (http://hdl.handle.net/11299/124271).
  !
  ! ksoildry - thermal conductivity for dry soil
  real(kind=dp), dimension(11)     :: ksoildry
  data ksoildry /0.30,0.28,0.26,0.20,0.16,0.20,0.16,0.16,0.20,0.16,0.16/
  ! ksoilwet - thermal conductivity for wet soil
  real(kind=dp), dimension(11)     :: ksoilwet
  data ksoilwet /2.20,2.10,2.04,1.66,1.28,1.70,1.28,1.28,1.74,1.28,1.28/

  ! q - emissions rate of all species at the current time (molecules/cm3-s)
  real(kind=dp), dimension(npts, ninteg) :: q
  ! vd - dry deposition exchange coefficient (cm/s)
  real(kind=dp), dimension(npts, ninteg) :: vd
  ! rb - leaf boundary layer resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rb
  ! rm - leaf mesophyll resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rm
  ! rc - leaf cuticular resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rc
  ! rs - leaf stomatal resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rs

  ! gp - dry deposition compensation point :: (molecules/cm3)
  real(kind=dp), dimension(npts, ninteg) :: gp

  ! rakv - aerodynamic resistance calculated from Kv profile (s/cm)
  real(kind=dp), allocatable :: rakv(:)
  ! ra - aerodynamic resistance calculated from Viney et al. (s/cm)
  real(kind=dp), allocatable :: ra(:)
  ! zolout - Monin-Obukhov stability parameter
  real(kind=dp), allocatable :: zolout(:)
  ! ribout - local Richardson number stability parameter
  real(kind=dp), allocatable :: ribout(:)
  ! gaero - aerodynamic conductance (mol/m2-s)
  real(kind=dp), allocatable :: gaero(:)

  ! ctot - total solution array for current integration timestep (molecules/cm3)
  real(kind=dp), dimension(npts, ntotal) :: ctot
  ! ctotp - total solution array for previous integration timestep (molecules/cm3)
  real(kind=dp), dimension(npts, ntotal) :: ctotp

  ! cint - integrated species solution array for current integration timestep (??)
  real(kind=dp), dimension(npts, ninteg) :: cint
  ! cintp - integrated species solution array for previous integration timestep (??)
  real(kind=dp), dimension(npts, ninteg) :: cintp

  ! cfxed - fixed species array (molecules/cm3)
  real(kind=dp), dimension(npts, nfixed) :: cfxed

  ! css - steady state species solution array for current integration time step (molecules/cm3)
  real(kind=dp), dimension(npts, nsstate) :: css
  ! cssp - steady state species solution array for previous integration time step (molecules/cm3)
  real(kind=dp), dimension(npts, nsstate) :: cssp

  ! cout - saved total solution array over entire simulation (molecules/cm3)
  real(kind=dp), allocatable :: cout(:,:,:)
  ! vdout - saved deposition velocities array over entire simulation (cm/s)
  real(kind=dp), allocatable :: vdout(:,:,:)
  ! vsout - saved surface exchange velocities over entire simulation (cm/s)
  real(kind=dp), allocatable :: vsout(:,:)

  ! rbout - saved leaf boundary layer resistance (s/cm)
  real(kind=dp), allocatable :: rbout(:,:,:)
  ! rcout - saved leaf cuticular resistance (s/cm)
  real(kind=dp), allocatable :: rcout(:,:,:)
  ! rmout - saved leaf mesophyll resistance (s/cm)
  real(kind=dp), allocatable :: rmout(:,:,:)
  ! rsout - saved leaf stomatal resistance (s/cm)
  real(kind=dp), allocatable :: rsout(:,:,:)
  ! rsoillout - saved soil resistance of chemical species (s/cm)
  real(kind=dp), allocatable :: rsoillout(:,:)

  ! qout - saved emissions array over entire simulation (molecules/cm3-s)
  real(kind=dp), allocatable :: qout(:,:,:)
  ! vfout - saved vertical fluxes over entire simulation (nmol/m2-s))
  real(kind=dp), allocatable :: vfout(:,:,:)
  ! emtout - saved budget emissions rates for entire simulation (ppbv/s)
  real(kind=dp), allocatable :: emtout(:,:,:)
  ! cntout - saved constrained species rates for entire simulation (ppbv/s)
  real(kind=dp), allocatable :: cntout(:,:,:)
  ! dptout - saved budget deposition rates for entire simulation (ppbv/s)
  real(kind=dp), allocatable :: dptout(:,:,:)
  ! vttout - saved budget vertical transport rates for entire simulation (ppbv/s)
  real(kind=dp), allocatable :: vttout(:,:,:)
  ! chtout - saved budget chemical production/loss rates for entire simulation (ppbv/s)
  real(kind=dp), allocatable :: chtout(:,:,:)
  ! rxnout - saved reaction rates for entire simulation (molecules/cm3-s)
  real(kind=dp), allocatable :: rxnout(:,:,:)
  ! ksout - saved reaction rate coefficients for entire simulation (molecules-cm-s units)
  real(kind=dp), allocatable :: ksout(:,:,:)
  ! timeout - saved t values corresponding to output saved in cout
  real(kind=dp), allocatable :: timeout(:)
  ! sdtout - saved sdatetime values corresponding to output saved in cout
  character(len=19), allocatable :: sdtout(:)

  ! Canopy physics data
  ! ppfd_direct - direct beam component of incoming PPFD (umol/m2-s)
  real(kind=dp)                  :: ppfd_direct
  ! ppfd_diffus - diffuse component of incoming PPFD (umol/m2-s)
  real(kind=dp)                  :: ppfd_diffus
  ! nir_direct - direct beam component of incoming NIR (W/m2)
  real(kind=dp)                  :: nir_direct
  ! nir_diffus - diffuse component of incoming NIR (W/m2)
  real(kind=dp)                  :: nir_diffus
  ! ppfd_sun - sunlit fraction PPFD (umol/m2-s)
  real(kind=dp), dimension(npts) :: ppfd_sun
  ! ppfd_shd - shaded fraction PPFD (umol/m2-s)
  real(kind=dp), dimension(npts) :: ppfd_shd
  ! ppfd_wgt - sun/shade weighted PPFD (umol/m2-s)
  real(kind=dp), dimension(npts) :: ppfd_wgt
  ! nir_sun - sunlit fraction NIR (W/m2)
  real(kind=dp), dimension(npts) :: nir_sun
  ! nir_shd - shaded fraction NIR (W/m2)
  real(kind=dp), dimension(npts) :: nir_shd
  ! nir_wgt - sun/shade weighted NIR (W/m2)
  real(kind=dp), dimension(npts) :: nir_wgt
  ! lw_up - upwelling long-wave radiation (W/m2)
  real(kind=dp), dimension(npts) :: lw_up
  ! lw_dn - downwelling long-wave radiation (W/m2)
  real(kind=dp), dimension(npts) :: lw_dn
  ! rt_sun - total absorbed PPFD and NIR radiation in sunlit fraction (W/m2)
  real(kind=dp), dimension(npts) :: rt_sun
  ! rt_shd - total absorbed PPFD and NIR radiation in shaded fraction (W/m2)
  real(kind=dp), dimension(npts) :: rt_shd
  ! rt_wgt - sun/shade weighted total absorbed PPFD and NIR radiation (W/m2)
  real(kind=dp), dimension(npts) :: rt_wgt
  ! rabs_sun - total absorbed radiation (PPFD+NIR+IR) in sunlit fraction (W/m2)
  real(kind=dp), dimension(npts) :: rabs_sun
  ! rabs_shd - total absorbed radiation (PPFD+NIR+IR) in shaded fraction (W/m2)
  real(kind=dp), dimension(npts) :: rabs_shd
  ! rabs_wgt - sun/shade weighted total absorbed radiation (PPFD+NIR+IR) (W/m2)
  real(kind=dp), dimension(npts) :: rabs_wgt
  ! fsun - sunlit fraction of canopy ()
  real(kind=dp), dimension(npts) :: fsun
  ! fshd - shaded fraction of canopy ()
  real(kind=dp), dimension(npts) :: fshd
  ! tl_sun - leaf temperature in sunlit fraction (K)
  real(kind=dp), dimension(npts) :: tl_sun
  ! tl_shd - leaf temperature in shaded fraction (K)
  real(kind=dp), dimension(npts) :: tl_shd
  ! tlp_sun - previous iteration leaf temperature in sunlit fraction (K)
  real(kind=dp), dimension(npts) :: tlp_sun
  ! tlp_shd - previous iteration leaf temperature in shaded fraction (K)
  real(kind=dp), dimension(npts) :: tlp_shd
  ! tl_wgt - sun/shade weighted leaf temperature (K)
  real(kind=dp), dimension(npts) :: tl_wgt
  ! gs_sun - leaf stomatal conductance in sunlit fraction (mol/m2-s)
  real(kind=dp), dimension(npts) :: gs_sun
  ! gs_shd - leaf stomatal conductance in shaded fraction (mol/m2-s)
  real(kind=dp), dimension(npts) :: gs_shd
  ! gs_wgt - sun/shade weighted leaf stomatal conductance (mol/m2-s)
  real(kind=dp), dimension(npts) :: gs_wgt

  ! gl_sun - leaf water vapor conductance in sunlit fraction (mol/m2-s)
  real(kind=dp), dimension(npts) :: gl_sun
  ! gl_shd - leaf water vapor conductance in shaded fraction (mol/m2-s)
  real(kind=dp), dimension(npts) :: gl_shd
  ! gb - leaf boundary layer conductance (mol/m2-s)
  real(kind=dp), dimension(npts) :: gb
  ! esun - evapotranspiration flux for the sunlit canopy fraction (mol/m2-s)
  real(kind=dp), dimension(npts) :: esun
  ! eshd - evapotranspiration flux for the shaded canopy fraction (mol/m2-s)
  real(kind=dp), dimension(npts) :: eshd
  ! hsun - sensible heat flux for the sunlit canopy fraction (W/m2)
  real(kind=dp), dimension(npts) :: hsun
  ! hshd - sensible heat flux for the shaded canopy fraction (W/m2)
  real(kind=dp), dimension(npts) :: hshd

  ! rs_sun - leaf stomatal resistance in sunlit fraction (s/m)
  real(kind=dp), dimension(npts) :: rs_sun
  ! rs_shd - leaf stomatal resistance in shaded fraction (s/m)
  real(kind=dp), dimension(npts) :: rs_shd
  ! rs_wgt - sun/shade weighted leaf stomatal resistance (s/m)
  real(kind=dp), dimension(npts) :: rs_wgt
  ! anet_sun - net photosynthetic assimilation in sunlit fraction (umol/m2-s)
  real(kind=dp), dimension(npts) :: anet_sun
  ! anet_shd - net photosynthetic assimilation in shaded fraction (umol/m2-s)
  real(kind=dp), dimension(npts) :: anet_shd
  ! anet_wgt - sun/shade weighted net photosynthetic assimilation (umol/m2-s)
  real(kind=dp), dimension(npts) :: anet_wgt
  ! epsleaf - leaf emissivity ()
  real(kind=dp), parameter       :: epsleaf=0.97
  ! epsgrnd - soil/litter surface emissivity ()
  real(kind=dp), parameter       :: epsgrnd=0.95
  ! sbsig - Stefan-Boltzmann constant (W/m2-K4)
  real(kind=dp), parameter       :: sbsig=5.67D-08
  ! cpair - heat capacity of air, J/mol-K
  real(kind=dp), parameter       :: cpair=29.3

  ! x - leaf angle distribution parameter in ellipsoidal kb function
  real(kind=dp)                  :: x
  ! kd - diffuse radiation extinction coefficient
  real(kind=dp)                  :: kd
  ! g0 - Medlyn et al. stomatal conductance parameter
  real(kind=dp)                  :: g0
  ! g1 - Medlyn et al. stomatal conductance parameter
  real(kind=dp)                  :: g1
  ! gsmax - maximum stomatal conductance (mol/m2-s)
  real(kind=dp), parameter       :: gsmax=10.0
  ! dleaf - characteristic leaf dimension = 0.75*(leaf width) (m)
  real(kind=dp)                  :: dleaf
  ! cca - atmospheric CO2 mixing ratio (umol/mol)
  real(kind=dp), parameter       :: cca=405.0
  ! ltol - iteration tolerance used in calculating leaf temperature
  real(kind=dp), parameter       :: ltol=0.001
  ! ctol - iteration tolerance used in calculating overall canopy physics
  real(kind=dp), parameter       :: ctol=0.01
 
  ! saved canopy physics data
  ! ppfddirout - direct beam component of incoming PPFD (umol/m2-s)
  real(kind=dp), allocatable     :: ppfddirout(:)
  ! ppfddifout - diffuse component of incoming PPFD (umol/m2-s)
  real(kind=dp), allocatable     :: ppfddifout(:)
  ! nirdirout - direct beam component of incoming NIR (W/m2)
  real(kind=dp), allocatable     :: nirdirout(:)
  ! nirdifout - diffuse component of incoming NIR (W/m2)
  real(kind=dp), allocatable     :: nirdifout(:)
  ! ppfdsunout - sunlit fraction PPFD (umol/m2-s)
  real(kind=dp), allocatable     :: ppfdsunout(:,:)
  ! ppfdshdout - shaded fraction PPFD (umol/m2-s)
  real(kind=dp), allocatable     :: ppfdshdout(:,:)
  ! ppfdwgtout - sun/shade weighted PPFD (umol/m2-s)
  real(kind=dp), allocatable     :: ppfdwgtout(:,:)
  ! nirsunout - sunlit fraction NIR (W/m2)
  real(kind=dp), allocatable     :: nirsunout(:,:)
  ! nirshdout - shaded fraction NIR (W/m2)
  real(kind=dp), allocatable     :: nirshdout(:,:)
  ! nirwgtout - sun/shade weighted NIR (W/m2)
  real(kind=dp), allocatable     :: nirwgtout(:,:)
  ! lwupout - upwelling long-wave radiation (W/m2)
  real(kind=dp), allocatable     :: lwupout(:,:)
  ! lwdnout - downwelling long-wave radiation (W/m2)
  real(kind=dp), allocatable     :: lwdnout(:,:)
  ! rtsunout - total absorbed PPFD and NIR radiation in sunlit fraction (W/m2)
  real(kind=dp), allocatable     :: rtsunout(:,:)
  ! rtshdout - total absorbed PPFD and NIR radiation in shaded fraction (W/m2)
  real(kind=dp), allocatable     :: rtshdout(:,:)
  ! rtwgtout - sun/shade weighted total absorbed PPFD and NIR radiation (W/m2)
  real(kind=dp), allocatable     :: rtwgtout(:,:)
  ! rabssunout - total absorbed radiation (PPFD+NIR+IR) in sunlit fraction (W/m2)
  real(kind=dp), allocatable     :: rabssunout(:,:)
  ! rabsshdout - total absorbed radiation (PPFD+NIR+IR) in shaded fraction (W/m2)
  real(kind=dp), allocatable     :: rabsshdout(:,:)
  ! rabswgtout - sun/shade weighted total absorbed radiation (PPFD+NIR+IR) (W/m2)
  real(kind=dp), allocatable     :: rabswgtout(:,:)
  ! fsunout - sunlit fraction of canopy ()
  real(kind=dp), allocatable     :: fsunout(:,:)
  ! fshdout - shaded fraction of canopy ()
  real(kind=dp), allocatable     :: fshdout(:,:)
  ! tlsunout - leaf temperature in sunlit fraction (K)
  real(kind=dp), allocatable     :: tlsunout(:,:)
  ! tlshdout - leaf temperature in shaded fraction (K)
  real(kind=dp), allocatable     :: tlshdout(:,:)
  ! tlwgtout - sun/shade weighted leaf temperature (K)
  real(kind=dp), allocatable     :: tlwgtout(:,:)
  ! gssunout - leaf stomatal conductance in sunlit fraction (mol/m2-s)
  real(kind=dp), allocatable     :: gssunout(:,:)
  ! gsshdout - leaf stomatal conductance in shaded fraction (mol/m2-s)
  real(kind=dp), allocatable     :: gsshdout(:,:)
  ! gswgtout - sun/shade weighted leaf stomatal conductance (mol/m2-s)
  real(kind=dp), allocatable     :: gswgtout(:,:)
  ! rssunout - leaf stomatal resistance in sunlit fraction (s/m)
  real(kind=dp), allocatable     :: rssunout(:,:)
  ! rsshdout - leaf stomatal resistance in shaded fraction (s/m)
  real(kind=dp), allocatable     :: rsshdout(:,:)
  ! rswgtout - sun/shade weighted leaf stomatal resistance (s/m)
  real(kind=dp), allocatable     :: rswgtout(:,:)
  ! anetsunout - net photosynthetic assimilation in sunlit fraction (umol/m2-s)
  real(kind=dp), allocatable     :: anetsunout(:,:)
  ! anetshdout - net photosynthetic assimilation in shaded fraction (umol/m2-s)
  real(kind=dp), allocatable     :: anetshdout(:,:)
  ! anetwgtout - sun/shade weighted net photosynthetic assimilation (umol/m2-s)
  real(kind=dp), allocatable     :: anetwgtout(:,:)

  ! Calculated met data saved for output
  ! tkout - saved temperature profiles (K)
  real(kind=dp), allocatable :: tkout(:,:)
  ! pmbout - saved pressure profiles (mb)
  real(kind=dp), allocatable :: pmbout(:,:)
  ! qhout - saved specific humidity profiles (g/kg)
  real(kind=dp), allocatable :: qhout(:,:)
  ! rhout - saved relative humidity profiles (%)
  real(kind=dp), allocatable :: rhout(:,:)
  ! ubarout - saved mean wind speed profiles (cm/s)
  real(kind=dp), allocatable :: ubarout(:,:)
  ! kvout - saved eddy diffusivity profiles (cm2/s)
  real(kind=dp), allocatable :: kvout(:,:)
  ! ppfdout - saved PPFD profiles (umol/m2-s)
  real(kind=dp), allocatable :: ppfdout(:,:)
  ! fjout - saved attentuated radiative flux profiles ()
  real(kind=dp), allocatable :: fjout(:,:)
  ! cairout - saved cair profiles (molecules/cm3)
  real(kind=dp), allocatable :: cairout(:,:)
  ! h2oout - saved h2o profiles (molecules/cm3)
  real(kind=dp), allocatable :: h2oout(:,:)

  ! Meteorological data

  ! Met data at reference height
  ! pmbzref - air pressure (mb)
  real(kind=dp)                  :: pmbzref
  ! umszref - wind speed (m/s)
  real(kind=dp)                  :: umszref
  ! ubzref - wind speed (cm/s)
  real(kind=dp)                  :: ubzref
  ! rhzref - relative humidity (%)
  real(kind=dp)                  :: rhzref
  ! qzref - molar water vapor concentration (moles/m3)
  real(kind=dp)                  :: qzref
  ! tkzref - air temperature (K)
  real(kind=dp)                  :: tkzref
  ! ppfdzref - PPFD (umol/m2-s)
  real(kind=dp)                  :: ppfdzref
  ! sradzref - solar radiation (W/m2)
  real(kind=dp)                  :: sradzref
  ! razref - aerodynamic resistance (s/cm)
  real(kind=dp)                  :: razref
  ! eatm - water vapor pressure (kPa)
  real(kind=dp)                  :: eatm

  ! tk - vertical profile of temperature at current simulation time (K)
  real(kind=dp), dimension(npts) :: tk
  ! tkn - vertical profile of temperature at next dtenv time step (K)
  real(kind=dp), dimension(npts) :: tkn
  ! tkp - vertical profile of temperature at previous dtenv time step (K)
  real(kind=dp), dimension(npts) :: tkp
  ! qh - vertical profile of specific humidity at current simulation time (g/kg)
  real(kind=dp), dimension(npts) :: qh 
  ! qhn - vertical profile of specific humidity at next dtenv time step (g/kg)
  real(kind=dp), dimension(npts) :: qhn
  ! qhp - vertical profile of specific humidity at previous dtenv time step (g/kg)
  real(kind=dp), dimension(npts) :: qhp
  ! pmb - vertical profile of pressure at current simulation time  (mb)
  real(kind=dp), dimension(npts) :: pmb
  ! pmbn - vertical profile of pressure at next dtenv time step (mb)
  real(kind=dp), dimension(npts) :: pmbn
  ! pmbp - vertical profile of pressure at previous dtenv time step (mb)
  real(kind=dp), dimension(npts) :: pmbp
  ! kv - eddy diffusivity at current simulation time (cm2/s)
  real(kind=dp), dimension(npts) :: kv
  ! kvn - eddy diffusivity at next dtenv time step (cm2/s)
  real(kind=dp), dimension(npts) :: kvn
  ! kvp - eddy diffusivity at previous dtenv time step (cm2/s)
  real(kind=dp), dimension(npts) :: kvp
  ! fj - vertical profile of actinic flux at current simulation time (??)
  ! TODO: initially this will only be the actinic flux normalized to the
  ! canopy top (i.e., dimensionless factor applied to Jphoto rates)
  real(kind=dp), dimension(npts) :: fj
  ! fjn - vertical profile of actinic flux at next dtenv time step (??)
  real(kind=dp), dimension(npts) :: fjn
  ! fjp - vertical profile of actinic flux at previous dtenv time step (??)
  real(kind=dp), dimension(npts) :: fjp
  ! ppfd - vertical profile of photosynthetic photon flux (umol/m2-s) at current simulation time
  real(kind=dp), dimension(npts) :: ppfd
  ! ppfdn - vertical profile of ppfd at next dtenv time step
  real(kind=dp), dimension(npts) :: ppfdn
  ! ppfdp - vertical profile of ppfd at previous dtenv time step
  real(kind=dp), dimension(npts) :: ppfdp
  ! ubar - vertical profile of mean wind speed at current simulation time (cm/s)
  real(kind=dp), dimension(npts) :: ubar
  ! ubarn - vertical profile of mean wind speed at next dtenv time step (cm/s)
  real(kind=dp), dimension(npts) :: ubarn
  ! ubarp - vertical profile of mean wind speed at previous dtenv time step (cm/s)
  real(kind=dp), dimension(npts) :: ubarp

  ! Calculated meteorological data
  ! cair - vertical air density at current simulation time (molec/cm3)
  real(kind=dp), dimension(npts) :: cair
  ! h2o - vertical water vapor concentration at current simulation time (molec/cm3)
  real(kind=dp), dimension(npts) :: h2o

  ! Constrained species data
  ! kcns - mixing rate for constrained species (1/s)
  real(kind=dp)    :: kcns
  ! ncns - number of constrained species
  integer(kind=i4) :: ncns
  ! cnssp - indices of constrained species
  integer(kind=i4), allocatable :: cnssp(:)
  ! cnsu - flag for ppbv to molec/cm3 conversion
  integer(kind=i4), allocatable :: cnsu(:)
  ! ccns - constrained species concentrations (molecules/cm3)
  real(kind=dp), allocatable    :: ccns(:,:)

  ! Emissions fluxes data
  ! tlk - vertical profile of leaf temperature at current simulation time (K)
  real(kind=dp), dimension(npts) :: tlk

  ! TODO:  Figure out how to make these dynamically set at runtime
  ! fevergreen - fraction of forest which is composed of evergreen trees
  real(kind=dp), parameter    :: fevergreen=0.1
  ! scale_bvoc - scaling factor for bvoc emissions
  real(kind=dp), parameter    :: scale_bvoc=0.2

  ! Entrainment exchange data
  ! naloft - number of species with non-zero concentrations aloft
  integer(kind=i4) :: naloft
  ! caloft - above ABL species concentrations (molecules/cm3)
  real(kind=dp), dimension(ninteg) :: caloft
  ! alfsp - indices of aloft species
  integer(kind=i4), allocatable :: alfsp(:)
  ! alfu - flag for ppbv to molec/cm3 conversion
  integer(kind=i4), allocatable :: alfu(:)

  ! Physical-Chemical data
  ! mdiffstp - molecular diffusivities of species in air at 0degC and 1 atm (cm2/s)
  real(kind=dp), dimension(ntotal) :: mdiffstp
  ! mdiffstp_default - default value of mdiffstp (cm2/s) for species with no reliable data
  real(kind=dp), parameter         :: mdiffstp_default=0.100
  ! hstar - effective Henry's Law coefficients (M/atm)
  real(kind=dp), dimension(ntotal) :: hstar
  ! hstar_default - default value of hstar (M/atm)
  real(kind=dp), parameter         :: hstar_default=1.000
  ! f0 - reactivity parameter used in resistance calculation (dimensionless)
  real(kind=dp), dimension(ntotal) :: f0
  ! f0_default - default value of f0
  real(kind=dp), parameter         :: f0_default=0.0
  ! molecmass - molecular mass of simulated species
  real(kind=dp), dimension(ntotal) :: molecmass

  ! simulation control file unit number
  integer(kind=i4), parameter :: UCTRL=11
  ! simulation control file name
  character(len=14), parameter :: filectrl='accessCTRL.dat'
  ! initial conditions file unit number
  integer(kind=i4), parameter :: UICS=12
  ! initial conditions file name
  character(len=16) :: icfile
  ! met data file unit number
  integer(kind=i4), parameter :: UENV=13
  ! met data file name
  character(len=16) :: envfile
  ! canopy morphology file unit number
  integer(kind=i4), parameter :: UCNPY=14
  ! canopy morphology data file name
  character(len=16) :: cnpyfile
  ! concentrations aloft file unit number
  integer(kind=i4), parameter :: UALFT=15
  ! concentrations aloft file name
  character(len=16) :: cafile
  ! soil concentrations file unit number
  integer(kind=i4), parameter :: USOIL=16
  ! soil concentrations file name
  character(len=16) :: csfile
  ! constrained species file unit number
  integer(kind=i4), parameter :: UCNS=17
  ! constrained species file name
  character(len=16) :: cnsfile
  ! grid definition file unit number
  integer(kind=i4), parameter :: UGRID=18
  ! grid definition file name
  character(len=16) :: grdfile

  ! simdescr - descriptive simulation identification
  character(len=100) :: simdescr
  ! simname - unique simulation name
  character(len=16) :: simname

  ! slat - latitude (deg)
  real(kind=dp) :: slat
  ! slon - longitude (deg)
  real(kind=dp) :: slon

  ! year - year of simulation
  integer(kind=i4) :: year
  ! month - month of simulation (at start)
  integer(kind=i4) :: month
  ! daymonth - day of the month (at start)
  integer(kind=i4) :: daymonth
  ! doy - day of year at beginning of simulation
  integer(kind=i4) :: doy
  ! hz, mz, sz - local time at beginning of simulation
  integer(kind=i4) :: hz, mz, sz
  ! tzone - time zone differential
  integer(kind=i4) :: tzone

  ! unit number for output CTRL simulation summary
  integer(kind=i4), parameter :: USUMM=20
  ! name of CTRL simulation summary file
  character(len=32) :: simsummary
  ! unit number for runtime output file
  integer(kind=i4), parameter :: URUN=21
  ! name of runtime output file
  character(len=32) :: simrunfile
  ! begining output file unit number for species files
  integer(kind=i4), parameter :: UOUT=30

  ! stdsp - indices of species printed to STDOUT
  integer(kind=i4) :: nstdsp
  integer(kind=i4), allocatable :: stdsp(:)

  ! options
  ! OPT_SIMTYPE - dynamic or steady-state simulation
  integer(kind=i4)            :: OPT_SIMTYPE
  integer(kind=i4), parameter :: DYNAMIC=1      ! solar time changes with t
  integer(kind=i4), parameter :: SSTATE=2       ! solar time fixed at tstart
  ! AQPHASE - option for a aqueous chemistry
  logical :: AQPHASE
  ! CHEMISTRY - option for chemistry
  logical :: CHEMISTRY
  ! CONSTRAIN - option for constrained species
  logical :: CONSTRAIN 
  ! VERTTRANS - option for vertical transport
  logical :: VERTTRANS
  ! EMISSION - option for emissions
  logical :: EMISSION
  ! DRYDEPOS - option for dry deposition
  logical :: DRYDEPOS
  ! METINTERP - option for interpolating met data to current time during integration
  logical :: METINTERP
  ! INTGWVAP - option for water vapor integration
  logical :: INTGWVAP
  ! INTGTAIR - option for air temperature integration
  logical :: INTGTAIR
  !
  ! KVSELECT - Kv parameterization selection
  integer(kind=i4)            :: KVSELECT
  integer(kind=i4), parameter :: KVRA=1        ! Aerodynamic resistance parameterization
  integer(kind=i4), parameter :: KVSTULL=2     ! Stull-based
  integer(kind=i4), parameter :: KVULKE=3      ! Ulke-based

  ! rtol, atol - ode relative and absolute tolerances
  real(kind=dp), parameter :: rtol=1.0D-02
  real(kind=dp), parameter :: atol=1.0D+01

  ! Avogadro's number
  real(kind=dp), parameter :: navo=6.022D+23

  ! pi
  real(kind=dp), parameter :: pi=4.0*atan(1.0)

end module GlobalData
!======================================================================================================================!
