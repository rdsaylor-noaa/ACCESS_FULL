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
!     Module:       Emissions                                                                                          !
!                                                                                                                      !
!     Description:  Contains emissions algorithms                                                                      !
!                   ACCESS 2.x.x algorithms adapted from MEGANv2.1                                                     !
!                                                                                                                      !
!     Guenther et al. (2012) The Model of Emissions and Gases and Aerosols from                                        !
!     Nature version 2.1 (MEGAN2.1): an extended and updated framework for modeling                                    !
!     biogenic emissions                                                                                               !
!                                                                                                                      !
!======================================================================================================================!
module Emissions
  use GlobalData
  implicit none

  private CalcRawEmissions, MapRaw2Mech, gam_t, gam_t_lif, gam_t_ldf, gam_l, gam_l_ldf, q_no_soil
  public GetEmissions

!**********************************************************************************************************************!
! Define data used only in this module
!**********************************************************************************************************************!

  ! number of explict species
  integer(kind=i4), parameter :: nbemis=28

  ! species indices
  integer(kind=i4), parameter :: bisop        = 1
  integer(kind=i4), parameter :: bmyrc        = 2
  integer(kind=i4), parameter :: bsabi        = 3
  integer(kind=i4), parameter :: blimo        = 4
  integer(kind=i4), parameter :: bcare        = 5
  integer(kind=i4), parameter :: bocim        = 6
  integer(kind=i4), parameter :: bapin        = 7
  integer(kind=i4), parameter :: bbpin        = 8
  integer(kind=i4), parameter :: bpcym        = 9
  integer(kind=i4), parameter :: bocym        = 10
  integer(kind=i4), parameter :: batrp        = 11
  integer(kind=i4), parameter :: bgtrp        = 12
  integer(kind=i4), parameter :: bbphe        = 13
  integer(kind=i4), parameter :: bafrn        = 14
  integer(kind=i4), parameter :: bbcro        = 15
  integer(kind=i4), parameter :: bmbo         = 16
  integer(kind=i4), parameter :: bmeoh        = 17
  integer(kind=i4), parameter :: bacet        = 18
  integer(kind=i4), parameter :: bethe        = 19
  integer(kind=i4), parameter :: bbute        = 20
  integer(kind=i4), parameter :: betha        = 21
  integer(kind=i4), parameter :: betoh        = 22
  integer(kind=i4), parameter :: baald        = 23
  integer(kind=i4), parameter :: bhcho        = 24
  integer(kind=i4), parameter :: bacta        = 25
  integer(kind=i4), parameter :: bform        = 26
  integer(kind=i4), parameter :: bno          = 27
  integer(kind=i4), parameter :: bco          = 28

  ! short species names
  character(len=6), parameter, dimension(nbemis) :: ssname = (/           &
                                                        'isopre', &
                                                        'myrcen', &
                                                        'sabine', &
                                                        'limone', &
                                                        'carene', &
                                                        'ocimen', &
                                                        'a-pine', &
                                                        'b-pine', &
                                                        'p-cyme', &
                                                        'o-cyme', &
                                                        'a-terp', &
                                                        'g-terp', &
                                                        'b-phel', &
                                                        'a-farn', &
                                                        'b-cary', &
                                                        'MBO   ', &
                                                        'methan', &
                                                        'aceton', &
                                                        'ethene', &
                                                        'butene', &
                                                        'ethane', &
                                                        'ethano', &
                                                        'acetal', &
                                                        'formal', &
                                                        'acetic', &
                                                        'formic', &
                                                        'NO    ', &
                                                        'CO    ' /)
  ! full species names
  character(len=16), parameter, dimension(nbemis) :: sename = (/           &
                                                        'isoprene        ', &
                                                        'myrcene         ', &
                                                        'sabinene        ', &
                                                        'limonene        ', &
                                                        '3-carene        ', &
                                                        'ocimene         ', &
                                                        'a-pinene        ', &
                                                        'b-pinene        ', &
                                                        'p-cymene        ', &
                                                        'o-cymene        ', &
                                                        'a-terpinene     ', &
                                                        'g-terpinene     ', &
                                                        'b-phellandrene  ', &
                                                        'a-farnesene     ', &
                                                        'b-caryophyllene ', &
                                                        'MBO             ', &
                                                        'methanol        ', &
                                                        'acetone         ', &
                                                        'ethene          ', &
                                                        'butene          ', &
                                                        'ethane          ', &
                                                        'ethanol         ', &
                                                        'acetaldehyde    ', &
                                                        'formaldehyde    ', &
                                                        'acetic acid     ', &
                                                        'formic acid     ', &
                                                        'nitric oxide    ', &
                                                        'carbon monoxide ' /)

  ! emissions of each explicit species prior to mapping to a specific 
  ! mechanism (molecules/cm3-s)
  real(kind=dp), dimension(npts,nbemis) :: qraw

  ! molecular weight (g/gmole) of each explicit species
  real(kind=dp), parameter, dimension(nbemis) :: mwbe = (/            &
                                                               68.12, &
                                                              136.20, &
                                                              136.20, &
                                                              136.20, &
                                                              136.20, &
                                                              136.20, &
                                                              136.20, &
                                                              136.20, &
                                                              134.20, &
                                                              134.20, &
                                                              136.20, &
                                                              136.20, &
                                                              136.20, &
                                                              204.35, &
                                                              204.35, &
                                                               86.13, &
                                                               32.04, &
                                                               58.08, &
                                                               28.05, &
                                                               56.11, &
                                                               30.07, &
                                                               46.07, &
                                                               44.05, &
                                                               30.03, &
                                                               60.05, &
                                                               46.03, &
                                                               30.01, &
                                                               28.01 /)

  ! emission factors for explict species
  integer(kind=i4), parameter :: EVERGREEN=1   ! temperate needleleaf evergreen
  integer(kind=i4), parameter :: DECIDUOUS=2   ! temperate broadleaf deciduous
  real(kind=dp), parameter, dimension(2,nbemis) :: efs = reshape( (/  &
                                                      600.0, 10000.0, &
                                                       70.0,    30.0, &
                                                       70.0,    50.0, &
                                                      100.0,    80.0, &
                                                      160.0,    30.0, &
                                                       70.0,   120.0, &
                                                      500.0,   400.0, &
                                                      300.0,   130.0, &
                                                        9.9,    10.3, &
                                                        3.1,     6.1, &
                                                        9.9,    10.3, &
                                                        9.9,    10.3, &
                                                       28.8,    10.3, &
                                                       40.0,    40.0, &
                                                       80.0,    40.0, &
                                                      200.0,    0.01, &
                                                      800.0,   800.0, &
                                                      240.0,   240.0, &
                                                      165.0,   165.0, &
                                                       30.8,    30.8, &
                                                        1.4,     1.4, &
                                                      200.0,   200.0, &
                                                      200.0,   200.0, &
                                                       40.0,    40.0, &
                                                       30.0,    30.0, &
                                                       30.0,    30.0, &
                                                        2.0,     2.0, &
                                                      600.0,   600.0 /), (/ 2, nbemis /) )
  
  ! various model parameters (adapted from MEGANv2.1)
  real(kind=dp), dimension(nbemis) :: betai, ldfi, ct1i, ceoi

  data betai/0.13,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.17, &
             0.17,0.13,0.08,0.10,0.10,0.10,0.10,0.13,0.13,0.13,0.13,0.13,0.08,0.08/ 
  data ldfi /1.00,0.60,0.60,0.20,0.20,0.80,0.60,0.20,0.40,0.40,0.40,0.40,0.40,0.50, &
             0.50,1.00,0.80,0.20,0.80,0.20,0.20,0.80,0.80,0.80,0.80,0.80,1.00,1.00/ 
  data ct1i /95000.,80000.,80000.,80000.,80000.,80000.,80000.,80000.,80000.,80000., &
             80000.,80000.,80000.,130000.,130000.,95000.,60000.,80000.,80000.,80000., &
             80000.,95000.,95000.,95000.,95000.,95000.,60000.,60000./
  data ceoi /2.00,1.83,1.83,1.83,1.83,1.83,1.83,1.83,1.83,1.83,1.83,1.83,1.83,2.37, &
             2.37,2.00,1.60,1.83,1.83,1.83,1.83,2.00,2.00,2.00,2.00,2.00,1.60,1.60/ 

contains

!**********************************************************************************************************************!
! subroutine GetEmissions - calculate emission rates for all emitted species
!                           q - molecules species/cm3-s
!**********************************************************************************************************************!
subroutine GetEmissions()
  integer(kind=i4) :: i, l

  ! zero emissions array
  do l=1,ninteg
    do i=1,npts
      q(i,l)=0.0_dp
    end do
  end do

  ! get explicit species emissions
  call CalcRawEmissions()

  ! map explicit emissions to a specific mechanism
  call MapRaw2Mech()

 800 format(a,8x,a9,1x,a9,3x,a,9x,a,9x,a,9x,a)
 900 format(a, 2f8.1, 4(2x,e10.3))
1000 format(1x,a,28a10)
1001 format(f5.1,28e10.3)
  return
end subroutine GetEmissions

!**********************************************************************************************************************!
! subroutine CalcRawEmissions - calculate emission rates for all explicit species
!**********************************************************************************************************************!
subroutine CalcRawEmissions()
  integer(kind=i4)                 :: iz, lraw
  real(kind=dp)                    :: efinit
  real(kind=dp)                    :: qbraw
  real(kind=dp)                    :: gamt, gaml
  real(kind=dp), dimension(nbemis) :: ef

  ! using the fraction evergreen parameter, calculate appropriate emission factors
  ! for each explict species
  do lraw=1,nbemis
    efinit = (fevergreen*efs(EVERGREEN,lraw) + &
             (1.0-fevergreen)*efs(DECIDUOUS,lraw))*0.278/mwbe(lraw)
    ef(lraw)=efinit*scale_bvoc
  end do 

  ! for each model level with vegetation, calculate explict emissions in
  ! molecules/cm3-s
  do lraw=1,nbemis
  do iz=1,npts
    if(lai(iz) > 0.0) then
      qbraw=ef(lraw)*lai(iz)*navo*1.D-13/dzhc
      qraw(iz,lraw)=qbraw*gam_t(lraw,tl_wgt(iz))*gam_l(lraw,ppfd_wgt(iz))
      gamt = gam_t(lraw,tl_wgt(iz))
      gaml = gam_l(lraw,ppfd_wgt(iz))
    else
      qraw(iz,lraw)=0.0
    end if
  end do
  end do

  ! soil NO (emitted into the lowest model layer)
  qraw(1,bno)=q_no_soil()

  return
end subroutine CalcRawEmissions

!**********************************************************************************************************************!
! function gam_t - emission activity factor for temperature
!**********************************************************************************************************************!
function gam_t(lraw,tleaf)
  integer(kind=i4), intent(in) :: lraw
  real(kind=dp)   , intent(in) :: tleaf
  real(kind=dp)                :: gam_t

  gam_t=(1.0-ldfi(lraw))*gam_t_lif(lraw,tleaf)+ldfi(lraw)*gam_t_ldf(lraw,tleaf)

end function gam_t

!**********************************************************************************************************************!
! function gam_t_lif - temperature emission activity factor - light
!                      independent fraction
!**********************************************************************************************************************!
function gam_t_lif(lraw,tleaf)
  integer(kind=i4), intent(in) :: lraw
  real(kind=dp)   , intent(in) :: tleaf
  real(kind=dp), parameter     :: ts=297.0_dp
  real(kind=dp)                :: gam_t_lif

  gam_t_lif = exp(betai(lraw)*(tleaf-ts))

end function gam_t_lif

!**********************************************************************************************************************!
! function gam_t_ldf - temperature emission activity factor - light
!                      dependent fraction
!**********************************************************************************************************************!
function gam_t_ldf(lraw,tleaf)
  integer(kind=i4), intent(in) :: lraw
  real(kind=dp)   , intent(in) :: tleaf
  real(kind=dp), parameter :: ct2=230000.0_dp
  real(kind=dp), parameter :: topt=313.0_dp
  real(kind=dp), parameter :: rgas=8.314_dp
  real(kind=dp)            :: ect1, ect2
  real(kind=dp)            :: gam_t_ldf

  ect1 = exp(ct1i(lraw)*(tleaf-topt)/(rgas*tleaf*topt))
  ect2 = exp(ct2*(tleaf-topt)/(rgas*tleaf*topt))

  gam_t_ldf = ceoi(lraw)*ct2*ect1/(ct2-ct1i(lraw)*(1.0-ect2))

end function gam_t_ldf

!**********************************************************************************************************************!
! function gam_l - emission activity factor for light
!**********************************************************************************************************************!
function gam_l(lraw,ppfdiz)
  integer(kind=i4), intent(in) :: lraw
  real(kind=dp)   , intent(in) :: ppfdiz
  real(kind=dp)                :: gam_l

  gam_l = (1.0-ldfi(lraw)) + ldfi(lraw)*gam_l_ldf(ppfdiz)

end function gam_l

!**********************************************************************************************************************!
! function gam_l_ldf - light emission activity factor - light dependent
!                      fraction
!**********************************************************************************************************************!
function gam_l_ldf(ppfdiz)
  real(kind=dp), intent(in) :: ppfdiz
  real(kind=dp), parameter  :: alpha=0.0027_dp
  real(kind=dp), parameter  :: cp=1.066_dp
  real(kind=dp)             :: gam_l_ldf

  gam_l_ldf = alpha*cp*ppfdiz/sqrt(1.0+alpha*alpha*ppfdiz*ppfdiz)

end function gam_l_ldf

!**********************************************************************************************************************!
! function q_no_soil - soil NO emissions
!**********************************************************************************************************************!
function q_no_soil()
  real(kind=dp) :: q_no_soil
  real(kind=dp) :: tsoil                      ! estimated soil temperature (deg C)
  real(kind=dp) :: e_no                       ! temperature corrected soil NO emissions 
                                              ! (nmol/m2-s) 
  real(kind=dp) :: be_no                      ! basal soil NO emissions (nmol/m2-s)

  be_no = 0.0_dp   ! no NO emissions from soil in this version

  tsoil=tsoilk-273.15
  if (tsoil < 30.0) then
    e_no=be_no*(tsoil/30.0)
  else
    e_no=be_no
  end if

  ! convert from nmol/m2-s to molec/cm3-s
  q_no_soil=e_no*navo/(1.0D+13*(z(2)-z(1)))

end function q_no_soil

!$INSERT MapRaw2Mech

end module Emissions
!======================================================================================================================!
