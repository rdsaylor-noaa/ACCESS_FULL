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
!     Module:       DryDep                                                                                             !
!                                                                                                                      !
!     Description:  contains dry deposition velocities algorithms                                                      !
!                                                                                                                      !
!======================================================================================================================!
module DryDep
  use GlobalData
  use CanopyPhysics
  use PhysChemData
  use Utils
  implicit none

  private  vdl
  public   GetDryDepExCoeffs, GetSoilDryDepExCoeffs

contains

!**********************************************************************************************************************!
! subroutine GetDryDepExCoeffs - calculate leaf-scale dry deposition velocities
!                                and compensation points
!**********************************************************************************************************************!
subroutine GetDryDepExCoeffs()
  integer(kind=i4)  :: i, l
  real(kind=dp)     :: mdiffl     ! molecular diffusivity of species l in air (cm2/s)
  real(kind=dp)     :: relhumi    ! relative humidity (%)
  real(kind=dp)     :: hstarl     ! effective Henry's Law coefficient (M/atm)
  real(kind=dp)     :: f0l        ! Wesley's reactivity parameter (dimensionless)
  real(kind=dp)     :: srad       ! solar radiation at canopy top - W/m2
  character(len=10) :: rs_select  ! selection of stomatal resistance algorithm

  rs_select = 'gs_medlyn'         ! TODO: make this selection via input file

  srad = ppfd(int(hc)+1)/4.57_dp   ! canopy top PPFD converted from umol/m2-s to W/m2
  do l=1,ninteg
    do i=1,npts
      if(lai(i) .gt. 0.0) then
        ! within canopy

        ! molecular diffusivity (cm2/s)
        mdiffl=MolecDiff(l, tk(i), pmb(i))

        ! relative humidity at level i
        relhumi=RelativeHumidity(tk(i), pmb(i), qh(i))

        ! leaf boundary layer resistance (s/cm)
        rb(i,l)=rbl(mdiffl, ubar(i))

        ! leaf cuticular resistance (s/cm)
        hstarl=EffHenrysLawCoeff(l)
        rc(i,l)=rcl(hstarl, f0l)

        ! leaf mesophyll resistance (s/cm)
        f0l=ReactivityParam(l)
        rm(i,l)=rml(hstarl, f0l)

        ! leaf stomatal resistance (s/cm)
        select case (rs_select)
          case('zhang_df')
            ! Zhang et al. (2002, 2003) with hardcodes deciduous forest parameters
            rs(i,l)=rs_zhang_df(mdiffl,tk(i),pmb(i),ppfd(i),srad,relhumi)
          case('gs_medlyn')
            ! rs from Medlyn et al. (2011) gs calculated in CanopyPhysics
            rs(i,l)=(mdiffh2o(tk(i),pmb(i))/mdiffl)*rs_wgt(i)*0.01   ! convert rs_wgt from s/m to s/cm
          case default
            ! Zhang et al. (2002, 2003) with hardcodes deciduous forest parameters
            rs(i,l)=rs_zhang_df(mdiffl,tk(i),pmb(i),ppfd(i),srad,relhumi)
        end select

        vd(i,l)=vdl(rb(i,l), rm(i,l), rc(i,l), rs(i,l))   ! deposition velocity (cm/s)
        gp(i,l)=gpla(i,l)                                 ! compensation point concentration (molec/cm3)
      else
        ! out of the canopy
        rb(i,l)=0.0_dp
        rc(i,l)=0.0_dp
        rm(i,l)=0.0_dp
        rs(i,l)=0.0_dp
        vd(i,l)=0.0_dp
        gp(i,l)=0.0_dp
      end if
    end do
  end do

  return
end subroutine GetDryDepExCoeffs

!**********************************************************************************************************************!
! subroutine GetSoilDryDepExCoeffs - calculate dry deposition velocities for
!                                    deposition to the ground surface
!
! Uses formulation as derived from ...
!    Sakaguchi & Zeng (2009) JGR, D01107, doi: 10.1029/2008JD010834.
!    Schuepp (1977) BLM, 12, 171-186.
!    Philip (1957) J. of Met., 14, 354-366.
!
! With data from ...
!    Rawls et al. (1982) Trans. ASAE, 25, 1316-1320.
!    Clapp and Hornberger (1978) Water Resources Res., 14, 601-604.
!
!**********************************************************************************************************************!
subroutine GetSoilDryDepExCoeffs()
  integer(kind=i4)          :: l
  real(kind=dp)             :: mdiffl         ! molecular diffusivity (cm2/s)

  ! soil exchange coeffs
  do l=1,ninteg

    ! ground boundary layer resistance (Rbg, s/cm) calculated in CanopyPhysics
    ! and assumed to be invariant over species

    ! soil diffusion resistance (s/cm)
    mdiffl=MolecDiff(l, tsoilk, pmb(1))
    rsoill(l) = SoilResist(mdiffl)

    ! deposition velocity to ground surface (cm/s)
    vs(l)=1.0/(rbg+rsoill(l))       
  end do

  return
end subroutine GetSoilDryDepExCoeffs

!**********************************************************************************************************************!
! function vdl - calculate dry deposition velocity
!
! Uses formulation suggested by Wolfe (2012) personal communication
!
!**********************************************************************************************************************!
function vdl(rb, rm, rc, rs)
  real(kind=dp), intent(in) :: rb    ! leaf boundary layer resistance (s/cm)
  real(kind=dp), intent(in) :: rm    ! mesophyll resistance (s/cm)
  real(kind=dp), intent(in) :: rc    ! cuticle resistance (s/cm)
  real(kind=dp), intent(in) :: rs    ! stomatal resistance (s/cm)
  real(kind=dp)             :: vdl   ! deposition velocity (cm/s)
  real(kind=dp)             :: rnum, rden, rl
  rnum = rc*(rs+rm)
  rden = rc+2.0*(rs+rm)
  rl = rb + (rnum/rden)
  vdl = 1.0/rl
  return 
end function vdl

end module DryDep
!======================================================================================================================!
