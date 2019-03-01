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
!     Module:       Integrate                                                                                          !
!                                                                                                                      !
!     Description:  Contains all top-level routines to integrate the model                                             !
!                                                                                                                      !
!======================================================================================================================!
module Integrate
  use GlobalData
  use EnvironData
  use Emissions
  use DryDep
  use Chem
  use GasChem
  use VertTransport
  use Utils
  implicit none

  private IntegrateEmissions, &
          IntegrateConstrained, ReadConstrainedData, &
          IntegrateVertExchange, &
          IntegrateChemistry, &
          IntegrateDryDeposition, &
          UpdateSolutionArrays, InterpMetData
  public IntegrateModel

contains

!**********************************************************************************************************************!
! subroutine IntegrateModel - handles overall model integration for one output
!                             time step
!**********************************************************************************************************************!
subroutine IntegrateModel(tin, tout)
  real(kind=dp), intent(in) :: tin, tout
  real(kind=dp) :: ts             ! current integration time [tin, tout]

  ts = tin

  zadeg = SolarZenithAngle(ts/(24.*3600.))
  zarad = zadeg*pi/180.

  ! read the next time slice of met data
  ! if interpolation is not enabled, met data at tin is used over the
  ! entire dtenv time step
  call ReadEnvironData()

  ! loop through integration sequence until tout
  ! dt is the integration time step
  ! Strang operator splitting is used
  do

    if (METINTERP .eqv. .true.) then
      call InterpMetData(ts, dthalf, tin, tout)        ! interpolate met data to ts+dthalf
    end if
      
    if (EMISSION  .eqv. .true.) then
      call IntegrateEmissions(ts, dt)                  ! emissions 
    end if

    if (CONSTRAIN .eqv. .true.) then
      call IntegrateConstrained(ts, dthalf)            ! constrained species
    end if

    if (VERTTRANS .eqv. .true.) then
      call IntegrateVertExchange(ts, dthalf)           ! vertical turbulent diffusion
    end if

    if (CHEMISTRY .eqv. .true.) then
      call IntegrateChemistry(ts, dt)                  ! chemistry of all phases
    end if

    if (CONSTRAIN .eqv. .true.) then
      call IntegrateConstrained(ts+dthalf, dthalf)     ! constrained species
    end if

    if (DRYDEPOS  .eqv. .true.) then
      call IntegrateDryDeposition(ts, dt)              ! dry deposition
    end if

    if (VERTTRANS .eqv. .true.) then
      call IntegrateVertExchange(ts+dthalf, dthalf)    ! vertical turbulent diffusion
    end if

    ts = ts+dt
    write(*,1001) (ts-tstart)/3600.          ! write elapsed simulation time to STDOUT
    write(URUN,1001) (ts-tstart)/3600.       ! write elapsed simulation time to URUN
    call flush(URUN)

    call UpdateSolutionArrays()

    if (ts >= tout) exit
 
  end do

  ! read the next time slice of constrained species data
  call ReadConstrainedData()

1001 format('t-t0 = ', f8.3, ' hrs')
  return
end subroutine IntegrateModel

!**********************************************************************************************************************!
! subroutine IntegrateEmissions - adds emissions of all emitted species over tdstep
!                                 emissions are allowed at all vertical grid points
!**********************************************************************************************************************!
subroutine IntegrateEmissions(ts, dtstep)
  real(kind=dp), intent(in) :: ts, dtstep
  real(kind=dp)    :: dzi, dz1, dzn
  integer(kind=i4) :: i, l
   
  ! obtains emissions, q (molecules/cm3-s)
  call GetEmissions()

  ! adds emissions for all integrated species
  do l=1,ninteg
    do i=1,npts
      ! integration over dtstep
      cint(i,l) = cint(i,l) + q(i,l)*dtstep

      ! accumulating budget emissions (ppbv)
      emtout(i,l,nt)=emtout(i,l,nt) + q(i,l)*dtstep*1.0D+9/cair(i)
    end do
  end do

  return
end subroutine IntegrateEmissions

!**********************************************************************************************************************!
! subroutine IntegrateConstrained - integrates constrained species term
!**********************************************************************************************************************!
subroutine IntegrateConstrained(ts, dtstep)
  real(kind=dp), intent(in) :: ts, dtstep
  integer(kind=i4) :: i, l

  ! integration for dtstep
  do l=1,ncns
    do i=1,npts
      ! integration over dtstep
      cint(i,cnssp(l)) = cint(i,cnssp(l)) - kcns*dtstep*(cint(i,cnssp(l))-ccns(i,l))
      cint(i,cnssp(l)) = max(cint(i,cnssp(l)), 0.0)
      
      ! accumulating budget constrained species process (ppbv)
      cntout(i,cnssp(l),nt)=cntout(i,cnssp(l),nt) &
                           - kcns*dtstep*(cint(i,cnssp(l))-ccns(i,l))*1.0D+09/cair(i)
    end do
  end do

  return
end subroutine IntegrateConstrained

!**********************************************************************************************************************!
! subroutine ReadConstrainedData - read next time slice of constrained species data
!**********************************************************************************************************************!
subroutine ReadConstrainedData()
  real(kind=dp)     :: lin_slope
  integer(kind=i4)  :: i, l
  character(len=10) :: cn_date
  character(len=8)  :: cn_time

  ! read the next time slice of constrained species data
  read(UCNS,*) cn_date, cn_time, (ccns(1,l),l=1,ncns) 

  ! convert the ppbv species (where cnsul)=1) to molec/cm3
  do l=1,ncns
    if (cnsu(l) .eq. 1) then
      ccns(1,l) = ConvertPPBToMolecCC(ccns(1,l),1)
    end if
  end do

  ! both ccns and caloft are now in molec/cm3
  do l=1,ncns

    ! calculate slope: caloft contains all species, ccns only contains constrained species
    lin_slope = (caloft(cnssp(l)) - ccns(1,l))/(z(npts)-zcnsref)

    ! below or at z=zcnsref, constant value
    ! above z=zcnsref, linearly interpolate to caloft value
    do i=2,npts
      if (z(i) .le. zcnsref) then
        ccns(i,l)=ccns(1,l)
      else
        ccns(i,l)=ccns(i-1,l)+lin_slope*(z(i)-z(i-1))
      end if
    end do

  end do

  return

end subroutine ReadConstrainedData

!**********************************************************************************************************************!
! subroutine IntegrateVertExchange - integrates the diffusion equation from tstart to 
!                                    tstart+dtstep; uses Crank-Nicolson time 
!                                    discretization as implemented in Patankar (1980)
!                                    where the diffusivity is a function of z and t
!**********************************************************************************************************************!
subroutine IntegrateVertExchange(tstart, dtstep)
  ! tstart  = simulation time at the begining of the timestep 
  ! dtstep  = timestep (seconds)
  ! 
  real(kind=dp), intent(in) :: tstart, dtstep
  integer(kind=i4) :: i, l
  real(kind=dp), dimension(npts) :: phi, phip

  ! calculate soil surface exchange velocities (vs)
  call GetSoilDryDepExCoeffs()

  ! integrate each species in turn
  do l=1,ninteg

    ! fill local species arrays for species l
    phi = 0.0_dp
    do i=1,npts
      phip(i) = cint(i,l)

      ! accumulating budget vertical transport process (ppbv)
      vttout(i,l,nt)=vttout(i,l,nt) - phip(i)*1.0D+9/cair(i)
    end do

    ! perform integration over dtstep
    call SubIntegrateVertTransportSplitBC(phi, phip, kv, caloft(l), csoil(l), vs(l), cair, dtstep)

    ! update integrated species solution array with result
    do i=1,npts
      cint(i,l) = phi(i)

      ! accumulating budget vertical transport process (ppbv)
      ! equivalent to adding increment (phi-phip) of integration
      vttout(i,l,nt)=vttout(i,l,nt) + phi(i)*1.0D+9/cair(i)
    end do

  end do   ! end species integration loop

  return
end subroutine IntegrateVertExchange

!**********************************************************************************************************************!
! subroutine IntegrateChemistry - loop over each vertical grid point and integrate
!                                  chemistry in each layer from ts to ts+dtstep
!**********************************************************************************************************************!
subroutine IntegrateChemistry(ts, dtstep)
  real(kind=dp), intent(in) :: ts, dtstep
  integer(kind=i4) :: ierror
  integer(kind=i4) :: i, l, j
  real(kind=dp), parameter         :: rmin = 0.0D0
! real(kind=dp), parameter         :: rmin = 1.0D-20
  real(kind=dp), dimension(ninteg) :: y
  real(kind=dp), dimension(nrxn)   :: k
  real(kind=dp), dimension(nrxn)   :: r
   
  ! loop over each vertical grid point
  do i=1,npts

    ! set current z value (in meters) and current vertical grid point
    zcurrent = z(i)*0.01
    zpt = i
   
    ! at grid point "i", map species concentrations to array "y"
    ! Units:  molecules/cm3
    do l=1,ninteg
      y(l)=cint(i,l)

      ! accumulating budget chemistry process (ppbv)
      chtout(i,l,nt)=chtout(i,l,nt) - y(l)*1.0D+9/cair(i)
    end do

    ! integrate chemistry at grid point "i" from "ts" to "ts+dtstep"
    call IntegChemOneLevel(ninteg, y, ts, ts+dtstep, ierror)

    ! get reaction rates for output later
    call GasRateCoeffs(ts, k)
    call RxnRates(ninteg, k, y, r)
    do j=1,nrxn
      ksout(i,j,nt)=k(j)                ! molecules-cm-s units
      rxnout(i,j,nt)=max(r(j), rmin)    ! molecules/cm3-s
    end do

    ! map "y" back to "cint" 
    do l=1,ninteg
      cint(i,l) = y(l)

      ! accumulating budget chemistry process (ppbv)
      ! equivalent to adding increment (yafter-ybefore)
      chtout(i,l,nt)=chtout(i,l,nt) + y(l)*1.0D+9/cair(i)
    end do

  end do   ! end of vertical grid loop

  return
end subroutine IntegrateChemistry

!**********************************************************************************************************************!
! subroutine IntegrateDryDeposition - for each species removes mass via dry deposition
!**********************************************************************************************************************!
subroutine IntegrateDryDeposition(ts, dtstep)
  real(kind=dp), intent(in) :: ts, dtstep
  integer(kind=i4) :: i, l
  real(kind=dp) :: cnew

  ! get deposition velocities
  call GetDryDepExCoeffs()

  ! loop over all species and levels, removing mass via dry deposition
  ! do not allow concentrations to go < 0 (can only remove what's there to start)
  do l=1,ninteg
    do i=2,npts-1
      cnew = cint(i,l) - (cint(i,l)-gp(i,l))*vd(i,l)*lai(i)*dtstep/dzhc
      cnew = max(cnew, 0.0)

      ! accumulating budget dry deposition process (ppbv)
      dptout(i,l,nt)=dptout(i,l,nt) + (cnew - cint(i,l))*1.0D+09/cair(i)
      cint(i,l) = cnew
    end do
  end do

  return
end subroutine IntegrateDryDeposition

!**********************************************************************************************************************!
! subroutine UpdateSolutionArrays - updates all storage arrays after successful t-step
!**********************************************************************************************************************!
subroutine UpdateSolutionArrays()
  integer(kind=i4) :: i, l

  do l=1,ninteg
    do i=1,npts
      ctot(i,l) = cint(i,l)
    end do
  end do

  return
end subroutine UpdateSolutionArrays

!**********************************************************************************************************************!
! subroutine InterpMetData - interpolates met data variables between tin and tout
!                            at the dthalf point of the current integration time
!                            step, i.e., at ts+dthalf
!**********************************************************************************************************************!
subroutine InterpMetData(tm0, dtmhalf, tin, tout)
  integer(kind=i4) :: n
  real(kind=dp) :: tm0, dtmhalf, tin, tout
  real(kind=dp) :: alpha, a, a1

  ! interpolate met data to tm0+dtmhalf
  alpha = (tout - (tm0+dtmhalf))/(tout - tin) 
  a = alpha
  a1 = 1.0 - alpha

  do n=1,npts
    tk(n) = a1*tkn(n) + a*tkp(n)
    qh(n) = a1*qhn(n) + a*qhp(n)
    pmb(n) = a1*pmbn(n) + a*pmbp(n)
    kv(n) = a1*kvn(n) + a*kvp(n)
    fj(n) = a1*fjn(n) + a*fjp(n)
    ppfd(n) = a1*ppfdn(n) + a*ppfdp(n)
    ubar(n) = a1*ubarn(n) + a*ubarp(n)
    ! calculate derived met data variables
    cair(n) = pmb(n)*7.2428D+18/tk(n)   ! molec/cm3 
    h2o(n)  = 0.001611*cair(n)          ! molec/cm3
  end do

  return
end subroutine InterpMetData

end module Integrate
!======================================================================================================================!
