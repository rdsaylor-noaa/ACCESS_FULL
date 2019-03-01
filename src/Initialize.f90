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
!     Module:       Initialize                                                                                         !
!                                                                                                                      !
!     Description:  Contains routines related to model initialization                                                  !
!                                                                                                                      !
!======================================================================================================================!
module Initialize
  use GlobalData
  use EnvironData
  use PhysChemData
  use DryDep
  use Utils
  use Output
  implicit none

  private ReadVerticalGridData, SetInitialConditions, SetSimulationData, GetCanopyData, GetConstrainedData, &
          GetAloftData, GetSoilData
  public InitializeModel

contains

!**********************************************************************************************************************!
! InitializeModel - performs all the model initialization steps 
!                   Called from Main.f90
!**********************************************************************************************************************!
subroutine InitializeModel()
  real(kind=dp)    :: delta, leap, mjd, hrlocal, hrutc, time21
  integer(kind=i4) :: m
  integer(kind=i4), dimension(12) :: months

  data months /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

  ! read CTRL file
  call SetSimulationData()

  ! calculate doy (day-of-year)
  doy = 0.0
  do m=1,month-1
    doy = doy + dble(months(m))
  end do
  doy = doy + dble(daymonth)

  ! calculate decimal UTC hour
  hrlocal = dble(hz) + (dble(mz)/60.) + (dble(sz)/3600.)
  hrutc = hrlocal - tzone

  ! years from 1949
  delta = dble(year)-1949.

  ! leap days over this period (note correction below)
  leap = aint(delta/4.)

  ! modified Julian Day
  mjd=32916.5+delta*365.+leap+dble(doy)+hrutc/24.

  ! the last year of century is not a leap year unless divisible by 400
  if (dmod(dble(year),100.0D+00) .eq. 0.0 .and. dmod(dble(year),400.0D+00) .ne. 0.0) then
    mjd=mjd-1.
  end if

  ! time in days since Jan 1, 2000 Greenwich Noon
  time21 = mjd-51545.0 
 
  ! time in seconds since Jan 1, 2000 Greenwich Noon
  ! simulation time is kept in seconds
  tstart=time21*24.*3600.

  zadeg = SolarZenithAngle(time21)
  zarad = zadeg*pi/180.

  ! ending time
  tend=tstart+dtenv*dble(ntenv)
  ntout = int((tend-tstart)/dtout)
  print *, 'ntout=',ntout

  dthalf = 0.5_dp*dt                   ! half time step
  t=tstart                             ! simulation time start
  nt=0                                 ! output number start
  nte=0                                ! met data number start

  ! allocate storage for final output arrays
  allocate(cout(npts,ntotal,0:ntout))
  allocate(vdout(npts,ninteg,0:ntout))
  allocate(vsout(ninteg,0:ntout))
  allocate(qout(npts,ninteg,0:ntout))
  allocate(vfout(npts,ninteg,0:ntout))
  allocate(emtout(npts,ninteg,0:ntout))
  allocate(cntout(npts,ninteg,0:ntout))
  allocate(dptout(npts,ninteg,0:ntout))
  allocate(vttout(npts,ninteg,0:ntout))
  allocate(chtout(npts,ninteg,0:ntout))
  allocate(ksout(npts,nrxn,0:ntout))
  allocate(rxnout(npts,nrxn,0:ntout))
  allocate(timeout(0:ntout))

  allocate(rbout(npts,ninteg,0:ntout))
  allocate(rmout(npts,ninteg,0:ntout))
  allocate(rcout(npts,ninteg,0:ntout))
  allocate(rsout(npts,ninteg,0:ntout))
  allocate(rsoillout(ninteg,0:ntout))

  allocate(sdtout(0:ntenv))
  allocate(ra(0:ntenv))
  allocate(rakv(0:ntenv))
  allocate(zolout(0:ntenv))
  allocate(gaero(0:ntenv))
  allocate(ribout(0:ntenv))
  allocate(ppfddirout(0:ntenv))
  allocate(ppfddifout(0:ntenv))
  allocate(nirdirout(0:ntenv))
  allocate(nirdifout(0:ntenv))

  allocate(tkout(npts,0:ntenv))
  allocate(pmbout(npts,0:ntenv))
  allocate(qhout(npts,0:ntenv))
  allocate(ubarout(npts,0:ntenv))
  allocate(kvout(npts,0:ntenv))
  allocate(ppfdout(npts,0:ntenv))
  allocate(fjout(npts,0:ntenv))
  allocate(cairout(npts,0:ntenv))
  allocate(h2oout(npts,0:ntenv))
  allocate(rhout(npts,0:ntenv))

  allocate(ppfdsunout(npts,0:ntenv))
  allocate(ppfdshdout(npts,0:ntenv))
  allocate(ppfdwgtout(npts,0:ntenv))
  allocate(nirsunout(npts,0:ntenv))
  allocate(nirshdout(npts,0:ntenv))
  allocate(nirwgtout(npts,0:ntenv))
  allocate(lwupout(npts,0:ntenv))
  allocate(lwdnout(npts,0:ntenv))
  allocate(rtsunout(npts,0:ntenv))
  allocate(rtshdout(npts,0:ntenv))
  allocate(rtwgtout(npts,0:ntenv))
  allocate(rabssunout(npts,0:ntenv))
  allocate(rabsshdout(npts,0:ntenv))
  allocate(rabswgtout(npts,0:ntenv))
  allocate(fsunout(npts,0:ntenv))
  allocate(fshdout(npts,0:ntenv))
  allocate(tlsunout(npts,0:ntenv))
  allocate(tlshdout(npts,0:ntenv))
  allocate(tlwgtout(npts,0:ntenv))
  allocate(gssunout(npts,0:ntenv))
  allocate(gsshdout(npts,0:ntenv))
  allocate(gswgtout(npts,0:ntenv))
  allocate(rssunout(npts,0:ntenv))
  allocate(rsshdout(npts,0:ntenv))
  allocate(rswgtout(npts,0:ntenv))
  allocate(anetsunout(npts,0:ntenv))
  allocate(anetshdout(npts,0:ntenv))
  allocate(anetwgtout(npts,0:ntenv))

  allocate(vsh2oout(0:ntenv))
  allocate(qsoilout(0:ntenv))
  allocate(effrhsoilout(0:ntenv))
  allocate(rbgout(0:ntenv))
  allocate(gbgout(0:ntenv))
  allocate(rsoilout(0:ntenv))
  allocate(tsoilkout(0:ntenv))
  allocate(tk0out(0:ntenv))

  ! make sure tendency arrays are zeroed initially
  emtout=0.0_dp
  cntout=0.0_dp
  dptout=0.0_dp
  vttout=0.0_dp
  chtout=0.0_dp
  ksout=0.0_dp
  rxnout=0.0_dp
  ra=0.0_dp
  rakv=0.0_dp
  tkout=0.0_dp
  pmbout=0.0_dp
  qhout=0.0_dp
  ubarout=0.0_dp
  kvout=0.0_dp
  ppfdout=0.0_dp
  fjout=0.0_dp
  cairout=0.0_dp
  h2oout=0.0_dp
  rhout=0.0_dp

  rbout=0.0_dp
  rmout=0.0_dp
  rcout=0.0_dp
  rsout=0.0_dp
  rsoillout=0.0_dp

  ppfddirout=0.0_dp
  ppfddifout=0.0_dp
  nirdirout=0.0_dp
  nirdifout=0.0_dp

  ppfdsunout=0.0_dp
  ppfdshdout=0.0_dp
  ppfdwgtout=0.0_dp
  nirsunout=0.0_dp
  nirshdout=0.0_dp
  nirwgtout=0.0_dp
  lwupout=0.0_dp
  lwdnout=0.0_dp
  rtsunout=0.0_dp
  rtshdout=0.0_dp
  rtwgtout=0.0_dp
  rabssunout=0.0_dp
  rabsshdout=0.0_dp
  rabswgtout=0.0_dp
  fsunout=0.0_dp
  fshdout=0.0_dp
  tlsunout=0.0_dp
  tlshdout=0.0_dp
  tlwgtout=0.0_dp
  gssunout=0.0_dp
  gsshdout=0.0_dp
  gswgtout=0.0_dp
  rssunout=0.0_dp
  rsshdout=0.0_dp
  rswgtout=0.0_dp
  anetsunout=0.0_dp
  anetshdout=0.0_dp
  anetwgtout=0.0_dp
  vsh2oout=0.0_dp
  qsoilout=0.0_dp
  effrhsoilout=0.0_dp
  rbgout=0.0_dp
  gbgout=0.0_dp
  rsoilout=0.0_dp
  tsoilkout=0.0_dp
  tk0out=0.0_dp

  ! set all physical-chemical data
  call SetPhysChemData()

  ! read z values of vertical grid
  call ReadVerticalGridData()

  ! read canopy morphology data
  call GetCanopyData()

  ! read soil compensation concs file
  call GetSoilData()

  ! read first time slice of environmental data file
  call ReadEnvironData()

  ! read initial conditions file
  call SetInitialConditions()

  ! get concentrations aloft
  call GetAloftData()

  ! read constrained concentration file
  call GetConstrainedData()

  return
end subroutine InitializeModel

!**********************************************************************************************************************!
! ReadVerticalGridData - read vertical grid definition from an input file
!                  
!**********************************************************************************************************************!
subroutine ReadVerticalGridData()
  integer(kind=i4) :: i, k, nptsdef
  real(kind=dp)    :: zm, dzm, dzhcm
 
  ! read GridDef file
  open(unit=UGRID,file=('./data/' // grdfile))
  read(UGRID,*) nptsdef, hc, ncnpy, alfa
  
  ! check file consistency
  if (nptsdef .ne. npts) then
    write(*,101) grdfile
    write(*,102) nptsdef, npts
    close(UGRID)
    stop
  end if

  ! ok, so read data
  do i=1,npts
    read(UGRID,*) k, zm, dzm
    z(k) = zm*100.0
    if (k == 2) dzhcm=dzm
  end do

  dzhc = dzhcm*100.0

  ! domain boundaries
  z0 = z(1)
  zi = z(npts)

  ! canopy height in cm
  hccm = hc*100.0

  ! set aerodynamic parameters
  !  From the recommendations of ...
  !  Monteith and Unsworth (2013) Principles of Environmental Physics
  d   = 0.667*hccm
  z0m = 0.097*hccm
  z0h = z0m/7.3

  close(UGRID) 

101 format('***GridDef file ', a, ' is inconsistent with current ACCESS configuration!')
102 format('***GridDef npts = ', i4 /'***GlobalData npts = ', i4)
  return
end subroutine ReadVerticalGridData

!**********************************************************************************************************************!
! SetSimulationData - reads input data from CTRL file for simulation
!                     and writes to STDOUT and the simulation summary file
!**********************************************************************************************************************!
subroutine SetSimulationData()
  integer(kind=i4)  :: j, isp
  character(len=35) :: outdir
  character(len=35) :: strmkdir
  logical           :: eoutdir
  integer(kind=i4), parameter   :: nhlines=7

  ! read CTRL file
  open(unit=UCTRL,file=('./ctrl/' // filectrl))
  ! skip header lines
  do j=1,nhlines
    read(UCTRL,*)
  end do
  read(UCTRL,101) simdescr
  read(UCTRL,*) 
  read(UCTRL,102) simname 
  read(UCTRL,*) OPT_SIMTYPE
  read(UCTRL,*) slat
  read(UCTRL,*) slon
  read(UCTRL,*) year
  read(UCTRL,*) month
  read(UCTRL,*) daymonth
  read(UCTRL,*) hz, mz, sz
  read(UCTRL,*) tzone
  read(UCTRL,*) dt
  read(UCTRL,*) dtenv
  read(UCTRL,*) ntenv
  read(UCTRL,*) dtout
  read(UCTRL,*) nstdsp
  allocate(stdsp(nstdsp))
  read(UCTRL,*) (stdsp(isp),isp=1,nstdsp)
  read(UCTRL,*) CHEMISTRY
  read(UCTRL,*) AQPHASE
  read(UCTRL,*) CONSTRAIN
  read(UCTRL,*) VERTTRANS
  read(UCTRL,*) EMISSION
  read(UCTRL,*) DRYDEPOS
  read(UCTRL,*) METINTERP
  read(UCTRL,*) INTGWVAP
  read(UCTRL,*) INTGTAIR
  read(UCTRL,*)
  read(UCTRL,*) grdfile
  read(UCTRL,*)
  read(UCTRL,*) icfile
  read(UCTRL,*)
  read(UCTRL,*) envfile
  read(UCTRL,*)
  read(UCTRL,*) cnpyfile
  read(UCTRL,*)
  read(UCTRL,*) cafile
  read(UCTRL,*)
  read(UCTRL,*) csfile
  if (CONSTRAIN .eqv. .TRUE.) then
    read(UCTRL,*)
    read(UCTRL,*) cnsfile
  end if
  close(UCTRL)

100 format(a6)
101 format(a100)
102 format(a16)

  ! write CTRL data to STDOUT
  write(6,300)
  write(6,*) simdescr
  write(6,2010) simname 
  if (OPT_SIMTYPE .eq. DYNAMIC) then
    write(6,2015)
  else
    write(6,2016)
  end if
  write(6,2020) slat
  write(6,2030) slon
  write(6,2031) year
  write(6,2040) month
  write(6,2041) daymonth
  write(6,2050) hz, mz, sz
  write(6,2051) tzone
  write(6,2060) dt
  write(6,2070) dtenv
  write(6,2071) ntenv
  write(6,2072) dtout
  write(6,2080) nstdsp
  write(6,2090) (trim(sspc(stdsp(isp))),isp=1,nstdsp)
  write(6,2200) 
  write(6,2110) CHEMISTRY
  write(6,2120) AQPHASE
  write(6,2140) CONSTRAIN
  write(6,2150) VERTTRANS
  write(6,2160) EMISSION
  write(6,2170) DRYDEPOS
  write(6,2171) METINTERP
  write(6,2172) INTGWVAP
  write(6,2173) INTGTAIR
  write(6,2210) 
  write(6,2179) grdfile
  write(6,2180) icfile
  write(6,2190) envfile
  write(6,2300) cnpyfile
  write(6,2302) cafile
  write(6,2303) csfile
  if (CONSTRAIN .eqv. .TRUE.) then
    write(6,2310) cnsfile
  end if
  write(6,300)

  ! if output directory in 'out' does not exist, create it
  ! and all necessary subdirectories
  outdir='./out/' // trim(simname) // '/.'
  inquire(file=outdir, exist=eoutdir)
  if(eoutdir .eqv. .FALSE.) then
    strmkdir = 'mkdir ./out/' // trim(simname)
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/grid'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/sp'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/emis'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/vd'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/vflux'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/budget'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/ks'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/rates'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/r'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/met'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/canopy'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/soil'
    call system(strmkdir)
  end if

  ! open simulation runtime file to capture STDOUT
  ! runtime file is normally closed by Utils:CleanUp
  ! at the end of the simulation
  simrunfile='./out/' // trim(simname) // '/ACCESS_STD.out'
  simrunfile=trim(simrunfile)
  open(unit=URUN,file=simrunfile)

  ! write header for runtime file
  write(URUN,300)
  write(URUN,*) simdescr
  write(URUN,2010) simname 
  write(URUN,300)

  ! write CTRL data to USUMM 
  simsummary='./out/' // trim(simname) // '/ACCESS_SUMM.out'
  simsummary=trim(simsummary)
  open(unit=USUMM,file=simsummary)
  write(USUMM,300)
  write(USUMM,*) simdescr
  write(USUMM,2010) simname 
  if (OPT_SIMTYPE .eq. DYNAMIC) then
    write(USUMM,2015)
  else
    write(USUMM,2016)
  end if
  write(USUMM,2020) slat
  write(USUMM,2030) slon
  write(USUMM,2031) year
  write(USUMM,2040) month
  write(USUMM,2041) daymonth
  write(USUMM,2050) hz, mz, sz
  write(USUMM,2051) tzone
  write(USUMM,2060) dt
  write(USUMM,2070) dtenv
  write(USUMM,2071) ntenv
  write(USUMM,2072) dtout
  write(USUMM,2080) nstdsp
  write(USUMM,2090) (trim(sspc(stdsp(isp))),isp=1,nstdsp)
  write(USUMM,2200) 
  write(USUMM,2110) CHEMISTRY
  write(USUMM,2120) AQPHASE
  write(USUMM,2140) CONSTRAIN
  write(USUMM,2150) VERTTRANS
  write(USUMM,2160) EMISSION
  write(USUMM,2170) DRYDEPOS
  write(USUMM,2171) METINTERP
  write(USUMM,2172) INTGWVAP
  write(USUMM,2173) INTGTAIR
  write(USUMM,2210) 
  write(USUMM,2179) grdfile
  write(USUMM,2180) icfile
  write(USUMM,2190) envfile
  write(USUMM,2300) cnpyfile
  write(USUMM,2302) cafile
  write(USUMM,2303) csfile
  if (CONSTRAIN .eqv. .TRUE.) then
    write(USUMM,2310) cnsfile
  end if
  write(USUMM,300)
  close(USUMM)

2000 format(/'Starting ACCESS v', a, ' simulation...'/)
2010 format(' Short sim name = ', a)
2015 format(' OPT_SIMTYPE    = DYNAMIC')
2016 format(' OPT_SIMTYPE    = SSTATE')
2020 format(' Latitude       = ', f7.2, ' deg')
2030 format(' Longitude      = ', f7.2, ' deg')
2031 format(' Year           = ', i4)
2040 format(' Month           = ', i4)
2041 format(' Day             = ', i4)
2050 format(' Start time     = ', i2.2, ':', i2.2, ':', i2.2, ' LT')
2051 format(' Time zone diff = ', i3)
2060 format(' Integration dt = ', f6.1, ' s')
2070 format(' Met data dt    = ', f7.1, ' s')
2071 format(' # metdata stps = ', i4)
2072 format(' Output dt      = ', f7.1, ' s')
2080 format(/' Species to STDOUT = ', i4)
2090 format(100(1x,a))
2110 format(' CHEMISTRY      = ', l2)
2120 format(' AQPHASE        = ', l2)
2140 format(' CONSTRAIN      = ', l2)
2150 format(' VERTTRANS      = ', l2)
2160 format(' EMISSION       = ', l2)
2170 format(' DRYDEPOS       = ', l2)
2171 format(' METINTERP      = ', l2)
2172 format(' INTGWVAP       = ', l2)
2173 format(' INTGTAIR       = ', l2)
2179 format(' GRD file name  = ', a)
2180 format(' IC file name   = ', a)
2190 format(' ENVDATA file name  = ', a)
2200 format(/' Model Options:')
2210 format(/' Input Files:')
2300 format(' CNPY file name  = ', a)
2302 format(' CALFT file name = ', a)
2303 format(' CSOIL file name = ', a)
2310 format(' CNS file name   = ', a)
300 format(80('='))
301 format(' ACCESS v',a)
8000 format('  !!!!ERROR!!!!!!'/'  ACCESS version ',6a /'  CTRL version   ',6a &
           /'  Stopping ...')  

  return
end subroutine SetSimulationData

!**********************************************************************************************************************!
! GetSoilData - reads soil data and soil compensation point concentrations
!**********************************************************************************************************************!
subroutine GetSoilData()
  real(kind=dp)     :: s_z0, s_zi
  integer(kind=i4)  :: i, l, j
  integer(kind=i4)  :: s_npts
  character(len=16) :: s_gasmech
  integer(kind=i4), parameter  :: nhlines=8
  integer(kind=i4)  :: ssp, slu
  real(kind=dp)     :: cs

  ! read soil data file
  open(unit=USOIL,file=('./data/' // csfile))
  ! skip header lines
  do j=1,nhlines
    read(USOIL,*)
  end do
  ! read consistency data and check
  read(USOIL,101) s_gasmech
  s_gasmech=trim(s_gasmech)
  read(USOIL,*)   s_npts
  read(USOIL,*)   s_z0
  read(USOIL,*)   s_zi 
  if (s_gasmech .ne. gasmech) then
    write(*,201) csfile
    write(*,202) s_gasmech, gasmech
    close(USOIL)
    stop
  else if (s_npts .ne. npts) then
    write(*,201) csfile
    write(*,203) s_npts, npts
    close(USOIL)
    stop
  else if (s_z0*100. .ne. z0) then
    write(*,201) csfile
    write(*,204) s_z0, z0
    close(USOIL)
    stop
  else if (s_zi*100. .ne. zi) then
    write(*,201) csfile
    write(*,205) s_zi, zi
    close(USOIL)
    stop
  end if

  ! made it this far, so read the data

  ! read soil type
  do j=1,14
    read(USOIL,*)
  end do
  read(USOIL,*) isoiltype

  ! read volumetric soil water content
  read(USOIL,*)
  read(USOIL,*)
  read(USOIL,*)
  read(USOIL,*) stheta

  ! read depth of topsoil
  read(USOIL,*)
  read(USOIL,*)
  read(USOIL,*)
  read(USOIL,*) dsoil

  ! based on soil type, set soil data parameters
  sattheta = xsattheta(isoiltype)
  rtheta   = xrtheta(isoiltype)
  sbcoef   = xsbcoef(isoiltype)
  satphi   = xsatphi(isoiltype)

  ! read compensation points
  read(USOIL,*)
  read(USOIL,*)
  read(USOIL,*) nsoil
  read(USOIL,*)
  read(USOIL,*)

  do l=1,ninteg
    csoil(l)=0.0_dp
  end do

  do l=1,nsoil
    read(USOIL,*) ssp, slu, cs
    if (slu .eq. 1) then
      csoil(ssp) = ConvertPPBToMolecCC(cs,npts)
    else
      csoil(ssp) = cs
    end if
  end do
  close(USOIL)

101 format(a16)
201 format('***Soil concentrations file ', a, ' is inconsistent with current ACCESS configuration!')
202 format('***SOIL gasmech = ', a /'***ACCESS gasmech = ', a)
203 format('***SOIL npts = ', i4 /'***ACCESS npts = ', i4)
204 format('***SOIL z0 = ', e12.4 /'***ACCESS z0 = ', e12.4)
205 format('***SOIL zi = ', e12.4 /'***ACCESS zi = ', e12.4)
  return
end subroutine GetSoilData

!*************************************************************************************!
! GetAloftData - reads concentrations aloft file
!*************************************************************************************!
subroutine GetAloftData()
  real(kind=dp)     :: al_z0, al_zi
  integer(kind=i4)  :: i, l, j
  integer(kind=i4)  :: al_npts
  character(len=16) :: al_gasmech
  integer(kind=i4), parameter  :: nhlines=8
  integer(kind=i4)  :: alsp, alu
  real(kind=dp)     :: ca

  ! read concentrations aloft file
  open(unit=UALFT,file=('./data/' // cafile))
  ! skip header lines
  do j=1,nhlines
    read(UALFT,*)
  end do
  ! read consistency data and check
  read(UALFT,101) al_gasmech
  al_gasmech=trim(al_gasmech)
  read(UALFT,*)   al_npts
  read(UALFT,*)   al_z0
  read(UALFT,*)   al_zi 
  if (al_gasmech .ne. gasmech) then
    write(*,201) cafile
    write(*,202) al_gasmech, gasmech
    close(UALFT)
    stop
  else if (al_npts .ne. npts) then
    write(*,201) cafile
    write(*,203) al_npts, npts
    close(UALFT)
    stop
  else if (al_z0 .ne. z0) then
    write(*,201) cafile
    write(*,204) al_z0, z0
    close(UALFT)
    stop
  else if (al_zi .ne. zi) then
    write(*,201) cafile
    write(*,205) al_zi, zi
    close(UALFT)
    stop
  end if

  ! made it this far, so read the data
  read(UALFT,*)
  read(UALFT,*)
  read(UALFT,*) naloft
  read(UALFT,*)
  read(UALFT,*)

  do l=1,ninteg
    caloft(l)=0.0_dp
  end do

  allocate(alfsp(naloft))
  allocate(alfu(naloft))
  do l=1,naloft
    read(UALFT,*) alsp, alu, ca
    alfsp(l) = alsp
    alfu(l)  = alu
    if (alu .eq. 1) then
      caloft(alsp) = ConvertPPBToMolecCC(ca,npts)
    else
      caloft(alsp) = ca
    end if
  end do
  close(UALFT)

101 format(a16)
201 format('***Aloft concentrations file ', a, ' is inconsistent with current ACCESS configuration!')
202 format('***ALFT gasmech = ', a /'***ACCESS gasmech = ', a)
203 format('***ALFT npts = ', i4 /'***ACCESS npts = ', i4)
204 format('***ALFT z0 = ', e12.4 /'***ACCESS z0 = ', e12.4)
205 format('***ALFT zi = ', e12.4 /'***ACCESS zi = ', e12.4)
  return
end subroutine GetAloftData

!**********************************************************************************************************************!
! GetConstrainedData - reads constrained species file
!**********************************************************************************************************************!
subroutine GetConstrainedData()
  real(kind=dp)     :: cn_z0, cn_zi, lin_slope
  integer(kind=i4)  :: i, l, j
  integer(kind=i4)  :: cn_npts
  character(len=16) :: cn_gasmech
  character(len=10) :: cn_date
  character(len=8)  :: cn_time
  integer(kind=i4), parameter  :: nhlines=8

  ! read constrained species file
  open(unit=UCNS,file=('./data/' // cnsfile))
  ! skip header lines
  do j=1,nhlines
    read(UCNS,*)
  end do
  ! read consistency data and check
  read(UCNS,101) cn_gasmech
  cn_gasmech=trim(cn_gasmech)
  if (cn_gasmech .ne. gasmech) then
    write(*,201) cnsfile
    write(*,202) cn_gasmech, gasmech
    close(UCNS)
    stop
  end if

  ! made it this far, so read the data
  read(UCNS,*)
  read(UCNS,*) zcnsref
  read(UCNS,*)
  read(UCNS,*) kcns
  read(UCNS,*)
  read(UCNS,*) ncns
  read(UCNS,*)
  read(UCNS,*)
  read(UCNS,*)

  zcnsref=zcnsref*100.   ! convert from m to cm

  allocate(cnssp(ncns))
  allocate(cnsu(ncns))
  do l=1,ncns
    read(UCNS,*) cnssp(l), cnsu(l)
  end do
  read(UCNS,*)
  read(UCNS,*)
  ! now read the concentrations of each species
  allocate(ccns(npts,ncns))

  ! read concentrations for initial time
  read(UCNS,*) cn_date, cn_time, (ccns(1,l),l=1,ncns) 

  ! convert the ppbv species (where cnsu(l)=1) to molec/cm3
  print *, 'ccns(1,iO3), cair(1) = ', ccns(1,1), cair(1)
  do l=1,ncns
    if (cnsu(l) .eq. 1) then
      ccns(1,l) = ConvertPPBToMolecCC(ccns(1,l),1)
    end if
  end do

  print *, 'ccns(1,iO3) = ', ccns(1,1), cair(1)

  ! both ccns and caloft are now in molec/cm3
  do l=1,ncns

    ! calculate slope: caloft contains all species, ccns only contains
    !  constrained species
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

101 format(a16)
201 format('***Constrained species concentrations file ', a, ' is inconsistent with current ACCESS configuration!')
202 format('***CNS gasmech = ', a /'***ACCESS gasmech = ', a)
203 format('***CNS npts = ', i4 /'***ACCESS npts = ', i4)
204 format('***CNS z0 = ', e12.4 /'***ACCESS z0 = ', e12.4)
205 format('***CNS zi = ', e12.4 /'***ACCESS zi = ', e12.4)
  return
end subroutine GetConstrainedData

!**********************************************************************************************************************!
! SetInitialConditions - reads initial conditions file and sets ICs
!**********************************************************************************************************************!
subroutine SetInitialConditions()
  real(kind=dp)     :: ic_z0, ic_zi
  integer(kind=i4)  :: i, l, j
  integer(kind=i4)  :: ic_npts
  integer(kind=i4)  :: nics
  character(len=16) :: ic_gasmech
  integer(kind=i4), parameter   :: niclines=8
  integer(kind=i4), allocatable :: icsp(:), icu(:)

  ! initially set all concentrations to zero
  ! only those specified in IC file are non-zero
  do i=1,npts
    do l=1,ninteg
      cint(i,l) = 0.0D+00
    end do
  end do

  ! read IC file
  open(unit=UICS,file=('./data/' // icfile))
  ! skip header lines
  do j=1,niclines
    read(UICS,*)
  end do
  ! read consistency data and check
  read(UICS,101) ic_gasmech
  ic_gasmech=trim(ic_gasmech)
  read(UICS,*)   ic_npts
  read(UICS,*)   ic_z0
  read(UICS,*)   ic_zi 
  if (ic_gasmech .ne. gasmech) then
    write(*,201) icfile
    write(*,202) ic_gasmech, gasmech
    close(UICS)
    stop
  else if (ic_npts .ne. npts) then
    write(*,201) icfile
    write(*,203) ic_npts, npts
    close(UICS)
    stop
  else if (ic_z0*100. .ne. z0) then
    write(*,201) icfile
    write(*,204) ic_z0, z0
    close(UICS)
    stop
  else if (ic_zi*100. .ne. zi) then
    write(*,201) icfile
    write(*,205) ic_zi, zi
    close(UICS)
    stop
  end if

  ! made it this far, so read the species and units
  read(UICS,*)
  read(UICS,*)
  read(UICS,*)
  read(UICS,*) nics
  read(UICS,*)
  read(UICS,*)
  allocate(icsp(nics))
  allocate(icu(nics))
  do l=1,nics
    read(UICS,*) icsp(l), icu(l)
  end do
  read(UICS,*)
  read(UICS,*)
  ! now read the vertical profiles of each species
  do i=1,npts
    read(UICS,*) (cint(i,icsp(l)),l=1,nics) 
  end do
  do i=1,npts
    print *, 'Init: O3 = ', cint(i, iO3)
  end do
  ! convert the ppbv species (where icu(l)=1) to molec/cm3
  do l=1,nics
    if (icu(l) .eq. 1) then
      do i=1,npts
        cint(i,icsp(l)) = ConvertPPBToMolecCC(cint(i,icsp(l)),i)
      end do
    end if
  end do

  cintp = cint

  deallocate(icsp)
  deallocate(icu)

101 format(a16)
201 format('***Initial conditions file ', a, ' is inconsistent with current ACCESS configuration!')
202 format('***IC gasmech = ', a /'***ACCESS gasmech = ', a)
203 format('***IC npts = ', i4 /'***ACCESS npts = ', i4)
204 format('***IC z0 = ', e12.4 /'***ACCESS z0 = ', e12.4)
205 format('***IC zi = ', e12.4 /'***ACCESS zi = ', e12.4)

  return
end subroutine SetInitialConditions

!**********************************************************************************************************************!
! GetCanopyData - reads canopy morphology data from input canopy file
!**********************************************************************************************************************!
subroutine GetCanopyData()
  integer(kind=i4) :: i, j
  integer(kind=i4), parameter :: nchlines=15
  integer(kind=i4) :: cm_npts
  real(kind=dp)    :: cm_z0, cm_zi, cm_hc, cm_alfa

  open(unit=UCNPY,file=('./data/' // cnpyfile))
  do j=1,nchlines
    read(UCNPY,*)
  end do
  ! read grid consistency data and check
  read(UCNPY,*)   cm_npts
  read(UCNPY,*)   cm_z0
  read(UCNPY,*)   cm_zi 
  read(UCNPY,*)   cm_hc
  read(UCNPY,*)   cm_alfa
  if (cm_npts .ne. npts) then
    write(*,201) cnpyfile
    write(*,203) cm_npts, npts
    close(UCNPY)
    stop
  else if (cm_z0 .ne. z0) then
    write(*,201) cnpyfile
    write(*,204) cm_z0, z0
    close(UCNPY)
    stop
  else if (cm_zi .ne. zi) then
    write(*,201) cnpyfile
    write(*,205) cm_zi, zi
    close(UCNPY)
    stop
  else if (cm_hc .ne. hc) then
    write(*,201) cnpyfile
    write(*,206) cm_hc, hc
    close(UCNPY)
    stop
  else if (cm_alfa .ne. alfa) then
    write(*,201) cnpyfile
    write(*,207) cm_alfa, alfa
    close(UCNPY)
    stop
  end if 

  ! passed checks, now read canopy data

  ! leaf angle distribution parameter
  read(UCNPY,*)
  read(UCNPY,*)
  read(UCNPY,*)   x
  read(UCNPY,*)
  read(UCNPY,*)

  ! diffuse radiation extinction coefficient
  read(UCNPY,*)   kd 
  read(UCNPY,*)
  read(UCNPY,*)

  ! g0, Medlyn et al. stomatal conductance parameter
  read(UCNPY,*)   g0 
  read(UCNPY,*)
  read(UCNPY,*)
 
  ! g1, Medlyn et al. stomatal conductance parameter
  read(UCNPY,*)   g1 
  read(UCNPY,*)
  read(UCNPY,*)
 
  ! dleaf, characteristic leaf dimension
  read(UCNPY,*)   dleaf
  read(UCNPY,*)
  read(UCNPY,*)
  read(UCNPY,*)

  ! zero arrays
  do i=1,npts
    lad(i) =0.0
    lai(i) =0.0
    clai(i)=0.0
  end do

  ! read canopy lad profile
  laitot=0.0
  do i=ncnpy,1,-1
    read(UCNPY,*) lad(i)
    lai(i)=lad(i)*dzhc             ! within canopy grid resolution is always constant!
    laitot=laitot+lad(i)*dzhc
    clai(i)=laitot
  end do
  close(UCNPY)

201 format('***Canopy data file ', a, ' is inconsistent with current ACCESS configuration!')
203 format('***Canopy npts = ', i4 /'***ACCESS npts = ', i4)
204 format('***Canopy z0 = ', e12.4 /'***ACCESS z0 = ', e12.4)
205 format('***Canopy zi = ', e12.4 /'***ACCESS zi = ', e12.4)
206 format('***Canopy hc = ', e12.4 /'***ACCESS hc = ', e12.4)
207 format('***Canopy alfa = ', e12.4 /'***ACCESS alfa = ', e12.4)
  return

end subroutine GetCanopyData


end module Initialize
!======================================================================================================================!
