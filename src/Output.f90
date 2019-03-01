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
!     Module:       Output                                                                                             !
!                                                                                                                      !
!     Description:  contains routines related to output of results                                                     !
!                                                                                                                      !
!======================================================================================================================!
module Output
  use GlobalData
  use Utils
  implicit none

  private PrinttoSTDOUT, SaveResults, CalculateVertFluxes, CalculateBudgetRates
  public OutputResult, PrintFinaltoFile, PrintCPUtime
contains

!**********************************************************************************************************************!
! OutputResult - prints results to STDOUT and stores for later
!**********************************************************************************************************************!
subroutine OutputResult()

  ! print defined results to STDOUT and simulation runtime file
  call PrinttoSTDOUT()

  ! save results for printing at end of simulation
  call SaveResults()

  nt = nt + 1

  return
end subroutine OutputResult

!**********************************************************************************************************************!
! PrinttoSTDOUT - print selected results to STDOUT
!                 and simultaneously to simulation runtime file
!**********************************************************************************************************************!
subroutine PrinttoSTDOUT()
  integer(kind=i4)  :: i
  integer(kind=i4)  :: l
  character(len=2)  :: ncols
  character(len=35) :: f0
  character(len=30) :: f1
  character(len=30) :: f2
  character(len=4)  :: col1

  ! dynamically create format string
  write(ncols,'(i2)') nstdsp
  f0='(1x, a, 2x, a, e12.3, 2x, a, e12.3)'
  f1='(4x,a4,1x,' // ncols // '(5x,a4,5x))'
  f2='(f8.1,' // ncols // 'e14.6)'
  col1='z(m)'

  ! array stdsp maps species indices for STDOUT output
  write(*,1001) (t-tstart)/3600.0_dp
  write(*,fmt=f0) sdtout(nte-1), 'Rib= ', ribout(nte-1), 'z/L= ', zolout(nte-1)
  write(*,fmt=f1) col1, (sspc(stdsp(l)),l=1,nstdsp)

  write(URUN,1001) (t-tstart)/3600.0_dp
  write(URUN,fmt=f0) sdtout(nte-1), 'Rib= ', ribout(nte-1), 'z/L= ', zolout(nte-1)
  write(URUN,fmt=f1) col1, (sspc(stdsp(l)),l=1,nstdsp)

  ! -ConvertOneToOutputFormat determines the appropriate output units
  !  from data in CTRL file
  do i=npts,1,-1
    write(*,f2) 0.01*z(i), (ConvertOneToOutputFormat(cint(i,stdsp(l)),stdsp(l),i),l=1,nstdsp) 

    write(URUN,f2) 0.01*z(i), (ConvertOneToOutputFormat(cint(i,stdsp(l)),stdsp(l),i),l=1,nstdsp) 
  end do

1001 format(/' t-t0 = ', f6.2, ' hrs')

  return
end subroutine PrinttoSTDOUT

!**********************************************************************************************************************!
! SaveResults - save all results to storage array
!**********************************************************************************************************************!
subroutine SaveResults()
  integer(kind=i4) :: l, i

  do l=1,ninteg
    ! soil exchange
    vsout(l,nt) = vs(l)
    rsoillout(l,nt) = rsoill(l)
 
    ! concentrations, emissions and deposition
    do i=1,npts
      cout(i,l,nt) = cint(i,l)
      rsout(i,l,nt) = rs(i,l)
      rbout(i,l,nt) = rb(i,l)
      rmout(i,l,nt) = rm(i,l)
      rcout(i,l,nt) = rc(i,l)
      vdout(i,l,nt) = vd(i,l)
      qout(i,l,nt)  = q(i,l)
    end do
  end do

  do i=1,npts
    ! met data
    tkout(i,nt)    = tk(i)
    pmbout(i,nt)   = pmb(i)
    qhout(i,nt)    = qh(i)
    ubarout(i,nt)  = ubar(i)
    kvout(i,nt)    = kv(i)
    ppfdout(i,nt)  = ppfd(i)
    cairout(i,nt)  = cair(i)
    h2oout(i,nt)   = h2o(i) 
    rhout(i,nt)    = RelativeHumidity(tk(i), pmb(i), qh(i))

    ! canopy physics data
    ppfdsunout(i,nt) = ppfd_sun(i)
    ppfdshdout(i,nt) = ppfd_shd(i)
    ppfdwgtout(i,nt) = ppfd_wgt(i)
    nirsunout(i,nt)  = nir_sun(i)
    nirshdout(i,nt)  = nir_shd(i)
    nirwgtout(i,nt)  = nir_wgt(i)
    rabssunout(i,nt) = rabs_sun(i)
    rabsshdout(i,nt) = rabs_shd(i)
    rabswgtout(i,nt) = rabs_wgt(i)
    rtsunout(i,nt)   = rt_sun(i)
    rtshdout(i,nt)   = rt_shd(i)
    rtwgtout(i,nt)   = rt_wgt(i)
    lwupout(i,nt)    = lw_up(i)
    lwdnout(i,nt)    = lw_dn(i)
    tlsunout(i,nt)   = tl_sun(i)
    tlshdout(i,nt)   = tl_shd(i)
    tlwgtout(i,nt)   = tl_wgt(i)
    gssunout(i,nt)   = gs_sun(i)
    gsshdout(i,nt)   = gs_shd(i)
    gswgtout(i,nt)   = gs_wgt(i)
    rssunout(i,nt)   = rs_sun(i)
    rsshdout(i,nt)   = rs_shd(i)
    rswgtout(i,nt)   = rs_wgt(i)
    anetsunout(i,nt) = anet_sun(i)
    anetshdout(i,nt) = anet_shd(i)
    anetwgtout(i,nt) = anet_wgt(i)
    fsunout(i,nt)    = fsun(i)
    fshdout(i,nt)    = fshd(i)
  end do
  ppfddirout(nt) = ppfd_direct
  ppfddifout(nt) = ppfd_diffus
  nirdirout(nt)  = nir_direct
  nirdifout(nt)  = nir_diffus

  ! soil physics data
  vsh2oout(nt)  = vsh2o
  qsoilout(nt)  = qsoil
  effrhsoilout(nt) = effrhsoil
  rbgout(nt)    = rbg
  gbgout(nt)    = gbg
  rsoilout(nt)  = rsoil
  tsoilkout(nt) = tsoilk
  tk0out(nt)    = tk0

  call CalculateVertFluxes()
  call CalculateBudgetRates()

  timeout(nt)=(t-tstart)/3600.0_dp     ! convert from seconds to hrs

  return
end subroutine SaveResults

!**********************************************************************************************************************!
! PrintFinaltoFile - prints selected results to individual species files
!                    output species are defined in CTRL file
!**********************************************************************************************************************!
subroutine PrintFinaltoFile()
  integer(kind=i4) :: i, l, m, nstp, me, j
  character(len=3)  :: ncols
  character(len=30) :: f1, f2, f3, f4, f5, f6
  character(len=4)  :: col1
  character(len=8)  :: srxn
  character(len=45) :: ofname
  character(len=40) :: hdr

  ! dynamically create format strings for output
  if( (ntout+1) < 100) then
    write(ncols,'(i2)') ntout+1
  else
    write(ncols,'(i3)') ntout+1
  end if
  f1='(6x,a4,3x,' // trim(ncols) // '(f8.1,5x))'
  f2='(f10.1,' // trim(ncols) // '(e13.5))'
  f3='(a)'
  f4='(a,i5.5)'
  f5='(6a10)'
  f6='(2i10,4f10.4)'
  col1='z(m)'

  ! output grid info
  ofname='./out/' // trim(simname) // '/grid/grid.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f5) 'npts', 'ncnpy', 'dz(m)', 'z0(m)', 'H(m)', 'hc(m)'
  write(UOUT,fmt=f6) npts, ncnpy, dzhc*0.01, z0*0.01, zi*0.01, hc    ! all heights in meters
  close(UOUT)

  ! number of met data time steps per output time step
  nstp=int(dtout/dtenv)

  ! output all species profiles
  ofname='./out/' // trim(simname) // '/sp/sp.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ConvertOneToOutputFormat(cout(i,l,m),l,i),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output deposition velocities
  ofname='./out/' // trim(simname) // '/vd/vd.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(vdout(i,l,m),m=0,ntout)
    end do
  end do
  close(UOUT) 

  ! output leaf boundary layer resistances over simulation
  ofname='./out/' // trim(simname) // '/r/rb.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(rbout(i,l,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output stomatal resistances over simulation
  ofname='./out/' // trim(simname) // '/r/rs.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(rsout(i,l,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output mesophyllic resistances over simulation
  ofname='./out/' // trim(simname) // '/r/rm.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(rmout(i,l,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output cuticular resistances over simulation
  ofname='./out/' // trim(simname) // '/r/rc.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(rcout(i,l,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output surface exchange velocities
  ofname='./out/' // trim(simname) // '/r/vs.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    write(UOUT,fmt=f2) 0.0,(vsout(l,m),m=0,ntout)
  end do
  close(UOUT)

  ! output soil resistances over simulation
  ofname='./out/' // trim(simname) // '/r/rsoill.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    write(UOUT,fmt=f2) 0.0,(rsoillout(l,m),m=0,ntout)
  enddo
  close(UOUT)

  ! output emissions
  ofname='./out/' // trim(simname) // '/emis/emis.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(qout(i,l,m),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output vertical fluxes
  ofname='./out/' // trim(simname) // '/vflux/vflux.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ConvertFluxToStdFormat(vfout(i,l,m)),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output budget emission rates
  ofname='./out/' // trim(simname) // '/budget/bem.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(emtout(i,l,m),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output budget constrained species rates
  ofname='./out/' // trim(simname) // '/budget/bcn.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(cntout(i,l,m),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output budget deposition rates
  ofname='./out/' // trim(simname) // '/budget/bdp.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(dptout(i,l,m),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output budget vertical transport rates
  ofname='./out/' // trim(simname) // '/budget/bvt.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(vttout(i,l,m),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output budget chemical transformation  rates
  ofname='./out/' // trim(simname) // '/budget/bch.out'
  open(UOUT,file=ofname)
  do l=1,ninteg
    write(UOUT,fmt=f3) sspc(l)
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)  
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(chtout(i,l,m),m=0,ntout)
    end do
  enddo
  close(UOUT)

  ! output instantaneous reactions rate coefficients over simulation
  ofname='./out/' // trim(simname) // '/ks/ks.out'
  open(UOUT,file=ofname)
  do j=1,nrxn
    write(srxn,fmt=f4) 'rxn', j
    write(UOUT,fmt=f3) srxn
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(ksout(i,j,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output instantaneous reactions rates over simulation
  ofname='./out/' // trim(simname) // '/rates/rates.out'
  open(UOUT,file=ofname)
  do j=1,nrxn
    write(srxn,fmt=f4) 'rxn', j
    write(UOUT,fmt=f3) srxn
    write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
    do i=1,npts
      write(UOUT,fmt=f2) 0.01*z(i),(rxnout(i,j,m),m=0,ntout)
    end do
  end do
  close(UOUT)

  ! output temperature profiles over simulation
  ofname='./out/' // trim(simname) // '/met/tk.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(tkout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output pressure profiles over simulation
  ofname='./out/' // trim(simname) // '/met/pmb.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(pmbout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output specific humidity profiles over simulation
  ofname='./out/' // trim(simname) // '/met/qh.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(qhout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output mean wind speed profiles over simulation
  ofname='./out/' // trim(simname) // '/met/ubar.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(ubarout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output eddy diffusivity profiles over simulation
  ofname='./out/' // trim(simname) // '/met/kv.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(kvout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output PPFD profiles over simulation
  ofname='./out/' // trim(simname) // '/met/ppfd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(ppfdout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output radiation attenuation profiles over simulation
  ofname='./out/' // trim(simname) // '/met/fj.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(fjout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output air density profiles over simulation
  ofname='./out/' // trim(simname) // '/met/cair.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(cairout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output h2o profiles over simulation
  ofname='./out/' // trim(simname) // '/met/h2o.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(h2oout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output rh profiles over simulation
  ofname='./out/' // trim(simname) // '/met/rh.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rhout(i,m*nstp),m=0,ntout)
  end do
  close(UOUT)

  ! output (z-d)/L over simulation
  ofname='./out/' // trim(simname) // '/met/zol.dat'
  open(UOUT,file=ofname)
  hdr = ' hr          zol()'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), zolout(m)
  end do
  close(UOUT)

  ! output aerodynamic conductances over simulation
  ofname='./out/' // trim(simname) // '/met/gaero.dat'
  open(UOUT,file=ofname)
  hdr = ' hr          gaero(mol/m2-s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), gaero(m)
  end do
  close(UOUT)

  ! output aerodynamic resistances over simulation
  ofname='./out/' // trim(simname) // '/met/ra.dat'
  open(UOUT,file=ofname)
  hdr = ' hr          ra(s/cm)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), ra(m)
  end do
  close(UOUT)

  ! output ppfd_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/ppfdsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(ppfdsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output ppfd_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/ppfdshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(ppfdshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output ppfd_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/ppfdwgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(ppfdwgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output nir_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/nirsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(nirsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output nir_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/nirshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(nirshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output nir_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/nirwgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(nirwgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output lw_up profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/lwup.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(lwupout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output lw_dn profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/lwdn.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(lwdnout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rt_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rtsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rtsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rt_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rtshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rtshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rt_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rtwgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rtwgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rabs_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rabssun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rabssunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rabs_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rabsshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rabsshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rabs_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rabswgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rabswgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output tl_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/tlsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(tlsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output tl_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/tlshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(tlshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output tl_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/tlwgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(tlwgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output gs_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/gssun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(gssunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output gs_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/gsshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(gsshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output gs_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/gswgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(gswgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rs_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rssun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rssunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rs_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rsshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rsshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output rs_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/rswgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(rswgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output anet_sun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/anetsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(anetsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output anet_shd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/anetshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(anetshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output anet_wgt profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/anetwgt.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(anetwgtout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output fsun profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/fsun.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(fsunout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output fshd profiles over simulation
  ofname='./out/' // trim(simname) // '/canopy/fshd.dat'
  open(UOUT,file=ofname)
  write(UOUT,fmt=f1) col1,(timeout(m),m=0,ntout)
  do i=1,npts
    write(UOUT,fmt=f2) 0.01*z(i),(fshdout(i,m),m=0,ntout)
  end do
  close(UOUT)

  ! output vsh2o exchange coeffients over simulation
  ofname='./out/' // trim(simname) // '/soil/vsh2o.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         vsh2o(cm/s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), vsh2oout(m)
  end do
  close(UOUT)

  ! output qsoil values over simulation
  ofname='./out/' // trim(simname) // '/soil/qsoil.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         qsoil(mol/cm3)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), qsoilout(m)
  end do
  close(UOUT)

  ! output effrhsoil over simulation
  ofname='./out/' // trim(simname) // '/soil/effrhsoil.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         effrhsoil'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), effrhsoilout(m)
  end do
  close(UOUT)

  ! output rbg over simulation
  ofname='./out/' // trim(simname) // '/soil/rbg.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         rbg(s/cm)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), rbgout(m)
  end do
  close(UOUT)

  ! output gbg over simulation
  ofname='./out/' // trim(simname) // '/soil/gbg.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         gbg(mol/m2-s)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), gbgout(m)
  end do
  close(UOUT)

  ! output rsoil of water vapor over simulation
  ofname='./out/' // trim(simname) // '/soil/rsoil.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         rsoil(s/cm)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), rsoilout(m)
  end do
  close(UOUT)

  ! output tsoilk over simulation
  ofname='./out/' // trim(simname) // '/soil/tsoilk.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         tsoilk(K)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), tsoilkout(m)
  end do
  close(UOUT)

  ! output tk0 over simulation
  ofname='./out/' // trim(simname) // '/soil/tk0.dat'
  open(UOUT,file=ofname)
  hdr = ' hr         tk0(K)'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    write(UOUT,1001) timeout(m), tk0out(m)
  end do
  close(UOUT)


  ! output elapsed hour/datetime key file
  ofname='./out/' // trim(simname) // '/ACCESS_timekey.dat'
  open(UOUT,file=ofname)
  hdr = ' hr          datetime'
  write(UOUT,1000) trim(hdr)
  do m=0,ntout
    me = m*nstp
    write(UOUT,1002) timeout(m), sdtout(me)
  end do
  close(UOUT)

  ! output species units key file
  ofname='./out/' // trim(simname) // '/ACCESS_ppbv.dat'
  open(UOUT,file=ofname)
  hdr = 'ispec        sspc'           
  write(UOUT,1000) trim(hdr)
  do i=1,noutppb
    write(UOUT,1003) outppb(i), sspc(outppb(i))
  end do
  close(UOUT)

1000 format(4x, a)
1001 format(f8.1, 5x, e13.5)
1002 format(f8.1, 5x, a)
1003 format(3x, i4, 10x, a)

  return
end subroutine PrintFinaltoFile

!**********************************************************************************************************************!
! subroutine PrintCPUtime - prints simulation CPU time to STDOUT and summary
!**********************************************************************************************************************!
subroutine PrintCPUtime(tcpu)
   real(kind=4) :: tcpu

   ! to STDOUT
   write(*,1001) tcpu/3600.,tcpu

   ! to USUMM for archive
   open(USUMM,file=simsummary,status='old',action='write',position='append')
   write(USUMM,1001) tcpu/3600.,tcpu
   close(USUMM)

1001 format('Total CPU time: ', f12.3,' hrs'/'                ',f12.3,' sec')

end subroutine PrintCPUtime

!**********************************************************************************************************************!
! subroutine CalculateVertFluxes - calculate vertical fluxes of all species 
!                                  at all grid points
!**********************************************************************************************************************!
subroutine CalculateVertFluxes()
  integer(kind=i4) :: i, l
  real(kind=dp)    :: kvm, grad

  do l=1,ninteg
    do i=1,npts-1
      kvm = 0.5*(kv(i+1)+kv(i))
      grad  = (cint(i+1,l)-cint(i,l))/(z(i+1)-z(i))
      vfout(i,l,nt) = - kvm*grad
    end do
    kvm = kv(npts)
    grad = (caloft(l)-cint(npts,l))/(z(npts)-z(npts-1))
    vfout(npts,l,nt) = - kvm*grad
  end do
    
  return
end subroutine CalculateVertFluxes

!**********************************************************************************************************************!
! subroutine CalculateBudgetRates - calculate budget rates for all processes 
!                                   (ppbv/hr)
!**********************************************************************************************************************!
subroutine CalculateBudgetRates()
  integer(kind=i4) :: i, l

  ! xxtout accumulations are in ppbv
  ! dividing by dtout gives the rate in ppbv/s
  ! multiplying by 3600.0 converts to ppbv/hr
  do l=1,ninteg
    do i=1,npts
      emtout(i,l,nt)=emtout(i,l,nt)*3600.0/dtout
      cntout(i,l,nt)=cntout(i,l,nt)*3600.0/dtout
      dptout(i,l,nt)=dptout(i,l,nt)*3600.0/dtout
      vttout(i,l,nt)=vttout(i,l,nt)*3600.0/dtout
      chtout(i,l,nt)=chtout(i,l,nt)*3600.0/dtout
    end do
  end do 

  return
end subroutine CalculateBudgetRates

end module Output
!======================================================================================================================!
