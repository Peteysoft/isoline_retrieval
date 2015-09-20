PROGRAM rttovscatt_test
  !
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  !
  ! Description:
  ! 1- read ECMWF profiles on model levels
  ! 2- interpolate the T, q profiles to the 43-level rttov grid
  ! 4- run rttovscatt direct model
  ! 5- test TL, AD and K codes of the rttovscatt package
  !
  ! Method:
  ! see comments in program
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date        Comment
  ! -------   ----        -------
  !           10/2002     Initial version (F. Chevallier)
  !           05/2003     New F90 code with structures (F. Chevallier P. Brunel)
  !           10/2003     Add K (F. Chevallier)
  !           10/3/04     Merged in polarimetry code (R Saunders)
  !           23/3/04     Added new precision test (R Saunders)
  !
  ! Code Description:
  !   Language:          Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  !
  Use rttov_const, only :   &
       & errorstatus_fatal,   &
       & errorstatus_success, &
       & default_err_unit,    &
       & sensor_id_mw,        &
       & npolar_return,       &
       & npolar_compute

  Use rttov_types, only : &
       & geometry_type        ,&
       & rttov_coef           ,&
       & rttov_scatt_coef     ,&
       & profile_type         ,&
       & profile_cloud_type   ,&
       & transmission_type    ,&
       & radiance_cloud_type

  Use parkind1, Only : jpim     ,jprb
  IMPLICIT NONE

#include "rttov_errorhandling.interface"
#include "rttov_readcoeffs.interface"
#include "rttov_initcoeffs.interface"
#include "rttov_intex.interface"
#include "rttov_readscattcoeffs.interface"
#include "rttov_scatt.interface"
#include "rttov_scatt_tl.interface"
#include "rttov_scatt_ad.interface"
#include "rttov_scatt_k.interface"

  ! Program arguments:

  ! Local parameters:
  Integer(Kind=jpim), parameter :: idim=1000
  Integer(Kind=jpim), parameter :: nwp_levels=60

  type( rttov_coef )                    :: coef
  type( rttov_coef )                    :: coef_rttov    ! (Only one instrument)
  type( rttov_scatt_coef )              :: coef_scatt
  type(profile_type), allocatable       :: profiles(:)
  type(profile_cloud_type), allocatable :: cld_profiles(:)
  type(radiance_cloud_type)             :: radiance
  Real(Kind=jprb),    Allocatable                  :: emissivity (:)

  ! Taylor test
  type(profile_type), allocatable       :: profiles2(:)
  type(profile_cloud_type), allocatable :: cld_profiles2(:)
  type(radiance_cloud_type)             :: radiance2
  Real(Kind=jprb),    Allocatable                  :: emissivity2 (:)

  ! TL arrays
  type(profile_type), allocatable       :: prof_inc(:)
  type(profile_cloud_type), allocatable :: cld_prof_inc(:)
  type(radiance_cloud_type)             :: radiance_tl
  Real(Kind=jprb),    Allocatable                  :: emissivity_inc (:)

  type(profile_type), allocatable       :: prof_inc2(:)
  type(profile_cloud_type), allocatable :: cld_prof_inc2(:)
  type(radiance_cloud_type)             :: radiance_tl2
  Real(Kind=jprb),    Allocatable                  :: emissivity_inc2 (:)

  ! AD arrays
  type(profile_type), allocatable       :: profiles_ad(:)
  type(profile_cloud_type), allocatable :: cld_profiles_ad(:)
  type(radiance_cloud_type)             :: radiance_inc
  Real(Kind=jprb),    Allocatable                  :: emissivity_ad (:)

  type(profile_type), allocatable       :: profiles_ad2(:)
  type(profile_cloud_type), allocatable :: cld_profiles_ad2(:)
  type(radiance_cloud_type)             :: radiance_inc2
  type(transmission_type)               :: transmission
  Real(Kind=jprb),    Allocatable                  :: emissivity_ad2 (:)

  ! K arrays
  type(profile_type), allocatable       :: profiles_k(:)
  type(profile_cloud_type), allocatable :: cld_profiles_k(:)
  type(radiance_cloud_type)             :: cld_radiance_k
  Real(Kind=jprb),    Allocatable                  :: emissivity_k (:)

  ! Local arrays:
  Real(Kind=jprb), allocatable      :: emis(:)
  Integer(Kind=jpim), allocatable   :: lchan(:)

  Integer(Kind=jpim)              :: coef_errorstatus      ! read coeffs error return code
  Integer(Kind=jpim), Allocatable :: rttov_errorstatus(:)  ! rttov error return code

  integer(Kind=jpim) :: nbtout
  integer(Kind=jpim) :: nfrequencies
  Integer(Kind=jpim) :: nchannels
  Integer(Kind=jpim) :: nprofiles
  Integer(Kind=jpim), Allocatable :: channels   (:)
  integer(Kind=jpim), Allocatable :: polarisations   (:,:)
  Integer(Kind=jpim), Allocatable :: lprofiles  (:)
  integer(Kind=jpim), Allocatable :: frequencies   (:)
  Real(Kind=jprb),    Allocatable :: radiance_total_ref (:)
  logical,            Allocatable :: calcemis  (:)


  Real(Kind=jprb), dimension(nwp_levels)     :: t, q, o3, co2, cc, clw, ciw, rain, sp
  !

  ! Local scalars:
  !Character (len=80) :: errMessage
  !Character (len=12) :: NameOfRoutine = 'rttovscatt_test '
  Integer(Kind=jpim) :: j
  Integer(Kind=jpim) :: ioff
  Integer(Kind=jpim) :: kinrad
  Real(Kind=jprb)       :: tbobs(7)
  Real(Kind=jprb)       :: st, t2m, q2m, u10, v10, psurf, zenangle, azangle
  Real(Kind=jprb)       :: rsatid, rlsm, rlon, rlat
  Real(Kind=jprb)       :: x, lambda
  Real(Kind=jprb)       :: ratio(4)
  Integer(Kind=jpim)    :: jdat
  Integer(Kind=jpim)    :: iatm, ichan
  Integer(Kind=jpim)    :: iexp
  Integer(Kind=jpim)    :: ioout, ioin
  Integer(Kind=jpim)    :: ich, isatid, nch, ich2
  Integer(Kind=jpim)    :: i, ii, nchan
  Integer(Kind=jpim)    :: lev , n, kpol
  Integer(Kind=jpim)    :: ichannels, ibtout, jch, pol_id, jchan
  Logical    :: switchrad  ! true if input is BT

  Integer(Kind=jpim) :: instrument(3) ! instrument triplet
  Real(Kind=jprb) ::  zdelta1, zdelta2
  Real(Kind=jprb) ::  z, eps, threshold, prec_factor

  Integer(Kind=jpim) :: Err_Unit        ! Logical error unit (<0 for default)
  Integer(Kind=jpim) :: verbosity_level ! (<0 for default)

  !!!!!  Some of my modifications (PM)   !!!!
  Integer(Kind=jpim) :: kice  		! ice crystal type
  Integer(Kind=jpim) :: kradip 		! ice effective size scheme
  Character(len=200) :: infile
  Character(len=200) :: resultsfile

  ! End of program arguments

  !-----End of header-----------------------------------------------------

  !Initialise error management with default value for
  ! the error unit number and
  ! Fatal error message output
  Err_unit = -1
  verbosity_level = 1
  call rttov_errorhandling(Err_unit, verbosity_level)

  ! Machine accuracy
  !eps = 1._JPRB
  !do while ((1+eps) > 1._JPRB)
  !  eps = eps /2._JPRB
  !enddo
  !   'threshold' is the maximum difference which is tolerated
  !   between two real numbers for them to be considered equal
  !   On some systems, 'threshold' can be set to a value as low as 10*eps
  !   Some other systems are not rigorous enough and larger values
  !   for prec_factor have to be used. It is set as default as 10.
  !prec_factor = 10000000._JPRB                         !  Edit for your machine
  !threshold = prec_factor * eps
  eps = 10._JPRB * epsilon( 1._JPRB )
  threshold =  eps

  write(*,*)
  write(*,*) 'Note that double precision should be specified in the Makefile'
  write(*,*) '   since lapack double precision routines are used'
  write(*,*)
  !
  write(*,*) 'Radiances(1) or Tbs(2)?'
  read(*,*) kinrad
  switchrad = kinrad == 2

  !!!  My modifications... (Petey)  !!!

  !
  ! Set satellite configuration
  !      only one satellite processed
  !instrument(1)=2  ! DMSP
  !write(*,*) ' DMSP sat id?'
  instrument(1)=1  ! DMSP
  write(*,*) ' NOAA sat id?'
  read(*,*) isatid
  instrument(2)=isatid
  write(*, *) ' Instrument id?'
  read(*,*) instrument(3)
  !instrument(3)=6  ! SSMI

  write(*,*) ' Zenith angle?'
  read(*,*) zenangle
  write(*,*) ' Ice crystal type (0=hex columns, 1=aggregate) ?'
  read(*,*) kice
  write(*,*) ' Ice effective size scheme'
  write(*,*) ' (0=Ou-Liou, 1=Wyser, 2=Boudala, 3=McFarquhar) ?'
  read(*,*) kradip

  write(*,*) ' Profile file name ?'
  read(*,*) infile

  write(*,*) ' Results file name ?'
  read(*,*) resultsfile

  !!!  End of my modifications...   !!!

  ! Read coef file for gaseous absorption
  call rttov_readcoeffs  (coef_errorstatus, coef_rttov, instrument)
  call rttov_initcoeffs  (coef_errorstatus, coef_rttov)
  nchan=coef_rttov%fmv_chn

  if(coef_errorstatus /= errorstatus_success ) then
     write ( ioout, * ) 'rttov_readcoeffs fatal error'
     stop
  endif

  if( any(coef_rttov%ff_val_chn( 1 : coef_rttov%fmv_chn ) /= 1 )) then
     WRITE(*,*) ' some requested channels have bad validity parameter'
     do i = 1, coef_rttov%fmv_chn
        write(*,*) i, coef_rttov%ff_val_chn(i)
     end do
  endif

  ! Read coef file for cloud/rain absorption/scattering
  Call rttov_readscattcoeffs  (coef_errorstatus, coef_rttov, coef_scatt)

  if(coef_errorstatus /= errorstatus_success ) then
     write ( ioout, * ) 'rttov_readscattcoeffs fatal error'
     stop
  endif

  !
  ! Open files
  !
  ioin = 1
!  open(ioin,file='profiles2_fmt',form='formatted',status='old')
  !!!! More modifications...   !!!!
  write (*, *) 'Opening file: ', infile, ' for input'
  open(ioin,file=infile,form='formatted',status='old')
  
  ioout = 2
!  open(ioout,file='outputscatt.ascii',form='formatted')
  write (*, *) 'Opening file: ', resultsfile, ' for output'
  open(ioout,file=resultsfile,form='formatted')

  !
  ! Count number of profiles
  !
  do iatm = 1,idim
     do i = 1,43
        read(ioin,*,end=50)
     enddo
  enddo
50 continue
   nprofiles = iatm - 1
  write(ioout,*) 'This dataset is made of ',nprofiles,' ECMWF model profiles'
  rewind(ioin)
 !
  ! Find out size of channel arrays summing all polarisation states required.
  nch = 0
  ichannels=0
  ibtout=0
  DO  J=1,nprofiles
        DO  JCH=1,NCHAN
           nch = nch +1
           If( coef_rttov%id_sensor /= sensor_id_mw) then
              ichannels=ichannels+1
              ibtout=ibtout+1
           End If
           If( coef_rttov % id_sensor == sensor_id_mw) then
              pol_id = coef_rttov % fastem_polar(jch) + 1
              ichannels=ichannels+npolar_compute(pol_id)
              ibtout=ibtout+npolar_return(pol_id)
           End If
        End Do
   End Do
   nchannels = ichannels
   nbtout = ibtout
  ! Set list of channels and corresponding emissivities
  !   (process all channels and let RTTOV compute the emissivity)
  !
  allocate(lchan (nchannels) )
  allocate(emis  (nchannels) )
  emis(:) = 0._JPRB
  do i = 1 , nchannels
    lchan(i) = i
  enddo
  !
  ! Initialisations and allocations
  ! NO allocation for CO2 profiles
  !
  nfrequencies = nprofiles * nchan
  Allocate( channels   ( nfrequencies ) )
  allocate( lprofiles  ( nfrequencies ) )
  allocate( emissivity ( nchannels ) )
  allocate( calcemis   ( nchannels ) )
  allocate( polarisations ( nchannels ,3) )
  allocate( frequencies ( nbtout ) )
  allocate( transmission % tau_surf      ( nchannels ) )
  allocate( transmission % tau_layer     ( coef_rttov % nlevels, nchannels ) )
  allocate( transmission % od_singlelayer( coef_rttov % nlevels, nchannels ) )

  allocate( rttov_errorstatus(nprofiles))

  ! Profiles on RTTOV pressure levels
  allocate( profiles(nprofiles))
  do j = 1, nprofiles
     ! allocate model profiles atmospheric arrays with model levels dimension
     profiles(j) % nlevels =  coef_rttov % nlevels
     allocate( profiles(j) % p  ( coef_rttov % nlevels ) )
     allocate( profiles(j) % t  ( coef_rttov % nlevels ) )
     allocate( profiles(j) % q  ( coef_rttov % nlevels ) )
     allocate( profiles(j) % o3 ( coef_rttov % nlevels ) )
     allocate( profiles(j) % clw( coef_rttov % nlevels ) )
     profiles(j) % p(:) = coef_rttov % ref_prfl_p(:)
  end do

  ! Cloud additional profiles
  allocate( cld_profiles(nprofiles))
  do j = 1, nprofiles
     ! allocate model profiles atmospheric arrays with model levels dimension
     cld_profiles(j) % nlevels =  nwp_levels
     allocate( cld_profiles(j) % p  ( nwp_levels ) )
     allocate( cld_profiles(j) % ph ( nwp_levels+1 ) )
     allocate( cld_profiles(j) % t  ( nwp_levels ) )
     allocate( cld_profiles(j) % cc ( nwp_levels ) )
     allocate( cld_profiles(j) % clw( nwp_levels ) )
     allocate( cld_profiles(j) % ciw( nwp_levels ) )
     allocate( cld_profiles(j) % rain( nwp_levels ) )
     allocate( cld_profiles(j) % sp  ( nwp_levels ) )
  end do

  ! allocate radiance results arrays with number of channels
  allocate( radiance % clear    ( nchannels ) )
  allocate( radiance % cloudy   ( nchannels ) )
  allocate( radiance % total    ( nchannels ) )
  allocate( radiance % bt       ( nchannels ) )
  allocate( radiance % bt_clear ( nchannels ) )
  allocate( radiance % upclear  ( nchannels ) )
  allocate( radiance % dnclear  ( nchannels ) )
  allocate( radiance % reflclear( nchannels ) )
  allocate( radiance % overcast ( nwp_levels, nchannels ) )
  allocate( radiance % freq_used( nchannels) )
  allocate( radiance_total_ref  ( nchannels ) )
  allocate( radiance % out  ( nbtout ) )
  allocate( radiance % out_clear( nbtout ) )
  allocate( radiance % total_out( nbtout ) )
  allocate( radiance % clear_out( nbtout ) )
  !
  ! Read profile dataset
  !

  iatmloop : do iatm = 1, nprofiles
     read(ioin,'(10e16.6)') rlon,         &! longitude (deg)
          & rlat,         &! latitude (deg)
          & rlsm,         &! land-sea mask (1=land)
          & st,           &! surface temperature (K)
          & psurf,        &! surface pressure (Pa)
          & t2m,          &! 2-meter temperature (K)
          & q2m,          &! 2-meter specific humidity (kg/kg)
          & u10,          &! 10-meter wind u (m/s)
          & v10            ! 10-meter wind u (m/s)
     read(ioin,'(10e16.6)') t   ! temperature (K)
     read(ioin,'(10e16.6)') q   ! specific humidity (kg/kg)
     read(ioin,'(10e16.6)') cc  ! cloud cover
     read(ioin,'(10e16.6)') clw ! liquid water (kg/kg)
     read(ioin,'(10e16.6)') ciw ! ice water (kg/kg)
     read(ioin,'(10e16.6)') rain  ! rain (kg/m2/s)
     read(ioin,'(10e16.6)') sp    ! solid precipitation (kg/m2/s)

     o3(:) = 1.e-7_JPRB
     !zenangle = 53.1_JPRB
     !zenangle = 0._JPRB	! Commented out (PM)--now read in
     azangle = 0._JPRB              ! azimuth angle
     q(:) = max(q(:),0._JPRB)
     clw(:) = max(clw(:),0._JPRB)
     ciw(:) = max(ciw(:),0._JPRB)
     rain(:) = max(rain(:),0._JPRB)
     sp(:) = max(sp(:),0._JPRB)

     !*get model vertical pressures from surface pressure (all Pa)
     call ec_p60l(                    &
           & psurf                    ,&
           & cld_profiles( iatm ) % p ,&
           & cld_profiles( iatm ) % ph )

     ! Convert to hPa
     cld_profiles( iatm ) % p(:)  = cld_profiles( iatm ) % p(:)  /100._JPRB
     cld_profiles( iatm ) % ph(:) = cld_profiles( iatm ) % ph(:) /100._JPRB

     ! Move to structures
     profiles( iatm ) % clw(:) = 0._JPRB
!     profiles( iatm ) % o3 (:) = 0._JPRB
     profiles( iatm ) % o3 (:) = 1.e-7_JPRB   !RWS
     profiles( iatm ) % s2m % p = psurf / 100._JPRB
     profiles( iatm ) % s2m % q = (q2m / (1-q2m)) * 1.60771704_JPRB*1e+06        ! ppmv
     profiles( iatm ) % s2m % o = 0._JPRB
     profiles( iatm ) % s2m % t = t2m
     profiles( iatm ) % s2m % u = u10
     profiles( iatm ) % s2m % v = v10
     profiles( iatm ) % skin % surftype = Int(1.0_JPRB - rlsm)
     profiles( iatm ) % skin % t        = st
     profiles( iatm ) % skin % fastem(:) = (/ 3.0_JPRB, 5.0_JPRB, 15.0_JPRB, 0.1_JPRB, 0.3_JPRB /)

     profiles( iatm ) % ozone_data = .false.  ! because no o3 available for test !!!!WARNING
     profiles( iatm ) % co2_data   = .false.
     profiles( iatm ) % clw_data   = .false.  ! switch for RTTOV only
     profiles( iatm ) % zenangle   = zenangle
     profiles( iatm ) % azangle    = azangle
     profiles( iatm ) % ctp        = 500._JPRB  ! default value
     profiles( iatm ) % cfraction  = 0._JPRB    ! default value


     cld_profiles( iatm ) % t(:)   = t(:)
     cld_profiles( iatm ) % cc(:)  = cc(:)
     cld_profiles( iatm ) % clw(:) = clw(:)
     cld_profiles( iatm ) % ciw(:) = ciw(:)
     cld_profiles( iatm ) % rain(:) = rain(:)
     cld_profiles( iatm ) % sp  (:) = sp  (:)

     ! ice size scheme:
     cld_profiles( iatm ) % kice = kice
     cld_profiles( iatm ) % kradip = kradip

     ! convert q to ppmv
     q(:)  = (q(:)  / (1- q(:)))  * 1.60771704_JPRB *1e+06

     ! convert input profile to RTTOV pressure levels
     call rttov_intex (                  &
         & nwp_levels,                     &
         & coef_rttov%nlevels,             &
         & cld_profiles(iatm) % p,         &
         & profiles(iatm) % p,             &
         & t(:),                           &
         & profiles(iatm) % t)
     call rttov_intex (                  &
         & nwp_levels,                     &
         & coef_rttov%nlevels,             &
         & cld_profiles(iatm) % p,         &
         & profiles(iatm) % p,             &
         & q(:),                           &
         & profiles(iatm) % q)

  enddo iatmloop

  close(ioin)

  ! Channel, profile list and emissivity arrays
  nch = 0
  ichannels=0
  ibtout=0
  DO iatm = 1, nprofiles
     ioff = (iatm - 1) * nchan
     channels(1+ioff:nchan+ioff)   = lchan(1:nchan)
     lprofiles(1+ioff:nchan+ioff)  = iatm
     !
     DO  JCH=1,nchan
       nch = nch +1
       polarisations(nch,1)=ichannels+1
       If( coef_rttov % id_sensor /= sensor_id_mw) then
          ! Note if input polarisation used this is only valid for a single polarisation option.
          emissivity( ichannels+1 ) = emis(jch)
          ichannels=ichannels+1
          polarisations(nch,2) = nch
          polarisations(nch,3) = 1
          frequencies(ibtout+1) = nch
          ibtout=ibtout+1
       End If
       If( coef_rttov% id_sensor == sensor_id_mw) then
          pol_id = coef_rttov % fastem_polar(jch) + 1
          Do ich2=1,npolar_compute(pol_id)
             emissivity(ichannels+ich2)=emis(jch)
          enddo
          Do n=ichannels+1,ichannels+npolar_compute(pol_id)
             polarisations(n,2)=nch
          End Do
          ichannels=ichannels+npolar_compute(pol_id)
          Do i=1, npolar_return(pol_id)
             frequencies(ibtout+i)=nch
          End Do
          polarisations(nch,3)=npolar_compute(pol_id)
          ! Note if input polarisation used this is only valid for a single polarisation option.
          ! We will need to know which frequency each element in the output array corresponds to
          ibtout=ibtout+npolar_return(pol_id)
       End If
     End Do
  End do
  nchannels = ichannels
  nbtout = ibtout
  calcemis(:)        = emissivity(:) < 0.01_JPRB

  !
  ! Call RTTOV_SCATT
  !

  write(ioout,*)
  write(ioout,*) 'Call to RTTOV_SCATT'
  write(ioout,*) '-------------------'
  write(ioout,*)
  write(6,*)' nfreq=',nfrequencies,' nchannels=',nchannels,' nbtout=',nbtout
  Call rttov_scatt( &
       & rttov_errorstatus,  &! out
       & nfrequencies,   &! in
       & nchannels,      &! in
       & nbtout,         &! in
       & nprofiles,      &! in
       & channels,       &! in
       & polarisations,  &! in
       & lprofiles,      &! in
       & profiles,       &! inout  (to invalid clw absorption)
       & cld_profiles,   &! in
       & coef_rttov,     &! in
       & coef_scatt,     &! in
       & calcemis,       &! in
       & emissivity,     &! inout
       & transmission,   &! inout
       & radiance )       ! inout

  If ( any( rttov_errorstatus(:) == errorstatus_fatal ) ) Then
     Do iatm = 1, nprofiles
        If ( rttov_errorstatus(iatm) == errorstatus_fatal ) Then
           write ( ioout, * ) 'rttov_scatt error for profile',iatm
        End If
     End Do
     Stop
  End If


  ! main output:
  ! radiance%total_out  = cloud-affected radiances
  ! radiance%clear_out  = clear-sky radiances
  ! radiance%out        = cloud-affected Tbs
  ! radiance%out_clear  = clear-sky Tbs
  !

  if (kinrad == 2) then
     write(ioout,*) 'Channel  cloudy Tb    clear Tb'
     do ichan = 1, nbtout
        write(ioout,'(i6,2x,30e23.16)')      &
              & ichan                       ,&
              & radiance%out(ichan)      ,&
              & radiance%out_clear(ichan)
     enddo
  else
     write(ioout,*) 'Channel cloudy Rad   clear Rad'
     do ichan = 1, nbtout
        write(ioout,'(i6,2x,30e23.16)')     &
              & ichan                      ,&
              & radiance%total_out(ichan)  ,&
              & radiance%clear_out(ichan)
     enddo
  endif

  !!!!  Important--more modifications--program exits here   !!!!
  close(ioout)

  stop

  !---------------------------------------------------------------------
  ! Test of TL
  !---------------------------------------------------------------------
  !
  write(ioout,*)
  write(ioout,*) 'Test TL'
  write(ioout,*) '-------'
  write(ioout,*)

  ! Set perturbations = 5% of initial profile
  allocate ( prof_inc(     nprofiles ))
  allocate ( cld_prof_inc( nprofiles ))
  Do j = 1, nprofiles
     prof_inc(j) % nlevels =  coef_rttov % nlevels
     allocate( prof_inc(j) % p  ( coef_rttov % nlevels ) )
     allocate( prof_inc(j) % t  ( coef_rttov % nlevels ) )
     allocate( prof_inc(j) % q  ( coef_rttov % nlevels ) )
     allocate( prof_inc(j) % o3 ( coef_rttov % nlevels ) )
     allocate( prof_inc(j) % clw( coef_rttov % nlevels ) )

     prof_inc(j) % ozone_Data = .False.  ! no meaning
     prof_inc(j) % co2_Data   = .False.  ! no meaning
     prof_inc(j) % clw_Data   = .False.  ! no meaning
     prof_inc(j) % zenangle   = -1  ! no meaning
     prof_inc(j) % azangle   = -1  ! no meaning

     ! increments for atmospheric variables
     prof_inc(j) % p(:)   = 0._JPRB    ! no tl on pressure levels
     prof_inc(j) % t(:)   = profiles(j) % t(:)  *0.05_JPRB
     prof_inc(j) % o3(:)  = profiles(j) % o3(:) *0.05_JPRB
     prof_inc(j) % clw(:) = profiles(j) % clw(:)*0.05_JPRB
     prof_inc(j) % q(:)   = profiles(j) % q(:)  *0.05_JPRB

     ! increments for air surface variables
     prof_inc(j) % s2m % t = profiles(j) % s2m % t  *0.05_JPRB
     prof_inc(j) % s2m % q = profiles(j) % s2m % q  *0.05_JPRB
     prof_inc(j) % s2m % o = profiles(j) % s2m % o  *0.05_JPRB
     prof_inc(j) % s2m % p = profiles(j) % s2m % p  *0.05_JPRB
     prof_inc(j) % s2m % u = profiles(j) % s2m % u  *0.05_JPRB
     prof_inc(j) % s2m % v = profiles(j) % s2m % v  *0.05_JPRB

     ! increments for skin variables
     prof_inc(j) % skin % surftype = -1  ! no meaning
     prof_inc(j) % skin % t        = profiles(j) % skin % t  *0.05_JPRB
     prof_inc(j) % skin % fastem(:)= profiles(j) % skin % fastem(:) *0.05_JPRB

     ! increments for cloud variables
     prof_inc(j) % ctp       = profiles(j) % ctp       *0.05_JPRB
     prof_inc(j) % cfraction = profiles(j) % cfraction *0.05_JPRB

     ! increments for cloud variables
     cld_prof_inc(j) % nlevels =  nwp_levels
     allocate( cld_prof_inc(j) % p  ( nwp_levels ) )
     allocate( cld_prof_inc(j) % ph ( nwp_levels+1 ) )
     allocate( cld_prof_inc(j) % t  ( nwp_levels ) )
     allocate( cld_prof_inc(j) % cc ( nwp_levels ) )
     allocate( cld_prof_inc(j) % clw( nwp_levels ) )
     allocate( cld_prof_inc(j) % ciw( nwp_levels ) )
     allocate( cld_prof_inc(j) % rain( nwp_levels ) )
     allocate( cld_prof_inc(j) % sp  ( nwp_levels ) )
     cld_prof_inc(j) % p(:)   = cld_profiles(j) % p(:)   *0.05_JPRB
     cld_prof_inc(j) % ph(:)  = cld_profiles(j) % ph(:)  *0.05_JPRB
     cld_prof_inc(j) % t(:)   = cld_profiles(j) % t(:)   *0.05_JPRB
     cld_prof_inc(j) % cc(:)  = cld_profiles(j) % cc(:)  *0.05_JPRB
     cld_prof_inc(j) % clw(:) = cld_profiles(j) % clw(:) *0.05_JPRB
     cld_prof_inc(j) % ciw(:) = cld_profiles(j) % ciw(:) *0.05_JPRB
     cld_prof_inc(j) % rain(:) = cld_profiles(j) % rain(:) *0.05_JPRB
     cld_prof_inc(j) % sp  (:) = cld_profiles(j) % sp  (:) *0.05_JPRB
  End Do

  ! emissivity
  allocate( emissivity_inc( nchannels ))
  calcemis(:)      = .false.

  emissivity_inc(:) = emissivity(:) * 0.05_JPRB

  ! allocate radiance results arrays with number of channels
  allocate( radiance_tl % clear    ( nchannels ) )
  allocate( radiance_tl % cloudy   ( nchannels ) )
  allocate( radiance_tl % total    ( nchannels ) )
  allocate( radiance_tl % bt       ( nchannels ) )
  allocate( radiance_tl % bt_clear ( nchannels ) )
  allocate( radiance_tl % out       ( nbtout ) )
  allocate( radiance_tl % out_clear ( nbtout ) )
  allocate( radiance_tl % total_out ( nbtout ) )
  allocate( radiance_tl % clear_out ( nbtout ) )
  allocate( radiance_tl % upclear  ( nchannels ) )
  allocate( radiance_tl % reflclear( nchannels ) )
  allocate( radiance_tl % overcast ( nwp_levels, nchannels ) )

  !---------------------------
  Call rttov_scatt_tl ( &
     & rttov_errorstatus,  &! out
     & nfrequencies,    &! in
     & nchannels,       &! in
     & nbtout,          &! in
     & nprofiles,       &! in
     & channels,        &! in
     & polarisations,   &! in
     & lprofiles,       &! in
     & profiles,        &! in
     & cld_profiles,    &! in
     & coef_rttov,      &! in
     & coef_scatt,      &! in
     & calcemis,        &! in
     & emissivity,      &! inout
     & prof_inc,        &! in
     & cld_prof_inc,    &! in
     & emissivity_inc,  &! inout
     & radiance,        &! inout
     & radiance_tl     ) ! inout

  If ( any( rttov_errorstatus(:) == errorstatus_fatal ) ) Then
     Do iatm = 1, nprofiles
        If ( rttov_errorstatus(iatm) == errorstatus_fatal ) Then
           write ( ioout, * ) 'rttov_scatt_tl error for profile',iatm
        End If
     End Do
     Stop
  End If

  ! Save radiance as a reference for the trajectory
  ! TL is used instead of rttov_scatt because
  !   calcemis = F and reflectivities have not been saved
  radiance_total_ref(:) = radiance%total(:)

  !---------------------------
  ! second run of TL
  !---------------------------
  lambda = 0.5_JPRB
  allocate ( prof_inc2(     nprofiles ))
  allocate ( cld_prof_inc2( nprofiles ))
  Do j = 1, nprofiles
     prof_inc2(j) % nlevels =  coef_rttov % nlevels
     allocate( prof_inc2(j) % p  ( coef_rttov % nlevels ) )
     allocate( prof_inc2(j) % t  ( coef_rttov % nlevels ) )
     allocate( prof_inc2(j) % q  ( coef_rttov % nlevels ) )
     allocate( prof_inc2(j) % o3 ( coef_rttov % nlevels ) )
     allocate( prof_inc2(j) % clw( coef_rttov % nlevels ) )

     prof_inc2(j) % ozone_Data = .False.  ! no meaning
     prof_inc2(j) % co2_Data   = .False.  ! no meaning
     prof_inc2(j) % clw_Data   = .False.  ! no meaning
     prof_inc2(j) % zenangle   = -1  ! no meaning
     prof_inc2(j) % azangle    = -1  ! no meaning

     ! increments for atmospheric variables
     prof_inc2(j) % p(:)   = 0._JPRB    ! no tl on pressure levels
     prof_inc2(j) % t(:)   = prof_inc(j) % t(:)  *lambda
     prof_inc2(j) % o3(:)  = prof_inc(j) % o3(:) *lambda
     prof_inc2(j) % clw(:) = prof_inc(j) % clw(:)*lambda
     prof_inc2(j) % q(:)   = prof_inc(j) % q(:)  *lambda

     ! increments for air surface variables
     prof_inc2(j) % s2m % t = prof_inc(j) % s2m % t  *lambda
     prof_inc2(j) % s2m % q = prof_inc(j) % s2m % q  *lambda
     prof_inc2(j) % s2m % o = prof_inc(j) % s2m % o  *lambda
     prof_inc2(j) % s2m % p = prof_inc(j) % s2m % p  *lambda
     prof_inc2(j) % s2m % u = prof_inc(j) % s2m % u  *lambda
     prof_inc2(j) % s2m % v = prof_inc(j) % s2m % v  *lambda

     ! increments for skin variables
     prof_inc2(j) % skin % surftype = -1  ! no meaning
     prof_inc2(j) % skin % t        = prof_inc(j) % skin % t  *lambda
     prof_inc2(j) % skin % fastem(:)= prof_inc(j) % skin % fastem(:) *lambda

     ! increments for cloud variables
     prof_inc2(j) % ctp       = prof_inc(j) % ctp       *lambda
     prof_inc2(j) % cfraction = prof_inc(j) % cfraction *lambda

     ! increments for cloud variables
     cld_prof_inc2(j) % nlevels =  nwp_levels
     allocate( cld_prof_inc2(j) % p  ( nwp_levels ) )
     allocate( cld_prof_inc2(j) % ph ( nwp_levels+1 ) )
     allocate( cld_prof_inc2(j) % t  ( nwp_levels ) )
     allocate( cld_prof_inc2(j) % cc ( nwp_levels ) )
     allocate( cld_prof_inc2(j) % clw( nwp_levels ) )
     allocate( cld_prof_inc2(j) % ciw( nwp_levels ) )
     allocate( cld_prof_inc2(j) % rain( nwp_levels ) )
     allocate( cld_prof_inc2(j) % sp  ( nwp_levels ) )
     cld_prof_inc2(j) % p(:)   = cld_prof_inc(j) % p(:)   *lambda
     cld_prof_inc2(j) % ph(:)  = cld_prof_inc(j) % ph(:)  *lambda
     cld_prof_inc2(j) % t(:)   = cld_prof_inc(j) % t(:)   *lambda
     cld_prof_inc2(j) % cc(:)  = cld_prof_inc(j) % cc(:)  *lambda
     cld_prof_inc2(j) % clw(:) = cld_prof_inc(j) % clw(:) *lambda
     cld_prof_inc2(j) % ciw(:) = cld_prof_inc(j) % ciw(:) *lambda
     cld_prof_inc2(j) % rain(:) = cld_prof_inc(j) % rain(:) *lambda
     cld_prof_inc2(j) % sp  (:) = cld_prof_inc(j) % sp  (:) *lambda
  End Do

  ! emissivity
  allocate( emissivity_inc2( nchannels ))
  emissivity_inc2(:) = emissivity_inc(:) * lambda

  ! allocate radiance results arrays with number of channels
  allocate( radiance_tl2 % clear    ( nchannels ) )
  allocate( radiance_tl2 % cloudy   ( nchannels ) )
  allocate( radiance_tl2 % total    ( nchannels ) )
  allocate( radiance_tl2 % bt       ( nchannels ) )
  allocate( radiance_tl2 % bt_clear ( nchannels ) )
  allocate( radiance_tl2 % out       ( nbtout ) )
  allocate( radiance_tl2 % out_clear ( nbtout ) )
  allocate( radiance_tl2 % total_out ( nbtout ) )
  allocate( radiance_tl2 % clear_out ( nbtout ) )
  allocate( radiance_tl2 % upclear  ( nchannels ) )
  allocate( radiance_tl2 % reflclear( nchannels ) )
  allocate( radiance_tl2 % overcast ( nwp_levels, nchannels ) )

  !---------------------------
  Call rttov_scatt_tl ( &
     & rttov_errorstatus,  &! out
     & nfrequencies,    &! in
     & nchannels,       &! in
     & nbtout,          &! in
     & nprofiles,       &! in
     & channels,        &! in
     & polarisations,   &! in
     & lprofiles,       &! in
     & profiles,        &! in
     & cld_profiles,    &! in
     & coef_rttov,      &! in
     & coef_scatt,      &! in
     & calcemis,        &! in
     & emissivity,      &! inout
     & prof_inc2,       &! in
     & cld_prof_inc2,   &! in
     & emissivity_inc2, &! inout
     & radiance,        &! inout
     & radiance_tl2    ) ! inout

  If ( any( rttov_errorstatus(:) == errorstatus_fatal ) ) Then
     Do iatm = 1, nprofiles
        If ( rttov_errorstatus(iatm) == errorstatus_fatal ) Then
           write ( ioout, * ) 'rttov_scatt_tl error for profile',iatm
        End If
     End Do
     Stop
  End If

  !---------------------------

  do ichan = 1, nchannels
     if( abs(lambda * radiance_tl%clear(ichan) - radiance_tl2%clear(ichan)) > threshold ) then
        write(default_err_unit,*) 'TL test fails for radiance_tl%clear for channel ', ichan
        stop
     endif
     if( abs(lambda * radiance_tl%bt_clear(ichan) - radiance_tl2%bt_clear(ichan)) > threshold ) then
        write(default_err_unit,*) 'TL test fails for radiance_tl%bt_clear for channel ', ichan
        stop
     endif
     if( abs(lambda * radiance_tl%bt(ichan) - radiance_tl2%bt(ichan)) > threshold ) then
        write(default_err_unit,*) 'TL test fails for radiance_tl%bt for channel ', ichan
        stop
     endif
     if( abs(lambda * radiance_tl%total(ichan) - radiance_tl2%total(ichan)) > threshold ) then
        write(default_err_unit,*) 'TL test fails for radiance_tl%total for channel ', ichan
        stop
     endif
  end do


  ! Now run the Taylor test
  !-------------------------

  !Allocate new profiles for direct code
  ! Profiles on RTTOV pressure levels
  allocate( profiles2(nprofiles))
  do j = 1, nprofiles
     ! allocate model profiles atmospheric arrays with model levels dimension
     profiles2(j) % nlevels =  coef_rttov % nlevels
     allocate( profiles2(j) % p  ( coef_rttov % nlevels ) )
     allocate( profiles2(j) % t  ( coef_rttov % nlevels ) )
     allocate( profiles2(j) % q  ( coef_rttov % nlevels ) )
     allocate( profiles2(j) % o3 ( coef_rttov % nlevels ) )
     allocate( profiles2(j) % clw( coef_rttov % nlevels ) )
     profiles2(j) % p(:) = coef_rttov % ref_prfl_p(:)
  end do

  ! Cloud additional profiles
  allocate( cld_profiles2(nprofiles))
  do j = 1, nprofiles
     ! allocate model profiles atmospheric arrays with model levels dimension
     cld_profiles2(j) % nlevels =  nwp_levels
     allocate( cld_profiles2(j) % p  ( nwp_levels ) )
     allocate( cld_profiles2(j) % ph ( nwp_levels+1 ) )
     allocate( cld_profiles2(j) % t  ( nwp_levels ) )
     allocate( cld_profiles2(j) % cc ( nwp_levels ) )
     allocate( cld_profiles2(j) % clw( nwp_levels ) )
     allocate( cld_profiles2(j) % ciw( nwp_levels ) )
     allocate( cld_profiles2(j) % rain( nwp_levels ) )
     allocate( cld_profiles2(j) % sp  ( nwp_levels ) )
  end do

  allocate( emissivity2( nchannels ))

  ! allocate radiance results arrays with number of channels
  allocate( radiance2 % clear    ( nchannels ) )
  allocate( radiance2 % cloudy   ( nchannels ) )
  allocate( radiance2 % total    ( nchannels ) )
  allocate( radiance2 % bt       ( nchannels ) )
  allocate( radiance2 % bt_clear ( nchannels ) )
  allocate( radiance2 % out       ( nbtout ) )
  allocate( radiance2 % out_clear ( nbtout ) )
  allocate( radiance2 % total_out ( nbtout ) )
  allocate( radiance2 % clear_out ( nbtout ) )
  allocate( radiance2 % upclear  ( nchannels ) )
  allocate( radiance2 % dnclear  ( nchannels ) )
  allocate( radiance2 % reflclear( nchannels ) )
  allocate( radiance2 % overcast ( nwp_levels, nchannels ) )
  allocate( radiance2 % freq_used( nchannels ) )

! Goto 1000

  do ichan = 1, nchannels
     write(ioout,*)
     write(ioout,*) '(Profile x channel) no. ',ichan
     write(ioout,*) '         Lambda       Clear Rad      Cloudy Rad' &
                  & //'        Clear Tb       Cloudy Tb'
     do iexp = -10, 0
        lambda = 10**(real(iexp))

        do j = 1, nprofiles
           profiles2(j) % ozone_Data = profiles(j) % ozone_Data
           profiles2(j) % co2_Data   = profiles(j) % co2_Data
           profiles2(j) % clw_Data   = profiles(j) % clw_Data
           profiles2(j) % zenangle   = profiles(j) % zenangle
           profiles2(j) % azangle   = profiles(j) % azangle

           ! increments for atmospheric variables
           profiles2(j) % p(:)   = profiles(j) % p(:)
           profiles2(j) % t(:)   = profiles(j) % t(:)  + prof_inc(j) % t(:)  *lambda
           profiles2(j) % o3(:)  = profiles(j) % o3(:) + prof_inc(j) % o3(:) *lambda
           profiles2(j) % clw(:) = profiles(j) % clw(:)+ prof_inc(j) % clw(:)*lambda
           profiles2(j) % q(:)   = profiles(j) % q(:)  + prof_inc(j) % q(:)  *lambda

           ! increments for air surface variables
           profiles2(j) % s2m % t = profiles(j) % s2m % t + prof_inc(j) % s2m % t  *lambda
           profiles2(j) % s2m % q = profiles(j) % s2m % q + prof_inc(j) % s2m % q  *lambda
           profiles2(j) % s2m % o = profiles(j) % s2m % o + prof_inc(j) % s2m % o  *lambda
           profiles2(j) % s2m % p = profiles(j) % s2m % p + prof_inc(j) % s2m % p  *lambda
           profiles2(j) % s2m % u = profiles(j) % s2m % u + prof_inc(j) % s2m % u  *lambda
           profiles2(j) % s2m % v = profiles(j) % s2m % v + prof_inc(j) % s2m % v  *lambda

           ! increments for skin variables
           profiles2(j) % skin % surftype = profiles(j) % skin % surftype
           profiles2(j) % skin % t        = profiles(j) % skin % t  + prof_inc(j) % skin % t  *lambda
           profiles2(j) % skin % fastem(:)= profiles(j) % skin % fastem(:) + prof_inc(j) % skin % fastem(:) *lambda

           ! increments for cloud variables
           profiles2(j) % ctp       = profiles(j) % ctp       + prof_inc(j) % ctp       *lambda
           profiles2(j) % cfraction = profiles(j) % cfraction + prof_inc(j) % cfraction *lambda

           ! increments for cloud variables
           cld_profiles2(j) % nlevels =  nwp_levels
           cld_profiles2(j) % p(:)   = cld_profiles(j) % p(:)   + cld_prof_inc(j) % p(:)   *lambda
           cld_profiles2(j) % ph(:)  = cld_profiles(j) % ph(:)  + cld_prof_inc(j) % ph(:)  *lambda
           cld_profiles2(j) % t(:)   = cld_profiles(j) % t(:)   + cld_prof_inc(j) % t(:)   *lambda
           cld_profiles2(j) % cc(:)  = cld_profiles(j) % cc(:)  + cld_prof_inc(j) % cc(:)  *lambda
           cld_profiles2(j) % clw(:) = cld_profiles(j) % clw(:) + cld_prof_inc(j) % clw(:) *lambda
           cld_profiles2(j) % ciw(:) = cld_profiles(j) % ciw(:) + cld_prof_inc(j) % ciw(:) *lambda
           cld_profiles2(j) % rain(:) = cld_profiles(j) % rain(:) + cld_prof_inc(j) % rain(:) *lambda
           cld_profiles2(j) % sp  (:) = cld_profiles(j) % sp  (:) + cld_prof_inc(j) % sp  (:) *lambda
        end do
        emissivity2(:) = emissivity(:) + emissivity_inc(:) * lambda

        !---------------------------
        Call rttov_scatt( &
        & rttov_errorstatus,  &! out
        & nfrequencies,  &! in
        & nchannels,     &! in
        & nbtout,        &! in
        & nprofiles,     &! in
        & channels,      &! in
        & polarisations, &! in
        & lprofiles,     &! in
        & profiles2,     &! inout  (to invalid clw absorption)
        & cld_profiles2, &! in
        & coef_rttov,    &! in
        & coef_scatt,    &! in
        & calcemis,      &! in
        & emissivity2,   &! inout
        & transmission,  &! inout
        & radiance2     ) ! inout

        If ( any( rttov_errorstatus(:) == errorstatus_fatal ) ) Then
           Do iatm = 1, nprofiles
              If ( rttov_errorstatus(iatm) == errorstatus_fatal ) Then
                 write ( ioout, * ) 'rttov_scatt error for profile',iatm
              End If
           End Do
           Stop
        End If

         !---------------------------

        ratio(1) = (radiance2 % clear(ichan) - radiance % clear(ichan)) / (lambda * radiance_tl % clear(ichan))
        ratio(2) = (radiance2 % total(ichan) - radiance % total(ichan)) / (lambda * radiance_tl % total(ichan))
        ratio(3) = (radiance2 % bt_clear(ichan) - radiance % bt_clear(ichan)) / (lambda * radiance_tl % bt_clear(ichan))
        ratio(4) = (radiance2 % bt(ichan)    - radiance % bt(ichan))    / (lambda * radiance_tl % bt(ichan))
        write(ioout,'(5f16.10)') lambda, ratio

     End do

  End do

  ! End of TL tests

  !---------------------------------------------------------------------
  ! Test of AD
  !---------------------------------------------------------------------
  1000 continue
  write(ioout,*)
  write(ioout,*) 'Test AD'
  write(ioout,*) '-------'
  write(ioout,*)

  write(ioout,*) '1- Test linearity'
  write(ioout,*)
!
  !Allocate new profiles for AD code
  ! Profiles on RTTOV pressure levels
  allocate( profiles_ad(nprofiles))
  do j = 1, nprofiles
     ! allocate model profiles atmospheric arrays with model levels dimension
     profiles_ad(j) % nlevels =  coef_rttov % nlevels
     allocate( profiles_ad(j) % p  ( coef_rttov % nlevels ) )
     allocate( profiles_ad(j) % t  ( coef_rttov % nlevels ) )
     allocate( profiles_ad(j) % q  ( coef_rttov % nlevels ) )
     allocate( profiles_ad(j) % o3 ( coef_rttov % nlevels ) )
     allocate( profiles_ad(j) % clw( coef_rttov % nlevels ) )
     profiles_ad(j) % p(:) = coef_rttov % ref_prfl_p(:)
  end do

  ! Cloud additional profiles
  allocate( cld_profiles_ad(nprofiles))
  do j = 1, nprofiles
     ! allocate model profiles atmospheric arrays with model levels dimension
     cld_profiles_ad(j) % nlevels =  nwp_levels
     allocate( cld_profiles_ad(j) % p  ( nwp_levels ) )
     allocate( cld_profiles_ad(j) % ph ( nwp_levels+1 ) )
     allocate( cld_profiles_ad(j) % t  ( nwp_levels ) )
     allocate( cld_profiles_ad(j) % cc ( nwp_levels ) )
     allocate( cld_profiles_ad(j) % clw( nwp_levels ) )
     allocate( cld_profiles_ad(j) % ciw( nwp_levels ) )
     allocate( cld_profiles_ad(j) % rain( nwp_levels ) )
     allocate( cld_profiles_ad(j) % sp  ( nwp_levels ) )
  end do
  Do j = 1, nprofiles
     profiles_ad(j) % ozone_Data = .False.  ! no meaning
     profiles_ad(j) % co2_Data   = .False.  ! no meaning
     profiles_ad(j) % clw_Data   = .False.  ! no meaning
     profiles_ad(j) % zenangle   = -1       ! no meaning
     profiles_ad(j) % azangle    = -1       ! no meaning

     ! increments for atmospheric variables
     profiles_ad(j) % p(:)   = 0._JPRB ! no AD on pressure levels
     profiles_ad(j) % t(:)   = 0._JPRB ! temperarure
     profiles_ad(j) % o3(:)  = 0._JPRB ! O3 ppmv
     profiles_ad(j) % clw(:) = 0._JPRB ! clw
     profiles_ad(j) % q(:)   = 0._JPRB ! WV

     ! increments for air surface variables
     profiles_ad(j) % s2m % t = 0._JPRB!  temperarure
     profiles_ad(j) % s2m % q = 0 !  WV
     profiles_ad(j) % s2m % o = 0 !  O3
     profiles_ad(j) % s2m % p = 0._JPRB!  pressure
     profiles_ad(j) % s2m % u = 0._JPRB!  wind components
     profiles_ad(j) % s2m % v = 0._JPRB!  wind components

     ! increments for skin variables
     profiles_ad(j) % skin % surftype = -1  ! no meaning
     profiles_ad(j) % skin % t        = 0._JPRB  ! on temperarure
     profiles_ad(j) % skin % fastem   = 0._JPRB

     ! increments for cloud variables
     profiles_ad(j) % ctp       = 0._JPRB  ! pressure
     profiles_ad(j) % cfraction = 0._JPRB  ! cloud fraction

     ! Cloud profiles
     cld_profiles_ad(j) % p  (:) = 0._JPRB
     cld_profiles_ad(j) % ph (:) = 0._JPRB
     cld_profiles_ad(j) % t  (:) = 0._JPRB
     cld_profiles_ad(j) % cc (:) = 0._JPRB
     cld_profiles_ad(j) % clw(:) = 0._JPRB
     cld_profiles_ad(j) % ciw(:) = 0._JPRB
     cld_profiles_ad(j) % rain(:) = 0._JPRB
     cld_profiles_ad(j) % sp  (:) = 0._JPRB
  End Do

  allocate( emissivity_ad( nchannels ))
  emissivity_ad(:) = 0._JPRB

  ! Set perturbations
  !
  ! allocate radiance results arrays with number of channels
  allocate( radiance_inc % clear    ( nchannels ) )
  allocate( radiance_inc % cloudy   ( nchannels ) )
  allocate( radiance_inc % total    ( nchannels ) )
  allocate( radiance_inc % bt       ( nchannels ) )
  allocate( radiance_inc % bt_clear ( nchannels ) )
  allocate( radiance_inc % out       ( nbtout ) )
  allocate( radiance_inc % out_clear ( nbtout ) )
  allocate( radiance_inc % total_out ( nbtout ) )
  allocate( radiance_inc % clear_out ( nbtout ) )
  allocate( radiance_inc % upclear  ( nchannels ) )
  allocate( radiance_inc % reflclear( nchannels ) )
  allocate( radiance_inc % overcast ( nwp_levels, nchannels ) )
  if (kinrad == 2) then
     radiance_inc % clear_out(:)   = 0._JPRB
     radiance_inc % total_out(:)   = 0._JPRB
     radiance_inc % out_clear(:) = 0.05_JPRB * radiance % out_clear(:)
     radiance_inc % out(:)       = 0.05_JPRB * radiance % out(:)
  else
     radiance_inc % clear_out(:)    = 0.05_JPRB * radiance % clear_out(:)
     radiance_inc % total_out(:)    = 0.05_JPRB * radiance % total_out(:)
     radiance_inc % out_clear(:) = 0._JPRB
     radiance_inc % out(:)       = 0._JPRB
  endif
  radiance_inc % clear(:)    = 0._JPRB   ! AD does not work for radiance_inc % clear(:) because of switchrad in RTTOV
  radiance_inc % cloudy   (:)   = 0._JPRB
  radiance_inc % upclear  (:)   = 0._JPRB
  radiance_inc % reflclear(:)   = 0._JPRB
  radiance_inc % overcast (:,:) = 0._JPRB
  radiance_inc % bt       (:)   = 0._JPRB
  radiance_inc % bt_clear (:)   = 0._JPRB
  radiance_inc % total    (:)   = 0._JPRB


  !---------------------------
  Call Rttov_scatt_ad ( &
     & rttov_errorstatus,  &! out
     & nfrequencies,    &! in
     & nchannels,       &! in
     & nbtout,          &! in
     & nprofiles,       &! in
     & channels,        &! in
     & polarisations,   &! in
     & lprofiles,       &! in
     & profiles,        &! in
     & cld_profiles,    &! in
     & coef_rttov,      &! in
     & coef_scatt,      &! in
     & calcemis,        &! in
     & emissivity,      &! inout
     & profiles_ad,     &! inout
     & cld_profiles_ad, &! inout
     & emissivity_ad,   &! inout
     & radiance2,       &! inout
     & radiance_inc    ) ! inout

  If ( any( rttov_errorstatus(:) == errorstatus_fatal ) ) Then
     Do iatm = 1, nprofiles
        If ( rttov_errorstatus(iatm) == errorstatus_fatal ) Then
           write ( ioout, * ) 'rttov_scatt_ad error for profile',iatm
        End If
     End Do
     Stop
  End If

  If ( Any( abs(radiance_total_ref(:) - radiance2%total(:)) > eps * radiance_total_ref(:)  ))  Then
    write(default_err_unit,*) 'wrong forward model in AD'
    write(default_err_unit,*) radiance_total_ref(:)
    write(default_err_unit,*) abs(radiance_total_ref(:)-radiance2%total(:)) / ( eps * radiance_total_ref(:))
    Stop
  Endif


  !---------------------------
  ! Second run of AD

  !Allocate new profiles for AD code
  ! Profiles on RTTOV pressure levels
  allocate( profiles_ad2(nprofiles))
  do j = 1, nprofiles
     ! allocate model profiles atmospheric arrays with model levels dimension
     profiles_ad2(j) % nlevels =  coef_rttov % nlevels
     allocate( profiles_ad2(j) % p  ( coef_rttov % nlevels ) )
     allocate( profiles_ad2(j) % t  ( coef_rttov % nlevels ) )
     allocate( profiles_ad2(j) % q  ( coef_rttov % nlevels ) )
     allocate( profiles_ad2(j) % o3 ( coef_rttov % nlevels ) )
     allocate( profiles_ad2(j) % clw( coef_rttov % nlevels ) )
     profiles_ad2(j) % p(:) = coef_rttov % ref_prfl_p(:)
  end do

  ! Cloud additional profiles
  allocate( cld_profiles_ad2(nprofiles))
  do j = 1, nprofiles
     ! allocate model profiles atmospheric arrays with model levels dimension
     cld_profiles_ad2(j) % nlevels =  nwp_levels
     allocate( cld_profiles_ad2(j) % p  ( nwp_levels ) )
     allocate( cld_profiles_ad2(j) % ph ( nwp_levels+1 ) )
     allocate( cld_profiles_ad2(j) % t  ( nwp_levels ) )
     allocate( cld_profiles_ad2(j) % cc ( nwp_levels ) )
     allocate( cld_profiles_ad2(j) % clw( nwp_levels ) )
     allocate( cld_profiles_ad2(j) % ciw( nwp_levels ) )
     allocate( cld_profiles_ad2(j) % rain( nwp_levels ) )
     allocate( cld_profiles_ad2(j) % sp  ( nwp_levels ) )
  end do

  Do j = 1, nprofiles
     profiles_ad2(j) % ozone_Data = .False.  ! no meaning
     profiles_ad2(j) % co2_Data   = .False.  ! no meaning
     profiles_ad2(j) % clw_Data   = .False.  ! no meaning
     profiles_ad2(j) % zenangle   = -1       ! no meaning
     profiles_ad2(j) % azangle    = -1       ! no meaning

     ! increments for atmospheric variables
     profiles_ad2(j) % p(:)   = 0._JPRB ! no AD on pressure levels
     profiles_ad2(j) % t(:)   = 0._JPRB ! temperarure
     profiles_ad2(j) % o3(:)  = 0._JPRB ! O3 ppmv
     profiles_ad2(j) % clw(:) = 0._JPRB ! clw
     profiles_ad2(j) % q(:)   = 0._JPRB ! WV

     ! increments for air surface variables
     profiles_ad2(j) % s2m % t = 0._JPRB!  temperarure
     profiles_ad2(j) % s2m % q = 0 !  WV
     profiles_ad2(j) % s2m % o = 0 !  O3
     profiles_ad2(j) % s2m % p = 0._JPRB!  pressure
     profiles_ad2(j) % s2m % u = 0._JPRB!  wind components
     profiles_ad2(j) % s2m % v = 0._JPRB!  wind components

     ! increments for skin variables
     profiles_ad2(j) % skin % surftype = -1  ! no meaning
     profiles_ad2(j) % skin % t        = 0._JPRB  ! on temperarure
     profiles_ad2(j) % skin % fastem   = 0._JPRB

     ! increments for cloud variables
     profiles_ad2(j) % ctp       = 0._JPRB  ! pressure
     profiles_ad2(j) % cfraction = 0._JPRB  ! cloud fraction

     ! Cloud profiles
     cld_profiles_ad2(j) % p  (:) = 0._JPRB
     cld_profiles_ad2(j) % ph (:) = 0._JPRB
     cld_profiles_ad2(j) % t  (:) = 0._JPRB
     cld_profiles_ad2(j) % cc (:) = 0._JPRB
     cld_profiles_ad2(j) % clw(:) = 0._JPRB
     cld_profiles_ad2(j) % ciw(:) = 0._JPRB
     cld_profiles_ad2(j) % rain(:) = 0._JPRB
     cld_profiles_ad2(j) % sp  (:) = 0._JPRB
  End Do

  allocate( emissivity_ad2( nchannels ))
  emissivity_ad2(:) = 0._JPRB

  ! allocate radiance results arrays with number of channels
  allocate( radiance_inc2 % clear    ( nchannels ) )
  allocate( radiance_inc2 % cloudy   ( nchannels ) )
  allocate( radiance_inc2 % total    ( nchannels ) )
  allocate( radiance_inc2 % bt       ( nchannels ) )
  allocate( radiance_inc2 % bt_clear ( nchannels ) )
  allocate( radiance_inc2 % out      ( nbtout ) )
  allocate( radiance_inc2 % out_clear( nbtout ) )
  allocate( radiance_inc2 % total_out( nbtout ) )
  allocate( radiance_inc2 % clear_out( nbtout ) )
  allocate( radiance_inc2 % upclear  ( nchannels ) )
  allocate( radiance_inc2 % reflclear( nchannels ) )
  allocate( radiance_inc2 % overcast ( nwp_levels, nchannels ) )

  lambda = 0.5_JPRB
  if (kinrad == 2) then
     radiance_inc2 % clear_out(:)   = 0._JPRB
     radiance_inc2 % total_out(:)   = 0._JPRB
     radiance_inc2 % out_clear(:) = 0.05_JPRB * radiance % out_clear(:)* lambda
     radiance_inc2 % out(:)       = 0.05_JPRB * radiance % out(:)* lambda
  else
     radiance_inc2 % clear_out(:)    = 0.05_JPRB * radiance % clear_out(:) * lambda
     radiance_inc2 % total_out(:)    = 0.05_JPRB * radiance % total_out(:) * lambda
     radiance_inc2 % out_clear(:) = 0._JPRB
     radiance_inc2 % out(:)       = 0._JPRB
  endif
  radiance_inc2 % clear(:)    = 0._JPRB   ! AD does not work for radiance_inc % clear(:) because of switchrad in RTTOV
  radiance_inc2 % cloudy   (:)   = 0._JPRB
  radiance_inc2 % upclear  (:)   = 0._JPRB
  radiance_inc2 % reflclear(:)   = 0._JPRB
  radiance_inc2 % overcast (:,:) = 0._JPRB
  radiance_inc2 % bt       (:)   = 0._JPRB
  radiance_inc2 % bt_clear (:)   = 0._JPRB
  radiance_inc2 % total    (:)   = 0._JPRB

  !---------------------------
  Call Rttov_scatt_ad ( &
     & rttov_errorstatus,  &! out
     & nfrequencies,     &! in
     & nchannels,        &! in
     & nbtout,           &! in
     & nprofiles,        &! in
     & channels,         &! in
     & polarisations,    &! in
     & lprofiles,        &! in
     & profiles,         &! in
     & cld_profiles,     &! in
     & coef_rttov,       &! in
     & coef_scatt,       &! in
     & calcemis,         &! in
     & emissivity,       &! inout
     & profiles_ad2,     &! inout
     & cld_profiles_ad2, &! inout
     & emissivity_ad2,   &! inout
     & radiance2,        &! inout
     & radiance_inc2    ) ! inout

  If ( any( rttov_errorstatus(:) == errorstatus_fatal ) ) Then
     Do iatm = 1, nprofiles
        If ( rttov_errorstatus(iatm) == errorstatus_fatal ) Then
           write ( ioout, * ) 'rttov_scatt_ad error for profile',iatm
        End If
     End Do
     Stop
  End If


  do j = 1, nprofiles
     do lev = 1, profiles_ad(j) % nlevels
        if ( abs(lambda * profiles_ad(j)%t(lev) - profiles_ad2(j)%t(lev)) > threshold ) then
           write(default_err_unit,*) 'test AD 1 fails',lev
           stop
        End If
        if ( abs(lambda * profiles_ad(j)%q(lev) - profiles_ad2(j)%q(lev)) > threshold ) Then
           write(default_err_unit,*) 'test AD 2 fails',lev
           stop
        End If
        if ( abs(lambda * profiles_ad(j)%o3(lev) - profiles_ad2(j)%o3(lev)) > threshold ) Then
           write(default_err_unit,*) 'test AD 3 fails',lev
           stop
        End If
     enddo
  enddo

  do j = 1, nprofiles
     do lev = 1, cld_profiles_ad(j) % nlevels
        if ( abs(lambda * cld_profiles_ad(j)%p(lev)   - cld_profiles_ad2(j)%p(lev)) > threshold ) Then
           write(default_err_unit,*) 'test AD 4 fails',lev
           stop
        End If
        if ( abs(lambda * cld_profiles_ad(j)%ph(lev)  - cld_profiles_ad2(j)%ph(lev)) > threshold ) Then
           write(default_err_unit,*) 'test AD 5 fails',lev
           stop
        End If
        if ( abs(lambda * cld_profiles_ad(j)%t(lev)   - cld_profiles_ad2(j)%t(lev)) > threshold ) Then
           write(default_err_unit,*) 'test AD 6 fails',lev
           stop
        End If
        if ( abs(lambda * cld_profiles_ad(j)%cc(lev)  - cld_profiles_ad2(j)%cc(lev)) > threshold ) Then
           write(default_err_unit,*) 'test AD 7 fails',lev
           stop
        End If
        if ( abs(lambda * cld_profiles_ad(j)%clw(lev) - cld_profiles_ad2(j)%clw(lev)) > threshold ) Then
           write(default_err_unit,*) 'test AD 8 fails',lev
           stop
        End If
        if ( abs(lambda * cld_profiles_ad(j)%ciw(lev) - cld_profiles_ad2(j)%ciw(lev)) > threshold ) Then
           write(default_err_unit,*) 'test AD 9 fails',lev
           stop
        End If
        if ( abs(lambda * cld_profiles_ad(j)%rain(lev) - cld_profiles_ad2(j)%rain(lev)) > threshold ) Then
           write(default_err_unit,*) 'test AD 10 fails',lev
           stop
        End If
        if ( abs(lambda * cld_profiles_ad(j)%sp  (lev) - cld_profiles_ad2(j)%sp  (lev)) > threshold ) Then
           write(default_err_unit,*) 'test AD 11 fails',lev
           stop
        End If
     enddo
     lev = cld_profiles_ad(j) % nlevels+1
     if ( abs(lambda * cld_profiles_ad(j)%ph(lev)  - cld_profiles_ad2(j)%ph(lev)) > threshold ) Then
        write(default_err_unit,*) 'test AD 12 fails',lev
        stop
     End If
  enddo

  do j = 1, nprofiles
     if ( abs(lambda * profiles_ad(j)%s2m%t - profiles_ad2(j)%s2m%t) > threshold ) Then
        write(default_err_unit,*) 'test AD 13 fails',j
        stop
     End If
     if ( abs(lambda * profiles_ad(j)%s2m%q - profiles_ad2(j)%s2m%q) > threshold ) Then
        write(default_err_unit,*) 'test AD 14 fails',j
        stop
     End If
     if ( abs(lambda * profiles_ad(j)%s2m%p - profiles_ad2(j)%s2m%p) > threshold ) Then
        write(default_err_unit,*) 'test AD 15 fails',j
        stop
     End If
     if ( abs(lambda * profiles_ad(j)%s2m%u - profiles_ad2(j)%s2m%u) > threshold ) Then
        write(default_err_unit,*) 'test AD 16 fails',j
        stop
     End If
     if ( abs(lambda * profiles_ad(j)%s2m%v - profiles_ad2(j)%s2m%v) > threshold ) Then
        write(default_err_unit,*) 'test AD 17 fails',j
        stop
     End If

     if ( abs(lambda * profiles_ad(j)%skin%t - profiles_ad2(j)%skin%t) > threshold ) Then
        write(default_err_unit,*) 'test AD 18 fails',j
        stop
     End If
  enddo

  do j = 1, nchannels
     if ( abs(lambda * emissivity_ad(j) - emissivity_ad2(j)) > threshold ) Then
        write(default_err_unit,*) 'test AD 19 fails',j
        stop
     End If
  enddo

  write(ioout,*) '2- Test equality of norms'
  write(ioout,*)
!
! Set perturbations
!

  ! Set perturbations = 5% of initial profile
  Do j = 1, nprofiles
     prof_inc(j) % nlevels =  coef_rttov % nlevels

     prof_inc(j) % ozone_Data = .False.  ! no meaning
     prof_inc(j) % co2_Data   = .False.  ! no meaning
     prof_inc(j) % clw_Data   = .False.  ! no meaning
     prof_inc(j) % zenangle   = -1  ! no meaning
     prof_inc(j) % azangle    = -1  ! no meaning

     ! increments for atmospheric variables
     prof_inc(j) % p(:)   = 0._JPRB    ! no tl on pressure levels
     prof_inc(j) % t(:)   = profiles(j) % t(:)  *0.05_JPRB
     prof_inc(j) % o3(:)  = profiles(j) % o3(:) *0.05_JPRB
     prof_inc(j) % clw(:) = profiles(j) % clw(:)*0.05_JPRB
     prof_inc(j) % q(:)   = profiles(j) % q(:)  *0.05_JPRB

     ! increments for air surface variables
     prof_inc(j) % s2m % t = profiles(j) % s2m % t  *0.05_JPRB
     prof_inc(j) % s2m % q = profiles(j) % s2m % q  *0.05_JPRB
     prof_inc(j) % s2m % o = profiles(j) % s2m % o  *0.05_JPRB
     prof_inc(j) % s2m % p = profiles(j) % s2m % p  *0.05_JPRB
     prof_inc(j) % s2m % u = profiles(j) % s2m % u  *0.05_JPRB
     prof_inc(j) % s2m % v = profiles(j) % s2m % v  *0.05_JPRB

     ! increments for skin variables
     prof_inc(j) % skin % surftype = -1  ! no meaning
     prof_inc(j) % skin % t        = profiles(j) % skin % t  *0.05_JPRB
     prof_inc(j) % skin % fastem(:)= profiles(j) % skin % fastem(:) *0.05_JPRB

     ! increments for cloud variables
     prof_inc(j) % ctp       = profiles(j) % ctp       *0.05_JPRB
     prof_inc(j) % cfraction = profiles(j) % cfraction *0.05_JPRB

     ! increments for cloud variables
     cld_prof_inc(j) % nlevels =  nwp_levels
     cld_prof_inc(j) % p(:)   = cld_profiles(j) % p(:)   *0.05_JPRB
     cld_prof_inc(j) % ph(:)  = cld_profiles(j) % ph(:)  *0.05_JPRB
     cld_prof_inc(j) % t(:)   = cld_profiles(j) % t(:)   *0.05_JPRB
     cld_prof_inc(j) % cc(:)  = cld_profiles(j) % cc(:)  *0.05_JPRB
     cld_prof_inc(j) % clw(:) = cld_profiles(j) % clw(:) *0.05_JPRB
     cld_prof_inc(j) % ciw(:) = cld_profiles(j) % ciw(:) *0.05_JPRB
     cld_prof_inc(j) % rain(:) = cld_profiles(j) % rain(:) *0.05_JPRB
     cld_prof_inc(j) % sp  (:) = cld_profiles(j) % sp  (:) *0.05_JPRB
  End Do

  ! emissivity
  emissivity_inc(:) = emissivity(:) * 0.05_JPRB
  emissivity_inc(:) = 0._JPRB

  radiance_tl % total(:)    = 0._JPRB
  radiance_tl % bt_clear(:) = 0._JPRB
  radiance_tl % bt(:)       = 0._JPRB

  !---------------------------
  Call rttov_scatt_tl ( &
     & rttov_errorstatus,  &! out
     & nfrequencies,    &! in
     & nchannels,       &! in
     & nbtout,          &! in
     & nprofiles,       &! in
     & channels,        &! in
     & polarisations,   &! in
     & lprofiles,       &! in
     & profiles,        &! in
     & cld_profiles,    &! in
     & coef_rttov,      &! in
     & coef_scatt,      &! in
     & calcemis,        &! in
     & emissivity,      &! inout
     & prof_inc,        &! in
     & cld_prof_inc,    &! in
     & emissivity_inc,  &! inout
     & radiance,        &! inout
     & radiance_tl     ) ! inout


  If ( any( rttov_errorstatus(:) == errorstatus_fatal ) ) Then
     Do iatm = 1, nprofiles
        If ( rttov_errorstatus(iatm) == errorstatus_fatal ) Then
           write ( ioout, * ) 'rttov_scatt_tl error for profile',iatm
        End If
     End Do
     Stop
  End If
  if (kinrad == 2) then
     radiance_tl % clear_out(:)    = 0._JPRB
     radiance_tl % total_out(:)    = 0._JPRB
  else
     radiance_tl % out_clear(:) = 0._JPRB
     radiance_tl % out(:)       = 0._JPRB
  endif
!  radiance_tl % clear_out(:)    = 0._JPRB   ! AD does not work for radiance_inc % clear(:) because of switchrad in RTTOV

  !* compute <subtl(delta_x),delta_z>
  zdelta1 = 0._JPRB
  do j = 1, nbtout
    zdelta1 = zdelta1 + radiance_tl % total_out(j)**2 + radiance_tl % clear_out(j)**2
    zdelta1 = zdelta1 + radiance_tl % out(j)**2 + radiance_tl % out_clear(j)**2
  enddo
  write(ioout,fmt='('' delta1 = '',2e24.17)') zdelta1

  !---------------------------
  ! Now run AD code with TL radiances in input
  Do j = 1, nprofiles
     profiles_ad(j) % ozone_Data = .False.  ! no meaning
     profiles_ad(j) % co2_Data   = .False.  ! no meaning
     profiles_ad(j) % clw_Data   = .False.  ! no meaning
     profiles_ad(j) % zenangle   = -1       ! no meaning
     profiles_ad(j) % azangle    = -1       ! no meaning

     ! increments for atmospheric variables
     profiles_ad(j) % p(:)   = 0._JPRB ! no AD on pressure levels
     profiles_ad(j) % t(:)   = 0._JPRB ! temperarure
     profiles_ad(j) % o3(:)  = 0._JPRB ! O3 ppmv
     profiles_ad(j) % clw(:) = 0._JPRB ! clw
     profiles_ad(j) % q(:)   = 0._JPRB ! WV

     ! increments for air surface variables
     profiles_ad(j) % s2m % t = 0._JPRB!  temperarure
     profiles_ad(j) % s2m % q = 0._JPRB !  WV
     profiles_ad(j) % s2m % o = 0._JPRB !  O3
     profiles_ad(j) % s2m % p = 0._JPRB!  pressure
     profiles_ad(j) % s2m % u = 0._JPRB!  wind components
     profiles_ad(j) % s2m % v = 0._JPRB!  wind components

     ! increments for skin variables
     profiles_ad(j) % skin % surftype = -1  ! no meaning
     profiles_ad(j) % skin % t        = 0._JPRB  ! on temperarure
     profiles_ad(j) % skin % fastem   = 0._JPRB

     ! increments for cloud variables
     profiles_ad(j) % ctp       = 0._JPRB  ! pressure
     profiles_ad(j) % cfraction = 0._JPRB  ! cloud fraction

     ! Cloud profiles
     cld_profiles_ad(j) % p  (:) = 0._JPRB
     cld_profiles_ad(j) % ph (:) = 0._JPRB
     cld_profiles_ad(j) % t  (:) = 0._JPRB
     cld_profiles_ad(j) % cc (:) = 0._JPRB
     cld_profiles_ad(j) % clw(:) = 0._JPRB
     cld_profiles_ad(j) % ciw(:) = 0._JPRB
     cld_profiles_ad(j) % rain(:) = 0._JPRB
     cld_profiles_ad(j) % sp  (:) = 0._JPRB
  End Do

  emissivity_ad(:) = 0._JPRB

  ! move TL results to AD radiance increments
  if (kinrad == 2) then
     radiance_inc % clear_out(:)   = 0._JPRB
     radiance_inc % total_out(:)   = 0._JPRB
     radiance_inc % out_clear(:) = radiance_tl % out_clear(:)
     radiance_inc % out(:)       = radiance_tl % out(:)
  else
     radiance_inc % clear_out(:)    = radiance_tl % clear_out(:)
     radiance_inc % total_out(:)    = radiance_tl % total_out(:)
     radiance_inc % out_clear(:) = 0._JPRB
     radiance_inc % out(:)       = 0._JPRB
  endif
  radiance_inc % clear(:)    = 0._JPRB   ! AD does not work for radiance_inc % clear(:) because of switchrad in RTTOV
  radiance_inc % cloudy   (:)   = 0._JPRB
  radiance_inc % upclear  (:)   = 0._JPRB
  radiance_inc % reflclear(:)   = 0._JPRB
  radiance_inc % overcast (:,:) = 0._JPRB
  radiance_inc % bt       (:)   = 0._JPRB
  radiance_inc % bt_clear (:)   = 0._JPRB
  radiance_inc % total    (:)   = 0._JPRB

  Call Rttov_scatt_ad ( &
     & rttov_errorstatus,  &! out
     & nfrequencies,    &! in
     & nchannels,       &! in
     & nbtout,          &! in
     & nprofiles,       &! in
     & channels,        &! in
     & polarisations,   &! in
     & lprofiles,       &! in
     & profiles,        &! in
     & cld_profiles,    &! in
     & coef_rttov,      &! in
     & coef_scatt,      &! in
     & calcemis,        &! in
     & emissivity,      &! inout
     & profiles_ad,     &! inout
     & cld_profiles_ad, &! inout
     & emissivity_ad,   &! inout
     & radiance2,       &! inout
     & radiance_inc    ) ! inout

  If ( any( rttov_errorstatus(:) == errorstatus_fatal ) ) Then
     Do iatm = 1, nprofiles
        If ( rttov_errorstatus(iatm) == errorstatus_fatal ) Then
           write ( ioout, * ) 'rttov_scatt_ad error for profile',iatm
        End If
     End Do
     Stop
  End If

  !* compute <delta_x,subad(delta_z)>
  zdelta2 = 0._JPRB
  do j = 1, nprofiles
     do lev = 1, prof_inc(j) % nlevels
        zdelta2 = zdelta2 + &
              & prof_inc(j)%t(lev)  * profiles_ad(j)%t(lev)   +&
              & prof_inc(j)%q(lev)  * profiles_ad(j)%q(lev)   +&
              & prof_inc(j)%o3(lev) * profiles_ad(j)%o3(lev)
    enddo
  enddo

  do j = 1, nprofiles
     do lev = 1, cld_prof_inc(j) % nlevels
        zdelta2 = zdelta2 + &
              & cld_prof_inc(j)%p(lev)* cld_profiles_ad(j)%p(lev)   +&
              & cld_prof_inc(j)%ph(lev)* cld_profiles_ad(j)%ph(lev)   +&
              & cld_prof_inc(j)%t(lev)* cld_profiles_ad(j)%t(lev)   +&
              & cld_prof_inc(j)%cc(lev)* cld_profiles_ad(j)%cc(lev)   +&
              & cld_prof_inc(j)%clw(lev)* cld_profiles_ad(j)%clw(lev)   +&
              & cld_prof_inc(j)%ciw(lev)* cld_profiles_ad(j)%ciw(lev)   +&
              & cld_prof_inc(j)%rain(lev)* cld_profiles_ad(j)%rain(lev)   +&
              & cld_prof_inc(j)%sp  (lev)* cld_profiles_ad(j)%sp  (lev)
    enddo
    lev = cld_prof_inc(j) % nlevels+1
        zdelta2 = zdelta2 + &
              & cld_prof_inc(j)%ph(lev) * cld_profiles_ad(j)%ph(lev)
  enddo

  do j = 1, nprofiles
     zdelta2 = zdelta2 + &
           & prof_inc(j)%s2m%t * profiles_ad(j)%s2m%t  + &
           & prof_inc(j)%s2m%q * profiles_ad(j)%s2m%q + &
           & prof_inc(j)%s2m%o * profiles_ad(j)%s2m%o + &
           & prof_inc(j)%s2m%p * profiles_ad(j)%s2m%p  + &
           & prof_inc(j)%s2m%u * profiles_ad(j)%s2m%u  + &
           & prof_inc(j)%s2m%v * profiles_ad(j)%s2m%v  + &
           & prof_inc(j)%skin%t * profiles_ad(j)%skin%t
  enddo

  do j = 1, nchannels
     zdelta2 = zdelta2 + &
           & emissivity_inc(j) * emissivity_ad(j)
  enddo
  write(ioout,fmt='('' delta2 = '',2e24.17)') zdelta2

  if (zdelta2 == 0._JPRB) then
    z = 1._JPRB
  else
    z = zdelta2
  endif

  write (ioout,fmt= &
   & '('' The difference is '',f9.3, '' times the zero of the machine '')') &
 & abs(zdelta2-zdelta1)/eps/z

  !---------------------------------------------------------------------
  ! Test of K
  !---------------------------------------------------------------------

  write(ioout,*)
  write(ioout,*) 'Test K'
  write(ioout,*) '------'
  write(ioout,*)

  !Allocate new profiles for K code
  ! Profiles on RTTOV pressure levels
  allocate( profiles_k(nchannels))
  do j = 1, nchannels
     ! allocate model profiles atmospheric arrays with model levels dimension
     profiles_k(j) % nlevels =  coef_rttov % nlevels
     allocate( profiles_k(j) % p  ( coef_rttov % nlevels ) )
     allocate( profiles_k(j) % t  ( coef_rttov % nlevels ) )
     allocate( profiles_k(j) % q  ( coef_rttov % nlevels ) )
     allocate( profiles_k(j) % o3 ( coef_rttov % nlevels ) )
     allocate( profiles_k(j) % clw( coef_rttov % nlevels ) )
     profiles_k(j) % p(:) = coef_rttov % ref_prfl_p(:)
  end do

  ! Cloud additional profiles
  allocate( cld_profiles_k(nchannels))
  do j = 1, nchannels
     ! allocate model profiles atmospheric arrays with model levels dimension
     cld_profiles_k(j) % nlevels =  nwp_levels
     allocate( cld_profiles_k(j) % p  ( nwp_levels ) )
     allocate( cld_profiles_k(j) % ph ( nwp_levels+1 ) )
     allocate( cld_profiles_k(j) % t  ( nwp_levels ) )
     allocate( cld_profiles_k(j) % cc ( nwp_levels ) )
     allocate( cld_profiles_k(j) % clw( nwp_levels ) )
     allocate( cld_profiles_k(j) % ciw( nwp_levels ) )
     allocate( cld_profiles_k(j) %rain( nwp_levels ) )
     allocate( cld_profiles_k(j) % sp ( nwp_levels ) )
  end do
  allocate( emissivity_k( nchannels ))

  Call Rttov_scatt_k    ( &
     & rttov_errorstatus,  &! out
     & nfrequencies,       &! in
     & nchannels,          &! in
     & nbtout,             &! in
     & nprofiles,          &! in
     & channels,           &! in
     & polarisations,      &! in
     & lprofiles,          &! in
     & profiles,           &! in
     & cld_profiles,       &! in
     & coef_rttov,         &! in
     & coef_scatt,         &! in
     & switchrad,          &! in
     & calcemis,           &! in
     & emissivity,         &! inout
     & profiles_k ,        &! inout
     & cld_profiles_k,     &! inout
     & emissivity_k,       &! inout
     & radiance)            ! inout

  If ( Any( abs(radiance_total_ref(:) - radiance%total(:)) > eps * radiance_total_ref(:)  ))  Then
    write(default_err_unit,*) 'wrong forward model in K'
    write(default_err_unit,*) radiance_total_ref(:)
    write(default_err_unit,*) abs(radiance_total_ref(:)-radiance%total(:)) / ( eps * radiance_total_ref(:))
    Stop
  Endif

  !---------------------------
  ! Compares K to AD
  ! Actually Rttov_scatt_k uses AD, but in an economical way
  !   the test here checks that values are correctly located in the matrix
  !   using the simplest (expensive) approach
  jchan = 0
  Do ichan = 1, nfrequencies

  Do j = 1, nprofiles
     ! increments for atmospheric variables
     profiles_ad(j) % p(:)   = 0._JPRB ! no AD on pressure levels
     profiles_ad(j) % t(:)   = 0._JPRB ! temperarure
     profiles_ad(j) % o3(:)  = 0._JPRB ! O3 ppmv
     profiles_ad(j) % clw(:) = 0._JPRB ! clw
     profiles_ad(j) % q(:)   = 0._JPRB ! WV

     ! increments for air surface variables
     profiles_ad(j) % s2m % t = 0._JPRB!  temperarure
     profiles_ad(j) % s2m % q = 0 !  WV
     profiles_ad(j) % s2m % o = 0 !  WV
     profiles_ad(j) % s2m % p = 0._JPRB!  pressure
     profiles_ad(j) % s2m % u = 0._JPRB!  wind components
     profiles_ad(j) % s2m % v = 0._JPRB!  wind components

     ! increments for skin variables
     profiles_ad(j) % skin % surftype = -1  ! no meaning
     profiles_ad(j) % skin % t        = 0._JPRB  ! on temperarure
     profiles_ad(j) % skin % fastem   = 0._JPRB

     ! increments for cloud variables
     profiles_ad(j) % ctp       = 0._JPRB  ! pressure
     profiles_ad(j) % cfraction = 0._JPRB  ! cloud fraction

     ! Cloud profiles
     cld_profiles_ad(j) % p  (:) = 0._JPRB
     cld_profiles_ad(j) % ph (:) = 0._JPRB
     cld_profiles_ad(j) % t  (:) = 0._JPRB
     cld_profiles_ad(j) % cc (:) = 0._JPRB
     cld_profiles_ad(j) % clw(:) = 0._JPRB
     cld_profiles_ad(j) % ciw(:) = 0._JPRB
     cld_profiles_ad(j) %rain(:) = 0._JPRB
     cld_profiles_ad(j) % sp (:) = 0._JPRB
  End Do

  emissivity_ad(:) = 0._JPRB

  radiance_inc % cloudy   (:)   = 0._JPRB
  radiance_inc % upclear  (:)   = 0._JPRB
  radiance_inc % reflclear(:)   = 0._JPRB
  radiance_inc % overcast (:,:) = 0._JPRB
  radiance_inc % clear(:)    = 0._JPRB
  radiance_inc % total(:)    = 0._JPRB
  radiance_inc % bt_clear(:) = 0._JPRB
  radiance_inc % bt(:)       = 0._JPRB
  radiance_inc % clear_out(:)    = 0._JPRB
  radiance_inc % total_out(:)    = 0._JPRB
  radiance_inc % out_clear(:) = 0._JPRB
  radiance_inc % out(:)       = 0._JPRB

  if (kinrad == 2) then
     radiance_inc % out(ichan)       = 1._JPRB
  else
     radiance_inc % total_out(ichan)    = 1._JPRB
  endif

  Call Rttov_scatt_ad ( &
     & rttov_errorstatus,  &! out
     & nfrequencies,    &! in
     & nchannels,       &! in
     & nbtout,          &! in
     & nprofiles,       &! in
     & channels,        &! in
     & polarisations,   &! in
     & lprofiles,       &! in
     & profiles,        &! in
     & cld_profiles,    &! in
     & coef_rttov,      &! in
     & coef_scatt,      &! in
     & calcemis,        &! in
     & emissivity,      &! inout
     & profiles_ad,     &! inout
     & cld_profiles_ad, &! inout
     & emissivity_ad,   &! inout
     & radiance2,       &! inout
     & radiance_inc    ) ! inout

  If ( any( rttov_errorstatus(:) == errorstatus_fatal ) ) Then
     Do iatm = 1, nprofiles
        If ( rttov_errorstatus(iatm) == errorstatus_fatal ) Then
           write ( ioout, * ) 'rttov_scatt_ad error for profile',iatm
        End If
     End Do
     Stop
  End If

  do j = lprofiles(ichan), lprofiles(ichan)
     do lev = 1, profiles_ad(j) % nlevels
        if ( abs(profiles_ad(j)%t(lev) - profiles_k(ichan)%t(lev)) > threshold ) then
           write(default_err_unit,*) 'test K 1 fails',lev
           !stop
        End If
        if ( abs(profiles_ad(j)%q(lev) - profiles_k(ichan)%q(lev)) > threshold ) Then
           write(default_err_unit,*) 'test K 2 fails',lev
           !stop
        End If
        if ( abs(profiles_ad(j)%o3(lev) - profiles_k(ichan)%o3(lev)) > threshold ) Then
           write(default_err_unit,*) 'test K 3 fails',lev
           !stop
        End If
     enddo
  enddo

  do j = lprofiles(ichan), lprofiles(ichan)
     do lev = 1, cld_profiles_ad(j) % nlevels
        if ( abs(cld_profiles_ad(j)%p(lev)   - cld_profiles_k(ichan)%p(lev)) > threshold ) Then
           write(default_err_unit,*) 'test K 4 fails',lev
           stop
        End If
        if ( abs(cld_profiles_ad(j)%ph(lev)  - cld_profiles_k(ichan)%ph(lev)) > threshold ) Then
           write(default_err_unit,*) 'test K 5 fails',lev
           stop
        End If
        if ( abs(cld_profiles_ad(j)%t(lev)   - cld_profiles_k(ichan)%t(lev)) > threshold ) Then
           write(default_err_unit,*) 'test K 6 fails',lev
           stop
        End If
        if ( abs(cld_profiles_ad(j)%cc(lev)  - cld_profiles_k(ichan)%cc(lev)) > threshold ) Then
           write(default_err_unit,*) 'test K 7 fails',lev
           stop
        End If
        if ( abs(cld_profiles_ad(j)%clw(lev) - cld_profiles_k(ichan)%clw(lev)) > threshold ) Then
           write(default_err_unit,*) 'test K 8 fails',lev
           stop
        End If
        if ( abs(cld_profiles_ad(j)%ciw(lev) - cld_profiles_k(ichan)%ciw(lev)) > threshold ) Then
           write(default_err_unit,*) 'test K 9 fails',lev
           stop
        End If
        if ( abs(cld_profiles_ad(j)%rain(lev) - cld_profiles_k(ichan)%rain(lev)) > threshold ) Then
           write(default_err_unit,*) 'test K 18 fails',lev
           stop
        End If
        if ( abs(cld_profiles_ad(j)%sp(lev) - cld_profiles_k(ichan)%sp(lev)) > threshold ) Then
           write(default_err_unit,*) 'test K 19 fails',lev
           stop
        End If
     enddo
     lev = cld_profiles_ad(j) % nlevels+1
     if ( abs(cld_profiles_ad(j)%ph(lev)  - cld_profiles_k(ichan)%ph(lev)) > threshold ) Then
        write(default_err_unit,*) 'test K 10 fails',lev
        stop
     End If
  Enddo

  do j = lprofiles(ichan), lprofiles(ichan)
     if ( abs(profiles_ad(j)%s2m%t - profiles_k(ichan)%s2m%t) > threshold ) Then
        write(default_err_unit,*) 'test K 11 fails',j
        stop
     End If
     if ( abs(profiles_ad(j)%s2m%q - profiles_k(ichan)%s2m%q) > threshold ) Then
        write(default_err_unit,*) 'test K 12 fails',j
        stop
     End If
     if ( abs(profiles_ad(j)%s2m%p - profiles_k(ichan)%s2m%p) > threshold ) Then
        write(default_err_unit,*) 'test K 13 fails',j
        stop
     End If
     if ( abs(profiles_ad(j)%s2m%u - profiles_k(ichan)%s2m%u) > threshold ) Then
        write(default_err_unit,*) 'test K 14 fails',j
        stop
     End If
     if ( abs(profiles_ad(j)%s2m%v - profiles_k(ichan)%s2m%v) > threshold ) Then
        write(default_err_unit,*) 'test K 15 fails',j
        stop
     End If

     if ( abs(profiles_ad(j)%skin%t - profiles_k(ichan)%skin%t) > threshold ) Then
        write(default_err_unit,*) 'test K 16 fails',j
        stop
     End If
  enddo


  Do kpol = 1 , polarisations(ichan,3)
    jchan = jchan + 1
    if ( abs(emissivity_ad(jchan) - emissivity_k(jchan)) > threshold ) Then
       write(default_err_unit,*) 'test K 17 fails',j
       stop
    End If
  Enddo

  Enddo

  write(ioout,*) 'K is ok'
  write(ioout,*)
  write(ioout,*) 'End of RTTOVSCATT tests'

  Stop

CONTAINS
!
! ----------------------------------------
!
      SUBROUTINE EC_P60l(spres,pap,paph)
!
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
!     Description:
!     Computes the 60-level vertical pressure grid
!       associated to the input surface pressure
!     All pressures are in Pa

  Use parkind1, Only : jpim     ,jprb
      implicit none
  Integer(Kind=jpim), parameter    :: nlev=60
  Integer(Kind=jpim)               :: jk
  Real(Kind=jprb)          :: spres
  Real(Kind=jprb)          :: aam(nlev+1), bbm(nlev+1)
  Real(Kind=jprb)          :: pap(nlev), paph(nlev+1)

      data aam / &
          & 0.000000_JPRB,    20.000000_JPRB,    38.425343_JPRB, &
         & 63.647804_JPRB,    95.636963_JPRB,   134.483307_JPRB, &
        & 180.584351_JPRB,   234.779053_JPRB,   298.495789_JPRB, &
        & 373.971924_JPRB,   464.618134_JPRB,   575.651001_JPRB, &
        & 713.218079_JPRB,   883.660522_JPRB,  1094.834717_JPRB, &
       & 1356.474609_JPRB,  1680.640259_JPRB,  2082.273926_JPRB, &
       & 2579.888672_JPRB,  3196.421631_JPRB,  3960.291504_JPRB, &
       & 4906.708496_JPRB,  6018.019531_JPRB,  7306.631348_JPRB, &
       & 8765.053711_JPRB, 10376.126953_JPRB, 12077.446289_JPRB, &
      & 13775.325195_JPRB, 15379.805664_JPRB, 16819.474609_JPRB, &
      & 18045.183594_JPRB, 19027.695313_JPRB, 19755.109375_JPRB, &
      & 20222.205078_JPRB, 20429.863281_JPRB, 20384.480469_JPRB, &
      & 20097.402344_JPRB, 19584.330078_JPRB, 18864.750000_JPRB, &
      & 17961.357422_JPRB, 16899.468750_JPRB, 15706.447266_JPRB, &
      & 14411.124023_JPRB, 13043.218750_JPRB, 11632.758789_JPRB, &
      & 10209.500977_JPRB,  8802.356445_JPRB,  7438.803223_JPRB, &
       & 6144.314941_JPRB,  4941.778320_JPRB,  3850.913330_JPRB, &
       & 2887.696533_JPRB,  2063.779785_JPRB,  1385.912598_JPRB, &
        & 855.361755_JPRB,   467.333588_JPRB,   210.393890_JPRB, &
         & 65.889244_JPRB,     7.367743_JPRB,     0.000000_JPRB, &
          & 0.000000_JPRB &
              & /
      data bbm / &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000758235_JPRB, 0.0004613950_JPRB, 0.0018151561_JPRB, &
      & 0.0050811190_JPRB, 0.0111429105_JPRB, 0.0206778757_JPRB, &
      & 0.0341211632_JPRB, 0.0516904071_JPRB, 0.0735338330_JPRB, &
      & 0.0996746942_JPRB, 0.1300225109_JPRB, 0.1643843204_JPRB, &
      & 0.2024759352_JPRB, 0.2439331412_JPRB, 0.2883229554_JPRB, &
      & 0.3351548910_JPRB, 0.3838921487_JPRB, 0.4339629412_JPRB, &
      & 0.4847715795_JPRB, 0.5357099175_JPRB, 0.5861684084_JPRB, &
      & 0.6355474591_JPRB, 0.6832686067_JPRB, 0.7287858129_JPRB, &
      & 0.7715966105_JPRB, 0.8112534285_JPRB, 0.8473749161_JPRB, &
      & 0.8796569109_JPRB, 0.9078838825_JPRB, 0.9319403172_JPRB, &
      & 0.9518215060_JPRB, 0.9676452279_JPRB, 0.9796627164_JPRB, &
      & 0.9882701039_JPRB, 0.9940194488_JPRB, 0.9976301193_JPRB, &
      & 1.0000000000_JPRB &
              & /

      do jk=1,nlev+1
      paph(jk)=aam(jk)+bbm(jk)*spres
      end do
      do jk=1,nlev
      pap(jk)=0.5_JPRB*(paph(jk)+paph(jk+1))
      end do

      end subroutine ec_p60l

END PROGRAM RTTOVSCATT_TEST
