! Copyright (c) 2026, MPAS-PartMC.
!
! Wind nudging (grid nudging of the horizontal wind toward analyses).
!
! The normal velocity u is relaxed toward analyses that are linearly interpolated in time:
!
!    d(ru)/dt += rho_edge * W(z) * (u_a(t) - u) / tau
!
! added to the physics tendency of ru once per time step, like the tendencies of the PBL and
! convection schemes. Only the wind is nudged; temperature, humidity, clouds and precipitation
! remain the model's own, as recommended for aerosol studies by Zhang et al. (2014, ACP 14,
! 8631) and as done in specified-dynamics CAM-chem and E3SM.
!
! The analyses come in through the input stream "nudging" (defined in the run's
! streams.atmosphere; run_120km.sbatch adds it with NUDGE=dir). The file read at time T holds
! the analysis valid at T in nudge_u_a and the analysis valid at T + interval in nudge_u_b, so
! that the target at any time t in [T, T + interval) is
!
!    u_a(t) = (1 - w) nudge_u_a + w nudge_u_b,   w = (t - T) / interval,
!
! with T the last multiple of the interval since 00 UTC. Nothing depends on earlier steps, so
! restarts give the same result as a continuous run. Files are made with tools/make_nudge.py
! from init_atmosphere output (config_init_case = 7 on the run's static file and vertical
! grid), so the analysis u lives on the model's own edges and levels.
!
! Environment:
!   MPAS_NUDGE=1           switch on (the stream must exist)
!   MPAS_NUDGE_TAU         relaxation time (s), default 21600
!   MPAS_NUDGE_INTERVAL    time between analyses (s), default 21600
!   MPAS_NUDGE_ZBOT        height above the ground (m) below which no nudging is applied,
!                          default 0; the weight rises linearly to 1 over MPAS_NUDGE_ZRAMP
!                          (m, default 500) above it
!   MPAS_NUDGE_ZTOP        height above sea level (m) above which the weight falls linearly to 0
!                          over MPAS_NUDGE_ZRAMP, default 1e9 (no limit)
module mpas_atm_nudging

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log, only : mpas_log_write
   use mpas_timekeeping

   implicit none
   private
   public :: atm_nudge_add_tend

   logical, save :: initialized = .false.
   logical, save :: nudge_on = .false.
   real (kind=RKIND), save :: tau = 21600.0_RKIND, interval = 21600.0_RKIND
   real (kind=RKIND), save :: zbot = 0.0_RKIND, ztop = 1.0e9_RKIND, zramp = 500.0_RKIND
   real (kind=RKIND), allocatable, save :: wgt(:,:)
   real (kind=RKIND), save :: sums(4) = 0.0_RKIND

contains

   subroutine read_env()
      character(len=64) :: v
      integer :: stat
      call get_environment_variable('MPAS_NUDGE', v, status=stat)
      if (stat == 0) nudge_on = (trim(v) == '1')
      call get_environment_variable('MPAS_NUDGE_TAU', v, status=stat)
      if (stat == 0) read(v, *) tau
      call get_environment_variable('MPAS_NUDGE_INTERVAL', v, status=stat)
      if (stat == 0) read(v, *) interval
      call get_environment_variable('MPAS_NUDGE_ZBOT', v, status=stat)
      if (stat == 0) read(v, *) zbot
      call get_environment_variable('MPAS_NUDGE_ZTOP', v, status=stat)
      if (stat == 0) read(v, *) ztop
      call get_environment_variable('MPAS_NUDGE_ZRAMP', v, status=stat)
      if (stat == 0) read(v, *) zramp
   end subroutine read_env


   !> Add the nudging tendency to tend_ru (the coupled physics tendency of ru).
   subroutine atm_nudge_add_tend(domain, block, itimestep, tend_ru)

      type (domain_type), intent(in) :: domain
      type (block_type), intent(in) :: block
      integer, intent(in) :: itimestep
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend_ru

      type (mpas_pool_type), pointer :: mesh, state, diag, nudge
      real (kind=RKIND), dimension(:,:), pointer :: u, rho_edge, ua, ub, zgrid
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: nEdgesSolve, nVertLevels, nCells
      integer :: iEdge, k, c1, c2, yr, mo, dy, h, mi, sc, ierr
      real (kind=RKIND) :: w, sod, z, zs, du, target
      type (MPAS_Time_type) :: now

      if (.not. initialized) then
         call read_env()
         initialized = .true.
         if (nudge_on) call mpas_log_write('Nudging of the horizontal wind: tau = $r s, analyses every $r s, ' &
                                           // 'from $r m above the ground (ramp $r m), up to $r m', &
                                           realArgs=(/tau, interval, zbot, zramp, ztop/))
      end if
      if (.not. nudge_on) return

      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_subpool(block % structs, 'nudge', nudge)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(state, 'u', u, 1)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)
      call mpas_pool_get_array(nudge, 'nudge_u_a', ua)
      call mpas_pool_get_array(nudge, 'nudge_u_b', ub)

      ! vertical weights, once: height of the layer middle at the edge (mean of its two cells)
      if (.not. allocated(wgt)) then
         allocate(wgt(nVertLevels, nEdgesSolve))
         do iEdge = 1, nEdgesSolve
            c1 = cellsOnEdge(1,iEdge)
            c2 = cellsOnEdge(2,iEdge)
            if (c1 > nCells) c1 = c2
            if (c2 > nCells) c2 = c1
            zs = 0.5_RKIND * (zgrid(1,c1) + zgrid(1,c2))
            do k = 1, nVertLevels
               z = 0.25_RKIND * (zgrid(k,c1) + zgrid(k+1,c1) + zgrid(k,c2) + zgrid(k+1,c2))
               wgt(k,iEdge) = min(1.0_RKIND, max(0.0_RKIND, (ztop - z) / max(zramp, 1.0_RKIND) + 1.0_RKIND))
               if (zbot > 0.0_RKIND) wgt(k,iEdge) = wgt(k,iEdge) &
                    * min(1.0_RKIND, max(0.0_RKIND, (z - zs - zbot) / max(zramp, 1.0_RKIND)))
            end do
         end do
      end if

      ! time weight from the time of day (files at multiples of the interval from 00 UTC)
      now = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(now, YYYY=yr, MM=mo, DD=dy, H=h, M=mi, S=sc, ierr=ierr)
      sod = real(h*3600 + mi*60 + sc, kind=RKIND)
      w = modulo(sod, interval) / interval
      ! a stream that was not read leaves zeros, which would pull the wind to rest
      if (w == 0.0_RKIND .or. itimestep == 1) then
         if (nEdgesSolve > 0) then
            if (maxval(abs(ua(:,1:nEdgesSolve))) == 0.0_RKIND .or. maxval(abs(ub(:,1:nEdgesSolve))) == 0.0_RKIND) &
               call mpas_log_write('Nudging: the analyses are zero; is the nudging stream defined and its file there?', &
                                   messageType=MPAS_LOG_CRIT)
         end if
      end if

      do iEdge = 1, nEdgesSolve
         do k = 1, nVertLevels
            if (wgt(k,iEdge) <= 0.0_RKIND) cycle
            target = (1.0_RKIND - w) * ua(k,iEdge) + w * ub(k,iEdge)
            du = target - u(k,iEdge)
            tend_ru(k,iEdge) = tend_ru(k,iEdge) + rho_edge(k,iEdge) * wgt(k,iEdge) * du / tau
            sums(1) = sums(1) + du**2
            sums(2) = sums(2) + 1.0_RKIND
         end do
      end do
      if (mod(itimestep, 12) == 0 .or. itimestep == 1) then
         call mpas_log_write('Nudging step $i: time weight $r, rms difference analysis - model on this ' &
                             // 'task since the last report $r m/s', intArgs=(/itimestep/), &
                             realArgs=(/w, sqrt(sums(1) / max(sums(2), 1.0_RKIND))/))
         sums = 0.0_RKIND
      end if

   end subroutine atm_nudge_add_tend

end module mpas_atm_nudging
