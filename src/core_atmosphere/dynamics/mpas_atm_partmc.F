! PartMC coupling for MPAS-Atmosphere (development).
!
! A number-mixing-ratio scalar (pmc_num) is advected by MPAS and its total
! time-integrated transport fluxes are captured in atm_advance_scalars_mono_work.
! PartMC particles held in every cell are moved stochastically with those fluxes.
!
! Run-time settings (environment variables):
!   PMC_N_PART        target computational particles per cell (0 = no particles)
!   PMC_TRANSPORT     flux (default): stochastic transfer with the captured fluxes;
!                     lagrangian: particles carry positions and follow the winds
!   PMC_LAGR_CONTROL  lagrangian mode: 1 = keep the particle count per cell near PMC_N_PART by
!                     splitting and thinning particles with power-of-two weights (default 0)
!   PMC_NUM_SCALE     factor applied to the initial number mixing ratios (default 1)
!   PMC_COAG          1 = PartMC Brownian coagulation in every cell each step (default 0)
!   PMC_KZ_MAX        peak of a prescribed vertical diffusivity profile (m^2/s, default 0 = no
!                     vertical mixing); applied to the Eulerian scalars and, in lagrangian
!                     mode, to the particles. PMC_KZ_TOP is the height where it reaches zero.
!   PMC_MOSAIC        1 = MOSAIC gas chemistry and gas-particle partitioning in every cell
!                     (needs a build with MOSAIC; default 0). PMC_GAS_DATA and PMC_GAS_INIT
!                     name the PartMC gas files (gas_data.dat, gas_init.dat); PMC_CHEM_TOP
!                     is the altitude (m) above which chemistry is skipped.
!   PMC_CHEM_BOX      "id:level,id:level": save these cells (global cell ID and level) before and
!                     after each of the first PMC_CHEM_BOX_STEPS (60) chemistry steps in
!                     PartMC's output format, in the directory box/
!   PMC_GAS_EMIT      PartMC gas emissions file (times, rate scale, mol m^-2 s^-1 per species),
!                     applied over land to the lowest PMC_EMIT_LEVELS layers and multiplied by
!                     PMC_GAS_EMIT_SCALE (default 1). Needs PMC_MOSAIC=1.
!   PMC_GAS_DRYDEP    file of gas dry deposition velocities (species name, m/s), applied to the
!                     lowest layer. Needs PMC_MOSAIC=1.
!   PMC_WET           1 = removal of particles and gases by grid-scale precipitation, in cloud
!                     (default 0), and with PMC_CONV=1 also in the convective updraft unless
!                     PMC_WET_CONV=0. 2 = only the diagnostics: nothing is removed. PMC_WET_SCRIT: only particles with a critical supersaturation
!                     up to this value (as a fraction, 0.002 = 0.2%) are removed (default: all).
!                     PMC_GAS_WET: file of the fraction of each gas that is in cloud water
!                     (species name, 0 to 1); needs PMC_MOSAIC=1.
!   PMC_EMIS_FILE     NetCDF file of emissions for every cell of the mesh (dimension nCells in
!                     global cell order): gas_<NAME> (mol m^-2 s^-1) for any gas species, bc and
!                     oc (kg m^-2 s^-1). The gases need PMC_MOSAIC=1; BC and OC are emitted as
!                     one carbonaceous mode (PMC_CARB_DG, default 6e-8 m, geometric standard
!                     deviation 1.8) with the cell's BC/OC ratio. PMC_EMIS_SCALE multiplies all.
!   PMC_SEASALT       1 = sea salt particles emitted over open water from the 10-m wind speed
!                     with the flux of Gong (2003), dry diameters 0.07-10 um; PMC_SEASALT_SCALE
!                     multiplies it
!   PMC_Q_BACKGROUND, PMC_Q_BLOB  initial number mixing ratio (per kg) of the background particles
!                     and of the blob (defaults 1e8 and 1e9)
!   PMC_EMIT_FLUX     number flux (m^-2 s^-1) of soot-like particles emitted over land (default 0)
!   PMC_EMIT_LEVELS   number of lowest layers the emission is spread over, uniformly in mixing
!                     ratio (default 1)
!   PMC_DRYDEP_V      dry deposition velocity (m/s) applied to the lowest layer (default 0)
!   PMC_DRYDEP_SCHEME 1 = lagrangian mode: a deposition velocity for each particle from its size
!                     and density with PartMC's scenario_loss_rate_drydep (Zhang et al. 2001
!                     form), using the wind speed of the lowest layer and the model's roughness
!                     length with a neutral surface layer; the reference scalars are not
!                     changed (default 0)
!   PMC_CONV          1 = convective transport of tracers and particles with the cu_ntiedtke
!                     scheme's mass fluxes (default 0)
!   PMC_CONV_CHECK    1 = with PMC_CONV, compare the change of moist static energy that the column
!                     transfer matrix gives with the convection scheme's own (its temperature
!                     and vapour tendencies), as a check of the matrix; diagnostic only
!   PMC_KZ_PHYS       1 = vertical mixing with the boundary-layer scheme's heat diffusivity
!                     (kzh) in place of the prescribed profile (default 0)
!   PMC_LAGR_SUBSTEPS lagrangian mode: trajectory sub-steps per model step (default 1)
!   PMC_LAGR_FLUXVEL  lagrangian mode: velocity built in each cell from its own edge and face
!                     mass fluxes instead of interpolating reconstructed cell-centre winds
!                     (default 0 = interpolate). 1: constant plus linear field with the model's
!                     discrete divergence. 2: also matches the mass flux through every edge
!                     (lowest-order Raviart-Thomas field on the triangles centre-vertex-vertex).
!   PMC_LAGR_SHUFFLE  lagrangian mode, with PMC_WET: 1 = before wet removal every particle is moved
!                     to a random point of its own cell and layer, which makes the cells well
!                     mixed as an Eulerian scheme assumes (a test; default 0)
!   PMC_LAGR_UCONST   lagrangian mode: 1 = horizontal velocity constant within a layer instead of
!                     linear between layer centres (default 0)
!   PMC_TEST_WIND     lagrangian mode: replace the model winds by an analytic flow with a
!                     12-day period: solid (tilted solid-body rotation) or deform
!                     (Nair and Lauritzen 2010 reversing deformational flow)
!   PMC_N_CLASS       flux mode: 1 (default) one weight class; 2 one class per source,
!                     each with its own advected number scalar
!   PMC_COUPLED       1 (default): the class scalars are reset from the particles every step
!   PMC_SCALE_HEIGHT  e-folding height (m) of the initial number mixing ratio
!   PMC_AERO_DATA     PartMC aerosol species file (default aero_data.dat)
!   PMC_RESTART_INTERVAL  write the complete particle state every this many steps (0 = never)
!                     to PMC_RESTART_DIR (default pmc_restart), one binary file per MPI task
!   PMC_RESTART_READ  step number of the particle state to read at start-up, for a model
!                     restart from the same time with the same number of tasks and partition
!   PMC_DIAG_INTERVAL write time-mean aerosol fields of every cell (number, CCN, species
!                     group masses, and the sums needed for the mixing state index chi) every
!                     this many steps to PMC_DIAG_DIR (default pmc_diag), one NetCDF file per
!                     task; the particles are sampled every PMC_DIAG_EVERY steps (default 5)
!   PMC_OUT_INTERVAL  write particle files every this many steps (0 = never)
!   PMC_OUT_FRAC      fraction of cells written each time (default 1)
!   PMC_OUT_DIR       directory for particle files (default pmc_out)

! Variable-length byte-buffer exchange between neighbouring MPI tasks, used to
! move PartMC particles across block boundaries. Kept in its own module so that
! the "mpi" module used here (and by PartMC) does not clash with "mpi_f08" in MPAS.
module mpas_atm_partmc_comm

#ifdef _MPI
   use mpi
#endif

   implicit none
   private

   public :: pmc_buf_t, pmc_buf_reserve, pmc_comm_exchange

   type pmc_buf_t
      character, allocatable :: data(:)
      integer :: pos = 0
   end type pmc_buf_t

contains

   !> Make sure the buffer can take n_extra more bytes.
   subroutine pmc_buf_reserve(buf, n_extra)
      type(pmc_buf_t), intent(inout) :: buf
      integer, intent(in) :: n_extra
      character, allocatable :: tmp(:)
      integer :: n_new

      if (.not. allocated(buf%data)) then
         allocate(buf%data(max(4096, 2*n_extra)))
         buf%pos = 0
      else if (buf%pos + n_extra > size(buf%data)) then
         n_new = max(2*size(buf%data), buf%pos + 2*n_extra)
         allocate(tmp(n_new))
         tmp(1:buf%pos) = buf%data(1:buf%pos)
         call move_alloc(tmp, buf%data)
      end if
   end subroutine pmc_buf_reserve


   !> Send send(i) to task nbr_rank(i) and receive one buffer from each neighbour.
   !> On return recv(i)%data has exactly the bytes received and recv(i)%pos = 0.
   subroutine pmc_comm_exchange(comm, nbr_rank, send, recv)
      integer, intent(in) :: comm
      integer, intent(in) :: nbr_rank(:)
      type(pmc_buf_t), intent(inout) :: send(:)
      type(pmc_buf_t), intent(inout) :: recv(:)

#ifdef _MPI
      integer, parameter :: tag = 7701
      integer :: i, n, ierr
      integer :: status(MPI_STATUS_SIZE)
      integer, allocatable :: req(:)

      n = size(nbr_rank)
      allocate(req(n))
      do i = 1, n
         if (.not. allocated(send(i)%data)) then
            allocate(send(i)%data(1))
            send(i)%pos = 0
         end if
         call MPI_Isend(send(i)%data, send(i)%pos, MPI_CHARACTER, nbr_rank(i), tag, comm, req(i), ierr)
      end do
      do i = 1, n
         call MPI_Probe(nbr_rank(i), tag, comm, status, ierr)
         if (allocated(recv(i)%data)) deallocate(recv(i)%data)
         call MPI_Get_count(status, MPI_CHARACTER, recv(i)%pos, ierr)
         allocate(recv(i)%data(max(1, recv(i)%pos)))
         call MPI_Recv(recv(i)%data, recv(i)%pos, MPI_CHARACTER, nbr_rank(i), tag, comm, status, ierr)
      end do
      call MPI_Waitall(n, req, MPI_STATUSES_IGNORE, ierr)
      do i = 1, n
         ! recv(i)%pos temporarily held the byte count; keep it in the array size
         if (recv(i)%pos < size(recv(i)%data)) then
            deallocate(recv(i)%data)
            allocate(recv(i)%data(0))
         end if
         recv(i)%pos = 0
         send(i)%pos = 0
      end do
#endif
   end subroutine pmc_comm_exchange

end module mpas_atm_partmc_comm


module mpas_atm_partmc

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log, only : mpas_log_write
   use mpas_constants, only : rvord, pii
   use mpas_vector_reconstruction, only : mpas_reconstruct
#ifdef MPAS_PARTMC
   use pmc_aero_data
   use pmc_aero_state
   use pmc_aero_particle
   use pmc_aero_weight
   use pmc_aero_weight_array
   use pmc_fractal
   use pmc_rand
   use pmc_mpi
   use pmc_spec_file
   use pmc_env_state
   use pmc_coagulation
   use pmc_scenario, only : scenario_t, scenario_loss_rate_drydep
#ifdef MPAS_PARTMC_MOSAIC
   use pmc_gas_data
   use pmc_gas_state
   use pmc_mosaic
   use pmc_output, only : output_state, OUTPUT_TYPE_DIST
#endif
   use pmc_coag_kernel
   use pmc_aero_sorted
   use pmc_aero_info_array
   use pmc_bin_grid
   use netcdf
   use mpas_atm_partmc_comm
#endif

   implicit none
   private

   public :: atm_partmc_init, atm_partmc_check, atm_partmc_particles_init, atm_partmc_particles_transport, &
             atm_partmc_particles_to_scalar, atm_partmc_output, atm_partmc_vertical_mixing, &
             atm_partmc_convection, atm_partmc_surface, atm_partmc_gas_surface, atm_partmc_chemistry, &
             atm_partmc_wet_removal, atm_partmc_restart_write, atm_partmc_diag, &
             atm_partmc_coagulation

   integer, parameter :: MODE_FLUX = 1, MODE_LAGR = 2
   integer, parameter :: MAX_CLASS = 2
   !> Reals stored per particle in Lagrangian mode: unit vector, zeta, initial unit vector.
   integer, parameter :: NPOS = 7
   integer, parameter :: WIND_MODEL = 0, WIND_SOLID = 1, WIND_DEFORM = 2
   !> Period of the analytic test winds (s).
   real (kind=RKIND), parameter :: test_period = 12.0_RKIND * 86400.0_RKIND

   !> Initial number mixing ratios (kg^-1) and shape of the blob.
   real (kind=RKIND) :: q_background = 1.0e8_RKIND
   real (kind=RKIND) :: q_blob = 1.0e9_RKIND
   real (kind=RKIND), parameter :: blob_lat = 0.6981317_RKIND   ! 40 N
   real (kind=RKIND), parameter :: blob_lon = 1.5707963_RKIND   ! 90 E
   real (kind=RKIND), parameter :: blob_width = 0.25_RKIND      ! radians of arc

   logical :: settings_read = .false.
   !> Target number of computational particles per cell (0 disables particles). Env: PMC_N_PART.
   integer :: pmc_n_part = 0
   !> Whether the class scalars are reset from the particle counts before every step (flux mode).
   logical :: pmc_coupled = .true.
   !> Transport scheme and number of weight classes (flux mode only).
   integer :: pmc_mode = MODE_FLUX
   integer :: n_class = 1
   real (kind=RKIND) :: scale_height = 0.0_RKIND
   character(len=256) :: aero_data_file = 'aero_data.dat'
   !> Lagrangian mode: winds used (model, or an analytic test flow) and sub-steps per model step.
   integer :: lagr_wind = WIND_MODEL
   integer :: lagr_substeps = 1
   !> Lagrangian mode: a particle in weight class g represents lagr_w * 2**(g - LAGR_G0) real
   !> particles. (PartMC weight classes, not weight groups: PartMC treats the groups of a class
   !> as one combined sample, so a group does not give a per-particle weight.) With population
   !> control, cells split heavy particles and thin light ones to stay near the target count.
   integer, parameter :: LAGR_GROUPS = 64, LAGR_G0 = 45
   logical :: lagr_control = .false.
   !> Particles removed by thinning and added by splitting since the start of the run.
   integer :: n_thinned = 0, n_split = 0
   !> Test vertical mixing: peak diffusivity (m^2/s, 0 = off) and the height where it reaches zero (m).
   real (kind=RKIND) :: pmc_kz_max = 0.0_RKIND, pmc_kz_top = 12000.0_RKIND
   !> Vertical mixing with the boundary-layer scheme's diffusivity (kzh) instead of the test profile.
   logical :: pmc_kz_phys = .false.
   !> Convective transport with the cumulus scheme's mass fluxes, and particle layer changes so far.
   logical :: pmc_conv = .false.
   !> MOSAIC chemistry: switch, gas files, top altitude (m), the model scalar index of each
   !> PartMC gas species, and the initial global sulfur (mol).
   logical :: pmc_chem = .false.
   character(len=256) :: gas_data_file = 'gas_data.dat', gas_init_file = 'gas_init.dat'
   real (kind=RKIND) :: pmc_chem_top = 1.0e9_RKIND, chem_s0 = 0.0_RKIND
   integer, allocatable :: gas_idx(:)
   !> Molar mass of dry air (kg/mol) used to turn the gas mixing ratios into moles. PartMC's
   !> MOSAIC interface turns mixing ratios into concentrations with p/(R T), R = 82.056 atm cm3
   !> mol-1 K-1 = 8.31432 J mol-1 K-1, and MPAS's pressure is rho_d * 287 * T * (1 + rvord qv),
   !> so the moles MOSAIC moves between gas and particles match the moles counted here only with
   !> M = 8.31432/287. The value 0.028966 used before differs by 1.5e-4, which showed as a
   !> steady loss of sulfur while sulfate condensed.
   real (kind=RKIND), parameter :: M_AIR = 82.056e-6_RKIND * 101325.0_RKIND / 287.0_RKIND
   !> Gas emission and dry deposition: files, scale, and the amounts (mol) emitted and
   !> deposited by this task so far for each species.
   character(len=256) :: gas_emit_file = '', gas_drydep_file = ''
   real (kind=RKIND) :: gas_emit_scale = 1.0_RKIND
   real (kind=RKIND), allocatable :: gas_emitted(:), gas_deposited(:), gas_content0(:), gas_wetdep(:)
   !> Wet removal: switch, largest critical supersaturation removed, gas file, and sums by this
   !> task so far: number removed (reference, particles), precipitation at the surface, column
   !> production of precipitation, and the part of the production counted as reaching the
   !> surface (kg); last, the number of particles expected to be removed given the particles
   !> that were in the layers (the sum of probability times number).
   logical :: pmc_wet = .false., pmc_wet_diag_only = .false., lagr_shuffle = .false.
   !> Removal in the convective updraft: switch, and sums by this task so far: number removed
   !> (reference, particles, expected from the particles present), convective precipitation at
   !> the surface, water condensed in the updraft, and the part of it counted as precipitation
   !> (kg).
   logical :: pmc_wet_conv = .true.
   !> Check of the convective matrix against the scheme's own tendencies: switch, and
   !> mass-weighted sums of matrix*scheme, scheme^2, matrix^2 by this task so far, for moist
   !> static energy, zonal wind, meridional wind and the scheme's passive test tracer
   !> exp(-z/3000 m) (second index), with the matrix result
   !> taken from the layer below, the same layer and the layer above (third index).
   logical :: pmc_conv_check = .false.
   real (kind=RKIND) :: conv_chk(3,4,3) = 0.0_RKIND
   !> The same sums for the changes added up over slabs of each column (below 2 km, 2-6 km,
   !> above 6 km; third index), which do not depend on how a change is spread within a slab.
   real (kind=RKIND) :: conv_slab(3,4,3) = 0.0_RKIND
   real (kind=RKIND) :: conv_wet(6) = 0.0_RKIND
   !> Fraction of each gas species that is in cloud water (from the PMC_GAS_WET file).
   real (kind=RKIND), allocatable :: gas_wet_frac(:)
   real (kind=RKIND) :: pmc_wet_scrit = 1.0e9_RKIND
   character(len=256) :: gas_wet_file = ''
   real (kind=RKIND) :: wet_sums(6) = 0.0_RKIND
   !> Sums weighted by the reference number removed, to describe where removal happens:
   !> probability, vertical Courant number of rising air, share in rising air, number of
   !> adjacent removing layers, and the share in runs of 1, 2-3, 4-7 and 8 or more layers.
   real (kind=RKIND) :: wet_where(8) = 0.0_RKIND
   !> Number removed from the reference (1) and expected from the particles present (2), by
   !> class: runs of 1, 2-3, 4-7, 8+ removing layers; rising, sinking air; lowest layer of its
   !> run, other layers; share of the horizontal inflow of air that comes from neighbours that
   !> are not removing in this layer: none, up to half, more than half.
   real (kind=RKIND) :: wet_cls(2,11) = 0.0_RKIND
   !> Sulfate (mol) in particles removed by dry and wet deposition, by this task so far.
   real (kind=RKIND) :: part_s_removed = 0.0_RKIND
   !> Net sulfate (mol) added to the particles by the thinning of the population control
   !> since the start of this run (not kept over restarts). Zero in expectation.
   real (kind=RKIND) :: part_s_thin = 0.0_RKIND
   !> Repeat sample cells' chemistry after another cell to check MOSAIC carries nothing over.
   logical :: pmc_chem_selftest = .false.
   !> Cells (global cell ID and level) whose state is saved in PartMC's own file format before
   !> and after each chemistry step, for replaying the step with the standalone PartMC program.
   integer :: n_box = 0, box_id(16), box_k(16), box_steps = 60
#ifdef MPAS_PARTMC_MOSAIC
   type(gas_data_t) :: gas_data
#endif
   !> Surface emission flux over land (m^-2 s^-1) and dry deposition velocity (m/s).
   real (kind=RKIND) :: pmc_emit_flux = 0.0_RKIND, pmc_drydep_v = 0.0_RKIND
   !> Real particles so far: emitted (reference, particles) and deposited (reference, particles).
   real (kind=RKIND) :: surf_sums(4) = 0.0_RKIND
   !> Size-dependent dry deposition: switch, and sums by this task so far: number expected to
   !> be removed (probability times number), number removed, and number-weighted sums of the
   !> deposition velocity and of one, for the mean velocity.
   logical :: pmc_drydep_scheme = .false.
   !> Emissions per cell from an inventory file: file, scale, the carbonaceous mode, fluxes of
   !> the owned cells (BC and OC in kg m^-2 s^-1, gases in mol m^-2 s^-1).
   character(len=256) :: emis_file = ''
   real (kind=RKIND) :: emis_scale = 1.0_RKIND, carb_dg = 6.0e-8_RKIND, carb_sigma = 1.8_RKIND
   real (kind=RKIND), allocatable :: emis_bc(:), emis_oc(:), emis_gas(:,:)
   logical :: emis_part_read = .false.
   !> Sea salt: switch, scale, and the size distribution of the flux at 1 m/s: bin edges of
   !> the dry diameter (m), cumulative share of the number, number flux (m^-2 s^-1).
   logical :: pmc_seasalt = .false.
   real (kind=RKIND) :: seasalt_scale = 1.0_RKIND, ss_unit = 0.0_RKIND
   integer, parameter :: N_SS = 80
   real (kind=RKIND) :: ss_edge(N_SS+1), ss_cum(N_SS)
   !> Sums by this task so far: carbonaceous number (reference, particles) and mass (inventory,
   !> particles, kg); sea salt number (expected, particles) and mass (particles, kg).
   real (kind=RKIND) :: emis_sums(7) = 0.0_RKIND
   real (kind=RKIND) :: drydep_sums(4) = 0.0_RKIND
   !> Number of lowest layers that receive the emission.
   integer :: pmc_emit_levels = 1
   integer :: n_conv = 0
   !> Lagrangian mode: horizontal velocity constant within each layer.
   logical :: lagr_u_const = .false.
   !> Lagrangian mode: velocity from each cell's own mass fluxes.
   integer :: lagr_flux_vel = 0
   !> For lagr_flux_vel = 2, per cell and per edge in counter-clockwise order: the edge, its
   !> length signed so that length*ruAvg is the outward mass flux, the two edge vertices relative
   !> to the cell centre in the tangent plane (m), and the area of the triangle they make with
   !> the centre. rt_n is the number of edges, 0 if the cell lacks an edge on this task.
   integer, allocatable :: rt_n(:), rt_edge(:,:)
   real (kind=RKIND), allocatable :: rt_len(:,:), rt_a(:,:,:), rt_b(:,:,:), rt_area(:,:)
   integer :: n_mixed = 0
   !> Brownian coagulation in every cell each step, and events so far.
   logical :: pmc_coag = .false.
   integer :: n_coag_total = 0
   !> Particle output: interval in steps (0 = off), fraction of cells written, directory.
   integer :: pmc_out_interval = 0
   !> Time-mean aerosol fields: write interval and sampling interval (steps), directory, the
   !> sums (variable, level, cell), the number of samples in them and the group of each species.
   integer :: pmc_diag_interval = 0, pmc_diag_every = 5, n_diag_samples = 0
   character(len=256) :: pmc_diag_dir = 'pmc_diag'
   integer, parameter :: N_DIAG = 21
   real (kind=RKIND), allocatable :: diag_acc(:,:,:)
   integer, allocatable :: diag_group(:)
   !> Restart files for the particle state: write interval (steps), step to read at start-up,
   !> directory.
   integer :: pmc_restart_interval = 0, pmc_restart_read = 0
   character(len=256) :: pmc_restart_dir = 'pmc_restart'
   real (kind=RKIND) :: pmc_out_frac = 1.0_RKIND
   character(len=256) :: pmc_out_dir = 'pmc_out'

   !> Captured transport fluxes and contents of one weight class.
   type flux_set_t
      real (kind=RKIND), dimension(:,:), pointer :: fh => null(), fv => null(), co => null(), cn => null()
   end type flux_set_t

   !> Positions of the particles of one cell (Lagrangian mode): unit vector, zeta, initial unit vector.
   type pos_t
      real (kind=RKIND), allocatable :: x(:,:)
      integer :: n = 0
   end type pos_t

#ifdef MPAS_PARTMC
   !> Aerosol species data.
   type(aero_data_t) :: aero_data
   !> Particle populations, (nVertLevels, nCellsSolve).
   type(aero_state_t), allocatable :: aero_states(:,:)
   !> Particles arriving in each cell during the current transport step.
   type(aero_state_t), allocatable :: aero_incoming(:,:)
   !> Send and receive buffers, one per neighbouring task.
   type(pmc_buf_t), allocatable :: send_buf(:), recv_buf(:)
#endif
   !> Particles of each class in each cell, and arrivals this step: (class, level, cell).
   integer, allocatable :: cnt(:,:,:), cnt_in(:,:,:)
   !> MPI tasks that own cells in this block's halo.
   integer, allocatable :: nbr_rank(:)
   !> For each local cell, position of its owner in nbr_rank (0 if owned here).
   integer, allocatable :: nbr_of_cell(:)
   !> Map from global cell ID to local owned-cell index (0 if not owned here).
   integer, allocatable :: global_to_local(:)
   !> Number of aerosol species.
   integer :: n_spec = 0
   !> Global mass (kg) of each species at initialization.
   real (kind=RKIND), allocatable :: mass_total_init(:)

   ! Lagrangian mode
   !> Particle positions, and positions of arrivals, (nVertLevels, nCellsSolve).
   type(pos_t), allocatable :: pos(:,:), pos_in(:,:)
   !> Real particles represented by a computational particle in weight group LAGR_G0.
   real (kind=RKIND) :: lagr_w = 0.0_RKIND
   !> Unit vectors of the cell centres (3, nCells), zeta of faces and of level centres.
   real (kind=RKIND), allocatable :: cxyz(:,:), zf(:), zc(:)
   !> Radius (on the unit sphere) of a disc around each cell centre that covers the cell.
   real (kind=RKIND), allocatable :: cell_rmax(:)
   real (kind=RKIND) :: radius = 1.0_RKIND
   !> Source-2 second moments at the first diagnostic step.
   real (kind=RKIND) :: m2_ref0 = 0.0_RKIND, m2_part0 = 0.0_RKIND

contains

   !> A name padded to the length PartMC reads from it.
   function pad_name(name)
      character(len=*), intent(in) :: name
      character(len=AERO_SOURCE_NAME_LEN) :: pad_name
      pad_name = name
   end function pad_name


   subroutine read_settings()
      character(len=256) :: envval
      integer :: stat, i
      real (kind=RKIND) :: num_scale

      if (settings_read) return
      settings_read = .true.
      call get_environment_variable('PMC_N_PART', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_n_part
      call get_environment_variable('PMC_COUPLED', envval, status=stat)
      if (stat == 0) pmc_coupled = (trim(envval) /= '0')
      call get_environment_variable('PMC_SCALE_HEIGHT', envval, status=stat)
      if (stat == 0) read(envval, *) scale_height
      call get_environment_variable('PMC_TRANSPORT', envval, status=stat)
      if (stat == 0) then
         if (trim(envval) == 'lagrangian') pmc_mode = MODE_LAGR
      end if
      call get_environment_variable('PMC_TEST_WIND', envval, status=stat)
      if (stat == 0) then
         if (trim(envval) == 'solid') lagr_wind = WIND_SOLID
         if (trim(envval) == 'deform') lagr_wind = WIND_DEFORM
      end if
      call get_environment_variable('PMC_LAGR_CONTROL', envval, status=stat)
      if (stat == 0) lagr_control = (trim(envval) /= '0')
      call get_environment_variable('PMC_NUM_SCALE', envval, status=stat)
      if (stat == 0) then
         read(envval, *) num_scale
         q_background = q_background * num_scale
         q_blob = q_blob * num_scale
      end if
      call get_environment_variable('PMC_COAG', envval, status=stat)
      if (stat == 0) pmc_coag = (trim(envval) /= '0')
      call get_environment_variable('PMC_KZ_MAX', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_kz_max
      call get_environment_variable('PMC_KZ_TOP', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_kz_top
      call get_environment_variable('PMC_LAGR_FLUXVEL', envval, status=stat)
      if (stat == 0) read(envval, *) lagr_flux_vel
      call get_environment_variable('PMC_LAGR_SHUFFLE', envval, status=stat)
      if (stat == 0) lagr_shuffle = (trim(envval) /= '0')
      call get_environment_variable('PMC_LAGR_UCONST', envval, status=stat)
      if (stat == 0) lagr_u_const = (trim(envval) /= '0')
      call get_environment_variable('PMC_MOSAIC', envval, status=stat)
      if (stat == 0) pmc_chem = (trim(envval) /= '0')
      call get_environment_variable('PMC_CHEM_SELFTEST', envval, status=stat)
      if (stat == 0) pmc_chem_selftest = (trim(envval) /= '0')
      ! PMC_CHEM_BOX="id:level,id:level,...", files go to the directory box/
      call get_environment_variable('PMC_CHEM_BOX', envval, status=stat)
      if (stat == 0 .and. len_trim(envval) > 0) then
         n_box = 1
         do i = 1, len_trim(envval)
            if (envval(i:i) == ',') n_box = n_box + 1
            if (envval(i:i) == ':') envval(i:i) = ','
         end do
         n_box = min(n_box, size(box_id))
         read(envval, *) (box_id(i), box_k(i), i = 1, n_box)
      end if
      call get_environment_variable('PMC_CHEM_BOX_STEPS', envval, status=stat)
      if (stat == 0) read(envval, *) box_steps
      call get_environment_variable('PMC_GAS_DATA', envval, status=stat)
      if (stat == 0) gas_data_file = envval
      call get_environment_variable('PMC_GAS_INIT', envval, status=stat)
      if (stat == 0) gas_init_file = envval
      call get_environment_variable('PMC_CHEM_TOP', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_chem_top
      call get_environment_variable('PMC_WET', envval, status=stat)
      if (stat == 0) pmc_wet = (trim(envval) /= '0')
      if (stat == 0) pmc_wet_diag_only = (trim(envval) == '2')
      call get_environment_variable('PMC_CONV_CHECK', envval, status=stat)
      if (stat == 0) pmc_conv_check = (trim(envval) /= '0')
      call get_environment_variable('PMC_WET_CONV', envval, status=stat)
      if (stat == 0) pmc_wet_conv = (trim(envval) /= '0')
      call get_environment_variable('PMC_WET_SCRIT', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_wet_scrit
      call get_environment_variable('PMC_GAS_WET', envval, status=stat)
      if (stat == 0) gas_wet_file = envval
      call get_environment_variable('PMC_GAS_EMIT', envval, status=stat)
      if (stat == 0) gas_emit_file = envval
      call get_environment_variable('PMC_GAS_EMIT_SCALE', envval, status=stat)
      if (stat == 0) read(envval, *) gas_emit_scale
      call get_environment_variable('PMC_GAS_DRYDEP', envval, status=stat)
      if (stat == 0) gas_drydep_file = envval
      call get_environment_variable('PMC_EMIS_FILE', envval, status=stat)
      if (stat == 0) emis_file = envval
      call get_environment_variable('PMC_EMIS_SCALE', envval, status=stat)
      if (stat == 0) read(envval, *) emis_scale
      call get_environment_variable('PMC_CARB_DG', envval, status=stat)
      if (stat == 0) read(envval, *) carb_dg
      call get_environment_variable('PMC_SEASALT', envval, status=stat)
      if (stat == 0) pmc_seasalt = (trim(envval) /= '0')
      call get_environment_variable('PMC_SEASALT_SCALE', envval, status=stat)
      if (stat == 0) read(envval, *) seasalt_scale
      call get_environment_variable('PMC_Q_BACKGROUND', envval, status=stat)
      if (stat == 0) read(envval, *) q_background
      call get_environment_variable('PMC_Q_BLOB', envval, status=stat)
      if (stat == 0) read(envval, *) q_blob
      call get_environment_variable('PMC_EMIT_FLUX', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_emit_flux
      call get_environment_variable('PMC_EMIT_LEVELS', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_emit_levels
      call get_environment_variable('PMC_DRYDEP_SCHEME', envval, status=stat)
      if (stat == 0) pmc_drydep_scheme = (trim(envval) /= '0')
      call get_environment_variable('PMC_DRYDEP_V', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_drydep_v
      call get_environment_variable('PMC_CONV', envval, status=stat)
      if (stat == 0) pmc_conv = (trim(envval) /= '0')
      call get_environment_variable('PMC_KZ_PHYS', envval, status=stat)
      if (stat == 0) pmc_kz_phys = (trim(envval) /= '0')
      call get_environment_variable('PMC_LAGR_SUBSTEPS', envval, status=stat)
      if (stat == 0) read(envval, *) lagr_substeps
      call get_environment_variable('PMC_N_CLASS', envval, status=stat)
      if (stat == 0) read(envval, *) n_class
      if (pmc_mode == MODE_LAGR) n_class = 1
      n_class = max(1, min(MAX_CLASS, n_class))
      call get_environment_variable('PMC_AERO_DATA', envval, status=stat)
      if (stat == 0) aero_data_file = envval
      call get_environment_variable('PMC_RESTART_INTERVAL', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_restart_interval
      call get_environment_variable('PMC_RESTART_READ', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_restart_read
      call get_environment_variable('PMC_RESTART_DIR', envval, status=stat)
      if (stat == 0) pmc_restart_dir = envval
      call get_environment_variable('PMC_DIAG_INTERVAL', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_diag_interval
      call get_environment_variable('PMC_DIAG_EVERY', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_diag_every
      call get_environment_variable('PMC_DIAG_DIR', envval, status=stat)
      if (stat == 0) pmc_diag_dir = envval
      call get_environment_variable('PMC_OUT_INTERVAL', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_out_interval
      call get_environment_variable('PMC_OUT_FRAC', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_out_frac
      call get_environment_variable('PMC_OUT_DIR', envval, status=stat)
      if (stat == 0) pmc_out_dir = envval
   end subroutine read_settings


   !> Scalar index of the number mixing ratio of weight class c.
   integer function class_index(state, c)
      type (mpas_pool_type), pointer :: state
      integer, intent(in) :: c
      integer, pointer :: idx
      if (c == 1) then
         call mpas_pool_get_dimension(state, 'index_pmc_num', idx)
      else
         call mpas_pool_get_dimension(state, 'index_pmc_num2', idx)
      end if
      class_index = idx
   end function class_index


   !> Captured fluxes of weight class c.
   subroutine class_fluxes(diag, c, fs)
      type (mpas_pool_type), pointer :: diag
      integer, intent(in) :: c
      type (flux_set_t), intent(out) :: fs
      if (c == 1) then
         call mpas_pool_get_array(diag, 'pmc_flux_h', fs % fh)
         call mpas_pool_get_array(diag, 'pmc_flux_v', fs % fv)
         call mpas_pool_get_array(diag, 'pmc_content_old', fs % co)
         call mpas_pool_get_array(diag, 'pmc_content_new', fs % cn)
      else
         call mpas_pool_get_array(diag, 'pmc2_flux_h', fs % fh)
         call mpas_pool_get_array(diag, 'pmc2_flux_v', fs % fv)
         call mpas_pool_get_array(diag, 'pmc2_content_old', fs % co)
         call mpas_pool_get_array(diag, 'pmc2_content_new', fs % cn)
      end if
   end subroutine class_fluxes


   !> Set the initial number fields (kg^-1): uniform background plus a Gaussian
   !> blob in the horizontal, unless the field was read in non-zero.
   !> pmc_ref is the total and pmc_ref_src2 the blob; they are never reset from particles.
   !> With one weight class pmc_num is the total; with two, pmc_num is the background
   !> and pmc_num2 the blob.
   subroutine atm_partmc_init(domain)

      type (domain_type), intent(inout) :: domain

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh, state
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      integer, pointer :: index_pmc_num, index_pmc_num2, index_pmc_ref, index_pmc_ref_src2, nCells, nVertLevels
      integer :: iCell, k
      real (kind=RKIND) :: d, cosd, z_mid, f_z, q1, q2
      real (kind=RKIND), dimension(:,:), pointer :: zgrid

      call read_settings()

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block % structs, 'state', state)
         nullify(index_pmc_num)
         call mpas_pool_get_dimension(state, 'index_pmc_num', index_pmc_num)
         if (.not. associated(index_pmc_num)) return
         call mpas_pool_get_dimension(state, 'index_pmc_num2', index_pmc_num2)
         call mpas_pool_get_dimension(state, 'index_pmc_ref', index_pmc_ref)
         call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', index_pmc_ref_src2)
         call mpas_pool_get_dimension(mesh, 'nCells', nCells)
         call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
         call mpas_pool_get_array(mesh, 'latCell', latCell)
         call mpas_pool_get_array(mesh, 'lonCell', lonCell)
         call mpas_pool_get_array(mesh, 'zgrid', zgrid)
         call mpas_pool_get_array(state, 'scalars', scalars, 1)

         if (maxval(scalars(index_pmc_ref,:,1:nCells)) <= 0.0_RKIND) then
            do iCell = 1, nCells
               cosd = sin(blob_lat)*sin(latCell(iCell)) + cos(blob_lat)*cos(latCell(iCell))*cos(lonCell(iCell)-blob_lon)
               d = acos(max(-1.0_RKIND, min(1.0_RKIND, cosd)))
               do k = 1, nVertLevels
                  f_z = 1.0_RKIND
                  if (scale_height > 0.0_RKIND) then
                     z_mid = 0.5_RKIND*(zgrid(k,iCell) + zgrid(k+1,iCell))
                     f_z = exp(-z_mid/scale_height)
                  end if
                  q1 = f_z * q_background
                  q2 = f_z * q_blob*exp(-(d/blob_width)**2)
                  scalars(index_pmc_ref_src2,k,iCell) = q2
                  scalars(index_pmc_ref,k,iCell) = q1 + q2
                  if (n_class == 1) then
                     scalars(index_pmc_num,k,iCell) = q1 + q2
                     scalars(index_pmc_num2,k,iCell) = 0.0_RKIND
                  else
                     scalars(index_pmc_num,k,iCell) = q1
                     scalars(index_pmc_num2,k,iCell) = q2
                  end if
               end do
            end do
            call mpas_log_write('PartMC: initialized number fields with background + blob')
         end if

         block => block % next
      end do

   end subroutine atm_partmc_init


   !> Check that the captured fluxes reproduce the transport update of pmc_num,
   !> and report the largest fraction of a cell's content leaving in one step.
   subroutine atm_partmc_check(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh, state, diag
      real (kind=RKIND), dimension(:,:), pointer :: flux_h, flux_v, c_old, c_new, edgesOnCell_sign
      real (kind=RKIND), dimension(:), pointer :: invAreaCell, areaCell, rdzw
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, pointer :: index_pmc_num, nCellsSolve, nVertLevels
      integer :: iCell, k, i, iEdge
      real (kind=RKIND) :: recon, out, s
      real (kind=RKIND) :: err_max, c_max, out_max, total_old, total_new
      real (kind=RKIND) :: g_err_max, g_c_max, g_out_max, g_total_old, g_total_new
      real (kind=RKIND) :: state_err, g_state_err
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars_new
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz_new

      err_max = 0.0_RKIND
      c_max = 0.0_RKIND
      out_max = 0.0_RKIND
      total_old = 0.0_RKIND
      total_new = 0.0_RKIND
      state_err = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'diag', diag)
         call mpas_pool_get_array(state, 'scalars', scalars_new, 2)
         call mpas_pool_get_array(state, 'rho_zz', rho_zz_new, 2)
         nullify(index_pmc_num)
         call mpas_pool_get_dimension(state, 'index_pmc_num', index_pmc_num)
         if (.not. associated(index_pmc_num)) return
         call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
         call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
         call mpas_pool_get_array(mesh, 'areaCell', areaCell)
         call mpas_pool_get_array(mesh, 'rdzw', rdzw)
         call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
         call mpas_pool_get_array(diag, 'pmc_flux_h', flux_h)
         call mpas_pool_get_array(diag, 'pmc_flux_v', flux_v)
         call mpas_pool_get_array(diag, 'pmc_content_old', c_old)
         call mpas_pool_get_array(diag, 'pmc_content_new', c_new)

         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               recon = c_old(k,iCell) - rdzw(k)*(flux_v(k+1,iCell) - flux_v(k,iCell))
               out = rdzw(k)*(max(0.0_RKIND, flux_v(k+1,iCell)) - min(0.0_RKIND, flux_v(k,iCell)))
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i,iCell)
                  s = edgesOnCell_sign(i,iCell)*flux_h(k,iEdge)*invAreaCell(iCell)
                  recon = recon - s
                  out = out + max(0.0_RKIND, s)
               end do
               err_max = max(err_max, abs(recon - c_new(k,iCell)))
               state_err = max(state_err, abs(rho_zz_new(k,iCell)*scalars_new(index_pmc_num,k,iCell) - c_new(k,iCell)))
               c_max = max(c_max, c_new(k,iCell))
               if (c_old(k,iCell) > 0.0_RKIND) out_max = max(out_max, out/c_old(k,iCell))
               total_old = total_old + c_old(k,iCell)*areaCell(iCell)/rdzw(k)
               total_new = total_new + c_new(k,iCell)*areaCell(iCell)/rdzw(k)
            end do
         end do

         block => block % next
      end do

      call mpas_dmpar_max_real(domain % dminfo, err_max, g_err_max)
      call mpas_dmpar_max_real(domain % dminfo, c_max, g_c_max)
      call mpas_dmpar_max_real(domain % dminfo, out_max, g_out_max)
      call mpas_dmpar_sum_real(domain % dminfo, total_old, g_total_old)
      call mpas_dmpar_sum_real(domain % dminfo, total_new, g_total_new)

      call mpas_dmpar_max_real(domain % dminfo, state_err, g_state_err)

      if (mod(itimestep, 12) == 0 .or. itimestep <= 2) then
         call mpas_log_write('PartMC check step $i: recon rel err = $r, max outflow frac = $r, total rel change = $r', &
                             intArgs=(/itimestep/), &
                             realArgs=(/g_err_max/g_c_max, g_out_max, (g_total_new - g_total_old)/g_total_old/))
         call mpas_log_write('PartMC check step $i: captured vs model state rel diff = $r', &
                             intArgs=(/itimestep/), realArgs=(/g_state_err/g_c_max/))
      end if

   end subroutine atm_partmc_check

   !> Zonal and meridional wind (m/s) of the analytic test flows at time t (s).
   subroutine test_wind(lat, lon, t, u, v)
      real (kind=RKIND), intent(in) :: lat, lon, t
      real (kind=RKIND), intent(out) :: u, v
      real (kind=RKIND) :: p(3), a(3), vel(3), lon_p, kk

      if (lagr_wind == WIND_SOLID) then
         ! rotation about an axis tilted 45 degrees from the pole
         a = (/ sqrt(0.5_RKIND), 0.0_RKIND, sqrt(0.5_RKIND) /)
         p = (/ cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat) /)
         vel(1) = a(2)*p(3) - a(3)*p(2)
         vel(2) = a(3)*p(1) - a(1)*p(3)
         vel(3) = a(1)*p(2) - a(2)*p(1)
         vel = vel * radius * 2.0_RKIND*pii/test_period
         u = -sin(lon)*vel(1) + cos(lon)*vel(2)
         v = -sin(lat)*cos(lon)*vel(1) - sin(lat)*sin(lon)*vel(2) + cos(lat)*vel(3)
      else
         ! Nair and Lauritzen (2010) case 4: deformation that reverses at half the period,
         ! on top of a zonal rotation
         kk = 10.0_RKIND * radius / test_period
         lon_p = lon - 2.0_RKIND*pii*t/test_period
         u = kk * sin(lon_p)**2 * sin(2.0_RKIND*lat) * cos(pii*t/test_period) &
             + 2.0_RKIND*pii*radius/test_period * cos(lat)
         v = kk * sin(2.0_RKIND*lon_p) * cos(lat) * cos(pii*t/test_period)
      end if
   end subroutine test_wind


   !> Unit vector p rotated by angle about the axis of the solid-body test flow.
   pure function solid_rotate(p, angle) result(q)
      real (kind=RKIND), intent(in) :: p(3), angle
      real (kind=RKIND) :: q(3), a(3), axp(3)

      a = (/ sqrt(0.5_RKIND), 0.0_RKIND, sqrt(0.5_RKIND) /)
      axp(1) = a(2)*p(3) - a(3)*p(2)
      axp(2) = a(3)*p(1) - a(1)*p(3)
      axp(3) = a(1)*p(2) - a(2)*p(1)
      q = p*cos(angle) + axp*sin(angle) + a*dot_product(a, p)*(1.0_RKIND - cos(angle))
   end function solid_rotate


   !> Physical volume (m^3) of cell (k,iCell).
   pure real (kind=RKIND) function cell_volume(areaCell, rdzw, zz)
      real (kind=RKIND), intent(in) :: areaCell, rdzw, zz
      cell_volume = areaCell / (rdzw * zz)
   end function cell_volume


   subroutine pos_append(p, x)
      type (pos_t), intent(inout) :: p
      real (kind=RKIND), intent(in) :: x(NPOS)
      real (kind=RKIND), allocatable :: tmp(:,:)

      if (.not. allocated(p % x)) then
         allocate(p % x(NPOS, 16))
      else if (p % n == size(p % x, 2)) then
         allocate(tmp(NPOS, 2*p % n))
         tmp(:, 1:p % n) = p % x(:, 1:p % n)
         call move_alloc(tmp, p % x)
      end if
      p % n = p % n + 1
      p % x(:, p % n) = x
   end subroutine pos_append


   !> Remove entry i by moving the last entry into its place, as PartMC does for particles.
   subroutine pos_remove(p, i)
      type (pos_t), intent(inout) :: p
      integer, intent(in) :: i
      if (i < p % n) p % x(:, i) = p % x(:, p % n)
      p % n = p % n - 1
   end subroutine pos_remove


   !> Local cell whose centre is nearest to unit vector p, searching from cell_hint.
   !> Geometry for the edge-flux-matching velocity (lagr_flux_vel = 2); see rt_n.
   subroutine lagr_rt_setup(mesh)
      type (mpas_pool_type), intent(in) :: mesh
      integer, pointer :: nCells, nEdges, nVertices
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell, verticesOnEdge, cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex, dvEdge
      integer :: c, j, i, n, e, iv, order(12), jj, tmp
      real (kind=RKIND) :: nrm(3), va(3,12), vb(3,12), ang(12), e1(3), e2(3), mid(3), t(3), sgn(12)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'xVertex', xVertex)
      call mpas_pool_get_array(mesh, 'yVertex', yVertex)
      call mpas_pool_get_array(mesh, 'zVertex', zVertex)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      allocate(rt_n(nCells), rt_edge(12,nCells), rt_len(12,nCells), rt_a(3,12,nCells), rt_b(3,12,nCells), &
               rt_area(12,nCells))
      rt_n = 0
      do c = 1, nCells
         n = nEdgesOnCell(c)
         nrm = cxyz(:,c)
         do j = 1, n
            e = edgesOnCell(j,c)
            if (e < 1 .or. e > nEdges) exit
            if (minval(verticesOnEdge(:,e)) < 1 .or. maxval(verticesOnEdge(:,e)) > nVertices) exit
            do i = 1, 2
               iv = verticesOnEdge(i,e)
               t = (/ xVertex(iv), yVertex(iv), zVertex(iv) /)
               t = radius * (t / norm2(t) - nrm)
               t = t - sum(t*nrm)*nrm
               if (i == 1) va(:,j) = t
               if (i == 2) vb(:,j) = t
            end do
            ! counter-clockwise seen from outside the sphere
            if (sum(nrm * cross3(va(:,j), vb(:,j))) < 0.0_RKIND) then
               t = va(:,j)
               va(:,j) = vb(:,j)
               vb(:,j) = t
            end if
            sgn(j) = 1.0_RKIND
            if (cellsOnEdge(2,e) == c) sgn(j) = -1.0_RKIND
         end do
         if (j <= n) cycle
         e1 = va(:,1) / norm2(va(:,1))
         e2 = cross3(nrm, e1)
         do j = 1, n
            mid = 0.5_RKIND*(va(:,j) + vb(:,j))
            ang(j) = atan2(sum(mid*e2), sum(mid*e1))
            order(j) = j
         end do
         do j = 2, n
            do jj = j, 2, -1
               if (ang(order(jj)) >= ang(order(jj-1))) exit
               tmp = order(jj)
               order(jj) = order(jj-1)
               order(jj-1) = tmp
            end do
         end do
         do i = 1, n
            j = order(i)
            rt_edge(i,c) = edgesOnCell(j,c)
            rt_len(i,c) = sgn(j) * dvEdge(edgesOnCell(j,c))
            rt_a(:,i,c) = va(:,j)
            rt_b(:,i,c) = vb(:,j)
            rt_area(i,c) = 0.5_RKIND * sum(nrm * cross3(va(:,j), vb(:,j)))
         end do
         rt_n(c) = n
      end do

   contains

      pure function cross3(a, b) result(c3)
         real (kind=RKIND), intent(in) :: a(3), b(3)
         real (kind=RKIND) :: c3(3)
         c3 = (/ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) /)
      end function cross3

   end subroutine lagr_rt_setup


   integer function find_cell(p, cell_hint, nCells, nEdgesOnCell, cellsOnCell)
      real (kind=RKIND), intent(in) :: p(3)
      integer, intent(in) :: cell_hint, nCells
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: cellsOnCell
      integer :: best, j, cj
      real (kind=RKIND) :: d_best, d
      logical :: moved

      best = cell_hint
      d_best = dot_product(p, cxyz(:,best))
      moved = .true.
      do while (moved)
         moved = .false.
         do j = 1, nEdgesOnCell(best)
            cj = cellsOnCell(j,best)
            if (cj < 1 .or. cj > nCells) cycle
            d = dot_product(p, cxyz(:,cj))
            if (d > d_best) then
               d_best = d
               best = cj
               moved = .true.
               exit
            end if
         end do
      end do
      find_cell = best
   end function find_cell


#ifdef MPAS_PARTMC
   !> Uniform random position (unit vector and zeta) in cell iCell, level k.
   subroutine random_in_cell(iCell, k, nCells, nEdgesOnCell, cellsOnCell, x)
      integer, intent(in) :: iCell, k, nCells
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: cellsOnCell
      real (kind=RKIND), intent(out) :: x(4)
      real (kind=RKIND) :: v(3), e1(3), e2(3), a, b

      ! tangent-plane basis at the cell centre
      v = (/ 0.0_RKIND, 0.0_RKIND, 1.0_RKIND /)
      if (abs(cxyz(3,iCell)) > 0.9_RKIND) v = (/ 1.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
      e1 = v - dot_product(v, cxyz(:,iCell)) * cxyz(:,iCell)
      e1 = e1 / norm2(e1)
      e2(1) = cxyz(2,iCell)*e1(3) - cxyz(3,iCell)*e1(2)
      e2(2) = cxyz(3,iCell)*e1(1) - cxyz(1,iCell)*e1(3)
      e2(3) = cxyz(1,iCell)*e1(2) - cxyz(2,iCell)*e1(1)
      do
         do
            a = 2.0_RKIND*pmc_random() - 1.0_RKIND
            b = 2.0_RKIND*pmc_random() - 1.0_RKIND
            if (a*a + b*b <= 1.0_RKIND) exit
         end do
         x(1:3) = cxyz(:,iCell) + cell_rmax(iCell)*(a*e1 + b*e2)
         x(1:3) = x(1:3) / norm2(x(1:3))
         if (find_cell(x(1:3), iCell, nCells, nEdgesOnCell, cellsOnCell) == iCell) exit
      end do
      x(4) = zf(k) + pmc_random()*(zf(k+1) - zf(k))
   end subroutine random_in_cell
#endif


   !> Level whose zeta range contains z.
   pure integer function find_level(z)
      real (kind=RKIND), intent(in) :: z
      integer :: k
      find_level = size(zc)
      do k = 1, size(zc)
         if (z < zf(k+1)) then
            find_level = k
            exit
         end if
      end do
   end function find_level


   !> Create a particle population in every owned cell, consistent with the number fields.
   subroutine atm_partmc_particles_init(domain)

      type (domain_type), intent(inout) :: domain

#ifdef MPAS_PARTMC
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh, state, diag
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, zz
      real (kind=RKIND), dimension(:), pointer :: areaCell, rdzw, xCell, yCell, zCell, xVertex, yVertex, zVertex
      real (kind=RKIND), pointer :: sphere_radius
      integer, dimension(:), pointer :: indexToCellID, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, verticesOnCell
      integer, pointer :: index_ref, index_src2, nCells, nCellsSolve, nVertLevels
      integer :: iCell, k, i, j, g, rank, n, n2, n_t, i_so4, i_nh4, i_oc, i_bc, i_src, g_n_cells
      type (field1DReal), pointer :: owner_field
      real (kind=RKIND), dimension(:), pointer :: owner
      type(spec_file_t) :: file
      real (kind=dp), allocatable :: vols(:)
      real (kind=RKIND) :: conc_tot, conc1, conc2, vol, f2, total, g_total, x(NPOS), v(3)
      type(aero_particle_t) :: particle

      call read_settings()
      if (pmc_n_part <= 0) then
         call mpas_log_write('PartMC: particles disabled (set PMC_N_PART > 0 to enable)')
         return
      end if
      if (pmc_out_interval > 0) call execute_command_line('mkdir -p ' // trim(pmc_out_dir))

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_dimension(state, 'index_pmc_ref', index_ref)
      call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', index_src2)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)

      ! species from a PartMC aero_data file; two sources with different composition
      call spec_file_open(trim(aero_data_file), file)
      call spec_file_read_aero_data(file, aero_data)
      call spec_file_close(file)
      call fractal_set_spherical(aero_data%fractal)
      i_src = aero_data_source_by_name(aero_data, pad_name('background'))
      i_src = aero_data_source_by_name(aero_data, pad_name('blob'))
      i_src = aero_data_source_by_name(aero_data, pad_name('carbon'))
      i_src = aero_data_source_by_name(aero_data, pad_name('seasalt'))
      if (n_class == 1) then
         i_src = aero_data_weight_class_by_name(aero_data, 'all')
      else
         i_src = aero_data_weight_class_by_name(aero_data, pad_name('background'))
         i_src = aero_data_weight_class_by_name(aero_data, pad_name('blob'))
      end if
      n_spec = aero_data_n_spec(aero_data)
      i_so4 = aero_data_spec_by_name(aero_data, 'SO4')
      i_nh4 = aero_data_spec_by_name(aero_data, 'NH4')
      i_oc = aero_data_spec_by_name(aero_data, 'OC')
      i_bc = aero_data_spec_by_name(aero_data, 'BC')
      if (min(i_so4, i_nh4, i_oc, i_bc) < 1) then
         call mpas_log_write('PartMC: aero_data must have species SO4, NH4, OC and BC', messageType=MPAS_LOG_CRIT)
      end if
      allocate(vols(n_spec))

      call pmc_srand(20261111, domain % dminfo % my_proc_id)

      ! find which task owns each halo cell by exchanging a field holding the owner's rank
      call mpas_pool_get_field(diag, 'pmc_owner', owner_field)
      call mpas_pool_get_array(diag, 'pmc_owner', owner)
      call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)
      owner(:) = -1.0_RKIND
      owner(1:nCellsSolve) = real(domain % dminfo % my_proc_id, kind=RKIND)
      call mpas_dmpar_exch_halo_field(owner_field)
      allocate(nbr_of_cell(nCells))
      allocate(nbr_rank(0))
      nbr_of_cell = 0
      do iCell = nCellsSolve + 1, nCells
         rank = nint(owner(iCell))
         if (rank < 0) call mpas_log_write('PartMC: halo cell $i has no owner', intArgs=(/iCell/), messageType=MPAS_LOG_CRIT)
         do i = 1, size(nbr_rank)
            if (nbr_rank(i) == rank) exit
         end do
         if (i > size(nbr_rank)) nbr_rank = [nbr_rank, rank]
         nbr_of_cell(iCell) = i
      end do
      allocate(send_buf(size(nbr_rank)), recv_buf(size(nbr_rank)))
      allocate(global_to_local(maxval(indexToCellID(1:nCells))))
      global_to_local = 0
      do iCell = 1, nCellsSolve
         global_to_local(indexToCellID(iCell)) = iCell
      end do

      allocate(aero_states(nVertLevels, nCellsSolve))
      allocate(aero_incoming(nVertLevels, nCellsSolve))
      allocate(cnt(n_class, nVertLevels, nCellsSolve), cnt_in(n_class, nVertLevels, nCellsSolve))
      cnt = 0
      cnt_in = 0

      if (pmc_mode == MODE_LAGR) then
         call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)
         call mpas_pool_get_array(mesh, 'xCell', xCell)
         call mpas_pool_get_array(mesh, 'yCell', yCell)
         call mpas_pool_get_array(mesh, 'zCell', zCell)
         call mpas_pool_get_array(mesh, 'xVertex', xVertex)
         call mpas_pool_get_array(mesh, 'yVertex', yVertex)
         call mpas_pool_get_array(mesh, 'zVertex', zVertex)
         radius = sphere_radius
         allocate(cxyz(3, nCells), zf(nVertLevels+1), zc(nVertLevels))
         do iCell = 1, nCells
            cxyz(:,iCell) = (/ xCell(iCell), yCell(iCell), zCell(iCell) /)
            cxyz(:,iCell) = cxyz(:,iCell) / norm2(cxyz(:,iCell))
         end do
         zf(1) = 0.0_RKIND
         do k = 1, nVertLevels
            zf(k+1) = zf(k) + 1.0_RKIND/rdzw(k)
            zc(k) = 0.5_RKIND*(zf(k) + zf(k+1))
         end do
         allocate(pos(nVertLevels, nCellsSolve), pos_in(nVertLevels, nCellsSolve))
         if (lagr_flux_vel == 2) call lagr_rt_setup(mesh)
         allocate(cell_rmax(nCells))
         cell_rmax = 0.0_RKIND
         do iCell = 1, nCellsSolve
            do j = 1, nEdgesOnCell(iCell)
               i = verticesOnCell(j,iCell)
               v = (/ xVertex(i), yVertex(i), zVertex(i) /)
               v = v / norm2(v)
               cell_rmax(iCell) = max(cell_rmax(iCell), 1.05_RKIND * norm2(v - cxyz(:,iCell)))
            end do
         end do

         ! every computational particle represents the same number of real particles
         total = 0.0_RKIND
         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               total = total + rho_zz(k,iCell) * zz(k,iCell) * scalars(index_ref,k,iCell) &
                               * cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))
            end do
         end do
         call mpas_dmpar_sum_real(domain % dminfo, total, g_total)
         call mpas_dmpar_sum_int(domain % dminfo, nCellsSolve*nVertLevels, g_n_cells)
         lagr_w = g_total / (real(pmc_n_part, kind=RKIND) * real(g_n_cells, kind=RKIND))
      end if

      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            call aero_state_zero(aero_states(k,iCell))
            call aero_state_set_weight(aero_states(k,iCell), aero_data, AERO_STATE_WEIGHT_FLAT_SPECIFIED)
            call aero_state_set_n_part_ideal(aero_states(k,iCell), real(pmc_n_part, kind=dp))
            call aero_state_zero(aero_incoming(k,iCell))
            call aero_state_set_weight(aero_incoming(k,iCell), aero_data, AERO_STATE_WEIGHT_FLAT_SPECIFIED)
            ! number concentration (m^-3) = rho * q = rho_zz * zz * q
            conc_tot = rho_zz(k,iCell) * zz(k,iCell) * scalars(index_ref,k,iCell)
            conc2 = rho_zz(k,iCell) * zz(k,iCell) * scalars(index_src2,k,iCell)
            conc1 = conc_tot - conc2
            vol = cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))

            if (pmc_mode == MODE_LAGR) then
               ! one weight per class, a factor of two apart
               call aero_weight_array_set_sizes(aero_states(k,iCell)%awa, 1, LAGR_GROUPS)
               aero_states(k,iCell)%awa%weight%type = AERO_WEIGHT_TYPE_NONE
               aero_states(k,iCell)%awa%weight%exponent = 0.0d0
               do g = 1, LAGR_GROUPS
                  aero_states(k,iCell)%awa%weight(1,g)%magnitude = lagr_w * 2.0d0**(g - LAGR_G0) / vol
               end do
               g = LAGR_G0
               if (lagr_control) then
                  ! the group that puts the count nearest the target in this cell
                  g = LAGR_G0 + nint(log(conc_tot * vol / (real(pmc_n_part, kind=RKIND) * lagr_w)) / log(2.0_RKIND))
                  g = max(1, min(LAGR_GROUPS, g))
               end if
               f2 = conc_tot * vol / (lagr_w * 2.0_RKIND**(g - LAGR_G0))
               n = int(f2)
               if (pmc_random() < f2 - real(n, kind=RKIND)) n = n + 1
               do i = 1, n
                  if (pmc_random() < conc2/conc_tot) then
                     call make_particle(2, g, 1)
                  else
                     call make_particle(1, g, 1)
                  end if
                  call random_in_cell(iCell, k, nCells, nEdgesOnCell, cellsOnCell, x(1:4))
                  x(5:7) = x(1:3)
                  call pos_append(pos(k,iCell), x)
               end do
               cnt(1,k,iCell) = n
            else if (n_class == 1) then
               aero_states(k,iCell)%awa%weight(1,1)%magnitude = conc_tot / real(pmc_n_part, kind=dp)
               ! number of particles from the second source, with randomized rounding
               f2 = real(pmc_n_part, kind=dp) * conc2 / conc_tot
               n2 = int(f2)
               if (pmc_random() < f2 - real(n2, kind=dp)) n2 = n2 + 1
               do i = 1, pmc_n_part
                  if (i <= n2) then
                     call make_particle(2, 1, 1)
                  else
                     call make_particle(1, 1, 1)
                  end if
               end do
               cnt(1,k,iCell) = pmc_n_part
            else
               ! one weight class per source, each with half of the particles
               n_t = pmc_n_part / n_class
               aero_states(k,iCell)%awa%weight(1,1)%magnitude = conc1 / real(n_t, kind=dp)
               do i = 1, n_t
                  call make_particle(1, 1, 1)
               end do
               cnt(1,k,iCell) = n_t
               if (conc2 > 1.0e-3_RKIND * conc1) then
                  aero_states(k,iCell)%awa%weight(1,2)%magnitude = conc2 / real(n_t, kind=dp)
                  do i = 1, n_t
                     call make_particle(2, 2, 1)
                  end do
                  cnt(2,k,iCell) = n_t
               else
                  aero_states(k,iCell)%awa%weight(1,2)%magnitude = 0.0d0
               end if
            end if
         end do
      end do

      call mpas_log_write('PartMC: target $i particles in each of $i cells', &
                          intArgs=(/pmc_n_part, nVertLevels*nCellsSolve/))
      call mpas_log_write('PartMC: $i species from ' // trim(aero_data_file), intArgs=(/n_spec/))
      allocate(mass_total_init(n_spec))
      call species_mass_totals(domain, mass_total_init)
      if (pmc_out_interval > 0) then
         call mpas_log_write('PartMC: particle output every $i steps for cell fraction $r to ' // trim(pmc_out_dir), &
                             intArgs=(/pmc_out_interval/), realArgs=(/pmc_out_frac/))
      end if
      if (pmc_mode == MODE_LAGR) then
         call mpas_log_write('PartMC: Lagrangian transport, each particle represents $r real particles', &
                             realArgs=(/lagr_w/))
      else if (pmc_coupled) then
         call mpas_log_write('PartMC: flux transport with $i weight classes, coupled', intArgs=(/n_class/))
      else
         call mpas_log_write('PartMC: flux transport with $i weight classes, one-way', intArgs=(/n_class/))
      end if

      if (pmc_restart_interval > 0) call execute_command_line('mkdir -p ' // trim(pmc_restart_dir))
      if (pmc_restart_read > 0) call restart_read(domain)

   contains

      !> Add a new particle from source i_source to cell (k,iCell) in weight class i_class.
      subroutine make_particle(i_source, i_class, i_group)
         integer, intent(in) :: i_source, i_class, i_group
         real (kind=dp) :: diam, vol_tot

         call aero_particle_zero(particle, aero_data)
         vols = 0.0d0
         if (i_source == 2) then
            ! blob: soot-like, lognormal with geometric mean diameter 50 nm
            diam = 0.05d-6 * exp(log(1.8d0) * rand_normal(0.0d0, 1.0d0))
            vol_tot = pii / 6.0d0 * diam**3
            vols(i_bc) = 0.6d0 * vol_tot
            vols(i_oc) = 0.4d0 * vol_tot
         else
            ! background: ammonium sulfate and organics, geometric mean diameter 100 nm
            diam = 0.1d-6 * exp(log(1.6d0) * rand_normal(0.0d0, 1.0d0))
            vol_tot = pii / 6.0d0 * diam**3
            vols(i_so4) = 0.5d0 * vol_tot
            vols(i_nh4) = 0.2d0 * vol_tot
            vols(i_oc) = 0.3d0 * vol_tot
         end if
         call aero_particle_set_component(particle, i_source, 0.0d0)
         particle%n_primary_parts = 1
         call aero_particle_set_vols(particle, vols)
         call aero_particle_new_id(particle)
         call aero_particle_set_weight(particle, i_group, i_class)
         call aero_state_add_particle(aero_states(k,iCell), particle, aero_data)
      end subroutine make_particle
#endif

   end subroutine atm_partmc_particles_init


   !> Name of this task's particle restart file for a step.
   function restart_name(domain, step)
      type (domain_type), intent(in) :: domain
      integer, intent(in) :: step
      character(len=512) :: restart_name
      write(restart_name, '(a,"/particles_s",i8.8,"_r",i5.5,".bin")') trim(pmc_restart_dir), step, &
            domain % dminfo % my_proc_id
   end function restart_name


   !> Write the complete particle state of this task: every cell's PartMC population (packed
   !> with PartMC's own routines, so nothing is left out), the particle positions, the state
   !> of the random number generator, the particle id counter, the weight unit and the
   !> running sums of the diagnostics. Together with the model's own restart file this is
   !> what a run needs to continue as if it had not stopped.
   subroutine atm_partmc_restart_write(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC
      type (mpas_pool_type), pointer :: mesh
      integer, dimension(:), pointer :: indexToCellID
      integer, pointer :: nCellsSolve, nVertLevels
      integer :: u, iCell, k, n, nb, p, n_seed, has_gas
      integer, allocatable :: seed(:)
      character, allocatable :: buf(:)
      integer (kind=I8KIND) :: n_tot

      if (pmc_n_part <= 0 .or. pmc_restart_interval <= 0) return
      if (mod(itimestep, pmc_restart_interval) /= 0) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)
      call random_seed(size=n_seed)
      allocate(seed(n_seed))
      call random_seed(get=seed)

      open(newunit=u, file=trim(restart_name(domain, itimestep)), access='stream', form='unformatted', status='replace')
      write(u) 20261113, itimestep, nCellsSolve, nVertLevels, NPOS, pmc_mode, n_seed
      write(u) indexToCellID(1:nCellsSolve)
      write(u) seed
      write(u) next_id
      write(u) lagr_w, chem_s0, part_s_removed, m2_ref0, m2_part0, surf_sums, wet_sums, conv_wet
      write(u) n_thinned, n_split, n_conv, n_mixed, n_coag_total
      has_gas = merge(size(gas_emitted), 0, allocated(gas_emitted))
      write(u) has_gas
      if (has_gas > 0) write(u) gas_emitted, gas_deposited, gas_wetdep, gas_content0
      n_tot = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            n = aero_state_n_part(aero_states(k,iCell))
            nb = pmc_mpi_pack_size_aero_state(aero_states(k,iCell))
            if (allocated(buf)) then
               if (size(buf) < nb) deallocate(buf)
            end if
            if (.not. allocated(buf)) allocate(buf(2*nb))
            p = 0
            call pmc_mpi_pack_aero_state(buf, p, aero_states(k,iCell))
            write(u) n, p
            write(u) buf(1:p)
            if (pmc_mode == MODE_LAGR .and. n > 0) write(u) pos(k,iCell) % x(:,1:n)
            n_tot = n_tot + n
         end do
      end do
      write(u) n_tot
      close(u)
      if (domain % dminfo % my_proc_id == 0) then
         call mpas_log_write('PartMC restart step $i: particle state written to ' // trim(pmc_restart_dir), &
                             intArgs=(/itimestep/))
      end if
#endif

   end subroutine atm_partmc_restart_write


   !> Read the particle state written by atm_partmc_restart_write, replacing the freshly
   !> made populations. The run must have the same number of tasks and the same partition.
   subroutine restart_read(domain)

      type (domain_type), intent(inout) :: domain

#ifdef MPAS_PARTMC
      type (mpas_pool_type), pointer :: mesh
      integer, dimension(:), pointer :: indexToCellID
      integer, pointer :: nCellsSolve, nVertLevels
      integer :: u, iCell, k, n, nb, p, n_seed, has_gas, hdr(7), ios, n_in
      integer, allocatable :: seed(:), ids(:)
      character, allocatable :: buf(:)
      integer (kind=I8KIND) :: n_tot, n_chk
      real (kind=RKIND) :: x(NPOS)
      real (kind=RKIND), allocatable :: xs(:,:)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)

      open(newunit=u, file=trim(restart_name(domain, pmc_restart_read)), access='stream', form='unformatted', &
           status='old', iostat=ios)
      if (ios /= 0) then
         call mpas_log_write('PartMC: cannot open ' // trim(restart_name(domain, pmc_restart_read)), &
                             messageType=MPAS_LOG_CRIT)
      end if
      read(u) hdr
      call random_seed(size=n_seed)
      if (hdr(1) /= 20261113 .or. hdr(3) /= nCellsSolve .or. hdr(4) /= nVertLevels .or. hdr(5) /= NPOS &
          .or. hdr(6) /= pmc_mode .or. hdr(7) /= n_seed) then
         call mpas_log_write('PartMC: the particle restart file does not fit this run (tasks, partition, ' &
                             // 'levels or transport mode differ)', messageType=MPAS_LOG_CRIT)
      end if
      allocate(ids(nCellsSolve), seed(n_seed))
      read(u) ids
      if (any(ids /= indexToCellID(1:nCellsSolve))) then
         call mpas_log_write('PartMC: the particle restart file is for a different partition', &
                             messageType=MPAS_LOG_CRIT)
      end if
      read(u) seed
      read(u) next_id
      read(u) lagr_w, chem_s0, part_s_removed, m2_ref0, m2_part0, surf_sums, wet_sums, conv_wet
      read(u) n_thinned, n_split, n_conv, n_mixed, n_coag_total
      read(u) has_gas
      if (has_gas > 0) then
         if (allocated(gas_emitted)) deallocate(gas_emitted, gas_deposited, gas_wetdep, gas_content0)
         allocate(gas_emitted(has_gas), gas_deposited(has_gas), gas_wetdep(has_gas), gas_content0(has_gas))
         read(u) gas_emitted, gas_deposited, gas_wetdep, gas_content0
      end if
      n_tot = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            read(u) n, nb
            if (allocated(buf)) then
               if (size(buf) < nb) deallocate(buf)
            end if
            if (.not. allocated(buf)) allocate(buf(2*nb))
            read(u) buf(1:nb)
            p = 0
            ! PartMC's unpack needs an empty population to fill
            call aero_state_zero(aero_states(k,iCell))
            call pmc_mpi_unpack_aero_state(buf, p, aero_states(k,iCell))
            if (aero_state_n_part(aero_states(k,iCell)) /= n) then
               call mpas_log_write('PartMC: particle restart file is inconsistent', messageType=MPAS_LOG_CRIT)
            end if
            if (pmc_mode == MODE_LAGR) then
               pos(k,iCell) % n = 0
               if (n > 0) then
                  allocate(xs(NPOS, n))
                  read(u) xs
                  do n_in = 1, n
                     x = xs(:,n_in)
                     call pos_append(pos(k,iCell), x)
                  end do
                  deallocate(xs)
               end if
            end if
            n_tot = n_tot + n
         end do
      end do
      read(u) n_chk
      close(u)
      if (n_chk /= n_tot) then
         call mpas_log_write('PartMC: particle restart file is truncated', messageType=MPAS_LOG_CRIT)
      end if
      ! last, so that nothing above has drawn random numbers after the state was set
      call random_seed(put=seed)
      if (domain % dminfo % my_proc_id == 0) then
         call mpas_log_write('PartMC: particle state read from ' // trim(pmc_restart_dir) // ', step $i', &
                             intArgs=(/pmc_restart_read/))
      end if
#endif

   end subroutine restart_read


   !> Coupled flux mode: set each class scalar (time level 1) from the particle counts, so that
   !> the host's scalar transport computes fluxes for the fields the particles represent.
   subroutine atm_partmc_particles_to_scalar(domain)

      type (domain_type), intent(inout) :: domain

#ifdef MPAS_PARTMC
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh, state
      type (field3DReal), pointer :: scalars_field
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, zz
      integer, pointer :: nCellsSolve, nVertLevels
      integer :: iCell, k, c, idx

      if (pmc_n_part <= 0 .or. .not. pmc_coupled .or. pmc_mode /= MODE_FLUX) return

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)

      do c = 1, n_class
         idx = class_index(state, c)
         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               ! number mixing ratio = number concentration / density, density = rho_zz * zz
               scalars(idx,k,iCell) = real(cnt(c,k,iCell), kind=RKIND) &
                    * aero_states(k,iCell)%awa%weight(1,c)%magnitude / (rho_zz(k,iCell) * zz(k,iCell))
            end do
         end do
      end do

      call mpas_pool_get_field(state, 'scalars', scalars_field, 1)
      call mpas_dmpar_exch_halo_field(scalars_field)
#endif

   end subroutine atm_partmc_particles_to_scalar


   !> Move particles between cells.
   !>
   !> Flux mode: each particle of weight class c in a cell leaves through face f with
   !> probability (outgoing flux of class c through f) / (class content of the cell before
   !> transport). A particle arriving in cell d becomes (W_c / W_d) particles on average,
   !> where W is the number of real particles represented by one computational particle.
   !>
   !> Lagrangian mode: each particle has a position, which is advanced with the
   !> step-averaged transport velocity (mass flux / density) interpolated to it.
   subroutine atm_partmc_particles_transport(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC
      integer, parameter :: MAX_FACES = 12
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh, diag, state
      type (flux_set_t) :: fs(MAX_CLASS)
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign, zz, rho_zz_new, rho_zz_old
      real (kind=RKIND), dimension(:,:), pointer :: ux, uy, uz, om, ruAvg, wwAvg, hdiv
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars_new
      real (kind=RKIND), dimension(:), pointer :: invAreaCell, areaCell, rdzw, fzm, fzp
      real (kind=RKIND), pointer :: config_dt
      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, indexToCellID
      integer, pointer :: nCells, nCellsSolve, nEdges, nVertLevels, index_ref, index_src2
      integer :: iCell, k, i, c, i_part, i_nbr, n_bytes, dest_global, comm, i_spec, k_worst
      integer :: n_part_total, n_part_min, n_part_max, n_moved, n_over, n_copy_max, n_clamp, n_t, n_thin, g_n_thin
      integer :: g_n_thinned, g_n_split
      integer :: g_n_part_total, g_n_part_min, g_n_part_max, g_n_moved, g_n_cells, g_n_over, g_n_copy_max, g_n_clamp
      integer :: counts(MAX_CLASS)
      real (kind=RKIND) :: w_src, vol, m, x(NPOS), xe(3), d
      real (kind=RKIND) :: terr(3), g_terr(3), terr_max, g_terr_max
      real (kind=RKIND), dimension(:), pointer :: angleEdge, latEdge, lonEdge
      real (kind=RKIND) :: c_part, c_ref, c_new, v_part, c1_part, v1_part, c1_ref, c2_part, v2_part, c2_ref
      real (kind=RKIND) :: p_sum_max, lost, ratio_max, ratio_min, lev_mean, lev_worst, mass_dev
      real (kind=RKIND) :: g_p_sum_max, g_lost, g_ratio_max, g_ratio_min
      ! sums: 1 err2, 2 norm2, 3 noise2 (total); 4-6 same for source 2; 7-9 same for source 1;
      ! 10 total from particles, 11 total reference, 12 one-step err2, 13 source-2 particles, 14 source-2 reference
      real (kind=RKIND) :: sums(16), g_sums(16), pruned, g_pruned
      real (kind=RKIND), allocatable :: err2_k(:), norm2_k(:), g_err2_k(:), g_norm2_k(:), mass_total(:)
      real (kind=RKIND), allocatable :: lev_sums(:,:), g_lev_sums(:,:)
      real (kind=RKIND) :: bias, bias_worst, bias_sigma, lev_all
      integer :: k_bias
      logical :: do_diag
      type(aero_particle_t) :: particle

      if (pmc_n_part <= 0) return

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_config(block % configs, 'config_dt', config_dt)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(state, 'scalars', scalars_new, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_new, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_old, 1)
      call mpas_pool_get_dimension(state, 'index_pmc_ref', index_ref)
      call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', index_src2)
      do c = 1, MAX_CLASS
         call class_fluxes(diag, c, fs(c))
      end do

      do_diag = (mod(itimestep, 12) == 0 .or. itimestep <= 2)
      n_moved = 0
      n_over = 0
      n_copy_max = 0
      n_clamp = 0
      n_thin = 0
      pruned = 0.0_RKIND
      p_sum_max = 0.0_RKIND
      lost = 0.0_RKIND
      ratio_max = 0.0_RKIND
      ratio_min = huge(1.0_RKIND)

      if (pmc_mode == MODE_LAGR) then
         call lagr_velocities()
         call lagr_move()
      else
         call flux_move()
      end if

      ! exchange particles that crossed into cells owned by other tasks
      if (size(nbr_rank) > 0) then
#ifdef MPAS_USE_MPI_F08
         comm = domain % dminfo % comm % MPI_VAL
#else
         comm = domain % dminfo % comm
#endif
         call pmc_comm_exchange(comm, nbr_rank, send_buf, recv_buf)
         do i_nbr = 1, size(nbr_rank)
            do while (recv_buf(i_nbr)%pos < size(recv_buf(i_nbr)%data))
               call pmc_mpi_unpack_integer(recv_buf(i_nbr)%data, recv_buf(i_nbr)%pos, dest_global)
               call pmc_mpi_unpack_integer(recv_buf(i_nbr)%data, recv_buf(i_nbr)%pos, k)
               if (pmc_mode == MODE_LAGR) then
                  do i = 1, NPOS
                     call pmc_mpi_unpack_real(recv_buf(i_nbr)%data, recv_buf(i_nbr)%pos, x(i))
                  end do
               else
                  call pmc_mpi_unpack_real(recv_buf(i_nbr)%data, recv_buf(i_nbr)%pos, w_src)
               end if
               ! pmc_mpi_unpack_aero_particle allocates component without checking
               if (allocated(particle%component)) deallocate(particle%component)
               call pmc_mpi_unpack_aero_particle(recv_buf(i_nbr)%data, recv_buf(i_nbr)%pos, particle)
               iCell = global_to_local(dest_global)
               if (iCell == 0) call mpas_log_write('PartMC: received particle for cell $i not owned here', &
                                                   intArgs=(/dest_global/), messageType=MPAS_LOG_CRIT)
               if (pmc_mode == MODE_LAGR) then
                  call deliver_lagr(particle, x, k, iCell)
               else
                  call deliver(particle, w_src, k, iCell)
               end if
            end do
         end do
      end if

      ! merge arrivals, rebalance, and compare with the Eulerian fields
      sums = 0.0_RKIND
      n_part_total = 0
      n_part_min = huge(1)
      n_part_max = 0
      allocate(err2_k(nVertLevels), norm2_k(nVertLevels), g_err2_k(nVertLevels), g_norm2_k(nVertLevels))
      err2_k = 0.0_RKIND
      norm2_k = 0.0_RKIND
      allocate(lev_sums(nVertLevels,3), g_lev_sums(nVertLevels,3))
      lev_sums = 0.0_RKIND
      n_t = pmc_n_part / n_class
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            if (pmc_mode == MODE_LAGR) then
               do i = 1, pos_in(k,iCell) % n
                  call aero_state_add_particle(aero_states(k,iCell), aero_incoming(k,iCell)%apa%particle(i), aero_data)
                  call pos_append(pos(k,iCell), pos_in(k,iCell) % x(:,i))
               end do
               pos_in(k,iCell) % n = 0
               if (lagr_control) call lagr_rebalance(k, iCell, nCells, nEdgesOnCell, cellsOnCell)
            else
               call aero_state_add_particles(aero_states(k,iCell), aero_incoming(k,iCell), aero_data)
            end if
            call aero_state_zero(aero_incoming(k,iCell))
            call aero_state_set_weight(aero_incoming(k,iCell), aero_data, AERO_STATE_WEIGHT_FLAT_SPECIFIED)

            counts = 0
            if (pmc_mode == MODE_FLUX) then
               call count_classes(aero_states(k,iCell), counts)
            else
               counts(1) = aero_state_n_part(aero_states(k,iCell))
            end if
            if (pmc_mode == MODE_FLUX .and. n_class > 1) then
               ! drop a class whose number here is negligible next to the largest class, so that
               ! particles are not spent on it
               m = 0.0_RKIND
               do c = 1, n_class
                  m = max(m, real(counts(c), kind=RKIND) * aero_states(k,iCell)%awa%weight(1,c)%magnitude)
               end do
               do c = 1, n_class
                  if (counts(c) == 0) cycle
                  w_src = real(counts(c), kind=RKIND) * aero_states(k,iCell)%awa%weight(1,c)%magnitude
                  if (w_src < 1.0e-4_RKIND * m) then
                     pruned = pruned + w_src * cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))
                     do i_part = aero_state_n_part(aero_states(k,iCell)), 1, -1
                        if (aero_states(k,iCell)%apa%particle(i_part)%weight_class == c) then
                           call aero_state_remove_particle_no_info(aero_states(k,iCell), i_part)
                        end if
                     end do
                     counts(c) = 0
                  end if
               end do
            end if
            if (pmc_mode == MODE_FLUX) then
               do c = 1, n_class
                  do while (counts(c) > 2*n_t)
                     call aero_state_halve(aero_states(k,iCell), 1, c)
                     call count_classes(aero_states(k,iCell), counts)
                  end do
                  do while (counts(c) < n_t/2 .and. counts(c) > 0)
                     call aero_state_double(aero_states(k,iCell), aero_data, 1, c)
                     call count_classes(aero_states(k,iCell), counts)
                  end do
               end do
            end if
            ! PartMC records every particle removed by halving; drop the records
            call aero_info_array_zero(aero_states(k,iCell)%aero_info_array)
            cnt(:,k,iCell) = counts(1:n_class)
            cnt_in(:,k,iCell) = 0

            if (.not. do_diag) cycle
            i = aero_state_n_part(aero_states(k,iCell))
            n_part_total = n_part_total + i
            n_part_min = min(n_part_min, i)
            n_part_max = max(n_part_max, i)
            vol = cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))
            ! particle contents in the same units as the captured contents (number conc / zz),
            ! and their variances if the counts were Poisson
            c_part = 0.0_RKIND
            v_part = 0.0_RKIND
            c1_part = 0.0_RKIND
            v1_part = 0.0_RKIND
            c2_part = 0.0_RKIND
            v2_part = 0.0_RKIND
            do i_part = 1, i
               c = aero_states(k,iCell)%apa%particle(i_part)%weight_class
               m = aero_states(k,iCell)%awa%weight(aero_states(k,iCell)%apa%particle(i_part)%weight_group,c)%magnitude &
                   / zz(k,iCell)
               c_part = c_part + m
               v_part = v_part + m*m
               if (aero_states(k,iCell)%apa%particle(i_part)%component(1)%source_id == 2) then
                  c2_part = c2_part + m
                  v2_part = v2_part + m*m
               else
                  c1_part = c1_part + m
                  v1_part = v1_part + m*m
               end if
            end do
            c_ref = rho_zz_new(k,iCell) * scalars_new(index_ref,k,iCell)
            c2_ref = rho_zz_new(k,iCell) * scalars_new(index_src2,k,iCell)
            ! a uniform mixing ratio stays uniform, so the background reference is exact
            c1_ref = rho_zz_new(k,iCell) * q_background
            c_new = 0.0_RKIND
            do c = 1, n_class
               c_new = c_new + fs(c) % cn(k,iCell)
            end do
            sums(1) = sums(1) + (c_part - c_ref)**2 * vol
            sums(2) = sums(2) + c_ref**2 * vol
            sums(3) = sums(3) + v_part * vol
            sums(4) = sums(4) + (c2_part - c2_ref)**2 * vol
            sums(5) = sums(5) + c2_ref**2 * vol
            sums(6) = sums(6) + v2_part * vol
            sums(7) = sums(7) + (c1_part - c1_ref)**2 * vol
            sums(8) = sums(8) + c1_ref**2 * vol
            sums(9) = sums(9) + v1_part * vol
            sums(10) = sums(10) + c_part * zz(k,iCell) * vol
            sums(11) = sums(11) + c_ref * zz(k,iCell) * vol
            sums(12) = sums(12) + (c_part - c_new)**2 * vol
            sums(13) = sums(13) + c2_part * zz(k,iCell) * vol
            sums(14) = sums(14) + c2_ref * zz(k,iCell) * vol
            ! integral of rho*q^2 for the source-2 mixing ratio q, which exact advection conserves;
            ! the particle value has its sampling variance removed
            sums(15) = sums(15) + c2_ref**2 / rho_zz_new(k,iCell) * zz(k,iCell) * vol
            sums(16) = sums(16) + (c2_part**2 - v2_part) / rho_zz_new(k,iCell) * zz(k,iCell) * vol
            err2_k(k) = err2_k(k) + (c_part - c_ref)**2 * vol
            norm2_k(k) = norm2_k(k) + c_ref**2 * vol
            ! level totals of the particles and the reference, and the sampling variance of the first
            lev_sums(k,1) = lev_sums(k,1) + c_part * zz(k,iCell) * vol
            lev_sums(k,2) = lev_sums(k,2) + c_ref * zz(k,iCell) * vol
            lev_sums(k,3) = lev_sums(k,3) + v_part * (zz(k,iCell) * vol)**2
         end do
      end do

      if (do_diag .and. pmc_mode == MODE_LAGR .and. lagr_wind /= WIND_MODEL) then
         ! distance of every particle from where it should be: the exact position for
         ! solid-body rotation, the starting position for the reversing flow
         terr = 0.0_RKIND
         terr_max = 0.0_RKIND
         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               do i = 1, pos(k,iCell) % n
                  if (lagr_wind == WIND_SOLID) then
                     xe = solid_rotate(pos(k,iCell) % x(5:7,i), &
                                       2.0_RKIND*pii/test_period * real(itimestep, kind=RKIND)*config_dt)
                  else
                     xe = pos(k,iCell) % x(5:7,i)
                  end if
                  d = 2.0_RKIND * radius * asin(min(1.0_RKIND, 0.5_RKIND*norm2(pos(k,iCell) % x(1:3,i) - xe)))
                  terr(1) = terr(1) + d
                  terr(2) = terr(2) + d*d
                  terr(3) = terr(3) + 1.0_RKIND
                  terr_max = max(terr_max, d)
               end do
            end do
         end do
         call mpas_dmpar_sum_real_array(domain % dminfo, 3, terr, g_terr)
         call mpas_dmpar_max_real(domain % dminfo, terr_max, g_terr_max)
         call mpas_log_write('PartMC trajectory test step $i: position error (km) mean = $r, rms = $r, max = $r', &
                             intArgs=(/itimestep/), realArgs=(/1.0e-3_RKIND*g_terr(1)/g_terr(3), &
                             1.0e-3_RKIND*sqrt(g_terr(2)/g_terr(3)), 1.0e-3_RKIND*g_terr_max/))
      end if

      if (do_diag) then
         call mpas_dmpar_sum_real_array(domain % dminfo, 16, sums, g_sums)
         if (m2_ref0 <= 0.0_RKIND) then
            m2_ref0 = g_sums(15)
            m2_part0 = g_sums(16)
         end if
         call mpas_dmpar_sum_int(domain % dminfo, n_part_total, g_n_part_total)
         call mpas_dmpar_sum_int(domain % dminfo, n_moved, g_n_moved)
         call mpas_dmpar_sum_int(domain % dminfo, nCellsSolve*nVertLevels, g_n_cells)
         call mpas_dmpar_min_int(domain % dminfo, n_part_min, g_n_part_min)
         call mpas_dmpar_max_int(domain % dminfo, n_part_max, g_n_part_max)
         call mpas_log_write('PartMC particles step $i: L2 rel err = $r, Poisson level = $r, total rel diff = $r, ' &
                             // 'moved frac = $r', intArgs=(/itimestep/), &
                             realArgs=(/sqrt(g_sums(1)/g_sums(2)), sqrt(g_sums(3)/g_sums(2)), &
                                        (g_sums(10) - g_sums(11))/g_sums(11), &
                                        real(g_n_moved, kind=RKIND)/real(max(1,g_n_part_total), kind=RKIND)/))
         call mpas_log_write('PartMC composition step $i: source-2 number L2 rel err = $r, Poisson level = $r, ' &
                             // 'total rel diff = $r', intArgs=(/itimestep/), &
                             realArgs=(/sqrt(g_sums(4)/g_sums(5)), sqrt(g_sums(6)/g_sums(5)), &
                                        (g_sums(13) - g_sums(14))/g_sums(14)/))
         call mpas_log_write('PartMC composition step $i: source-2 second moment relative to step 1: ' &
                             // 'reference = $r, particles = $r', intArgs=(/itimestep/), &
                             realArgs=(/g_sums(15)/m2_ref0, g_sums(16)/m2_part0/))
         if (scale_height <= 0.0_RKIND) then
            call mpas_log_write('PartMC uniform step $i: background number L2 rel err = $r, Poisson level = $r', &
                                intArgs=(/itimestep/), &
                                realArgs=(/sqrt(g_sums(7)/g_sums(8)), sqrt(g_sums(9)/g_sums(8))/))
         end if
         call mpas_log_write('PartMC particles step $i: n_part per cell mean/min/max = $i $i $i', &
                             intArgs=(/itimestep, g_n_part_total/g_n_cells, g_n_part_min, g_n_part_max/))
         call mpas_dmpar_sum_real_array(domain % dminfo, nVertLevels, err2_k, g_err2_k)
         call mpas_dmpar_sum_real_array(domain % dminfo, nVertLevels, norm2_k, g_norm2_k)
         lev_mean = 0.0_RKIND
         lev_worst = 0.0_RKIND
         k_worst = 0
         do k = 1, nVertLevels
            vol = sqrt(g_err2_k(k)/g_norm2_k(k))
            lev_mean = lev_mean + vol/real(nVertLevels, kind=RKIND)
            if (vol > lev_worst) then
               lev_worst = vol
               k_worst = k
            end if
         end do
         call mpas_log_write('PartMC particles step $i: per-level L2 rel err mean = $r, worst = $r at level $i', &
                             intArgs=(/itimestep, k_worst/), realArgs=(/lev_mean, lev_worst/))
         ! bias of each level's total against the reference, over levels holding at least 0.1% of
         ! the total, with the sampling standard deviation of the worst one
         call mpas_dmpar_sum_real_array(domain % dminfo, 3*nVertLevels, lev_sums, g_lev_sums)
         lev_all = sum(g_lev_sums(:,2))
         bias_worst = 0.0_RKIND
         bias_sigma = 0.0_RKIND
         k_bias = 0
         do k = 1, nVertLevels
            if (g_lev_sums(k,2) < 1.0e-3_RKIND*lev_all) cycle
            bias = g_lev_sums(k,1)/g_lev_sums(k,2) - 1.0_RKIND
            if (abs(bias) > abs(bias_worst)) then
               bias_worst = bias
               bias_sigma = sqrt(g_lev_sums(k,3))/g_lev_sums(k,2)
               k_bias = k
            end if
            if (mod(itimestep, 60) == 0) then
               call mpas_log_write('PartMC level step $i: level $i total rel diff = $r, share of total = $r', &
                                   intArgs=(/itimestep, k/), realArgs=(/bias, g_lev_sums(k,2)/lev_all/))
            end if
         end do
         call mpas_log_write('PartMC particles step $i: worst level total rel diff = $r at level $i, ' &
                             // 'sampling std = $r', intArgs=(/itimestep, k_bias/), realArgs=(/bias_worst, bias_sigma/))
         allocate(mass_total(n_spec))
         call species_mass_totals(domain, mass_total)
         mass_dev = 0.0_RKIND
         do i_spec = 1, n_spec
            if (mass_total_init(i_spec) > 0.0_RKIND) then
               mass_dev = max(mass_dev, abs(mass_total(i_spec)/mass_total_init(i_spec) - 1.0_RKIND))
            end if
         end do
         call mpas_log_write('PartMC composition step $i: largest species global mass rel change = $r', &
                             intArgs=(/itimestep/), realArgs=(/mass_dev/))
         if (pmc_mode == MODE_LAGR) then
            call mpas_dmpar_sum_int(domain % dminfo, n_clamp, g_n_clamp)
            call mpas_dmpar_sum_int(domain % dminfo, n_thinned, g_n_thinned)
            call mpas_dmpar_sum_int(domain % dminfo, n_split, g_n_split)
            call mpas_log_write('PartMC lagrangian step $i: particles held at the top or bottom = $i, ' &
                                // 'removed by thinning so far = $i, added by splitting so far = $i', &
                                intArgs=(/itimestep, g_n_clamp, g_n_thinned, g_n_split/))
         else
            call mpas_dmpar_max_real(domain % dminfo, p_sum_max, g_p_sum_max)
            call mpas_log_write('PartMC particles step $i: one-step L2 rel noise = $r, max outflow prob sum = $r', &
                                intArgs=(/itimestep/), realArgs=(/sqrt(g_sums(12)/g_sums(2)), g_p_sum_max/))
            call mpas_dmpar_max_real(domain % dminfo, ratio_max, g_ratio_max)
            call mpas_dmpar_min_real(domain % dminfo, ratio_min, g_ratio_min)
            call mpas_dmpar_max_int(domain % dminfo, n_copy_max, g_n_copy_max)
            call mpas_log_write('PartMC particles step $i: weight ratio src/dst min = $r, max = $r, ' &
                                // 'max copies of one particle = $i', &
                                intArgs=(/itimestep, g_n_copy_max/), realArgs=(/g_ratio_min, g_ratio_max/))
            call mpas_dmpar_sum_int(domain % dminfo, n_thin, g_n_thin)
            call mpas_dmpar_sum_real(domain % dminfo, pruned, g_pruned)
            call mpas_log_write('PartMC particles step $i: class coarsenings on arrival = $i, ' &
                                // 'number frac dropped as negligible = $r', &
                                intArgs=(/itimestep, g_n_thin/), realArgs=(/g_pruned/g_sums(11)/))
            call mpas_dmpar_sum_int(domain % dminfo, n_over, g_n_over)
            call mpas_dmpar_sum_real(domain % dminfo, lost, g_lost)
            call mpas_log_write('PartMC particles step $i: cells with outflow > content = $i of $i, ' &
                                // 'unmoved number frac = $r', &
                                intArgs=(/itimestep, g_n_over, g_n_cells/), realArgs=(/g_lost/g_sums(11)/))
         end if
      end if

   contains

      !> Number of particles in each weight class.
      subroutine count_classes(aero_state, n)
         type(aero_state_t), intent(in) :: aero_state
         integer, intent(out) :: n(MAX_CLASS)
         integer :: j, cc
         n = 0
         do j = 1, aero_state_n_part(aero_state)
            cc = aero_state%apa%particle(j)%weight_class
            n(cc) = n(cc) + 1
         end do
      end subroutine count_classes


      !> Pack a particle for the task that owns local halo cell cell_to.
      subroutine send_particle(part, k_to, cell_to, n_real, reals)
         type(aero_particle_t), intent(in) :: part
         integer, intent(in) :: k_to, cell_to, n_real
         real (kind=RKIND), intent(in) :: reals(n_real)
         integer :: nb, j, n_b

         nb = nbr_of_cell(cell_to)
         n_b = 2*pmc_mpi_pack_size_integer(1) + n_real*pmc_mpi_pack_size_real(reals(1)) &
               + pmc_mpi_pack_size_aero_particle(part)
         call pmc_buf_reserve(send_buf(nb), n_b)
         call pmc_mpi_pack_integer(send_buf(nb)%data, send_buf(nb)%pos, indexToCellID(cell_to))
         call pmc_mpi_pack_integer(send_buf(nb)%data, send_buf(nb)%pos, k_to)
         do j = 1, n_real
            call pmc_mpi_pack_real(send_buf(nb)%data, send_buf(nb)%pos, reals(j))
         end do
         call pmc_mpi_pack_aero_particle(send_buf(nb)%data, send_buf(nb)%pos, part)
      end subroutine send_particle


      !> Flux mode: choose which particles leave each cell and through which face.
      subroutine flux_move()
         integer :: dest_cell(MAX_FACES), dest_k(MAX_FACES)
         real (kind=RKIND) :: prob(MAX_FACES), p_rem, p_sum, w_class(MAX_CLASS), cell_vol
         integer, allocatable :: dest_face(:), idx(:)
         integer :: n_face, i_face, n_tot, n_c, n_rem, n_move, i_move, pick, j, cc, iEdge

         allocate(dest_face(4*pmc_n_part + 16), idx(4*pmc_n_part + 16))
         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               n_tot = aero_state_n_part(aero_states(k,iCell))
               if (n_tot == 0) cycle
               if (n_tot > size(dest_face)) then
                  deallocate(dest_face, idx)
                  allocate(dest_face(2*n_tot), idx(2*n_tot))
               end if
               dest_face(1:n_tot) = 0
               cell_vol = cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))

               ! faces: horizontal neighbours, then up, then down
               n_face = nEdgesOnCell(iCell)
               do i = 1, n_face
                  dest_cell(i) = cellsOnCell(i,iCell)
                  dest_k(i) = k
               end do
               if (k < nVertLevels) then
                  n_face = n_face + 1
                  dest_cell(n_face) = iCell
                  dest_k(n_face) = k + 1
               end if
               if (k > 1) then
                  n_face = n_face + 1
                  dest_cell(n_face) = iCell
                  dest_k(n_face) = k - 1
               end if

               do cc = 1, n_class
                  if (fs(cc) % co(k,iCell) <= 0.0_RKIND .or. cnt(cc,k,iCell) == 0) cycle

                  ! outgoing probability for each face
                  do i = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(i,iCell)
                     prob(i) = max(0.0_RKIND, edgesOnCell_sign(i,iCell)*fs(cc) % fh(k,iEdge)*invAreaCell(iCell)) &
                               / fs(cc) % co(k,iCell)
                  end do
                  i_face = nEdgesOnCell(iCell)
                  if (k < nVertLevels) then
                     i_face = i_face + 1
                     prob(i_face) = max(0.0_RKIND, rdzw(k)*fs(cc) % fv(k+1,iCell)) / fs(cc) % co(k,iCell)
                  end if
                  if (k > 1) then
                     i_face = i_face + 1
                     prob(i_face) = max(0.0_RKIND, -rdzw(k)*fs(cc) % fv(k,iCell)) / fs(cc) % co(k,iCell)
                  end if

                  ! The limiter bounds the net update, not the gross outflow, so in a noisy field a cell
                  ! with much less content than its neighbours can have outflow exceeding its content.
                  ! Scale the probabilities to sum to one (as WRF-PartMC does) and record the shortfall.
                  p_sum = sum(prob(1:n_face))
                  p_sum_max = max(p_sum_max, p_sum)
                  if (p_sum > 1.0_RKIND) then
                     n_over = n_over + 1
                     lost = lost + (p_sum - 1.0_RKIND) * fs(cc) % co(k,iCell) * zz(k,iCell) * cell_vol
                     prob(1:n_face) = prob(1:n_face) / p_sum
                  end if

                  ! particles of this class
                  n_c = 0
                  do j = 1, n_tot
                     if (aero_states(k,iCell)%apa%particle(j)%weight_class == cc) then
                        n_c = n_c + 1
                        idx(n_c) = j
                     end if
                  end do

                  ! multinomial sample of destinations by sequential binomials
                  n_rem = n_c
                  p_rem = 1.0_RKIND
                  do i_face = 1, n_face
                     if (n_rem <= 0) exit
                     if (prob(i_face) <= 0.0_RKIND) cycle
                     n_move = rand_binomial(n_rem, min(1.0_RKIND, prob(i_face)/p_rem))
                     p_rem = max(p_rem - prob(i_face), 1.0e-300_RKIND)
                     do i_move = 1, n_move
                        pick = pmc_rand_int(n_rem)
                        dest_face(idx(pick)) = i_face
                        idx(pick) = idx(n_rem)
                        n_rem = n_rem - 1
                     end do
                     n_moved = n_moved + n_move
                  end do
                  w_class(cc) = aero_states(k,iCell)%awa%weight(1,cc)%magnitude * cell_vol
               end do

               ! remove the chosen particles, highest index first so that PartMC's
               ! move-last-into-gap removal never disturbs one still to be handled
               do j = n_tot, 1, -1
                  i_face = dest_face(j)
                  if (i_face == 0) cycle
                  particle = aero_states(k,iCell)%apa%particle(j)
                  cc = particle%weight_class
                  call aero_state_remove_particle_no_info(aero_states(k,iCell), j)
                  cnt(cc,k,iCell) = cnt(cc,k,iCell) - 1
                  if (nbr_of_cell(dest_cell(i_face)) == 0) then
                     call deliver(particle, w_class(cc), dest_k(i_face), dest_cell(i_face))
                  else
                     call send_particle(particle, dest_k(i_face), dest_cell(i_face), 1, (/ w_class(cc) /))
                  end if
               end do
            end do
         end do
      end subroutine flux_move


      !> Put a particle that represented w_from real particles into owned cell (k_to, cell_to),
      !> duplicating or thinning it to match the weight of its class there.
      subroutine deliver(part, w_from, k_to, cell_to)
         type(aero_particle_t), intent(inout) :: part
         real (kind=RKIND), intent(in) :: w_from
         integer, intent(in) :: k_to, cell_to
         real (kind=RKIND) :: r, vol_to, r_max
         integer :: n, j, cc

         r_max = real(max(4, pmc_n_part / n_class), kind=RKIND)

         cc = part%weight_class
         vol_to = cell_volume(areaCell(cell_to), rdzw(k_to), zz(k_to,cell_to))
         do
            ! a class with no particles in the cell takes the weight of its first arrival
            if (cnt(cc,k_to,cell_to) + cnt_in(cc,k_to,cell_to) == 0) then
               aero_states(k_to,cell_to)%awa%weight(1,cc)%magnitude = w_from / vol_to
            end if
            r = w_from / (aero_states(k_to,cell_to)%awa%weight(1,cc)%magnitude * vol_to)
            if (r <= r_max) exit
            ! The arrival represents far more than a particle of this class here does (the class
            ! is nearly gone from the cell). Coarsen the class here instead of making many copies.
            call thin_class(aero_states(k_to,cell_to), cc, cnt(cc,k_to,cell_to))
            call thin_class(aero_incoming(k_to,cell_to), cc, cnt_in(cc,k_to,cell_to))
            aero_states(k_to,cell_to)%awa%weight(1,cc)%magnitude = &
                 2.0d0 * aero_states(k_to,cell_to)%awa%weight(1,cc)%magnitude
            n_thin = n_thin + 1
         end do
         n = int(r)
         if (pmc_random() < r - real(n, kind=RKIND)) n = n + 1
         ratio_max = max(ratio_max, r)
         ratio_min = min(ratio_min, r)
         n_copy_max = max(n_copy_max, n)
         do j = 1, n
            if (j > 1) call aero_particle_new_id(part)
            call aero_state_add_particle(aero_incoming(k_to,cell_to), part, aero_data)
         end do
         cnt_in(cc,k_to,cell_to) = cnt_in(cc,k_to,cell_to) + n
      end subroutine deliver


      !> Remove each particle of class cc with probability one half; n_left is how many remain.
      subroutine thin_class(aero_state, cc, n_left)
         type(aero_state_t), intent(inout) :: aero_state
         integer, intent(in) :: cc
         integer, intent(out) :: n_left
         integer :: j

         n_left = 0
         do j = aero_state_n_part(aero_state), 1, -1
            if (aero_state%apa%particle(j)%weight_class /= cc) cycle
            if (pmc_random() < 0.5d0) then
               call aero_state_remove_particle_no_info(aero_state, j)
            else
               n_left = n_left + 1
            end if
         end do
      end subroutine thin_class


      !> Lagrangian mode: put a particle with position xp into owned cell (k_to, cell_to).
      subroutine deliver_lagr(part, xp, k_to, cell_to)
         type(aero_particle_t), intent(in) :: part
         real (kind=RKIND), intent(in) :: xp(NPOS)
         integer, intent(in) :: k_to, cell_to

         call aero_state_add_particle(aero_incoming(k_to,cell_to), part, aero_data)
         call pos_append(pos_in(k_to,cell_to), xp)
      end subroutine deliver_lagr


      !> Step-averaged transport velocity: Cartesian horizontal velocity at cell centres
      !> and d(zeta)/dt at cell faces, from the mass fluxes used for scalar transport.
      subroutine lagr_velocities()
         type (field2DReal), pointer :: f_ux, f_uy, f_uz, f_om, f_div, f_ru
         real (kind=RKIND), allocatable :: u_e(:,:), u_zonal(:,:), u_merid(:,:)
         real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge, dvEdge, areaCell_l
         integer, dimension(:,:), pointer :: edgesOnCell_l
         real (kind=RKIND) :: rho_e, rho_f, u_test, v_test, xe(3), flux, sgn
         integer :: iEdge, c1, c2, kk, ic, je

         call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
         call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
         call mpas_pool_get_field(diag, 'pmc_ux', f_ux)
         call mpas_pool_get_field(diag, 'pmc_uy', f_uy)
         call mpas_pool_get_field(diag, 'pmc_uz', f_uz)
         call mpas_pool_get_field(diag, 'pmc_omega', f_om)
         ux => f_ux % array
         uy => f_uy % array
         uz => f_uz % array
         om => f_om % array

         if (lagr_flux_vel > 0 .and. lagr_wind == WIND_MODEL) then
            ! Mass-flux vector and half the mass divergence of every layer of every cell, from the
            ! cell's own edge fluxes: m(x) = M + (D/2)(x - x_c) has divergence D, and M is the
            ! area-weighted reconstruction sum over edges of (x_e - x_c) * outward flux.
            ! The vertical mass flux is kept at the faces.
            call mpas_pool_get_field(diag, 'pmc_div', f_div)
            hdiv => f_div % array
            call mpas_pool_get_array(mesh, 'xEdge', xEdge)
            call mpas_pool_get_array(mesh, 'yEdge', yEdge)
            call mpas_pool_get_array(mesh, 'zEdge', zEdge)
            call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
            call mpas_pool_get_array(mesh, 'areaCell', areaCell_l)
            call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell_l)
            do ic = 1, nCellsSolve
               ux(:,ic) = 0.0_RKIND
               uy(:,ic) = 0.0_RKIND
               uz(:,ic) = 0.0_RKIND
               hdiv(:,ic) = 0.0_RKIND
               do je = 1, nEdgesOnCell(ic)
                  iEdge = edgesOnCell_l(je,ic)
                  sgn = 1.0_RKIND
                  if (cellsOnEdge(2,iEdge) == ic) sgn = -1.0_RKIND
                  xe = (/ xEdge(iEdge), yEdge(iEdge), zEdge(iEdge) /)
                  xe = radius * (xe / norm2(xe) - cxyz(:,ic))
                  do kk = 1, nVertLevels
                     flux = sgn * ruAvg(kk,iEdge) * dvEdge(iEdge) / areaCell_l(ic)
                     ux(kk,ic) = ux(kk,ic) + xe(1)*flux
                     uy(kk,ic) = uy(kk,ic) + xe(2)*flux
                     uz(kk,ic) = uz(kk,ic) + xe(3)*flux
                     hdiv(kk,ic) = hdiv(kk,ic) + 0.5_RKIND*flux
                  end do
               end do
               om(1,ic) = 0.0_RKIND
               om(nVertLevels+1,ic) = 0.0_RKIND
               om(2:nVertLevels,ic) = wwAvg(2:nVertLevels,ic)
            end do
            call mpas_dmpar_exch_halo_field(f_ux)
            call mpas_dmpar_exch_halo_field(f_uy)
            call mpas_dmpar_exch_halo_field(f_uz)
            call mpas_dmpar_exch_halo_field(f_om)
            call mpas_dmpar_exch_halo_field(f_div)
            if (lagr_flux_vel == 2) then
               call mpas_pool_get_field(diag, 'ruAvg', f_ru)
               call mpas_dmpar_exch_halo_field(f_ru)
            end if
            return
         end if

         allocate(u_e(nVertLevels, nEdges+1), u_zonal(nVertLevels, nCells+1), u_merid(nVertLevels, nCells+1))
         u_e = 0.0_RKIND
         if (lagr_wind /= WIND_MODEL) then
            ! analytic horizontal wind at the middle of the step, as edge-normal components
            call mpas_pool_get_array(mesh, 'angleEdge', angleEdge)
            call mpas_pool_get_array(mesh, 'latEdge', latEdge)
            call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
            do iEdge = 1, nEdges
               call test_wind(latEdge(iEdge), lonEdge(iEdge), (real(itimestep, kind=RKIND) - 0.5_RKIND)*config_dt, &
                              u_test, v_test)
               u_e(:,iEdge) = u_test*cos(angleEdge(iEdge)) + v_test*sin(angleEdge(iEdge))
            end do
         else
            do iEdge = 1, nEdges
               c1 = cellsOnEdge(1,iEdge)
               c2 = cellsOnEdge(2,iEdge)
               if (c1 < 1 .or. c1 > nCells .or. c2 < 1 .or. c2 > nCells) cycle
               do kk = 1, nVertLevels
                  rho_e = 0.25_RKIND*(rho_zz_old(kk,c1) + rho_zz_new(kk,c1) + rho_zz_old(kk,c2) + rho_zz_new(kk,c2))
                  u_e(kk,iEdge) = ruAvg(kk,iEdge) / rho_e
               end do
            end do
         end if
         call mpas_reconstruct(mesh, u_e, ux, uy, uz, u_zonal, u_merid)
         do ic = 1, nCellsSolve
            om(1,ic) = 0.0_RKIND
            om(nVertLevels+1,ic) = 0.0_RKIND
            do kk = 2, nVertLevels
               rho_f = 0.5_RKIND*(fzm(kk)*(rho_zz_old(kk,ic) + rho_zz_new(kk,ic)) &
                                  + fzp(kk)*(rho_zz_old(kk-1,ic) + rho_zz_new(kk-1,ic)))
               om(kk,ic) = wwAvg(kk,ic) / rho_f
            end do
         end do
         if (lagr_wind /= WIND_MODEL) om = 0.0_RKIND
         call mpas_dmpar_exch_halo_field(f_ux)
         call mpas_dmpar_exch_halo_field(f_uy)
         call mpas_dmpar_exch_halo_field(f_uz)
         call mpas_dmpar_exch_halo_field(f_om)
      end subroutine lagr_velocities


      !> Velocity at unit vector p and height zp: linear in the triangle of cell centres
      !> that contains p, and linear in zeta.
      subroutine lagr_vel(p, zp, cell_hint, vel, omega)
         real (kind=RKIND), intent(in) :: p(3), zp
         integer, intent(in) :: cell_hint
         real (kind=RKIND), intent(out) :: vel(3), omega
         integer :: cell, j, iv, tri(3), best_tri(3), kk, k1, k2, kf, ct
         real (kind=RKIND) :: wt(3), best_wt(3), best_min, det, a(3), b(3), cc3(3), f, ff

         cell = find_cell(p, cell_hint, nCells, nEdgesOnCell, cellsOnCell)
         if (lagr_flux_vel > 0 .and. lagr_wind == WIND_MODEL) then
            kf = find_level(zp)
            ff = (zp - zf(kf)) / (zf(kf+1) - zf(kf))
            f = 2.0_RKIND / (rho_zz_old(kf,cell) + rho_zz_new(kf,cell))
            a = radius * (p - cxyz(:,cell))
            if (lagr_flux_vel == 2 .and. rt_n(cell) > 0) then
               call rt_mass_flux(cell, kf, a, vel)
               vel = f * vel
            else
               vel(1) = f * (ux(kf,cell) + hdiv(kf,cell)*a(1))
               vel(2) = f * (uy(kf,cell) + hdiv(kf,cell)*a(2))
               vel(3) = f * (uz(kf,cell) + hdiv(kf,cell)*a(3))
            end if
            vel = vel - sum(vel*p)*p
            omega = f * ((1.0_RKIND - ff)*om(kf,cell) + ff*om(kf+1,cell))
            return
         end if
         best_min = -huge(1.0_RKIND)
         best_tri = cell
         best_wt = (/ 1.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
         do j = 1, nEdgesOnCell(cell)
            iv = verticesOnCell(j,cell)
            tri = cellsOnVertex(:,iv)
            if (minval(tri) < 1 .or. maxval(tri) > nCells) cycle
            a = cxyz(:,tri(1))
            b = cxyz(:,tri(2))
            cc3 = cxyz(:,tri(3))
            det = triple(a, b, cc3)
            wt(1) = triple(p, b, cc3) / det
            wt(2) = triple(a, p, cc3) / det
            wt(3) = triple(a, b, p) / det
            if (minval(wt) > best_min) then
               best_min = minval(wt)
               best_tri = tri
               best_wt = wt
            end if
         end do
         best_wt = max(best_wt, 0.0_RKIND)
         best_wt = best_wt / sum(best_wt)

         ! level centres bracketing zp
         if (zp <= zc(1)) then
            k1 = 1
            k2 = 1
            f = 0.0_RKIND
         else if (zp >= zc(nVertLevels)) then
            k1 = nVertLevels
            k2 = nVertLevels
            f = 0.0_RKIND
         else
            k1 = 1
            do kk = 1, nVertLevels - 1
               if (zp < zc(kk+1)) then
                  k1 = kk
                  exit
               end if
            end do
            k2 = k1 + 1
            f = (zp - zc(k1)) / (zc(k2) - zc(k1))
         end if
         kf = find_level(zp)
         ff = (zp - zf(kf)) / (zf(kf+1) - zf(kf))
         if (lagr_u_const) then
            ! the layer's own velocity, whose divergence is what the layer's omega difference balances
            k1 = kf
            k2 = kf
            f = 0.0_RKIND
         end if

         vel = 0.0_RKIND
         omega = 0.0_RKIND
         do j = 1, 3
            ct = best_tri(j)
            vel(1) = vel(1) + best_wt(j) * ((1.0_RKIND - f)*ux(k1,ct) + f*ux(k2,ct))
            vel(2) = vel(2) + best_wt(j) * ((1.0_RKIND - f)*uy(k1,ct) + f*uy(k2,ct))
            vel(3) = vel(3) + best_wt(j) * ((1.0_RKIND - f)*uz(k1,ct) + f*uz(k2,ct))
            omega = omega + best_wt(j) * ((1.0_RKIND - ff)*om(kf,ct) + ff*om(kf+1,ct))
         end do
      end subroutine lagr_vel


      !> Horizontal mass flux per area at offset r (m) from the centre of a cell, in layer kk:
      !> the lowest-order Raviart-Thomas field on the triangles (centre, vertex, vertex), which
      !> has the model's mass flux through every edge and a uniform divergence. The fluxes
      !> through the spokes from the centre to the vertices follow from the uniform divergence,
      !> up to a rotation about the centre that is set to zero on average.
      subroutine rt_mass_flux(cell, kk, r, m)
         integer, intent(in) :: cell, kk
         real (kind=RKIND), intent(in) :: r(3)
         real (kind=RKIND), intent(out) :: m(3)
         integer :: n, i, ip, it
         real (kind=RKIND) :: fl(MAX_FACES), g(MAX_FACES), d, nrm(3), best, w1, w2, s1, s2

         n = rt_n(cell)
         nrm = cxyz(:,cell)
         do i = 1, n
            fl(i) = rt_len(i,cell) * ruAvg(kk, rt_edge(i,cell))
         end do
         d = sum(fl(1:n)) / sum(rt_area(1:n,cell))
         g(1) = 0.0_RKIND
         do i = 1, n - 1
            g(i+1) = g(i) + d*rt_area(i,cell) - fl(i)
         end do
         s1 = 0.0_RKIND
         s2 = 0.0_RKIND
         do i = 1, n
            w1 = 1.0_RKIND / norm2(rt_a(:,i,cell))
            s1 = s1 + g(i)*w1
            s2 = s2 + w1
         end do
         g(1:n) = g(1:n) - s1/s2
         ! the triangle that holds r
         it = 1
         best = -huge(1.0_RKIND)
         do i = 1, n
            w1 = triple(nrm, rt_a(:,i,cell), r)
            w2 = triple(nrm, r, rt_b(:,i,cell))
            if (min(w1, w2) > best) then
               best = min(w1, w2)
               it = i
            end if
         end do
         ip = it + 1
         if (ip > n) ip = 1
         m = (fl(it)*r + g(ip)*(r - rt_a(:,it,cell)) - g(it)*(r - rt_b(:,it,cell))) / (2.0_RKIND*rt_area(it,cell))
      end subroutine rt_mass_flux


      pure real (kind=RKIND) function triple(a, b, c3)
         real (kind=RKIND), intent(in) :: a(3), b(3), c3(3)
         triple = a(1)*(b(2)*c3(3) - b(3)*c3(2)) + a(2)*(b(3)*c3(1) - b(1)*c3(3)) + a(3)*(b(1)*c3(2) - b(2)*c3(1))
      end function triple


      !> Lagrangian mode: advance every particle with midpoint steps through the
      !> step-averaged velocity field.
      subroutine lagr_move()
         real (kind=RKIND) :: xo(NPOS), xm(4), xn(NPOS), vel(3), omega, z_lo, z_hi, dt
         integer :: j, cell_n, k_n, i_sub
         logical :: clamped

         dt = config_dt / real(lagr_substeps, kind=RKIND)
         z_lo = zf(1) + 1.0e-6_RKIND*(zf(2) - zf(1))
         z_hi = zf(nVertLevels+1) - 1.0e-6_RKIND*(zf(nVertLevels+1) - zf(nVertLevels))
         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               do j = pos(k,iCell) % n, 1, -1
                  xn = pos(k,iCell) % x(:,j)
                  cell_n = iCell
                  clamped = .false.
                  do i_sub = 1, lagr_substeps
                     xo = xn
                     call lagr_vel(xo(1:3), xo(4), cell_n, vel, omega)
                     xm(1:3) = xo(1:3) + 0.5_RKIND*dt*vel/radius
                     xm(1:3) = xm(1:3) / norm2(xm(1:3))
                     xm(4) = min(z_hi, max(z_lo, xo(4) + 0.5_RKIND*dt*omega))
                     call lagr_vel(xm(1:3), xm(4), cell_n, vel, omega)
                     xn(1:3) = xo(1:3) + dt*vel/radius
                     xn(1:3) = xn(1:3) / norm2(xn(1:3))
                     xn(4) = xo(4) + dt*omega
                     if (xn(4) < z_lo .or. xn(4) > z_hi) then
                        clamped = .true.
                        xn(4) = min(z_hi, max(z_lo, xn(4)))
                     end if
                     cell_n = find_cell(xn(1:3), cell_n, nCells, nEdgesOnCell, cellsOnCell)
                  end do
                  if (clamped) n_clamp = n_clamp + 1
                  k_n = find_level(xn(4))
                  if (cell_n == iCell .and. k_n == k) then
                     pos(k,iCell) % x(:,j) = xn
                     cycle
                  end if
                  particle = aero_states(k,iCell)%apa%particle(j)
                  call aero_state_remove_particle_no_info(aero_states(k,iCell), j)
                  call pos_remove(pos(k,iCell), j)
                  n_moved = n_moved + 1
                  if (cell_n <= nCellsSolve) then
                     call deliver_lagr(particle, xn, k_n, cell_n)
                  else
                     call send_particle(particle, k_n, cell_n, NPOS, xn)
                  end if
               end do
            end do
         end do
      end subroutine lagr_move
#endif

   end subroutine atm_partmc_particles_transport


#ifdef MPAS_PARTMC
   !> Lagrangian population control for cell (k,iCell). Too many particles: each particle of
   !> the lightest class present is either removed or moved up one class (twice the weight),
   !> with equal probability, which keeps the expected number and mass. Too few: each particle
   !> of the heaviest class present becomes two of the next class down (half the weight), which
   !> is exact; the new one is placed at a random point of the cell.
   subroutine lagr_rebalance(k, iCell, nCells, nEdgesOnCell, cellsOnCell)
      integer, intent(in) :: k, iCell, nCells
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: cellsOnCell
      integer :: n, j, g, g_sel
      real (kind=RKIND) :: xp(NPOS), s_w
      type(aero_particle_t) :: twin

      n = pos(k,iCell) % n
      do while (n > 2*pmc_n_part)
         g_sel = LAGR_GROUPS + 1
         do j = 1, n
            g_sel = min(g_sel, aero_states(k,iCell)%apa%particle(j)%weight_class)
         end do
         if (g_sel >= LAGR_GROUPS) exit
         do j = n, 1, -1
            if (aero_states(k,iCell)%apa%particle(j)%weight_class /= g_sel) cycle
            ! sulfate this particle stands for (mol): lost if removed, gained again if doubled
            s_w = particle_sulfate(aero_states(k,iCell)%apa%particle(j)) * lagr_w * 2.0_RKIND**(g_sel - LAGR_G0)
            if (pmc_random() < 0.5d0) then
               part_s_thin = part_s_thin - s_w
               call aero_state_remove_particle_no_info(aero_states(k,iCell), j)
               call pos_remove(pos(k,iCell), j)
               n_thinned = n_thinned + 1
            else
               part_s_thin = part_s_thin + s_w
               aero_states(k,iCell)%apa%particle(j)%weight_class = g_sel + 1
            end if
         end do
         n = pos(k,iCell) % n
      end do
      do while (n > 0 .and. n < pmc_n_part/2)
         g_sel = 0
         do j = 1, n
            g_sel = max(g_sel, aero_states(k,iCell)%apa%particle(j)%weight_class)
         end do
         if (g_sel <= 1) exit
         do j = 1, n
            if (aero_states(k,iCell)%apa%particle(j)%weight_class /= g_sel) cycle
            aero_states(k,iCell)%apa%particle(j)%weight_class = g_sel - 1
            twin = aero_states(k,iCell)%apa%particle(j)
            call aero_particle_new_id(twin)
            call aero_state_add_particle(aero_states(k,iCell), twin, aero_data)
            call random_in_cell(iCell, k, nCells, nEdgesOnCell, cellsOnCell, xp(1:4))
            xp(5:7) = pos(k,iCell) % x(5:7,j)
            call pos_append(pos(k,iCell), xp)
            n_split = n_split + 1
         end do
         n = pos(k,iCell) % n
      end do
   end subroutine lagr_rebalance


#endif


   !> Vertical sub-grid mixing with a prescribed diffusivity profile (a test stand-in for the
   !> host's boundary-layer and convective mixing).
   !>
   !> The implicit diffusion step A q_new = q_old for the mixing ratio in a column is applied to
   !> the Eulerian number scalars. For the particles the same step is written for layer contents,
   !> m_new = T m_old with T = D A^-1 D^-1 and D = air mass per area of each layer. T is a
   !> transition matrix (non-negative, columns sum to one): a particle in layer j moves to layer
   !> i with probability T(i,j), to a random height in that layer. A uniform mixing ratio is
   !> then preserved in expectation for any time step.
   subroutine atm_partmc_vertical_mixing(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC
      type (mpas_pool_type), pointer :: mesh, state, diag_physics
      type (field3DReal), pointer :: scalars_field
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, zz, zgrid, kzh
      real (kind=RKIND), dimension(:), pointer :: rdzw
      real (kind=RKIND), pointer :: config_dt
      integer, dimension(:,:), pointer :: cellsOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, pointer :: nCells, nCellsSolve, nVertLevels, index_pmc_num, index_pmc_num2, index_ref, index_src2
      integer :: iCell, k, kk, j, i_s, nz
      integer, allocatable :: idx(:)
      real (kind=RKIND), allocatable :: d_mass(:), e(:), lower(:), diag_a(:), upper(:), rhs(:), sol(:), tmat(:,:)
      real (kind=RKIND) :: s_h, kz, rho_f, dz_c, u, cum, xp(NPOS)
      type(aero_particle_t) :: particle

      if (pmc_kz_max <= 0.0_RKIND .and. .not. pmc_kz_phys) return

      nullify(kzh)
      if (pmc_kz_phys) then
         call mpas_pool_get_subpool(domain % blocklist % structs, 'diag_physics', diag_physics)
         call mpas_pool_get_array(diag_physics, 'kzh', kzh)
         if (.not. associated(kzh)) then
            call mpas_log_write('PartMC: PMC_KZ_PHYS needs a boundary-layer scheme that provides kzh', &
                                messageType=MPAS_LOG_CRIT)
         end if
      end if

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', config_dt)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_pmc_num', index_pmc_num)
      call mpas_pool_get_dimension(state, 'index_pmc_num2', index_pmc_num2)
      call mpas_pool_get_dimension(state, 'index_pmc_ref', index_ref)
      call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', index_src2)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mixed_scalars(state, idx)
      nz = nVertLevels
      allocate(d_mass(nz), e(nz+1), lower(nz), diag_a(nz), upper(nz), rhs(nz), sol(nz), tmat(nz,nz))

      do iCell = 1, nCellsSolve
         ! air mass per area of each layer, and the exchange coefficient at each interior face
         e = 0.0_RKIND
         do k = 1, nz
            d_mass(k) = rho_zz(k,iCell) / rdzw(k)
         end do
         do k = 2, nz
            s_h = zgrid(k,iCell) / pmc_kz_top
            kz = 0.0_RKIND
            if (s_h < 1.0_RKIND) kz = pmc_kz_max * 6.75_RKIND * s_h * (1.0_RKIND - s_h)**2
            ! the scheme's kzh(k) is at the face below layer k
            if (pmc_kz_phys) kz = kzh(k,iCell)
            rho_f = 0.5_RKIND*(rho_zz(k,iCell)*zz(k,iCell) + rho_zz(k-1,iCell)*zz(k-1,iCell))
            dz_c = 0.5_RKIND*(zgrid(k+1,iCell) - zgrid(k-1,iCell))
            e(k) = config_dt * rho_f * kz / dz_c
         end do
         do k = 1, nz
            lower(k) = -e(k) / d_mass(k)
            upper(k) = -e(k+1) / d_mass(k)
            diag_a(k) = 1.0_RKIND - lower(k) - upper(k)
         end do

         ! Eulerian scalars
         do i_s = 1, size(idx)
            rhs = scalars(idx(i_s),:,iCell)
            call tridiag_solve(nz, lower, diag_a, upper, rhs, sol)
            scalars(idx(i_s),:,iCell) = sol
         end do

         if (pmc_n_part <= 0 .or. pmc_mode /= MODE_LAGR) cycle

         ! transition matrix, one column per source layer
         do j = 1, nz
            rhs = 0.0_RKIND
            rhs(j) = 1.0_RKIND
            call tridiag_solve(nz, lower, diag_a, upper, rhs, sol)
            tmat(:,j) = max(0.0_RKIND, d_mass * sol / d_mass(j))
            tmat(:,j) = tmat(:,j) / sum(tmat(:,j))
         end do

         call column_jump(iCell, nz, tmat, n_mixed, nCells, nEdgesOnCell, cellsOnCell)
      end do

      call mpas_pool_get_field(state, 'scalars', scalars_field, 2)
      call mpas_dmpar_exch_halo_field(scalars_field)

      if (mod(itimestep, 12) == 0 .and. pmc_mode == MODE_LAGR .and. pmc_n_part > 0) then
         call mpas_dmpar_sum_int(domain % dminfo, n_mixed, j)
         call mpas_log_write('PartMC mixing step $i: particle layer changes by vertical mixing so far = $i', &
                             intArgs=(/itimestep, j/))
      end if
#endif

   end subroutine atm_partmc_vertical_mixing


   !> Indices of the scalars that vertical mixing and convection act on: the particle number
   !> scalars and, when chemistry is on, the gas tracers.
   subroutine mixed_scalars(state, idx)
      type (mpas_pool_type), intent(in) :: state
      integer, allocatable, intent(out) :: idx(:)
      integer, pointer :: i1, i2, i3, i4, gas_start, gas_end
      integer :: i, n_gas_t

      call mpas_pool_get_dimension(state, 'index_pmc_num', i1)
      call mpas_pool_get_dimension(state, 'index_pmc_num2', i2)
      call mpas_pool_get_dimension(state, 'index_pmc_ref', i3)
      call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', i4)
      n_gas_t = 0
      if (pmc_chem) then
         call mpas_pool_get_dimension(state, 'gas_start', gas_start)
         call mpas_pool_get_dimension(state, 'gas_end', gas_end)
         n_gas_t = gas_end - gas_start + 1
      end if
      allocate(idx(4 + n_gas_t))
      idx(1:4) = (/ i1, i2, i3, i4 /)
      do i = 1, n_gas_t
         idx(4+i) = gas_start + i - 1
      end do
   end subroutine mixed_scalars


   !> Move the particles of one column between layers: a particle in layer j goes to layer i
   !> with probability tmat(i,j), to a random height in that layer (Lagrangian mode only).
   !>
   !> With tmat_s, p_rem and temp (removal in the convective updraft): a particle in layer j
   !> that counts as activated is removed with probability p_rem(j) and otherwise goes to
   !> layer i with probability tmat_s(i,j); the others use tmat.
   subroutine column_jump(iCell, nz, tmat, n_jump, nCells, nEdgesOnCell, cellsOnCell, tmat_s, p_rem, temp)

      integer, intent(in) :: iCell, nz, nCells
      real (kind=RKIND), intent(in) :: tmat(nz,nz)
      integer, intent(inout) :: n_jump
      integer, intent(in) :: nEdgesOnCell(:), cellsOnCell(:,:)
      real (kind=RKIND), intent(in), optional :: tmat_s(nz,nz), p_rem(nz), temp(nz)

#ifdef MPAS_PARTMC
      integer :: k, kk, j, g
      real (kind=RKIND) :: u, cum, xp(NPOS), colp(nz), w
      type(aero_particle_t) :: particle
      type(env_state_t) :: env_state
      logical :: scav

      do k = 1, nz
         if (tmat(k,k) >= 1.0_RKIND) cycle
         if (present(tmat_s)) env_state%temp = temp(k)
         do j = pos(k,iCell) % n, 1, -1
            scav = present(tmat_s)
            if (scav .and. pmc_wet_scrit < 1.0e8_RKIND) then
               scav = (aero_particle_crit_rel_humid(aero_states(k,iCell)%apa%particle(j), aero_data, env_state) &
                       - 1.0d0 <= pmc_wet_scrit)
            end if
            if (scav) then
               g = aero_states(k,iCell)%apa%particle(j)%weight_class
               w = lagr_w * 2.0_RKIND**(g - LAGR_G0)
               conv_wet(3) = conv_wet(3) + p_rem(k) * w
               if (pmc_random() < p_rem(k)) then
                  conv_wet(2) = conv_wet(2) + w
                  part_s_removed = part_s_removed + particle_sulfate(aero_states(k,iCell)%apa%particle(j)) * w
                  call aero_state_remove_particle_no_info(aero_states(k,iCell), j)
                  call pos_remove(pos(k,iCell), j)
                  cycle
               end if
               colp = tmat_s(:,k)
            else
               colp = tmat(:,k)
            end if
            u = pmc_random()
            ! try the particle's own layer first, since that is the usual outcome
            if (u < colp(k)) cycle
            cum = colp(k)
            do kk = 1, nz
               if (kk == k) cycle
               cum = cum + colp(kk)
               if (u < cum) exit
            end do
            kk = min(kk, nz)
            if (kk == k) cycle
            particle = aero_states(k,iCell)%apa%particle(j)
            xp = pos(k,iCell) % x(:,j)
            xp(4) = zf(kk) + pmc_random()*(zf(kk+1) - zf(kk))
            call aero_state_remove_particle_no_info(aero_states(k,iCell), j)
            call pos_remove(pos(k,iCell), j)
            call aero_state_add_particle(aero_incoming(kk,iCell), particle, aero_data)
            call pos_append(pos_in(kk,iCell), xp)
            n_jump = n_jump + 1
         end do
      end do
      do k = 1, nz
         do j = 1, pos_in(k,iCell) % n
            call aero_state_add_particle(aero_states(k,iCell), aero_incoming(k,iCell)%apa%particle(j), aero_data)
            call pos_append(pos(k,iCell), pos_in(k,iCell) % x(:,j))
         end do
         if (pos_in(k,iCell) % n > 0) then
            pos_in(k,iCell) % n = 0
            call aero_state_zero(aero_incoming(k,iCell))
         end if
         if (lagr_control) call lagr_rebalance(k, iCell, nCells, nEdgesOnCell, cellsOnCell)
      end do
#endif

   end subroutine column_jump


   !> Convective transport of the tracers and particles with the cumulus scheme's mass fluxes.
   !>
   !> Each column gets a transfer matrix for the step from a bulk-plume picture: every layer
   !> gives air to the updraft and the downdraft at their entrainment rates; that air leaves
   !> the plume in the layers above (below) in proportion to the plume's detrainment there;
   !> and the environment moves between neighbouring layers to balance the plumes' mass flux
   !> (compensating subsidence). Entrainment is the scheme's detrainment plus the change of
   !> the mass flux across the layer, so every layer's air mass is unchanged and a uniform
   !> mixing ratio stays uniform. The matrix is applied to the Eulerian scalars and, in
   !> Lagrangian mode, sampled by the particles.
   subroutine atm_partmc_convection(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC
      type (mpas_pool_type), pointer :: mesh, state, diag_physics
      type (field3DReal), pointer :: scalars_field
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, cu_mfu, cu_mfd, cu_ude, cu_dde
      real (kind=RKIND), dimension(:,:), pointer :: theta_m, exner, pressure_p, pressure_base
      real (kind=RKIND), dimension(:), pointer :: rdzw, areaCell, raincv
      real (kind=RKIND), pointer :: config_dt
      integer, dimension(:,:), pointer :: cellsOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, pointer :: nCells, nCellsSolve, nVertLevels, index_pmc_num, index_pmc_num2, index_ref, index_src2
      integer, pointer :: index_qv
      integer :: iCell, k, j, i_s, nz, n_sub, i_sub, n_col, g_n_col, g_n_conv
      type (mpas_pool_type), pointer :: diag
      real (kind=RKIND), allocatable :: tsub_s(:,:), tmat_s(:,:), tk(:), qs(:), qv(:), p_rem(:), col_s(:)
      integer, allocatable :: gas_of(:)
      real (kind=RKIND) :: w_cond, pe, es, pres, lost, e_gas, sums(6), g_sums(6)
      logical :: scav
      type (mpas_pool_type), pointer :: tend_physics
      real (kind=RKIND), dimension(:,:), pointer :: rthcuten, rqvcuten, zgrid, rucuten, rvcuten, u_zon, u_mer, cu_trc
      real (kind=RKIND), allocatable :: hm(:), dh_m(:), dh_s(:)
      real (kind=RKIND) :: csum(36), g_csum(36), z_mid, wgt, slab_m(3), slab_s(3)
      integer :: i_sl
      integer :: i_f, i_sh, kk
      character(len=8), parameter :: chk_name(4) = (/ 'energy  ', 'u wind  ', 'v wind  ', 'tracer  ' /)
      integer, allocatable :: idx(:)
      real (kind=RKIND), allocatable :: m(:), mu_top(:), mu_bot(:), dn_top(:), dn_bot(:), e_u(:), d_u(:), &
                                        e_d(:), d_d(:), f_u(:), f_d(:), s_top(:), out(:), tsub(:,:), tmat(:,:), &
                                        tnew(:,:), col(:)
      real (kind=RKIND) :: dt, pass, o_max, g_o_max
      real (kind=RKIND), save :: o_max_seen = 0.0_RKIND

      if (.not. pmc_conv) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'diag_physics', diag_physics)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', config_dt)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_pmc_num', index_pmc_num)
      call mpas_pool_get_dimension(state, 'index_pmc_num2', index_pmc_num2)
      call mpas_pool_get_dimension(state, 'index_pmc_ref', index_ref)
      call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', index_src2)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      nullify(cu_mfu)
      call mpas_pool_get_array(diag_physics, 'cu_mfu', cu_mfu)
      call mpas_pool_get_array(diag_physics, 'cu_mfd', cu_mfd)
      call mpas_pool_get_array(diag_physics, 'cu_ude', cu_ude)
      call mpas_pool_get_array(diag_physics, 'cu_dde', cu_dde)
      if (.not. associated(cu_mfu)) then
         call mpas_log_write('PartMC: PMC_CONV needs the cu_ntiedtke convection scheme', messageType=MPAS_LOG_CRIT)
      end if
      call mixed_scalars(state, idx)
      nz = nVertLevels
      ! for removal in the updraft
      call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_base)
      nullify(raincv)
      call mpas_pool_get_array(diag_physics, 'raincv', raincv)
      if (pmc_conv_check) then
         call mpas_pool_get_subpool(domain % blocklist % structs, 'tend_physics', tend_physics)
         nullify(rthcuten)
         call mpas_pool_get_array(tend_physics, 'rthcuten', rthcuten)
         call mpas_pool_get_array(tend_physics, 'rqvcuten', rqvcuten)
         call mpas_pool_get_array(tend_physics, 'rucuten', rucuten)
         call mpas_pool_get_array(tend_physics, 'rvcuten', rvcuten)
         call mpas_pool_get_array(diag_physics, 'cu_trc', cu_trc)
         call mpas_pool_get_array(diag, 'uReconstructZonal', u_zon)
         call mpas_pool_get_array(diag, 'uReconstructMeridional', u_mer)
         call mpas_pool_get_array(mesh, 'zgrid', zgrid)
         if (.not. associated(rthcuten)) then
            call mpas_log_write('PartMC: PMC_CONV_CHECK needs rthcuten and rqvcuten', messageType=MPAS_LOG_CRIT)
         end if
         allocate(hm(nz), dh_m(nz), dh_s(nz))
      end if
      allocate(tsub_s(nz,nz), tmat_s(nz,nz), tk(nz), qs(nz), qv(nz), p_rem(nz), col_s(nz), gas_of(size(idx)))
      gas_of = 0
      if (pmc_wet .and. pmc_wet_conv) then
         call gas_wet_setup(domain)
         if (allocated(gas_wet_frac)) then
            ! the PartMC species of each mixed scalar
            do i_s = 5, size(idx)
               do j = 1, size(gas_idx)
                  if (gas_idx(j) == idx(i_s)) gas_of(i_s) = j
               end do
            end do
         end if
      end if
      allocate(m(nz), mu_top(nz), mu_bot(nz), dn_top(nz), dn_bot(nz), e_u(nz), d_u(nz), e_d(nz), d_d(nz), &
               f_u(nz), f_d(nz), s_top(nz), out(nz), tsub(nz,nz), tmat(nz,nz), tnew(nz,nz), col(nz))
      n_col = 0
      o_max = 0.0_RKIND

      do iCell = 1, nCellsSolve
         if (maxval(cu_mfu(:,iCell)) <= 1.0e-8_RKIND) cycle
         n_col = n_col + 1
         do k = 1, nz
            m(k) = rho_zz(k,iCell) / rdzw(k)
            mu_top(k) = max(0.0_RKIND, cu_mfu(k,iCell))
            dn_top(k) = max(0.0_RKIND, -cu_mfd(k,iCell))
         end do
         mu_top(nz) = 0.0_RKIND
         dn_top(nz) = 0.0_RKIND
         mu_bot(1) = 0.0_RKIND
         dn_bot(1) = 0.0_RKIND
         mu_bot(2:nz) = mu_top(1:nz-1)
         dn_bot(2:nz) = dn_top(1:nz-1)
         do k = 1, nz
            ! updraft: in at the bottom, out at the top
            d_u(k) = max(0.0_RKIND, cu_ude(k,iCell))
            e_u(k) = mu_top(k) - mu_bot(k) + d_u(k)
            if (e_u(k) < 0.0_RKIND) then
               d_u(k) = d_u(k) - e_u(k)
               e_u(k) = 0.0_RKIND
            end if
            pass = mu_bot(k) + e_u(k)
            f_u(k) = 1.0_RKIND
            if (pass > 0.0_RKIND) f_u(k) = min(1.0_RKIND, d_u(k)/pass)
            ! downdraft: in at the top, out at the bottom
            d_d(k) = max(0.0_RKIND, cu_dde(k,iCell))
            e_d(k) = dn_bot(k) - dn_top(k) + d_d(k)
            if (e_d(k) < 0.0_RKIND) then
               d_d(k) = d_d(k) - e_d(k)
               e_d(k) = 0.0_RKIND
            end if
            pass = dn_top(k) + e_d(k)
            f_d(k) = 1.0_RKIND
            if (pass > 0.0_RKIND) f_d(k) = min(1.0_RKIND, d_d(k)/pass)
            ! environment mass flux through the top of the layer, upward positive
            s_top(k) = dn_top(k) - mu_top(k)
         end do

         ! fraction of each layer that leaves it in the step sets the number of sub-steps
         do k = 1, nz
            out(k) = e_u(k) + e_d(k) + max(0.0_RKIND, s_top(k))
            if (k > 1) out(k) = out(k) + max(0.0_RKIND, -s_top(k-1))
            out(k) = out(k) * config_dt / m(k)
         end do
         o_max = max(o_max, maxval(out))
         n_sub = max(1, ceiling(maxval(out) / 0.5_RKIND))
         dt = config_dt / real(n_sub, kind=RKIND)

         ! Removal in the updraft. Air taken in at layer j with vapour qv(j) that leaves the
         ! updraft at layer k has condensed qv(j) - qs(k) on the way, if that is positive. The
         ! fraction pe of all the water condensed like this in the step falls as the scheme's
         ! convective precipitation at the surface; a tracer that is all in the condensate
         ! loses that fraction on every path with condensation.
         scav = .false.
         pe = 0.0_RKIND
         if (pmc_wet .and. pmc_wet_conv .and. associated(raincv)) then
            do k = 1, nz
               qv(k) = max(0.0_RKIND, scalars(index_qv,k,iCell))
               tk(k) = theta_m(k,iCell) / (1.0_RKIND + rvord*qv(k)) * exner(k,iCell)
               pres = pressure_p(k,iCell) + pressure_base(k,iCell)
               es = 611.2_RKIND * exp(17.67_RKIND * (tk(k) - 273.15_RKIND) / (tk(k) - 29.65_RKIND))
               qs(k) = 0.622_RKIND * es / max(pres - es, 1.0_RKIND)
            end do
            w_cond = 0.0_RKIND
            do j = 1, nz
               if (e_u(j) <= 0.0_RKIND) cycle
               pass = e_u(j) * config_dt
               do k = j, nz
                  w_cond = w_cond + pass*f_u(k) * max(0.0_RKIND, qv(j) - qs(k))
                  pass = pass*(1.0_RKIND - f_u(k))
                  if (pass <= 0.0_RKIND) exit
               end do
            end do
            conv_wet(4) = conv_wet(4) + max(0.0_RKIND, raincv(iCell)) * areaCell(iCell)
            conv_wet(5) = conv_wet(5) + w_cond * areaCell(iCell)
            if (w_cond > 0.0_RKIND .and. raincv(iCell) > 0.0_RKIND) then
               pe = min(1.0_RKIND, raincv(iCell) / w_cond)
               conv_wet(6) = conv_wet(6) + pe * w_cond * areaCell(iCell)
               scav = .not. pmc_wet_diag_only
            end if
         end if

         ! one sub-step: column j holds where the air of layer j goes
         tsub = 0.0_RKIND
         tsub_s = 0.0_RKIND
         do j = 1, nz
            tsub(j,j) = 1.0_RKIND - out(j)/real(n_sub, kind=RKIND)
            if (e_u(j) > 0.0_RKIND) then
               pass = e_u(j) * dt / m(j)
               do k = j, nz
                  tsub(k,j) = tsub(k,j) + pass*f_u(k)
                  if (scav) then
                     if (qv(j) > qs(k)) tsub_s(k,j) = tsub_s(k,j) - pe*pass*f_u(k)
                  end if
                  pass = pass*(1.0_RKIND - f_u(k))
                  if (pass <= 0.0_RKIND) exit
               end do
            end if
            if (e_d(j) > 0.0_RKIND) then
               pass = e_d(j) * dt / m(j)
               do k = j, 1, -1
                  tsub(k,j) = tsub(k,j) + pass*f_d(k)
                  pass = pass*(1.0_RKIND - f_d(k))
                  if (pass <= 0.0_RKIND) exit
               end do
            end if
            if (j < nz) then
               if (s_top(j) > 0.0_RKIND) tsub(j+1,j) = tsub(j+1,j) + s_top(j)*dt/m(j)
            end if
            if (j > 1) then
               if (s_top(j-1) < 0.0_RKIND) tsub(j-1,j) = tsub(j-1,j) - s_top(j-1)*dt/m(j)
            end if
         end do
         tmat = tsub
         do i_sub = 2, n_sub
            tnew = matmul(tsub, tmat)
            tmat = tnew
         end do
         ! Check against the scheme. Moist static energy h = cp T + g z + Lv qv is carried by
         ! the plumes almost without sources (condensation and evaporation move energy between
         ! its terms), and the scheme carries the horizontal wind as a passive quantity, so the
         ! matrix applied to them should give about the scheme's own changes.
         if (pmc_conv_check) then
            do i_f = 1, 4
               do k = 1, nz
                  if (i_f == 4) then
                     ! the scheme's own passive tracer: height above the ground
                     hm(k) = exp(-(0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell)) - zgrid(1,iCell)) / 3000.0_RKIND)
                     dh_s(k) = config_dt * cu_trc(k,iCell)
                  else if (i_f == 1) then
                     z_mid = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
                     hm(k) = 1004.5_RKIND * theta_m(k,iCell) &
                             / (1.0_RKIND + rvord*max(0.0_RKIND, scalars(index_qv,k,iCell))) * exner(k,iCell) &
                             + 9.80616_RKIND * z_mid + 2.5e6_RKIND * max(0.0_RKIND, scalars(index_qv,k,iCell))
                     dh_s(k) = config_dt * (1004.5_RKIND * rthcuten(k,iCell) * exner(k,iCell) &
                                            + 2.5e6_RKIND * rqvcuten(k,iCell))
                  else if (i_f == 2) then
                     hm(k) = u_zon(k,iCell)
                     dh_s(k) = config_dt * rucuten(k,iCell)
                  else
                     hm(k) = u_mer(k,iCell)
                     dh_s(k) = config_dt * rvcuten(k,iCell)
                  end if
               end do
               dh_m = matmul(tmat, hm * m) / m - hm
               do i_sh = 1, 3
                  do k = 1, nz
                     kk = k + i_sh - 2
                     if (kk < 1 .or. kk > nz) cycle
                     wgt = m(k) * areaCell(iCell)
                     conv_chk(1,i_f,i_sh) = conv_chk(1,i_f,i_sh) + wgt * dh_m(kk) * dh_s(k)
                     conv_chk(2,i_f,i_sh) = conv_chk(2,i_f,i_sh) + wgt * dh_s(k)**2
                     conv_chk(3,i_f,i_sh) = conv_chk(3,i_f,i_sh) + wgt * dh_m(kk)**2
                  end do
               end do
               slab_m = 0.0_RKIND
               slab_s = 0.0_RKIND
               do k = 1, nz
                  z_mid = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
                  i_sl = 1
                  if (z_mid > 2000.0_RKIND) i_sl = 2
                  if (z_mid > 6000.0_RKIND) i_sl = 3
                  slab_m(i_sl) = slab_m(i_sl) + m(k) * dh_m(k)
                  slab_s(i_sl) = slab_s(i_sl) + m(k) * dh_s(k)
               end do
               do i_sl = 1, 3
                  conv_slab(1,i_f,i_sl) = conv_slab(1,i_f,i_sl) + areaCell(iCell) * slab_m(i_sl) * slab_s(i_sl)
                  conv_slab(2,i_f,i_sl) = conv_slab(2,i_f,i_sl) + areaCell(iCell) * slab_s(i_sl)**2
                  conv_slab(3,i_f,i_sl) = conv_slab(3,i_f,i_sl) + areaCell(iCell) * slab_m(i_sl)**2
               end do
            end do
         end if

         ! the same with removal: tsub_s held the amounts removed
         if (scav) then
            tsub_s = tsub + tsub_s
            tmat_s = tsub_s
            do i_sub = 2, n_sub
               tnew = matmul(tsub_s, tmat_s)
               tmat_s = tnew
            end do
         end if

         ! Eulerian scalars: the mass of tracer in layer i is the sum over j of tmat(i,j) times
         ! the mass that was in layer j
         do i_s = 1, size(idx)
            col = scalars(idx(i_s),:,iCell) * m
            if (.not. scav) then
               scalars(idx(i_s),:,iCell) = matmul(tmat, col) / m
            else if (i_s <= 4) then
               ! number scalars: all in the condensate
               col_s = matmul(tmat_s, col)
               if (idx(i_s) == index_ref) conv_wet(1) = conv_wet(1) + (sum(col) - sum(col_s)) * areaCell(iCell)
               scalars(idx(i_s),:,iCell) = col_s / m
            else
               ! gases: the in-cloud fraction behaves as the number scalars, the rest is only moved
               e_gas = 0.0_RKIND
               if (gas_of(i_s) > 0) e_gas = gas_wet_frac(gas_of(i_s))
               col_s = matmul(tmat_s, col)
               lost = e_gas * (sum(col) - sum(col_s))
               scalars(idx(i_s),:,iCell) = (e_gas * col_s + (1.0_RKIND - e_gas) * matmul(tmat, col)) / m
               if (gas_of(i_s) > 0) then
                  gas_wetdep(gas_of(i_s)) = gas_wetdep(gas_of(i_s)) &
                                            + lost * 1.0e-9_RKIND / M_AIR * areaCell(iCell)
               end if
            end if
         end do

         if (pmc_n_part <= 0 .or. pmc_mode /= MODE_LAGR) cycle
         do j = 1, nz
            tmat(:,j) = max(0.0_RKIND, tmat(:,j))
            tmat(:,j) = tmat(:,j) / sum(tmat(:,j))
         end do
         if (scav) then
            do j = 1, nz
               tmat_s(:,j) = max(0.0_RKIND, tmat_s(:,j))
               p_rem(j) = min(1.0_RKIND, max(0.0_RKIND, 1.0_RKIND - sum(tmat_s(:,j))))
               if (p_rem(j) < 1.0_RKIND) tmat_s(:,j) = tmat_s(:,j) / sum(tmat_s(:,j))
            end do
            call column_jump(iCell, nz, tmat, n_conv, nCells, nEdgesOnCell, cellsOnCell, tmat_s, p_rem, tk)
         else
            call column_jump(iCell, nz, tmat, n_conv, nCells, nEdgesOnCell, cellsOnCell)
         end if
      end do

      call mpas_pool_get_field(state, 'scalars', scalars_field, 2)
      call mpas_dmpar_exch_halo_field(scalars_field)

      call mpas_dmpar_max_real(domain % dminfo, o_max, g_o_max)
      o_max_seen = max(o_max_seen, g_o_max)
      if (mod(itimestep, 12) == 0) then
         call mpas_dmpar_sum_int(domain % dminfo, n_col, g_n_col)
         call mpas_dmpar_sum_int(domain % dminfo, n_conv, g_n_conv)
         call mpas_log_write('PartMC convection step $i: convective columns = $i, particle layer changes so far = $i, ' &
                             // 'largest layer outflow fraction so far = $r', &
                             intArgs=(/itimestep, g_n_col, g_n_conv/), realArgs=(/o_max_seen/))
         if (pmc_conv_check) then
            csum = reshape(conv_chk, (/36/))
            call mpas_dmpar_sum_real_array(domain % dminfo, 36, csum, g_csum)
            do i_f = 1, 4
               ! element (a, i_f, i_sh) of conv_chk is at a + 3*(i_f-1) + 12*(i_sh-1)
               call mpas_log_write('PartMC convection check step $i: ' // trim(chk_name(i_f)) // ', matrix against ' &
                                   // 'scheme, slope and correlation with the matrix result taken from the layer ' &
                                   // 'below = $r $r, the same layer = $r $r, the layer above = $r $r', &
                                   intArgs=(/itimestep/), realArgs=(/ &
                  g_csum(1+3*(i_f-1))    / max(g_csum(2+3*(i_f-1)), 1.0e-30_RKIND), &
                  g_csum(1+3*(i_f-1))    / max(sqrt(g_csum(2+3*(i_f-1))*g_csum(3+3*(i_f-1))), 1.0e-30_RKIND), &
                  g_csum(1+3*(i_f-1)+12)  / max(g_csum(2+3*(i_f-1)+12), 1.0e-30_RKIND), &
                  g_csum(1+3*(i_f-1)+12)  / max(sqrt(g_csum(2+3*(i_f-1)+12)*g_csum(3+3*(i_f-1)+12)), 1.0e-30_RKIND), &
                  g_csum(1+3*(i_f-1)+24) / max(g_csum(2+3*(i_f-1)+24), 1.0e-30_RKIND), &
                  g_csum(1+3*(i_f-1)+24) / max(sqrt(g_csum(2+3*(i_f-1)+24)*g_csum(3+3*(i_f-1)+24)), 1.0e-30_RKIND)/))
            end do
            csum = reshape(conv_slab, (/36/))
            call mpas_dmpar_sum_real_array(domain % dminfo, 36, csum, g_csum)
            do i_f = 1, 4
               call mpas_log_write('PartMC convection check step $i: ' // trim(chk_name(i_f)) // ', column slabs, ' &
                                   // 'slope and correlation below 2 km = $r $r, 2-6 km = $r $r, above 6 km = $r $r', &
                                   intArgs=(/itimestep/), realArgs=(/ &
                  g_csum(1+3*(i_f-1))    / max(g_csum(2+3*(i_f-1)), 1.0e-30_RKIND), &
                  g_csum(1+3*(i_f-1))    / max(sqrt(g_csum(2+3*(i_f-1))*g_csum(3+3*(i_f-1))), 1.0e-30_RKIND), &
                  g_csum(1+3*(i_f-1)+12)  / max(g_csum(2+3*(i_f-1)+12), 1.0e-30_RKIND), &
                  g_csum(1+3*(i_f-1)+12)  / max(sqrt(g_csum(2+3*(i_f-1)+12)*g_csum(3+3*(i_f-1)+12)), 1.0e-30_RKIND), &
                  g_csum(1+3*(i_f-1)+24) / max(g_csum(2+3*(i_f-1)+24), 1.0e-30_RKIND), &
                  g_csum(1+3*(i_f-1)+24) / max(sqrt(g_csum(2+3*(i_f-1)+24)*g_csum(3+3*(i_f-1)+24)), 1.0e-30_RKIND)/))
            end do
         end if
         if (pmc_wet .and. pmc_wet_conv) then
            sums = conv_wet
            call mpas_dmpar_sum_real_array(domain % dminfo, 6, sums, g_sums)
            call mpas_log_write('PartMC convective removal step $i: number removed so far: reference = $r, ' &
                                // 'particles = $r, particles expected from the layers'' content = $r; so far (kg): ' &
                                // 'convective precipitation = $r, condensed in the updraft = $r, of that counted ' &
                                // 'as precipitation = $r', intArgs=(/itimestep/), realArgs=g_sums)
         end if
      end if
#endif

   end subroutine atm_partmc_convection


   !> Read one variable of the emissions file (a value for every cell of the mesh, in global
   !> cell order) for the owned cells. found is false if the file has no such variable.
   subroutine emis_read(domain, name, vals, found)

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: name
      real (kind=RKIND), intent(out) :: vals(:)
      logical, intent(out) :: found

#ifdef MPAS_PARTMC
      type (mpas_pool_type), pointer :: mesh
      integer, dimension(:), pointer :: indexToCellID
      integer :: ncid, varid, dimids(1), n_glob, iCell
      real (kind=RKIND), allocatable :: glob(:)

      found = .false.
      vals = 0.0_RKIND
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)
      call nc_check(nf90_open(trim(emis_file), NF90_NOWRITE, ncid), 'open ' // trim(emis_file))
      if (nf90_inq_varid(ncid, name, varid) == NF90_NOERR) then
         call nc_check(nf90_inquire_variable(ncid, varid, dimids=dimids), 'inquire ' // name)
         call nc_check(nf90_inquire_dimension(ncid, dimids(1), len=n_glob), 'dimension of ' // name)
         allocate(glob(n_glob))
         call nc_check(nf90_get_var(ncid, varid, glob), 'read ' // name)
         do iCell = 1, size(vals)
            vals(iCell) = glob(indexToCellID(iCell))
         end do
         found = .true.
      end if
      call nc_check(nf90_close(ncid), 'close ' // trim(emis_file))
#else
      found = .false.
      vals = 0.0_RKIND
#endif

   end subroutine emis_read


   !> Surface emission and dry deposition in the lowest layer, for the Eulerian number scalars
   !> and, in Lagrangian mode, the particles.
   !>
   !> Dry deposition: a constant deposition velocity v acts on the lowest layer of depth dz, so
   !> a fraction 1 - exp(-v dt / dz) of the layer's number is removed; each particle there is
   !> removed with that probability. Emission: a uniform number flux over land of soot-like
   !> particles (the composition and sizes of source 2). The reference scalars get the flux
   !> divided by the layer's air mass. New computational particles get the weight level that
   !> the cell's content calls for, and their number is rounded at random, so emission is
   !> exact in expectation.
   subroutine atm_partmc_surface(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC
      type (mpas_pool_type), pointer :: mesh, state, sfc_input
      type (field3DReal), pointer :: scalars_field
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, zz
      real (kind=RKIND), dimension(:), pointer :: rdzw, areaCell
      integer, dimension(:), pointer :: landmask
      real (kind=RKIND), pointer :: config_dt
      integer, dimension(:,:), pointer :: cellsOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, pointer :: nCells, nCellsSolve, index_pmc_num, index_pmc_num2, index_ref, index_src2
      integer :: iCell, j, k, g, n, i, i_oc, i_bc, idx(4), i_s
      real (kind=RKIND) :: dz, m1, mk, m_emit, p_dep, n_cell, w, f, x(NPOS), sums(4), g_sums(4)
      real (kind=dp) :: diam, vol_tot
      real (kind=dp), allocatable :: vols(:)
      type(aero_particle_t) :: particle
      type (mpas_pool_type), pointer :: diag, diag_physics
      real (kind=RKIND), dimension(:,:), pointer :: theta_m, exner, pressure_p, pressure_base, zgrid, u_zon, u_mer
      real (kind=RKIND), dimension(:), pointer :: znt
      integer, pointer :: index_qv
      type(scenario_t) :: scen
      type(env_state_t) :: env_state
      real (kind=dp) :: v_dep
      real (kind=RKIND) :: dsum(4), g_dsum(4)
      real (kind=dp), parameter :: test_diam(6) = (/ 1.0d-8, 3.0d-8, 1.0d-7, 3.0d-7, 1.0d-6, 1.0d-5 /)
      real (kind=RKIND) :: test_v(6)

      real (kind=RKIND), dimension(:), pointer :: u10, v10, xland
      type (mpas_pool_type), pointer :: diag_phys2
      real (kind=RKIND) :: f_bc, rho_mix, m_mean, n_flux, wind, r80, dfdr, esum(7), g_esum(7)
      logical :: found
      integer :: ib, i_na, i_cl

      if (pmc_emit_flux <= 0.0_RKIND .and. pmc_drydep_v <= 0.0_RKIND .and. .not. pmc_drydep_scheme &
          .and. len_trim(emis_file) == 0 .and. .not. pmc_seasalt) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'sfc_input', sfc_input)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', config_dt)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(state, 'index_pmc_num', index_pmc_num)
      call mpas_pool_get_dimension(state, 'index_pmc_num2', index_pmc_num2)
      call mpas_pool_get_dimension(state, 'index_pmc_ref', index_ref)
      call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', index_src2)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      ! landmask, not xland: the physics sets xland to land where there is sea ice
      nullify(landmask)
      call mpas_pool_get_array(sfc_input, 'landmask', landmask)
      idx = (/ index_pmc_num, index_pmc_num2, index_ref, index_src2 /)
      if (pmc_drydep_scheme) then
         call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
         call mpas_pool_get_subpool(domain % blocklist % structs, 'diag_physics', diag_physics)
         call mpas_pool_get_dimension(state, 'index_qv', index_qv)
         call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
         call mpas_pool_get_array(diag, 'exner', exner)
         call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
         call mpas_pool_get_array(diag, 'pressure_base', pressure_base)
         call mpas_pool_get_array(diag, 'uReconstructZonal', u_zon)
         call mpas_pool_get_array(diag, 'uReconstructMeridional', u_mer)
         call mpas_pool_get_array(mesh, 'zgrid', zgrid)
         nullify(znt)
         call mpas_pool_get_array(diag_physics, 'znt', znt)
         if (itimestep == 1 .and. domain % dminfo % my_proc_id == 0) then
            ! deposition velocity against diameter for fixed conditions, for a look at the curve
            env_state%temp = 288.0d0
            env_state%pressure = 1.0d5
            env_state%height = 1.0d0
            scen%drydep%z_ref = 25.0d0
            scen%drydep%u_mean = 5.0d0
            scen%drydep%z_rough = 0.01d0
            do i = 1, 6
               test_v(i) = scenario_loss_rate_drydep(pii/6.0d0*test_diam(i)**3, 1800.0d0, aero_data, env_state, scen)
            end do
            call mpas_log_write('PartMC dry deposition velocity (m/s) for diameters 10 nm, 30 nm, 100 nm, 300 nm, ' &
                                // '1 um, 10 um (density 1800 kg/m3, 5 m/s at 25 m, roughness 0.01 m): ' &
                                // '$r $r $r $r $r $r', realArgs=test_v)
         end if
      end if
      i_oc = aero_data_spec_by_name(aero_data, 'OC')
      i_bc = aero_data_spec_by_name(aero_data, 'BC')
      i_na = aero_data_spec_by_name(aero_data, 'Na')
      i_cl = aero_data_spec_by_name(aero_data, 'Cl')
      allocate(vols(n_spec))

      if (len_trim(emis_file) > 0 .and. .not. emis_part_read) then
         emis_part_read = .true.
         allocate(emis_bc(nCellsSolve), emis_oc(nCellsSolve))
         call emis_read(domain, 'bc', emis_bc, found)
         call emis_read(domain, 'oc', emis_oc, found)
         emis_bc = max(0.0_RKIND, emis_bc) * emis_scale
         emis_oc = max(0.0_RKIND, emis_oc) * emis_scale
         call mpas_log_write('PartMC: BC and OC emissions from ' // trim(emis_file) // ' times $r, one mode with ' &
                             // 'geometric mean diameter $r m', realArgs=(/emis_scale, carb_dg/))
      end if
      if (pmc_seasalt) then
         call mpas_pool_get_subpool(domain % blocklist % structs, 'diag_physics', diag_phys2)
         nullify(u10)
         call mpas_pool_get_array(diag_phys2, 'u10', u10)
         call mpas_pool_get_array(diag_phys2, 'v10', v10)
         nullify(xland)
         call mpas_pool_get_array(sfc_input, 'xland', xland)
         if (.not. associated(u10) .or. .not. associated(xland)) then
            call mpas_log_write('PartMC: PMC_SEASALT needs u10, v10 and xland', messageType=MPAS_LOG_CRIT)
         end if
         if (ss_unit <= 0.0_RKIND) then
            ! Gong (2003) flux dF/dr (m^-2 s^-1 um^-1) at 1 m/s against the radius at 80% humidity
            ! r (um), which is about the dry diameter; N_SS bins equal in log r from 0.07 to 10 um
            do ib = 1, N_SS + 1
               ss_edge(ib) = 0.07e-6_RKIND * (10.0_RKIND / 0.07_RKIND)**(real(ib - 1, kind=RKIND) / real(N_SS, kind=RKIND))
            end do
            do ib = 1, N_SS
               r80 = sqrt(ss_edge(ib) * ss_edge(ib+1)) * 1.0e6_RKIND
               dfdr = 1.373_RKIND * r80**(-4.7_RKIND * (1.0_RKIND + 30.0_RKIND*r80)**(-0.017_RKIND * r80**(-1.44_RKIND))) &
                      * (1.0_RKIND + 0.057_RKIND * r80**3.45_RKIND) &
                      * 10.0_RKIND**(1.607_RKIND * exp(-((0.433_RKIND - log10(r80)) / 0.433_RKIND)**2))
               ss_cum(ib) = dfdr * (ss_edge(ib+1) - ss_edge(ib)) * 1.0e6_RKIND
               if (ib > 1) ss_cum(ib) = ss_cum(ib) + ss_cum(ib-1)
            end do
            ss_unit = ss_cum(N_SS)
            ss_cum = ss_cum / ss_unit
            call mpas_log_write('PartMC: sea salt number flux at 10 m/s = $r m^-2 s^-1, times $r', &
                                realArgs=(/ss_unit * 10.0_RKIND**3.41_RKIND, seasalt_scale/))
         end if
      end if

      do iCell = 1, nCellsSolve
         dz = 1.0_RKIND / (rdzw(1) * zz(1,iCell))
         ! air mass per area of the lowest layer; number per area is m1 times the mixing ratio
         m1 = rho_zz(1,iCell) / rdzw(1)

         if (pmc_drydep_scheme .and. pmc_n_part > 0 .and. pmc_mode == MODE_LAGR) then
            ! a deposition velocity for every particle
            env_state%temp = theta_m(1,iCell) / (1.0_RKIND + rvord*max(0.0_RKIND, scalars(index_qv,1,iCell))) &
                             * exner(1,iCell)
            env_state%pressure = pressure_p(1,iCell) + pressure_base(1,iCell)
            env_state%height = 1.0d0
            scen%drydep%z_ref = 0.5d0 * (zgrid(2,iCell) - zgrid(1,iCell))
            scen%drydep%u_mean = max(0.5d0, sqrt(u_zon(1,iCell)**2 + u_mer(1,iCell)**2))
            scen%drydep%z_rough = 0.01d0
            if (associated(znt)) scen%drydep%z_rough = znt(iCell)
            scen%drydep%z_rough = min(0.3d0 * scen%drydep%z_ref, max(1.0d-5, scen%drydep%z_rough))
            do j = pos(1,iCell) % n, 1, -1
               v_dep = scenario_loss_rate_drydep(aero_particle_volume(aero_states(1,iCell)%apa%particle(j)), &
                          aero_particle_density(aero_states(1,iCell)%apa%particle(j), aero_data), aero_data, &
                          env_state, scen)
               p_dep = 1.0_RKIND - exp(-v_dep * config_dt / dz)
               g = aero_states(1,iCell)%apa%particle(j)%weight_class
               w = lagr_w * 2.0_RKIND**(g - LAGR_G0)
               drydep_sums(1) = drydep_sums(1) + p_dep * w
               drydep_sums(3) = drydep_sums(3) + v_dep * w
               drydep_sums(4) = drydep_sums(4) + w
               if (pmc_random() >= p_dep) cycle
               drydep_sums(2) = drydep_sums(2) + w
               part_s_removed = part_s_removed + particle_sulfate(aero_states(1,iCell)%apa%particle(j)) * w
               call aero_state_remove_particle_no_info(aero_states(1,iCell), j)
               call pos_remove(pos(1,iCell), j)
            end do
         else if (pmc_drydep_v > 0.0_RKIND) then
            p_dep = 1.0_RKIND - exp(-pmc_drydep_v * config_dt / dz)
            surf_sums(3) = surf_sums(3) + p_dep * m1 * scalars(index_ref,1,iCell) * areaCell(iCell)
            do i_s = 1, 4
               scalars(idx(i_s),1,iCell) = scalars(idx(i_s),1,iCell) * (1.0_RKIND - p_dep)
            end do
            if (pmc_n_part > 0 .and. pmc_mode == MODE_LAGR) then
               do j = pos(1,iCell) % n, 1, -1
                  if (pmc_random() >= p_dep) cycle
                  g = aero_states(1,iCell)%apa%particle(j)%weight_class
                  surf_sums(4) = surf_sums(4) + lagr_w * 2.0_RKIND**(g - LAGR_G0)
                  part_s_removed = part_s_removed + particle_sulfate(aero_states(1,iCell)%apa%particle(j)) &
                                                    * lagr_w * 2.0_RKIND**(g - LAGR_G0)
                  call aero_state_remove_particle_no_info(aero_states(1,iCell), j)
                  call pos_remove(pos(1,iCell), j)
               end do
            end if
         end if

         if (pmc_emit_flux > 0.0_RKIND .and. associated(landmask)) then
            if (landmask(iCell) == 1) then
               ! spread over the lowest layers in proportion to their air mass
               m_emit = 0.0_RKIND
               do k = 1, pmc_emit_levels
                  m_emit = m_emit + rho_zz(k,iCell) / rdzw(k)
               end do
               surf_sums(1) = surf_sums(1) + pmc_emit_flux * config_dt * areaCell(iCell)
               do k = 1, pmc_emit_levels
                  mk = rho_zz(k,iCell) / rdzw(k)
                  f = pmc_emit_flux * config_dt * mk / m_emit
                  scalars(index_ref,k,iCell) = scalars(index_ref,k,iCell) + f / mk
                  scalars(index_src2,k,iCell) = scalars(index_src2,k,iCell) + f / mk
                  if (pmc_n_part <= 0 .or. pmc_mode /= MODE_LAGR) cycle
                  n_cell = mk * scalars(index_ref,k,iCell) * areaCell(iCell)
                  g = LAGR_G0 + nint(log(n_cell / (real(pmc_n_part, kind=RKIND) * lagr_w)) / log(2.0_RKIND))
                  g = max(1, min(LAGR_GROUPS, g))
                  w = lagr_w * 2.0_RKIND**(g - LAGR_G0)
                  f = f * areaCell(iCell) / w
                  n = int(f)
                  if (pmc_random() < f - real(n, kind=RKIND)) n = n + 1
                  do i = 1, n
                     call aero_particle_zero(particle, aero_data)
                     diam = 0.05d-6 * exp(log(1.8d0) * rand_normal(0.0d0, 1.0d0))
                     vol_tot = pii / 6.0d0 * diam**3
                     vols = 0.0d0
                     vols(i_bc) = 0.6d0 * vol_tot
                     vols(i_oc) = 0.4d0 * vol_tot
                     call aero_particle_set_component(particle, 2, 0.0d0)
                     particle%n_primary_parts = 1
                     call aero_particle_set_vols(particle, vols)
                     call aero_particle_new_id(particle)
                     call aero_particle_set_weight(particle, 1, g)
                     call aero_state_add_particle(aero_states(k,iCell), particle, aero_data)
                     call random_in_cell(iCell, k, nCells, nEdgesOnCell, cellsOnCell, x(1:4))
                     x(5:7) = x(1:3)
                     call pos_append(pos(k,iCell), x)
                  end do
                  surf_sums(2) = surf_sums(2) + real(n, kind=RKIND) * w
                  if (lagr_control .and. k > 1) call lagr_rebalance(k, iCell, nCells, nEdgesOnCell, cellsOnCell)
               end do
            end if
         end if

         ! carbonaceous particles from the inventory
         if (allocated(emis_bc)) then
            if (emis_bc(iCell) + emis_oc(iCell) > 0.0_RKIND) then
               f_bc = emis_bc(iCell) / (emis_bc(iCell) + emis_oc(iCell))
               rho_mix = 1.0_RKIND / (f_bc / aero_data%density(i_bc) + (1.0_RKIND - f_bc) / aero_data%density(i_oc))
               m_mean = rho_mix * pii / 6.0_RKIND * carb_dg**3 * exp(4.5_RKIND * log(carb_sigma)**2)
               n_flux = (emis_bc(iCell) + emis_oc(iCell)) / m_mean
               emis_sums(3) = emis_sums(3) + (emis_bc(iCell) + emis_oc(iCell)) * config_dt * areaCell(iCell)
               m_emit = 0.0_RKIND
               do k = 1, pmc_emit_levels
                  m_emit = m_emit + rho_zz(k,iCell) / rdzw(k)
               end do
               do k = 1, pmc_emit_levels
                  call emit_kind(3, k, n_flux * config_dt * (rho_zz(k,iCell) / rdzw(k)) / m_emit)
               end do
            end if
         end if
         ! sea salt over open water
         if (pmc_seasalt) then
            if (nint(xland(iCell)) == 2) then
               wind = sqrt(u10(iCell)**2 + v10(iCell)**2)
               call emit_kind(4, 1, seasalt_scale * ss_unit * wind**3.41_RKIND * config_dt)
            end if
         end if

         if (pmc_n_part > 0 .and. pmc_mode == MODE_LAGR .and. lagr_control) then
            call lagr_rebalance(1, iCell, nCells, nEdgesOnCell, cellsOnCell)
         end if
      end do

      call mpas_pool_get_field(state, 'scalars', scalars_field, 2)
      call mpas_dmpar_exch_halo_field(scalars_field)

      if (mod(itimestep, 12) == 0 .and. (allocated(emis_bc) .or. pmc_seasalt)) then
         esum = emis_sums
         call mpas_dmpar_sum_real_array(domain % dminfo, 7, esum, g_esum)
         call mpas_log_write('PartMC emission step $i: carbonaceous number so far: reference = $r, particles = $r; ' &
                             // 'mass (kg): inventory = $r, particles = $r; sea salt number: expected = $r, ' &
                             // 'particles = $r, mass (kg) = $r', intArgs=(/itimestep/), realArgs=g_esum)
      end if
      if (mod(itimestep, 12) == 0) then
         sums = surf_sums
         call mpas_dmpar_sum_real_array(domain % dminfo, 4, sums, g_sums)
         call mpas_log_write('PartMC surface step $i: emitted so far: reference = $r, particles = $r; ' &
                             // 'deposited so far: reference = $r, particles = $r', &
                             intArgs=(/itimestep/), realArgs=g_sums)
         if (pmc_drydep_scheme) then
            dsum = drydep_sums
            call mpas_dmpar_sum_real_array(domain % dminfo, 4, dsum, g_dsum)
            call mpas_log_write('PartMC dry deposition step $i: size-dependent: number expected so far = $r, ' &
                                // 'removed = $r, number-weighted mean velocity in the lowest layer = $r m/s', &
                                intArgs=(/itimestep/), &
                                realArgs=(/g_dsum(1), g_dsum(2), g_dsum(3)/max(g_dsum(4), 1.0e-30_RKIND)/))
         end if
      end if

   contains

      !> Emit n_area real particles per m^2 of source i_kind (3 carbonaceous, 4 sea salt) into
      !> layer kk of the current cell: the reference scalar gets the number, and computational
      !> particles are made at the weight level that suits the cell, their number rounded at
      !> random.
      subroutine emit_kind(i_kind, kk, n_area)
         integer, intent(in) :: i_kind, kk
         real (kind=RKIND), intent(in) :: n_area
         real (kind=RKIND) :: m_kk, n_have, w_new, f_new, u
         integer :: g_new, n_new, i_new, jp, ibin
         real (kind=dp) :: d_new, v_new, m_new

         if (n_area <= 0.0_RKIND) return
         m_kk = rho_zz(kk,iCell) / rdzw(kk)
         scalars(index_ref,kk,iCell) = scalars(index_ref,kk,iCell) + n_area / m_kk
         if (i_kind == 3) emis_sums(1) = emis_sums(1) + n_area * areaCell(iCell)
         if (i_kind == 4) emis_sums(5) = emis_sums(5) + n_area * areaCell(iCell)
         if (pmc_n_part <= 0 .or. pmc_mode /= MODE_LAGR) return
         ! number of real particles the cell holds, from its particles if it has any
         n_have = 0.0_RKIND
         do jp = 1, pos(kk,iCell) % n
            n_have = n_have + lagr_w * 2.0_RKIND**(aero_states(kk,iCell)%apa%particle(jp)%weight_class - LAGR_G0)
         end do
         n_have = max(n_have, n_area * areaCell(iCell))
         g_new = LAGR_G0 + nint(log(n_have / (real(pmc_n_part, kind=RKIND) * lagr_w)) / log(2.0_RKIND))
         g_new = max(1, min(LAGR_GROUPS, g_new))
         w_new = lagr_w * 2.0_RKIND**(g_new - LAGR_G0)
         f_new = n_area * areaCell(iCell) / w_new
         n_new = int(f_new)
         if (pmc_random() < f_new - real(n_new, kind=RKIND)) n_new = n_new + 1
         do i_new = 1, n_new
            call aero_particle_zero(particle, aero_data)
            vols = 0.0d0
            if (i_kind == 3) then
               d_new = carb_dg * exp(log(carb_sigma) * rand_normal(0.0d0, 1.0d0))
               v_new = pii / 6.0d0 * d_new**3
               ! volumes from the mass fraction of BC
               vols(i_bc) = v_new * (f_bc / aero_data%density(i_bc)) * rho_mix
               vols(i_oc) = v_new - vols(i_bc)
               m_new = v_new * rho_mix
               emis_sums(2) = emis_sums(2) + w_new
               emis_sums(4) = emis_sums(4) + w_new * m_new
            else
               u = pmc_random()
               do ibin = 1, N_SS - 1
                  if (u <= ss_cum(ibin)) exit
               end do
               d_new = ss_edge(ibin) * (ss_edge(ibin+1) / ss_edge(ibin))**pmc_random()
               v_new = pii / 6.0d0 * d_new**3
               ! NaCl by mass: 39.34% Na, 60.66% Cl
               m_new = v_new / (0.3934d0 / aero_data%density(i_na) + 0.6066d0 / aero_data%density(i_cl))
               vols(i_na) = 0.3934d0 * m_new / aero_data%density(i_na)
               vols(i_cl) = 0.6066d0 * m_new / aero_data%density(i_cl)
               emis_sums(6) = emis_sums(6) + w_new
               emis_sums(7) = emis_sums(7) + w_new * m_new
            end if
            call aero_particle_set_component(particle, i_kind, 0.0d0)
            particle%n_primary_parts = 1
            call aero_particle_set_vols(particle, vols)
            call aero_particle_new_id(particle)
            call aero_particle_set_weight(particle, 1, g_new)
            call aero_state_add_particle(aero_states(kk,iCell), particle, aero_data)
            call random_in_cell(iCell, kk, nCells, nEdgesOnCell, cellsOnCell, x(1:4))
            x(5:7) = x(1:3)
            call pos_append(pos(kk,iCell), x)
         end do
         if (lagr_control .and. kk > 1) call lagr_rebalance(kk, iCell, nCells, nEdgesOnCell, cellsOnCell)
      end subroutine emit_kind
#endif

   end subroutine atm_partmc_surface


   !> Read the gas species and initial state, find each species' model scalar, set gases that
   !> are still zero everywhere to the initial state, and initialise MOSAIC. Does this once.
   subroutine gas_setup(domain)

      type (domain_type), intent(inout) :: domain

#ifdef MPAS_PARTMC_MOSAIC
      type (mpas_pool_type), pointer :: mesh, state
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), pointer :: config_dt
      integer, pointer :: nCellsSolve, nVertLevels, i_ptr
      integer :: iCell, k, j, n_gas
      type(spec_file_t) :: file
      type(env_state_t) :: env_state
      type(gas_state_t) :: gas_state
      real (kind=dp) :: del_t

      if (allocated(gas_idx)) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', config_dt)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      del_t = config_dt

      if (.not. allocated(gas_idx)) then
         call spec_file_open(trim(gas_data_file), file)
         call spec_file_read_gas_data(file, gas_data)
         call spec_file_close(file)
         call spec_file_open(trim(gas_init_file), file)
         call spec_file_read_gas_state(file, gas_data, gas_state)
         call spec_file_close(file)
         n_gas = gas_data_n_spec(gas_data)
         allocate(gas_idx(n_gas))
         do j = 1, n_gas
            nullify(i_ptr)
            call mpas_pool_get_dimension(state, 'index_gas_' // trim(gas_data%name(j)), i_ptr)
            if (.not. associated(i_ptr)) then
               call mpas_log_write('PartMC: no model scalar for gas ' // trim(gas_data%name(j)), &
                                   messageType=MPAS_LOG_CRIT)
            end if
            gas_idx(j) = i_ptr
         end do
         if (maxval(scalars(gas_idx,:,1:nCellsSolve)) <= 0.0_RKIND) then
            do iCell = 1, nCellsSolve
               do k = 1, nVertLevels
                  scalars(gas_idx,k,iCell) = gas_state%mix_rat
               end do
            end do
         end if
         env_state%temp = 288.0d0
         env_state%pressure = 1.0d5
         env_state%rel_humid = 0.5d0
         env_state%latitude = 0.0d0
         env_state%longitude = 0.0d0
         env_state%altitude = 0.0d0
         call mosaic_init(env_state, aero_data, del_t, .false.)
         call mpas_log_write('PartMC: MOSAIC with $i gas species from ' // trim(gas_data_file) // ', initial state ' &
                             // trim(gas_init_file) // ', below $r m', intArgs=(/n_gas/), realArgs=(/pmc_chem_top/))
      end if
#endif

   end subroutine gas_setup


   !> Global content (mol) of each gas species, summed over this task's cells.
   subroutine gas_content(domain, content)

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(out) :: content(:)

      type (mpas_pool_type), pointer :: mesh, state
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz
      real (kind=RKIND), dimension(:), pointer :: rdzw, areaCell
      integer, pointer :: nCellsSolve, nVertLevels
      integer :: iCell, k

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      content = 0.0_RKIND
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            ! mol of dry air in the cell times the mixing ratio
            content = content + scalars(gas_idx,k,iCell) * 1.0e-9_RKIND &
                                * rho_zz(k,iCell) / rdzw(k) * areaCell(iCell) / M_AIR
         end do
      end do

   end subroutine gas_content


   !> Surface emission and dry deposition of the gas tracers.
   !>
   !> Emission: the fluxes of a PartMC gas emissions file (mol m^-2 s^-1, with its time profile
   !> read against the time since the start of the run) over land, added to the lowest
   !> PMC_EMIT_LEVELS layers uniformly in mixing ratio. Dry deposition: a deposition velocity v
   !> for each species acts on the lowest layer of depth dz, removing the fraction
   !> 1 - exp(-v dt / dz). The amounts emitted and deposited are accumulated for a budget.
   subroutine atm_partmc_gas_surface(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC_MOSAIC
      type (mpas_pool_type), pointer :: mesh, state, sfc_input
      type (field3DReal), pointer :: scalars_field
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, zz
      real (kind=RKIND), dimension(:), pointer :: rdzw, areaCell
      integer, dimension(:), pointer :: landmask
      real (kind=RKIND), pointer :: config_dt
      integer, pointer :: nCellsSolve
      integer :: iCell, k, j, n_gas, j_e, j_d
      real (kind=RKIND) :: dz, air1, air_emit
      real (kind=RKIND), allocatable :: local(:), g_content(:), g_emit(:), g_dep(:), frac(:), amount(:)
      type(spec_file_t) :: file
      type(gas_state_t), save :: gas_drydep, emis
      type(gas_state_t), allocatable, save :: gas_emit_states(:)
      real (kind=dp), allocatable, save :: gas_emit_time(:), gas_emit_rate(:)
      real (kind=dp) :: rate
      logical, save :: gas_surface_read = .false.
      real (kind=RKIND), allocatable :: cell_vals(:)
      integer :: n_found
      logical :: found
      real (kind=RKIND), allocatable :: g_wet(:)
      integer :: j_w

      if (.not. pmc_chem .or. pmc_n_part <= 0) return
      if (len_trim(gas_emit_file) == 0 .and. len_trim(gas_drydep_file) == 0 .and. len_trim(emis_file) == 0 &
          .and. .not. (pmc_wet .and. len_trim(gas_wet_file) > 0)) return

      call gas_setup(domain)
      n_gas = gas_data_n_spec(gas_data)
      allocate(local(n_gas), g_content(n_gas), g_emit(n_gas), g_dep(n_gas), frac(n_gas), amount(n_gas), g_wet(n_gas))

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'sfc_input', sfc_input)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', config_dt)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      ! landmask, not xland: the physics sets xland to land where there is sea ice
      nullify(landmask)
      call mpas_pool_get_array(sfc_input, 'landmask', landmask)

      call gas_budget_setup(domain)
      if (.not. gas_surface_read) then
         gas_surface_read = .true.
         if (len_trim(emis_file) > 0) then
            allocate(emis_gas(n_gas, nCellsSolve), cell_vals(nCellsSolve))
            emis_gas = 0.0_RKIND
            n_found = 0
            do j = 1, n_gas
               call emis_read(domain, 'gas_' // trim(gas_data%name(j)), cell_vals, found)
               if (found) then
                  emis_gas(j,:) = max(0.0_RKIND, cell_vals) * emis_scale
                  n_found = n_found + 1
               end if
            end do
            call mpas_log_write('PartMC: gas emissions for $i species from ' // trim(emis_file) // ' times $r', &
                                intArgs=(/n_found/), realArgs=(/emis_scale/))
         end if
         if (len_trim(gas_emit_file) > 0) then
            call spec_file_open(trim(gas_emit_file), file)
            call spec_file_read_gas_states_times_rates(file, gas_data, gas_emit_time, gas_emit_rate, gas_emit_states)
            call spec_file_close(file)
            call mpas_log_write('PartMC: gas emissions from ' // trim(gas_emit_file) // ' times $r over land', &
                                realArgs=(/gas_emit_scale/))
         end if
         if (len_trim(gas_drydep_file) > 0) then
            call spec_file_open(trim(gas_drydep_file), file)
            call spec_file_read_gas_state(file, gas_data, gas_drydep)
            call spec_file_close(file)
            call mpas_log_write('PartMC: gas dry deposition velocities from ' // trim(gas_drydep_file) &
                                // ', largest $r m/s', realArgs=(/real(maxval(gas_drydep%mix_rat), kind=RKIND)/))
         end if
      end if

      if (allocated(gas_emit_states)) then
         ! mol m^-2 s^-1 in emis%mix_rat, times rate
         call gas_state_interp_1d(gas_emit_states, gas_emit_time, gas_emit_rate, &
                                  real(itimestep - 1, kind=dp) * real(config_dt, kind=dp), emis, rate)
      end if

      do iCell = 1, nCellsSolve
         ! mol of dry air per area in the lowest layer
         air1 = rho_zz(1,iCell) / rdzw(1) / M_AIR
         if (allocated(gas_drydep%mix_rat)) then
            dz = 1.0_RKIND / (rdzw(1) * zz(1,iCell))
            frac = 1.0_RKIND - exp(-gas_drydep%mix_rat * config_dt / dz)
            amount = frac * scalars(gas_idx,1,iCell)
            scalars(gas_idx,1,iCell) = scalars(gas_idx,1,iCell) - amount
            gas_deposited = gas_deposited + amount * 1.0e-9_RKIND * air1 * areaCell(iCell)
         end if
         if (allocated(gas_emit_states) .and. associated(landmask)) then
            if (landmask(iCell) == 1) then
               air_emit = 0.0_RKIND
               do k = 1, pmc_emit_levels
                  air_emit = air_emit + rho_zz(k,iCell) / rdzw(k) / M_AIR
               end do
               ! mol per area this step
               amount = emis%mix_rat * rate * gas_emit_scale * config_dt
               do k = 1, pmc_emit_levels
                  scalars(gas_idx,k,iCell) = scalars(gas_idx,k,iCell) + amount / air_emit * 1.0e9_RKIND
               end do
               gas_emitted = gas_emitted + amount * areaCell(iCell)
            end if
         end if
         if (allocated(emis_gas)) then
            ! the inventory's flux of this cell, mol per area this step
            air_emit = 0.0_RKIND
            do k = 1, pmc_emit_levels
               air_emit = air_emit + rho_zz(k,iCell) / rdzw(k) / M_AIR
            end do
            amount = emis_gas(:,iCell) * config_dt
            do k = 1, pmc_emit_levels
               scalars(gas_idx,k,iCell) = scalars(gas_idx,k,iCell) + amount / air_emit * 1.0e9_RKIND
            end do
            gas_emitted = gas_emitted + amount * areaCell(iCell)
         end if
      end do

      call mpas_pool_get_field(state, 'scalars', scalars_field, 2)
      call mpas_dmpar_exch_halo_field(scalars_field)

      if (mod(itimestep, 12) == 0 .or. itimestep == 1) then
         call gas_content(domain, local)
         call mpas_dmpar_sum_real_array(domain % dminfo, n_gas, local, g_content)
         call mpas_dmpar_sum_real_array(domain % dminfo, n_gas, gas_emitted, g_emit)
         call mpas_dmpar_sum_real_array(domain % dminfo, n_gas, gas_deposited, g_dep)
         call mpas_dmpar_sum_real_array(domain % dminfo, n_gas, gas_wetdep, g_wet)
         ! the species with the most emitted, dry deposited and wet deposited so far
         j_e = maxloc(g_emit, dim=1)
         j_d = maxloc(g_dep, dim=1)
         j_w = maxloc(g_wet, dim=1)
         do j = 1, n_gas
            if (j /= j_e .and. j /= j_d .and. j /= j_w) cycle
            call mpas_log_write('PartMC gas surface step $i: ' // trim(gas_data%name(j)) // ' (mol) content = $r, ' &
                                // 'at start = $r, emitted = $r, dry deposited = $r, wet deposited = $r, ' &
                                // 'budget residual / content = $r', &
                                intArgs=(/itimestep/), &
                                realArgs=(/g_content(j), gas_content0(j), g_emit(j), g_dep(j), g_wet(j), &
                                           (g_content(j) - gas_content0(j) - g_emit(j) + g_dep(j) + g_wet(j)) &
                                           / max(g_content(j), 1.0e-30_RKIND)/))
         end do
      end if
#endif

   end subroutine atm_partmc_gas_surface


   !> Sulfate (mol) in one real particle represented by a computational particle.
   function particle_sulfate(particle)
      type(aero_particle_t), intent(in) :: particle
      real (kind=RKIND) :: particle_sulfate
      integer, save :: i_so4 = -1
      if (i_so4 < 0) i_so4 = aero_data_spec_by_name(aero_data, 'SO4')
      particle_sulfate = 0.0_RKIND
      if (i_so4 > 0) particle_sulfate = particle%vol(i_so4) * aero_data%density(i_so4) / aero_data%molec_weight(i_so4)
   end function particle_sulfate


   !> Allocate the gas budget sums and record the content of each gas at the start. Does this
   !> once, at the first call.
   subroutine gas_budget_setup(domain)

      type (domain_type), intent(inout) :: domain

#ifdef MPAS_PARTMC_MOSAIC
      real (kind=RKIND), allocatable :: local(:)
      integer :: n_gas

      ! the sums may already be there from a restart file; the gas data are still needed
      call gas_setup(domain)
      if (allocated(gas_emitted)) return
      n_gas = gas_data_n_spec(gas_data)
      allocate(gas_emitted(n_gas), gas_deposited(n_gas), gas_wetdep(n_gas), gas_content0(n_gas), local(n_gas))
      gas_emitted = 0.0_RKIND
      gas_deposited = 0.0_RKIND
      gas_wetdep = 0.0_RKIND
      call gas_content(domain, local)
      call mpas_dmpar_sum_real_array(domain % dminfo, n_gas, local, gas_content0)
#endif

   end subroutine gas_budget_setup


   !> Read the in-cloud fractions of the gases for wet removal, once.
   subroutine gas_wet_setup(domain)

      type (domain_type), intent(inout) :: domain

#ifdef MPAS_PARTMC_MOSAIC
      type(spec_file_t) :: file
      type(gas_state_t) :: gas_wet

      if (allocated(gas_wet_frac)) return
      if (.not. pmc_chem .or. pmc_n_part <= 0 .or. len_trim(gas_wet_file) == 0) return
      call gas_budget_setup(domain)
      call spec_file_open(trim(gas_wet_file), file)
      call spec_file_read_gas_state(file, gas_data, gas_wet)
      call spec_file_close(file)
      allocate(gas_wet_frac(gas_data_n_spec(gas_data)))
      gas_wet_frac = min(1.0d0, max(0.0d0, gas_wet%mix_rat))
      call mpas_log_write('PartMC: wet removal of gases with the in-cloud fractions of ' // trim(gas_wet_file))
#endif

   end subroutine gas_wet_setup


   !> Removal of particles and gases in cloud by grid-scale precipitation.
   !>
   !> The microphysics gives, for each layer, the rate P at which cloud water and ice turn into
   !> precipitation (rainprod, kg/kg/s). Over one step the layer's condensate loses the fraction
   !>   f = P dt / (qc + qi + P dt)
   !> (qc and qi are the values after the microphysics, so the denominator is the condensate
   !> before the conversion). Precipitation that evaporates on the way down puts its content
   !> back in the air; that is not followed layer by layer here. Instead f is multiplied by the
   !> fraction r of the column's production that reaches the ground in this step,
   !>   r = min(1, surface precipitation / column sum of P dt),
   !> so that nothing is removed where no precipitation arrives. A tracer that is all in the
   !> cloud condensate loses the fraction r f. Particles: each one in the layer whose critical
   !> supersaturation is at most PMC_WET_SCRIT is removed with probability r f. The reference
   !> number scalars lose r f, so they agree with the particles only when all particles count
   !> as activated. Gases lose e r f, with e the fraction of the gas in cloud water from the
   !> PMC_GAS_WET file. The cell is taken as cloudy throughout wherever it has condensate,
   !> as the microphysics does. Convective precipitation and scavenging below cloud are not
   !> included.
   subroutine atm_partmc_wet_removal(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC
      type (mpas_pool_type), pointer :: mesh, state, diag, diag_physics
      type (field3DReal), pointer :: scalars_field
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, rainprod, theta_m, exner, wwAvg
      real (kind=RKIND), dimension(:), pointer :: rdzw, areaCell, rainncv
      real (kind=RKIND), pointer :: config_dt
      integer, dimension(:,:), pointer :: cellsOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, pointer :: nCells, nCellsSolve, nVertLevels, index_qc, index_qi, index_qv
      integer, pointer :: index_pmc_num, index_pmc_num2, index_ref, index_src2
      integer :: iCell, k, j, g, i_s, idx(4), n_gas
      real (kind=RKIND) :: prod, col_prod, r, p, mk, cond, w, sums(6), g_sums(6)
      real (kind=RKIND), allocatable :: amount(:)
      type(env_state_t) :: env_state
      logical :: removed
      integer :: k1, k2, n_run, i_cls(4), ic, iEdge, jCell
      real (kind=RKIND) :: amt, cn, wsum(8), g_wsum(8), csum(22), g_csum(22), amt_p, f_in, f_clean, fl
      real (kind=RKIND), dimension(:,:), pointer :: ruAvg, edgesOnCell_sign
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      integer, dimension(:,:), pointer :: edgesOnCell
      type (field2DReal), pointer :: rainprod_field
      real (kind=RKIND) :: xs(NPOS)

      if (.not. pmc_wet) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'diag_physics', diag_physics)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', config_dt)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qc', index_qc)
      call mpas_pool_get_dimension(state, 'index_qi', index_qi)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_dimension(state, 'index_pmc_num', index_pmc_num)
      call mpas_pool_get_dimension(state, 'index_pmc_num2', index_pmc_num2)
      call mpas_pool_get_dimension(state, 'index_pmc_ref', index_ref)
      call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', index_src2)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      nullify(rainprod)
      nullify(rainncv)
      call mpas_pool_get_array(diag_physics, 'rainprod', rainprod)
      call mpas_pool_get_array(diag_physics, 'rainncv', rainncv)
      if (.not. associated(rainprod) .or. .not. associated(rainncv)) then
         call mpas_log_write('PartMC: wet removal needs a microphysics scheme that gives rainprod', &
                             messageType=MPAS_LOG_CRIT)
      end if
      idx = (/ index_pmc_num, index_pmc_num2, index_ref, index_src2 /)
      ! neighbours' production is needed for the diagnostics
      call mpas_pool_get_field(diag_physics, 'rainprod', rainprod_field)
      call mpas_dmpar_exch_halo_field(rainprod_field)

      n_gas = 0
      call gas_wet_setup(domain)
      if (allocated(gas_wet_frac)) then
         n_gas = size(gas_wet_frac)
         allocate(amount(n_gas))
      end if

      if (lagr_shuffle .and. pmc_n_part > 0 .and. pmc_mode == MODE_LAGR) then
         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               do j = 1, pos(k,iCell) % n
                  call random_in_cell(iCell, k, nCells, nEdgesOnCell, cellsOnCell, xs(1:4))
                  xs(5:7) = xs(1:3)
                  pos(k,iCell) % x(:,j) = xs
               end do
            end do
         end do
      end if

      do iCell = 1, nCellsSolve
         ! precipitation produced in the column this step (kg m^-2) and the part reaching the ground
         col_prod = 0.0_RKIND
         do k = 1, nVertLevels
            col_prod = col_prod + max(0.0_RKIND, rainprod(k,iCell)) * config_dt * rho_zz(k,iCell) / rdzw(k)
         end do
         wet_sums(3) = wet_sums(3) + max(0.0_RKIND, rainncv(iCell)) * areaCell(iCell)
         wet_sums(4) = wet_sums(4) + col_prod * areaCell(iCell)
         if (col_prod <= 0.0_RKIND .or. rainncv(iCell) <= 0.0_RKIND) cycle
         r = min(1.0_RKIND, rainncv(iCell) / col_prod)
         wet_sums(5) = wet_sums(5) + r * col_prod * areaCell(iCell)

         do k = 1, nVertLevels
            prod = max(0.0_RKIND, rainprod(k,iCell)) * config_dt
            if (prod <= 0.0_RKIND) cycle
            cond = max(0.0_RKIND, scalars(index_qc,k,iCell)) + max(0.0_RKIND, scalars(index_qi,k,iCell))
            p = r * prod / (cond + prod)
            ! dry air mass per area; number per area is mk times the mixing ratio
            mk = rho_zz(k,iCell) / rdzw(k)

            amt = p * mk * scalars(index_ref,k,iCell) * areaCell(iCell)
            wet_sums(1) = wet_sums(1) + amt
            ! where this removal happens: run of adjacent layers with production, and the
            ! Courant number of the vertical mass flux through the layer
            k1 = k
            do while (k1 > 1)
               if (rainprod(k1-1,iCell) <= 0.0_RKIND) exit
               k1 = k1 - 1
            end do
            k2 = k
            do while (k2 < nVertLevels)
               if (rainprod(k2+1,iCell) <= 0.0_RKIND) exit
               k2 = k2 + 1
            end do
            n_run = k2 - k1 + 1
            cn = 0.5_RKIND * (wwAvg(k,iCell) + wwAvg(k+1,iCell)) * config_dt * rdzw(k) / rho_zz(k,iCell)
            wet_where(1) = wet_where(1) + amt * p
            if (cn > 0.0_RKIND) then
               wet_where(2) = wet_where(2) + amt * cn
               wet_where(3) = wet_where(3) + amt
            end if
            wet_where(4) = wet_where(4) + amt * real(n_run, kind=RKIND)
            if (n_run == 1) then
               i_cls(1) = 1
            else if (n_run <= 3) then
               i_cls(1) = 2
            else if (n_run <= 7) then
               i_cls(1) = 3
            else
               i_cls(1) = 4
            end if
            wet_where(4 + i_cls(1)) = wet_where(4 + i_cls(1)) + amt
            i_cls(2) = merge(5, 6, cn > 0.0_RKIND)
            i_cls(3) = merge(7, 8, k == k1)
            ! horizontal inflow of air to the cell in this layer, and the part from neighbours
            ! without production
            f_in = 0.0_RKIND
            f_clean = 0.0_RKIND
            do ic = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(ic,iCell)
               jCell = cellsOnCell(ic,iCell)
               fl = -edgesOnCell_sign(ic,iCell) * ruAvg(k,iEdge) * dvEdge(iEdge)
               if (fl <= 0.0_RKIND) cycle
               f_in = f_in + fl
               if (rainprod(k,jCell) <= 0.0_RKIND) f_clean = f_clean + fl
            end do
            if (f_clean <= 0.0_RKIND) then
               i_cls(4) = 9
            else if (f_clean <= 0.5_RKIND * f_in) then
               i_cls(4) = 10
            else
               i_cls(4) = 11
            end if
            amt_p = 0.0_RKIND
            if (pmc_n_part > 0 .and. pmc_mode == MODE_LAGR) then
               do j = 1, pos(k,iCell) % n
                  g = aero_states(k,iCell)%apa%particle(j)%weight_class
                  amt_p = amt_p + p * lagr_w * 2.0_RKIND**(g - LAGR_G0)
               end do
            end if
            do ic = 1, 4
               wet_cls(1,i_cls(ic)) = wet_cls(1,i_cls(ic)) + amt
               wet_cls(2,i_cls(ic)) = wet_cls(2,i_cls(ic)) + amt_p
            end do
            if (pmc_wet_diag_only) then
               wet_sums(6) = wet_sums(6) + amt_p
               cycle
            end if
            do i_s = 1, 4
               scalars(idx(i_s),k,iCell) = scalars(idx(i_s),k,iCell) * (1.0_RKIND - p)
            end do
#ifdef MPAS_PARTMC_MOSAIC
            if (n_gas > 0) then
               amount = gas_wet_frac * p * scalars(gas_idx,k,iCell)
               scalars(gas_idx,k,iCell) = scalars(gas_idx,k,iCell) - amount
               gas_wetdep = gas_wetdep + amount * 1.0e-9_RKIND * mk / M_AIR * areaCell(iCell)
            end if
#endif
            if (pmc_n_part <= 0 .or. pmc_mode /= MODE_LAGR) cycle
            removed = .false.
            if (pmc_wet_scrit < 1.0e8_RKIND) then
               env_state%temp = theta_m(k,iCell) / (1.0_RKIND + rvord*max(0.0_RKIND, scalars(index_qv,k,iCell))) &
                                * exner(k,iCell)
            end if
            do j = pos(k,iCell) % n, 1, -1
               g = aero_states(k,iCell)%apa%particle(j)%weight_class
               wet_sums(6) = wet_sums(6) + p * lagr_w * 2.0_RKIND**(g - LAGR_G0)
               if (pmc_random() >= p) cycle
               if (pmc_wet_scrit < 1.0e8_RKIND) then
                  if (aero_particle_crit_rel_humid(aero_states(k,iCell)%apa%particle(j), aero_data, env_state) &
                      - 1.0d0 > pmc_wet_scrit) cycle
               end if
               g = aero_states(k,iCell)%apa%particle(j)%weight_class
               w = lagr_w * 2.0_RKIND**(g - LAGR_G0)
               wet_sums(2) = wet_sums(2) + w
               part_s_removed = part_s_removed + particle_sulfate(aero_states(k,iCell)%apa%particle(j)) * w
               call aero_state_remove_particle_no_info(aero_states(k,iCell), j)
               call pos_remove(pos(k,iCell), j)
               removed = .true.
            end do
            if (removed .and. lagr_control) call lagr_rebalance(k, iCell, nCells, nEdgesOnCell, cellsOnCell)
         end do
      end do

      call mpas_pool_get_field(state, 'scalars', scalars_field, 2)
      call mpas_dmpar_exch_halo_field(scalars_field)

      if (mod(itimestep, 12) == 0 .or. itimestep <= 3) then
         sums = wet_sums
         call mpas_dmpar_sum_real_array(domain % dminfo, 6, sums, g_sums)
         call mpas_log_write('PartMC wet removal step $i: number removed so far: reference = $r, particles = $r; ' &
                             // 'so far (kg): surface precipitation = $r, column production = $r, ' &
                             // 'production counted as reaching the surface = $r; particles expected ' &
                             // 'from the layers'' content (all activated) = $r', &
                             intArgs=(/itimestep/), realArgs=g_sums)
         wsum = wet_where
         call mpas_dmpar_sum_real_array(domain % dminfo, 8, wsum, g_wsum)
         if (g_sums(1) > 0.0_RKIND) then
            call mpas_log_write('PartMC wet removal step $i: weighted by number removed: probability = $r, ' &
                                // 'Courant number of rising air = $r, share in rising air = $r, adjacent ' &
                                // 'removing layers = $r, share in runs of 1, 2-3, 4-7, 8+ layers = $r $r $r $r', &
                                intArgs=(/itimestep/), &
                                realArgs=(/g_wsum(1)/g_sums(1), g_wsum(2)/max(g_wsum(3), 1.0e-30_RKIND), &
                                           g_wsum(3)/g_sums(1), g_wsum(4)/g_sums(1), g_wsum(5:8)/g_sums(1)/))
            csum = reshape(wet_cls, (/22/))
            call mpas_dmpar_sum_real_array(domain % dminfo, 22, csum, g_csum)
            call mpas_log_write('PartMC wet removal step $i: particles expected / reference by class: runs of 1, ' &
                                // '2-3, 4-7, 8+ layers = $r $r $r $r; rising, sinking air = $r $r; lowest layer of ' &
                                // 'its run, other layers = $r $r', intArgs=(/itimestep/), &
                                realArgs=g_csum(2:16:2) / max(g_csum(1:15:2), 1.0e-30_RKIND))
            call mpas_log_write('PartMC wet removal step $i: share of the reference removal: rising air = $r, ' &
                                // 'lowest layer of its run = $r', intArgs=(/itimestep/), &
                                realArgs=(/g_csum(9)/g_sums(1), g_csum(13)/g_sums(1)/))
            call mpas_log_write('PartMC wet removal step $i: by share of horizontal inflow from neighbours that are ' &
                                // 'not removing (none, up to half, more than half): particles expected / reference ' &
                                // '= $r $r $r; share of the reference removal = $r $r $r', intArgs=(/itimestep/), &
                                realArgs=(/g_csum(18:22:2) / max(g_csum(17:21:2), 1.0e-30_RKIND), &
                                           g_csum(17:21:2) / g_sums(1)/))
         end if
      end if
#endif

   end subroutine atm_partmc_wet_removal


   !> Gas chemistry and gas-particle partitioning with MOSAIC (CBM-Z, ASTEM, MESA) in every
   !> owned cell, through PartMC's interface.
   !>
   !> The gases are model scalars gas_<name> (nmol per mol of dry air), so they are advected by
   !> the model and mixed and convected with the particle number scalars. They start from one
   !> PartMC gas_init file everywhere; atm_partmc_gas_surface does emission and deposition. MOSAIC
   !> works with mixing ratios relative to moist air, so values are converted on the way in and
   !> out. Temperature, pressure and relative humidity come from the model state; latitude,
   !> longitude, altitude and time are set for MOSAIC's solar zenith angle.
   subroutine atm_partmc_chemistry(domain, itimestep)

#ifdef MPAS_PARTMC_MOSAIC
      use module_data_mosaic_main, only : rlon, rlat, zalt_m
#endif

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC_MOSAIC
      type (mpas_pool_type), pointer :: mesh, state, diag
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: theta_m, exner, pressure_p, pressure_base, zgrid, zz
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, areaCell, rdzw
      real (kind=RKIND), pointer :: config_dt
      character (len=StrKIND), pointer :: config_start_time
      integer, pointer :: nCellsSolve, nVertLevels, index_qv
      integer :: iCell, k, j, n_before, n_cells, g_n_cells, n_gas, yr, mo, dy, hr, mi, sc, i_so2, i_h2so4, i_so4
      integer, parameter :: days_before(12) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
      type(spec_file_t) :: file
      type(env_state_t) :: env_state
      type(gas_state_t) :: gas_state
      character(len=PMC_UUID_LEN) :: uuid
      real (kind=RKIND) :: qv, vol, sums(4), g_sums(4), t0, t1, t_local, t_max, moist_fac, ch4_min, ch4_max
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz
      type (field3DReal), pointer :: scalars_field
      integer, pointer :: i_ptr
      real (kind=dp) :: del_t, e_s, doy, decl, hour_angle, cos_ref, rh_in
      real (kind=dp), allocatable :: gas_a(:), gas_pre(:)
      type(aero_state_t) :: state_a, state_b, state_pre
      type(env_state_t) :: env_b
      real (kind=RKIND) :: test_sums(3), g_test(3)
      integer :: kk, jc, n_test, g_n_test, i_box
      integer, dimension(:), pointer :: indexToCellID
      character(len=64) :: box_prefix
      type(aero_data_t) :: aero_data_box

      if (.not. pmc_chem .or. pmc_n_part <= 0) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', config_dt)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_start_time', config_start_time)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_base)

      del_t = config_dt
      uuid = 'mpas-partmc'
      read(config_start_time, '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') yr, mo, dy, hr, mi, sc
      env_state%start_day = days_before(mo) + dy
      if (mo > 2 .and. mod(yr, 4) == 0) env_state%start_day = env_state%start_day + 1
      env_state%start_time = real(hr*3600 + mi*60 + sc, kind=dp)
      env_state%elapsed_time = real(itimestep - 1, kind=dp) * del_t

      call gas_setup(domain)
      n_gas = gas_data_n_spec(gas_data)
      call gas_state_set_size(gas_state, n_gas)
      i_so2 = gas_data_spec_by_name(gas_data, 'SO2')
      i_h2so4 = gas_data_spec_by_name(gas_data, 'H2SO4')
      i_so4 = aero_data_spec_by_name(aero_data, 'SO4')

      n_cells = 0
      sums = 0.0_RKIND
      test_sums = 0.0_RKIND
      n_test = 0
      allocate(gas_a(n_gas), gas_pre(n_gas))
      call cpu_time(t0)
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            vol = cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))
            qv = max(0.0_RKIND, scalars(index_qv,k,iCell))
            env_state%temp = theta_m(k,iCell) / (1.0_RKIND + rvord*qv) * exner(k,iCell)
            env_state%pressure = pressure_p(k,iCell) + pressure_base(k,iCell)
            e_s = env_state_sat_vapor_pressure(env_state)
            env_state%rel_humid = min(0.99d0, max(1.0d-4, qv * env_state%pressure / (0.622d0 + qv) / e_s))
            env_state%latitude = latCell(iCell) * 180.0d0 / pii
            env_state%longitude = lonCell(iCell) * 180.0d0 / pii
            env_state%altitude = 0.5d0*(zgrid(k,iCell) + zgrid(k+1,iCell))
            ! moles of moist air per mole of dry air
            moist_fac = 1.0_RKIND + rvord*qv

            if (env_state%altitude < pmc_chem_top) then
               rlat = latCell(iCell)
               rlon = lonCell(iCell)
               zalt_m = env_state%altitude
               gas_state%mix_rat = scalars(gas_idx,k,iCell) / moist_fac
               n_before = aero_state_n_part(aero_states(k,iCell))
               if (pmc_chem_selftest .and. mod(iCell, 97) == 0 .and. mod(k, 9) == 1) then
                  state_pre = aero_states(k,iCell)
                  gas_pre = gas_state%mix_rat
                  rh_in = env_state%rel_humid
               end if
               i_box = 0
               if (itimestep <= box_steps) then
                  do j = 1, n_box
                     if (indexToCellID(iCell) == box_id(j) .and. k == box_k(j)) i_box = j
                  end do
               end if
               if (i_box > 0) then
                  write(box_prefix, '(a,i0,a,i0)') 'box/c', box_id(i_box), '_k', k
                  ! PartMC's file format wants a name for each weight class
                  aero_data_box = aero_data
                  if (allocated(aero_data_box%weight_class_name)) deallocate(aero_data_box%weight_class_name)
                  allocate(aero_data_box%weight_class_name(size(aero_states(k,iCell)%awa%weight, 2)))
                  do j = 1, size(aero_data_box%weight_class_name)
                     write(aero_data_box%weight_class_name(j), '(a,i0)') 'w', j
                  end do
                  call output_state(trim(box_prefix) // '_pre', OUTPUT_TYPE_DIST, aero_data_box, aero_states(k,iCell), &
                                    gas_data, gas_state, env_state, itimestep, env_state%elapsed_time, del_t, 1, &
                                    .false., .false., uuid)
               end if
               call mosaic_timestep(env_state, aero_data, aero_states(k,iCell), gas_data, gas_state, .false., uuid)
               if (i_box > 0) then
                  call output_state(trim(box_prefix) // '_post', OUTPUT_TYPE_DIST, aero_data_box, aero_states(k,iCell), &
                                    gas_data, gas_state, env_state, itimestep, env_state%elapsed_time, del_t, 1, &
                                    .false., .false., uuid)
               end if
               if (aero_state_n_part(aero_states(k,iCell)) /= n_before) then
                  call mpas_log_write('PartMC: MOSAIC changed the number of particles in a cell', &
                                      messageType=MPAS_LOG_CRIT)
               end if
               scalars(gas_idx,k,iCell) = gas_state%mix_rat * moist_fac
               n_cells = n_cells + 1

               if (pmc_chem_selftest .and. mod(iCell, 97) == 0 .and. mod(k, 9) == 1) then
                  ! Repeat this cell's step from the saved state after MOSAIC has worked on a
                  ! different cell: MOSAIC keeps its working state in module variables, so this
                  ! checks that nothing carries over from one cell to the next.
                  gas_a = gas_state%mix_rat
                  state_a = aero_states(k,iCell)
                  env_b = env_state
                  ! solar zenith angle against a direct calculation
                  doy = real(env_state%start_day, kind=dp) + (env_state%start_time + env_state%elapsed_time &
                        + 0.5d0*del_t) / 86400.0d0
                  decl = -23.44d0 * pii/180.0d0 * cos(2.0d0*pii*(doy + 9.0d0)/365.0d0)
                  hour_angle = 2.0d0*pii*(doy - floor(doy)) - pii + lonCell(iCell)
                  cos_ref = sin(latCell(iCell))*sin(decl) + cos(latCell(iCell))*cos(decl)*cos(hour_angle)
                  test_sums(3) = max(test_sums(3), abs(cos(env_state%solar_zenith_angle) - cos_ref))
                  ! a far-away cell, then this one again
                  kk = nVertLevels + 1 - k
                  jc = nCellsSolve + 1 - iCell
                  env_b%temp = theta_m(kk,jc) / (1.0_RKIND + rvord*max(0.0_RKIND, scalars(index_qv,kk,jc))) * exner(kk,jc)
                  env_b%pressure = pressure_p(kk,jc) + pressure_base(kk,jc)
                  env_b%rel_humid = 0.3d0
                  rlat = latCell(jc)
                  rlon = lonCell(jc)
                  state_b = aero_states(kk,jc)
                  gas_state%mix_rat = scalars(gas_idx,kk,jc)
                  call mosaic_timestep(env_b, aero_data, state_b, gas_data, gas_state, .false., uuid)
                  env_b = env_state
                  env_b%rel_humid = rh_in
                  rlat = latCell(iCell)
                  rlon = lonCell(iCell)
                  state_b = state_pre
                  gas_state%mix_rat = gas_pre
                  call mosaic_timestep(env_b, aero_data, state_b, gas_data, gas_state, .false., uuid)
                  do j = 1, n_gas
                     if (gas_a(j) > 1.0d-12) test_sums(1) = max(test_sums(1), abs(gas_state%mix_rat(j)/gas_a(j) - 1.0d0))
                  end do
                  do j = 1, aero_state_n_part(state_a)
                     test_sums(2) = max(test_sums(2), abs(sum(state_b%apa%particle(j)%vol) &
                                                          / sum(state_a%apa%particle(j)%vol) - 1.0d0))
                  end do
                  n_test = n_test + 1
               end if
            end if

            if (mod(itimestep, 12) == 0 .or. itimestep == 1) then
               ! sulfur (mol): gas, then particles
               sums(1) = sums(1) + (scalars(gas_idx(i_so2),k,iCell) + scalars(gas_idx(i_h2so4),k,iCell)) * 1.0d-9 &
                                   * rho_zz(k,iCell) * zz(k,iCell) * vol / M_AIR
               do j = 1, aero_state_n_part(aero_states(k,iCell))
                  sums(2) = sums(2) + aero_states(k,iCell)%apa%particle(j)%vol(i_so4) &
                                      * aero_data%density(i_so4) / aero_data%molec_weight(i_so4) &
                                      * aero_weight_array_num_conc(aero_states(k,iCell)%awa, &
                                           aero_states(k,iCell)%apa%particle(j), aero_data) * vol
                  sums(3) = sums(3) + aero_states(k,iCell)%apa%particle(j)%vol(aero_data%i_water) &
                                      * aero_data%density(aero_data%i_water) &
                                      * aero_weight_array_num_conc(aero_states(k,iCell)%awa, &
                                           aero_states(k,iCell)%apa%particle(j), aero_data) * vol
               end do
            end if
         end do
      end do
      call cpu_time(t1)
      call mpas_pool_get_field(state, 'scalars', scalars_field, 2)
      call mpas_dmpar_exch_halo_field(scalars_field)

      if (pmc_chem_selftest) then
         call mpas_dmpar_max_real_array(domain % dminfo, 3, test_sums, g_test)
         call mpas_dmpar_sum_int(domain % dminfo, n_test, g_n_test)
         call mpas_log_write('PartMC chemistry self-test step $i: $i cells repeated after another cell; largest ' &
                             // 'rel diff gas = $r, particle volume = $r; largest cos(zenith) diff from direct ' &
                             // 'calculation = $r', intArgs=(/itimestep, g_n_test/), realArgs=g_test)
      end if
      if (mod(itimestep, 12) == 0 .or. itimestep == 1) then
         t_local = t1 - t0
         ! net sulfur (mol) put in at the surface so far
         if (allocated(gas_emitted)) then
            sums(4) = gas_emitted(i_so2) + gas_emitted(i_h2so4) - gas_deposited(i_so2) - gas_deposited(i_h2so4) &
                      - gas_wetdep(i_so2) - gas_wetdep(i_h2so4)
         end if
         sums(4) = sums(4) - part_s_removed
         call mpas_dmpar_sum_real_array(domain % dminfo, 4, sums, g_sums)
         call mpas_dmpar_sum_int(domain % dminfo, n_cells, g_n_cells)
         call mpas_dmpar_max_real(domain % dminfo, t_local, t_max)
         if (chem_s0 <= 0.0_RKIND) chem_s0 = g_sums(1) + g_sums(2) - g_sums(4)
         ! range of a long-lived gas, which transport alone must leave uniform
         j = gas_idx(gas_data_spec_by_name(gas_data, 'CH4'))
         t_local = minval(scalars(j,:,1:nCellsSolve))
         call mpas_dmpar_min_real(domain % dminfo, t_local, ch4_min)
         t_local = maxval(scalars(j,:,1:nCellsSolve))
         call mpas_dmpar_max_real(domain % dminfo, t_local, ch4_max)
         call mpas_log_write('PartMC chemistry step $i: CH4 tracer min = $r, max = $r ppb', &
                             intArgs=(/itimestep/), realArgs=(/ch4_min, ch4_max/))
         t_local = part_s_thin
         call mpas_dmpar_sum_real(domain % dminfo, t_local, ch4_min)
         call mpas_log_write('PartMC chemistry step $i: sulfate added by the thinning of the population control ' &
                             // 'since the start of this run = $r mol, / total at start = $r; rel change of the ' &
                             // 'total less this = $r', intArgs=(/itimestep/), &
                             realArgs=(/ch4_min, ch4_min/chem_s0, &
                                        (g_sums(1) + g_sums(2) - g_sums(4) - ch4_min)/chem_s0 - 1.0_RKIND/))
         call mpas_log_write('PartMC chemistry step $i: cells = $i, sulfur (mol) gas = $r, particles = $r, ' &
                             // 'total less net surface input, rel change = $r, particle water (kg) = $r, ' &
                             // 'slowest task = $r s', &
                             intArgs=(/itimestep, g_n_cells/), &
                             realArgs=(/g_sums(1), g_sums(2), (g_sums(1) + g_sums(2) - g_sums(4))/chem_s0 - 1.0_RKIND, &
                                        g_sums(3), t_max/))
      end if
#endif

   end subroutine atm_partmc_chemistry


   !> PartMC Brownian coagulation in every owned cell for one model time step.
   !>
   !> In Lagrangian mode the weight classes present in the cell are renumbered 1..m for the
   !> call, so that PartMC does not loop over all the unused weight levels, and positions are
   !> carried through by particle id. A particle with an id that was not there before (the
   !> product of a coagulation) is placed at a random point of the cell.
   subroutine atm_partmc_coagulation(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC
      type (mpas_pool_type), pointer :: mesh, state, diag
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: theta_m, exner, pressure_p, pressure_base
      real (kind=RKIND), pointer :: config_dt
      integer, dimension(:,:), pointer :: cellsOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, pointer :: nCells, nCellsSolve, nVertLevels, index_qv
      integer :: iCell, k, j, jj, n0, n1, m, c, n_samp, n_coag, g_n_coag, n_new, g_n_new
      integer :: map(LAGR_GROUPS), present(LAGR_GROUPS)
      integer (kind=I8KIND), allocatable :: ids(:)
      real (kind=RKIND), allocatable :: xs(:,:)
      real (kind=RKIND) :: xp(NPOS)
      type(env_state_t) :: env_state
      type(aero_weight_array_t) :: awa_save
      type(aero_sorted_t) :: empty_sorted
      type(bin_grid_t) :: coag_grid
      integer, parameter :: n_bin = 36
      real (kind=dp), allocatable :: k_min_lev(:,:,:), k_max_lev(:,:,:), k_min(:,:), k_max(:,:)
      real (kind=RKIND) :: t_lo, t_hi, p_lo, p_hi, temp, pres
      integer :: i_env
      real (kind=dp) :: del_t

      if (.not. pmc_coag .or. pmc_n_part <= 0) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', config_dt)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_base)

      del_t = config_dt
      env_state%rel_humid = 0.0d0
      n_new = 0

      ! One size grid for all cells, and bounds on the kernel for each level that hold for the
      ! range of temperature and pressure on that level. PartMC otherwise rebuilds both for
      ! every cell, which costs far more than the coagulation itself.
      call bin_grid_make(coag_grid, BIN_GRID_TYPE_LOG, n_bin, 1.0d-10, 1.0d-4)
      allocate(k_min_lev(n_bin,n_bin,nVertLevels), k_max_lev(n_bin,n_bin,nVertLevels), &
               k_min(n_bin,n_bin), k_max(n_bin,n_bin))
      do k = 1, nVertLevels
         t_lo = huge(1.0_RKIND)
         t_hi = 0.0_RKIND
         p_lo = huge(1.0_RKIND)
         p_hi = 0.0_RKIND
         do iCell = 1, nCellsSolve
            temp = theta_m(k,iCell) / (1.0_RKIND + rvord*scalars(index_qv,k,iCell)) * exner(k,iCell)
            pres = pressure_p(k,iCell) + pressure_base(k,iCell)
            t_lo = min(t_lo, temp)
            t_hi = max(t_hi, temp)
            p_lo = min(p_lo, pres)
            p_hi = max(p_hi, pres)
         end do
         k_min_lev(:,:,k) = huge(1.0d0)
         k_max_lev(:,:,k) = 0.0d0
         if (nCellsSolve == 0) cycle
         do i_env = 1, 4
            env_state%temp = merge(t_lo, t_hi, i_env <= 2)
            env_state%pressure = merge(p_lo, p_hi, mod(i_env, 2) == 1)
            call est_k_minmax_binned_unweighted(coag_grid, COAG_KERNEL_TYPE_BROWN, aero_data, env_state, k_min, k_max)
            k_min_lev(:,:,k) = min(k_min_lev(:,:,k), k_min)
            k_max_lev(:,:,k) = max(k_max_lev(:,:,k), k_max)
         end do
         k_min_lev(:,:,k) = k_min_lev(:,:,k) / 1.05d0
         k_max_lev(:,:,k) = k_max_lev(:,:,k) * 1.05d0
      end do
      allocate(ids(4*pmc_n_part + 16), xs(NPOS, 4*pmc_n_part + 16))

      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            n0 = aero_state_n_part(aero_states(k,iCell))
            if (n0 < 2) cycle
            env_state%temp = theta_m(k,iCell) / (1.0_RKIND + rvord*scalars(index_qv,k,iCell)) * exner(k,iCell)
            env_state%pressure = pressure_p(k,iCell) + pressure_base(k,iCell)

            if (pmc_mode == MODE_LAGR) then
               if (n0 > size(ids)) then
                  deallocate(ids, xs)
                  allocate(ids(2*n0), xs(NPOS, 2*n0))
               end if
               map = 0
               m = 0
               do j = 1, n0
                  ids(j) = aero_states(k,iCell)%apa%particle(j)%id
                  xs(:,j) = pos(k,iCell) % x(:,j)
                  c = aero_states(k,iCell)%apa%particle(j)%weight_class
                  if (map(c) == 0) then
                     m = m + 1
                     map(c) = m
                     present(m) = c
                  end if
                  aero_states(k,iCell)%apa%particle(j)%weight_class = map(c)
               end do
               awa_save = aero_states(k,iCell)%awa
               call aero_weight_array_set_sizes(aero_states(k,iCell)%awa, 1, m)
               aero_states(k,iCell)%awa%weight%type = AERO_WEIGHT_TYPE_NONE
               aero_states(k,iCell)%awa%weight%exponent = 0.0d0
               do j = 1, m
                  aero_states(k,iCell)%awa%weight(1,j)%magnitude = awa_save%weight(1,present(j))%magnitude
               end do
            end if

            call aero_state_sort(aero_states(k,iCell), aero_data, bin_grid=coag_grid)
            aero_states(k,iCell)%aero_sorted%coag_kernel_min = k_min_lev(:,:,k)
            aero_states(k,iCell)%aero_sorted%coag_kernel_max = k_max_lev(:,:,k)
            aero_states(k,iCell)%aero_sorted%coag_kernel_bounds_valid = .true.
            call mc_coag(COAG_KERNEL_TYPE_BROWN, env_state, aero_data, aero_states(k,iCell), del_t, n_samp, n_coag)
            n_coag_total = n_coag_total + n_coag

            ! drop the sorting structures and removal records PartMC built, to save memory and so
            ! that particles can again be removed by index
            aero_states(k,iCell)%aero_sorted = empty_sorted
            aero_states(k,iCell)%valid_sort = .false.
            call aero_info_array_zero(aero_states(k,iCell)%aero_info_array)

            n1 = aero_state_n_part(aero_states(k,iCell))
            if (pmc_mode == MODE_LAGR) then
               aero_states(k,iCell)%awa = awa_save
               pos(k,iCell) % n = 0
               do j = 1, n1
                  aero_states(k,iCell)%apa%particle(j)%weight_class = &
                       present(aero_states(k,iCell)%apa%particle(j)%weight_class)
                  xp(1) = huge(1.0_RKIND)
                  ! most particles are untouched and keep their index
                  if (j <= n0) then
                     if (ids(j) == aero_states(k,iCell)%apa%particle(j)%id) xp = xs(:,j)
                  end if
                  if (xp(1) > 2.0_RKIND) then
                     do jj = 1, n0
                        if (ids(jj) == aero_states(k,iCell)%apa%particle(j)%id) then
                           xp = xs(:,jj)
                           exit
                        end if
                     end do
                  end if
                  if (xp(1) > 2.0_RKIND) then
                     call random_in_cell(iCell, k, nCells, nEdgesOnCell, cellsOnCell, xp(1:4))
                     xp(5:7) = xp(1:3)
                     n_new = n_new + 1
                  end if
                  call pos_append(pos(k,iCell), xp)
               end do
               if (lagr_control) call lagr_rebalance(k, iCell, nCells, nEdgesOnCell, cellsOnCell)
            else
               cnt(:,k,iCell) = 0
               do j = 1, n1
                  c = aero_states(k,iCell)%apa%particle(j)%weight_class
                  cnt(c,k,iCell) = cnt(c,k,iCell) + 1
               end do
            end if
         end do
      end do

      if (mod(itimestep, 12) == 0) then
         call mpas_dmpar_sum_int(domain % dminfo, n_coag_total, g_n_coag)
         call mpas_dmpar_sum_int(domain % dminfo, n_new, g_n_new)
         call mpas_log_write('PartMC coagulation step $i: events so far = $i, particles placed at random this step = $i', &
                             intArgs=(/itimestep, g_n_coag, g_n_new/))
      end if
#endif

   end subroutine atm_partmc_coagulation


   !> Solve a tridiagonal system (Thomas algorithm).
   pure subroutine tridiag_solve(n, lower, diag_a, upper, rhs, sol)
      integer, intent(in) :: n
      real (kind=RKIND), intent(in) :: lower(n), diag_a(n), upper(n), rhs(n)
      real (kind=RKIND), intent(out) :: sol(n)
      real (kind=RKIND) :: cp(n), dp_(n), denom
      integer :: i

      cp(1) = upper(1) / diag_a(1)
      dp_(1) = rhs(1) / diag_a(1)
      do i = 2, n
         denom = diag_a(i) - lower(i)*cp(i-1)
         cp(i) = upper(i) / denom
         dp_(i) = (rhs(i) - lower(i)*dp_(i-1)) / denom
      end do
      sol(n) = dp_(n)
      do i = n-1, 1, -1
         sol(i) = dp_(i) - cp(i)*sol(i+1)
      end do
   end subroutine tridiag_solve


#ifdef MPAS_PARTMC
   !> Global mass (kg) of each species, summed over all particles.
   subroutine species_mass_totals(domain, totals)

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(out) :: totals(:)

      type (mpas_pool_type), pointer :: mesh
      real (kind=RKIND), dimension(:,:), pointer :: zz
      real (kind=RKIND), dimension(:), pointer :: areaCell, rdzw
      integer, pointer :: nCellsSolve, nVertLevels
      real (kind=RKIND), allocatable :: local(:), vol_sum(:)
      integer :: iCell, k, i_part, c

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'zz', zz)

      allocate(local(n_spec), vol_sum(n_spec))
      local = 0.0_RKIND
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            vol_sum = 0.0_RKIND
            do i_part = 1, aero_state_n_part(aero_states(k,iCell))
               c = aero_states(k,iCell)%apa%particle(i_part)%weight_class
               vol_sum = vol_sum + aero_states(k,iCell)%apa%particle(i_part)%vol &
                                   * aero_states(k,iCell)%awa%weight( &
                                        aero_states(k,iCell)%apa%particle(i_part)%weight_group,c)%magnitude
            end do
            local = local + vol_sum * aero_data%density * cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))
         end do
      end do
      call mpas_dmpar_sum_real_array(domain % dminfo, n_spec, local, totals)

   end subroutine species_mass_totals


   !> Pseudo-random number in [0,1) that depends only on the cell, level and step, so that
   !> the cells chosen for output do not depend on the MPI decomposition.
   pure real (kind=RKIND) function cell_hash(cell_id, k, itimestep)
      integer, intent(in) :: cell_id, k, itimestep
      integer (kind=I8KIND), parameter :: m = 2147483647_I8KIND
      integer (kind=I8KIND) :: x

      x = mod(int(cell_id, kind=I8KIND)*2654435761_I8KIND + int(k, kind=I8KIND)*40503_I8KIND &
              + int(itimestep, kind=I8KIND)*104729_I8KIND, m)
      x = mod(x*48271_I8KIND + 11_I8KIND, m)
      x = mod(x*48271_I8KIND + 11_I8KIND, m)
      x = mod(x*48271_I8KIND + 11_I8KIND, m)
      cell_hash = real(x, kind=RKIND) / real(m, kind=RKIND)
   end function cell_hash


   subroutine nc_check(status, what)
      integer, intent(in) :: status
      character(len=*), intent(in) :: what
      if (status /= NF90_NOERR) then
         call mpas_log_write('PartMC output: ' // what // ': ' // trim(nf90_strerror(status)), messageType=MPAS_LOG_CRIT)
      end if
   end subroutine nc_check
#endif


   !> Time-mean aerosol fields of every cell, for maps.
   !>
   !> Every PMC_DIAG_EVERY steps the particles of each cell are added to sums; every
   !> PMC_DIAG_INTERVAL steps the means are written and the sums start again. Variables (per
   !> volume of air): number of all particles and of submicron ones (dry diameter below 1 um);
   !> number that activates at 0.1%, 0.3% and 1.0% supersaturation; dry mass of six species
   !> groups for all sizes and for submicron particles; and, for submicron particles, the
   !> mass-weighted sums of the per-particle mixing entropy for three choices of "species":
   !> the six groups (chi_a), black carbon against the rest (chi_o), and non-hygroscopic
   !> (BC, dust, primary organic) against hygroscopic (sea salt, secondary organic, secondary
   !> inorganic) (chi_h), as in Zheng et al. (2021). The mixing state index of the time and
   !> space aggregate follows from these sums: D_alpha = exp(sum(m H) / sum(m)), D_gamma from
   !> the group masses, chi = (D_alpha - 1) / (D_gamma - 1). Groups: BC; dust (OIN, Ca, CO3);
   !> sea salt (Na, Cl); primary organic (OC); secondary organic (ARO1, ARO2, ALK1, OLE1, API1,
   !> API2, LIM1, LIM2); secondary inorganic (SO4, NO3, NH4, MSA). Water is left out.
   subroutine atm_partmc_diag(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC
      type (mpas_pool_type), pointer :: mesh, state, diag
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: theta_m, exner, zgrid
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      integer, dimension(:), pointer :: indexToCellID
      integer, pointer :: nCellsSolve, nVertLevels, index_qv
      integer :: iCell, k, j, i_s, g, ncid, d_cell, d_lev, d_var, varid
      real (kind=dp) :: nc, d_dry, s_crit, mg(6), m_tot, h_a, h_o, h_h, p, m_nh, m_hy
      type(env_state_t) :: env_state
      real (kind=R4KIND), allocatable :: out(:,:,:), zm(:,:)
      character(len=512) :: fname
      character(len=*), parameter :: var_names = 'num,num_sub,ccn_0p1,ccn_0p3,ccn_1p0,' &
         // 'm_bc,m_dust,m_seasalt,m_pom,m_soa,m_sia,' &
         // 'msub_bc,msub_dust,msub_seasalt,msub_pom,msub_soa,msub_sia,' &
         // 'mh_a,mh_o,mh_h,m_water'

      if (pmc_n_part <= 0 .or. pmc_diag_interval <= 0) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

      if (.not. allocated(diag_acc)) then
         allocate(diag_acc(N_DIAG, nVertLevels, nCellsSolve), diag_group(n_spec))
         diag_acc = 0.0_RKIND
         diag_group = 0
         do i_s = 1, n_spec
            select case (trim(aero_data%name(i_s)))
            case ('BC')
               diag_group(i_s) = 1
            case ('OIN', 'Ca', 'CO3')
               diag_group(i_s) = 2
            case ('Na', 'Cl')
               diag_group(i_s) = 3
            case ('OC')
               diag_group(i_s) = 4
            case ('ARO1', 'ARO2', 'ALK1', 'OLE1', 'API1', 'API2', 'LIM1', 'LIM2')
               diag_group(i_s) = 5
            case ('SO4', 'NO3', 'NH4', 'MSA')
               diag_group(i_s) = 6
            end select
         end do
         call execute_command_line('mkdir -p ' // trim(pmc_diag_dir))
      end if

      if (mod(itimestep, pmc_diag_every) == 0) then
         call mpas_pool_get_dimension(state, 'index_qv', index_qv)
         call mpas_pool_get_array(state, 'scalars', scalars, 2)
         call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
         call mpas_pool_get_array(diag, 'exner', exner)
         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               env_state%temp = theta_m(k,iCell) / (1.0_RKIND + rvord*max(0.0_RKIND, scalars(index_qv,k,iCell))) &
                                * exner(k,iCell)
               do j = 1, aero_state_n_part(aero_states(k,iCell))
                  nc = aero_weight_array_num_conc(aero_states(k,iCell)%awa, aero_states(k,iCell)%apa%particle(j), &
                                                  aero_data)
                  mg = 0.0d0
                  do i_s = 1, n_spec
                     g = diag_group(i_s)
                     if (g > 0) mg(g) = mg(g) + aero_states(k,iCell)%apa%particle(j)%vol(i_s) * aero_data%density(i_s)
                  end do
                  m_tot = sum(mg)
                  d_dry = aero_particle_dry_diameter(aero_states(k,iCell)%apa%particle(j), aero_data)
                  s_crit = aero_particle_crit_rel_humid(aero_states(k,iCell)%apa%particle(j), aero_data, env_state) &
                           - 1.0d0
                  diag_acc(1,k,iCell) = diag_acc(1,k,iCell) + nc
                  if (s_crit <= 0.001d0) diag_acc(3,k,iCell) = diag_acc(3,k,iCell) + nc
                  if (s_crit <= 0.003d0) diag_acc(4,k,iCell) = diag_acc(4,k,iCell) + nc
                  if (s_crit <= 0.010d0) diag_acc(5,k,iCell) = diag_acc(5,k,iCell) + nc
                  diag_acc(6:11,k,iCell) = diag_acc(6:11,k,iCell) + nc * mg
                  diag_acc(21,k,iCell) = diag_acc(21,k,iCell) + nc &
                     * aero_states(k,iCell)%apa%particle(j)%vol(aero_data%i_water) * aero_data%density(aero_data%i_water)
                  if (d_dry >= 1.0d-6 .or. m_tot <= 0.0d0) cycle
                  diag_acc(2,k,iCell) = diag_acc(2,k,iCell) + nc
                  diag_acc(12:17,k,iCell) = diag_acc(12:17,k,iCell) + nc * mg
                  ! per-particle mixing entropies
                  h_a = 0.0d0
                  do g = 1, 6
                     p = mg(g) / m_tot
                     if (p > 0.0d0) h_a = h_a - p * log(p)
                  end do
                  h_o = 0.0d0
                  p = mg(1) / m_tot
                  if (p > 0.0d0 .and. p < 1.0d0) h_o = -p * log(p) - (1.0d0 - p) * log(1.0d0 - p)
                  h_h = 0.0d0
                  m_nh = mg(1) + mg(2) + mg(4)
                  m_hy = mg(3) + mg(5) + mg(6)
                  p = m_nh / m_tot
                  if (m_nh > 0.0d0 .and. m_hy > 0.0d0) h_h = -p * log(p) - (1.0d0 - p) * log(1.0d0 - p)
                  diag_acc(18,k,iCell) = diag_acc(18,k,iCell) + nc * m_tot * h_a
                  diag_acc(19,k,iCell) = diag_acc(19,k,iCell) + nc * m_tot * h_o
                  diag_acc(20,k,iCell) = diag_acc(20,k,iCell) + nc * m_tot * h_h
               end do
            end do
         end do
         n_diag_samples = n_diag_samples + 1
      end if

      if (mod(itimestep, pmc_diag_interval) /= 0 .or. n_diag_samples == 0) return

      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)
      allocate(out(N_DIAG, nVertLevels, nCellsSolve), zm(nVertLevels, nCellsSolve))
      out = real(diag_acc / real(n_diag_samples, kind=RKIND), kind=R4KIND)
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            zm(k,iCell) = real(0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell)), kind=R4KIND)
         end do
      end do
      write(fname, '(a,"/diag_s",i8.8,"_r",i5.5,".nc")') trim(pmc_diag_dir), itimestep, domain % dminfo % my_proc_id
      call nc_check(nf90_create(trim(fname), NF90_NETCDF4, ncid), 'create ' // trim(fname))
      call nc_check(nf90_def_dim(ncid, 'cell', nCellsSolve, d_cell), 'dim cell')
      call nc_check(nf90_def_dim(ncid, 'level', nVertLevels, d_lev), 'dim level')
      call nc_check(nf90_def_dim(ncid, 'variable', N_DIAG, d_var), 'dim variable')
      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'variables', var_names), 'att')
      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'units', 'number m-3; mass kg m-3; mh_* kg m-3 (mass times entropy)'), 'att')
      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'itimestep', itimestep), 'att')
      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'samples', n_diag_samples), 'att')
      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'steps_averaged', pmc_diag_interval), 'att')
      call nc_check(nf90_def_var(ncid, 'cell_id', NF90_INT, (/d_cell/), varid), 'def cell_id')
      call nc_check(nf90_put_var(ncid, varid, indexToCellID(1:nCellsSolve)), 'put cell_id')
      call nc_check(nf90_def_var(ncid, 'lat', NF90_DOUBLE, (/d_cell/), varid), 'def lat')
      call nc_check(nf90_put_var(ncid, varid, latCell(1:nCellsSolve)), 'put lat')
      call nc_check(nf90_def_var(ncid, 'lon', NF90_DOUBLE, (/d_cell/), varid), 'def lon')
      call nc_check(nf90_put_var(ncid, varid, lonCell(1:nCellsSolve)), 'put lon')
      call nc_check(nf90_def_var(ncid, 'height', NF90_FLOAT, (/d_lev, d_cell/), varid), 'def height')
      call nc_check(nf90_put_var(ncid, varid, zm), 'put height')
      call nc_check(nf90_def_var(ncid, 'mean', NF90_FLOAT, (/d_var, d_lev, d_cell/), varid, deflate_level=1), 'def mean')
      call nc_check(nf90_put_var(ncid, varid, out), 'put mean')
      call nc_check(nf90_close(ncid), 'close')
      if (domain % dminfo % my_proc_id == 0) then
         call mpas_log_write('PartMC diag step $i: time-mean aerosol fields written, $i samples', &
                             intArgs=(/itimestep, n_diag_samples/))
      end if
      diag_acc = 0.0_RKIND
      n_diag_samples = 0
#endif

   end subroutine atm_partmc_diag


   !> Write the particles of a sample of cells, one NetCDF file per MPI task.
   !>
   !> Layout: per-cell variables on dimension "cell", per-particle variables on
   !> "particle". The particles of cell c are part_start(c) .. part_start(c)+n_part(c)-1
   !> (zero-based). Each computational particle represents the number concentration
   !> part_num_conc. num_conc and mass_conc are exact totals over the whole population
   !> of the cell.
   subroutine atm_partmc_output(domain, itimestep, dt)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep
      real (kind=RKIND), intent(in) :: dt

#ifdef MPAS_PARTMC
      type (mpas_pool_type), pointer :: mesh, state, diag
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, zz, zgrid, theta_m, exner, pressure_p, pressure_base
      real (kind=RKIND), dimension(:), pointer :: areaCell, rdzw, latCell, lonCell
      integer, dimension(:), pointer :: indexToCellID
      integer, pointer :: nCellsSolve, nVertLevels, index_qv, index_pmc_ref
      integer :: iCell, k, i_part, n_cell_out, n_part_out, ic, ip, n, i_spec, c
      real (kind=RKIND), allocatable :: p_wt(:), p_lat(:), p_lon(:), p_zeta(:)
      real (kind=RKIND) :: wt
      integer :: g_n_cell_out, g_n_part_out, clock0, clock1, clock_rate
      integer :: ncid, d_cell, d_part, d_spec, varid, stat
      integer, allocatable :: c_id(:), c_lev(:), c_npart(:), c_start(:), p_source(:)
      integer (kind=I8KIND), allocatable :: p_id(:)
      real (kind=RKIND), allocatable :: c_lat(:), c_lon(:), c_z(:), c_vol(:), c_num(:), c_numref(:), &
                                        c_temp(:), c_pres(:), c_qv(:), c_rho(:), c_mass(:,:)
      real (kind=R4KIND), allocatable :: p_mass(:,:)
      logical, allocatable :: sel(:,:)
      character(len=512) :: fname, names
      real (kind=RKIND) :: wall, g_wall

      if (pmc_n_part <= 0 .or. pmc_out_interval <= 0) return
      if (mod(itimestep, pmc_out_interval) /= 0) return
      call system_clock(clock0, clock_rate)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_dimension(state, 'index_pmc_ref', index_pmc_ref)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_base)

      ! choose cells and count
      allocate(sel(nVertLevels, nCellsSolve))
      n_cell_out = 0
      n_part_out = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            sel(k,iCell) = (cell_hash(indexToCellID(iCell), k, itimestep) < pmc_out_frac)
            if (sel(k,iCell)) then
               n_cell_out = n_cell_out + 1
               n_part_out = n_part_out + aero_state_n_part(aero_states(k,iCell))
            end if
         end do
      end do

      if (n_cell_out > 0) then
         allocate(c_id(n_cell_out), c_lev(n_cell_out), c_npart(n_cell_out), c_start(n_cell_out), &
                  c_lat(n_cell_out), c_lon(n_cell_out), c_z(n_cell_out), c_vol(n_cell_out), &
                  c_num(n_cell_out), c_numref(n_cell_out), c_temp(n_cell_out), c_pres(n_cell_out), &
                  c_qv(n_cell_out), c_rho(n_cell_out), c_mass(n_spec, n_cell_out))
         allocate(p_mass(n_spec, n_part_out), p_source(n_part_out), p_id(n_part_out), p_wt(n_part_out))
         if (pmc_mode == MODE_LAGR) allocate(p_lat(n_part_out), p_lon(n_part_out), p_zeta(n_part_out))

         ic = 0
         ip = 0
         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               if (.not. sel(k,iCell)) cycle
               ic = ic + 1
               n = aero_state_n_part(aero_states(k,iCell))
               c_id(ic) = indexToCellID(iCell)
               c_lev(ic) = k
               c_npart(ic) = n
               c_start(ic) = ip
               c_lat(ic) = latCell(iCell)
               c_lon(ic) = lonCell(iCell)
               c_z(ic) = 0.5_RKIND*(zgrid(k,iCell) + zgrid(k+1,iCell))
               c_vol(ic) = cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))
               c_num(ic) = 0.0_RKIND
               c_rho(ic) = rho_zz(k,iCell) * zz(k,iCell)
               c_numref(ic) = c_rho(ic) * scalars(index_pmc_ref,k,iCell)
               c_qv(ic) = scalars(index_qv,k,iCell)
               c_temp(ic) = theta_m(k,iCell) / (1.0_RKIND + rvord*c_qv(ic)) * exner(k,iCell)
               c_pres(ic) = pressure_p(k,iCell) + pressure_base(k,iCell)
               c_mass(:,ic) = 0.0_RKIND
               do i_part = 1, n
                  ip = ip + 1
                  ! exact cell totals are accumulated in double precision before the particle
                  ! masses are stored in single precision
                  c = aero_states(k,iCell)%apa%particle(i_part)%weight_class
                  wt = aero_states(k,iCell)%awa%weight(aero_states(k,iCell)%apa%particle(i_part)%weight_group,c)%magnitude
                  p_wt(ip) = wt
                  if (pmc_mode == MODE_LAGR) then
                     p_lat(ip) = asin(max(-1.0_RKIND, min(1.0_RKIND, pos(k,iCell) % x(3,i_part))))
                     p_lon(ip) = atan2(pos(k,iCell) % x(2,i_part), pos(k,iCell) % x(1,i_part))
                     p_zeta(ip) = pos(k,iCell) % x(4,i_part)
                  end if
                  c_num(ic) = c_num(ic) + wt
                  c_mass(:,ic) = c_mass(:,ic) + aero_states(k,iCell)%apa%particle(i_part)%vol * aero_data%density * wt
                  p_mass(:,ip) = real(aero_states(k,iCell)%apa%particle(i_part)%vol * aero_data%density, kind=R4KIND)
                  p_source(ip) = aero_states(k,iCell)%apa%particle(i_part)%component(1)%source_id
                  p_id(ip) = int(aero_states(k,iCell)%apa%particle(i_part)%id, kind=I8KIND)
               end do
            end do
         end do

         write(fname, '(a,"/particles_s",i8.8,"_r",i5.5,".nc")') trim(pmc_out_dir), itimestep, domain % dminfo % my_proc_id
         call nc_check(nf90_create(trim(fname), NF90_NETCDF4, ncid), 'create ' // trim(fname))
         call nc_check(nf90_def_dim(ncid, 'cell', n_cell_out, d_cell), 'dim cell')
         call nc_check(nf90_def_dim(ncid, 'particle', max(1, n_part_out), d_part), 'dim particle')
         call nc_check(nf90_def_dim(ncid, 'species', n_spec, d_spec), 'dim species')

         names = ''
         do i_spec = 1, n_spec
            names = trim(names) // trim(aero_data%name(i_spec))
            if (i_spec < n_spec) names = trim(names) // ','
         end do
         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'format_version', 2), 'att')
         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'species', trim(names)), 'att')
         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'sources', 'background,blob'), 'att')
         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'itimestep', itimestep), 'att')
         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'time_s', real(itimestep, kind=RKIND)*dt), 'att')
         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'mpi_task', domain % dminfo % my_proc_id), 'att')
         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'n_particles', n_part_out), 'att')
         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'n_part_target', pmc_n_part), 'att')
         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'cell_fraction', pmc_out_frac), 'att')

         call def_put_int('cell_id', d_cell, 'global MPAS cell index', '1', c_id)
         call def_put_int('level', d_cell, 'model level, 1 = lowest', '1', c_lev)
         call def_put_int('n_part', d_cell, 'computational particles in the cell', '1', c_npart)
         call def_put_int('part_start', d_cell, 'zero-based index of the first particle of the cell', '1', c_start)
         call def_put_real('lat', d_cell, NF90_DOUBLE, 'latitude', 'rad', c_lat)
         call def_put_real('lon', d_cell, NF90_DOUBLE, 'longitude', 'rad', c_lon)
         call def_put_real('height', d_cell, NF90_FLOAT, 'height of the cell centre', 'm', c_z)
         call def_put_real('volume', d_cell, NF90_DOUBLE, 'cell volume', 'm3', c_vol)
         call def_put_real('temperature', d_cell, NF90_FLOAT, 'temperature', 'K', c_temp)
         call def_put_real('pressure', d_cell, NF90_FLOAT, 'pressure', 'Pa', c_pres)
         call def_put_real('qv', d_cell, NF90_FLOAT, 'water vapour mixing ratio', 'kg kg-1', c_qv)
         call def_put_real('dry_air_density', d_cell, NF90_FLOAT, 'dry air density', 'kg m-3', c_rho)
         call def_put_real('num_conc', d_cell, NF90_DOUBLE, 'total number concentration', 'm-3', c_num)
         call def_put_real('num_conc_ref', d_cell, NF90_DOUBLE, &
                           'number concentration of the Eulerian reference scalar', 'm-3', c_numref)

         call nc_check(nf90_def_var(ncid, 'mass_conc', NF90_DOUBLE, (/d_spec, d_cell/), varid), 'def mass_conc')
         call nc_check(nf90_put_att(ncid, varid, 'long_name', 'total mass concentration of each species'), 'att')
         call nc_check(nf90_put_att(ncid, varid, 'units', 'kg m-3'), 'att')
         call nc_check(nf90_put_var(ncid, varid, c_mass), 'put mass_conc')

         call nc_check(nf90_def_var(ncid, 'part_mass', NF90_FLOAT, (/d_spec, d_part/), varid, &
                                    deflate_level=1, shuffle=.true.), 'def part_mass')
         call nc_check(nf90_put_att(ncid, varid, 'long_name', 'mass of each species in one particle'), 'att')
         call nc_check(nf90_put_att(ncid, varid, 'units', 'kg'), 'att')
         if (n_part_out > 0) call nc_check(nf90_put_var(ncid, varid, p_mass), 'put part_mass')

         call nc_check(nf90_def_var(ncid, 'part_num_conc', NF90_DOUBLE, (/d_part/), varid), 'def part_num_conc')
         call nc_check(nf90_put_att(ncid, varid, 'long_name', &
                                    'number concentration represented by the computational particle'), 'att')
         call nc_check(nf90_put_att(ncid, varid, 'units', 'm-3'), 'att')
         if (n_part_out > 0) call nc_check(nf90_put_var(ncid, varid, p_wt), 'put part_num_conc')

         if (pmc_mode == MODE_LAGR .and. n_part_out > 0) then
            call nc_check(nf90_def_var(ncid, 'part_lat', NF90_DOUBLE, (/d_part/), varid), 'def part_lat')
            call nc_check(nf90_put_att(ncid, varid, 'units', 'rad'), 'att')
            call nc_check(nf90_put_var(ncid, varid, p_lat), 'put part_lat')
            call nc_check(nf90_def_var(ncid, 'part_lon', NF90_DOUBLE, (/d_part/), varid), 'def part_lon')
            call nc_check(nf90_put_att(ncid, varid, 'units', 'rad'), 'att')
            call nc_check(nf90_put_var(ncid, varid, p_lon), 'put part_lon')
            call nc_check(nf90_def_var(ncid, 'part_zeta', NF90_DOUBLE, (/d_part/), varid), 'def part_zeta')
            call nc_check(nf90_put_att(ncid, varid, 'long_name', 'height coordinate zeta of the particle'), 'att')
            call nc_check(nf90_put_att(ncid, varid, 'units', 'm'), 'att')
            call nc_check(nf90_put_var(ncid, varid, p_zeta), 'put part_zeta')
         end if

         call nc_check(nf90_def_var(ncid, 'part_source', NF90_BYTE, (/d_part/), varid), 'def part_source')
         call nc_check(nf90_put_att(ncid, varid, 'long_name', 'source of the particle, one-based index into sources'), 'att')
         if (n_part_out > 0) call nc_check(nf90_put_var(ncid, varid, p_source), 'put part_source')

         call nc_check(nf90_def_var(ncid, 'part_id', NF90_INT64, (/d_part/), varid), 'def part_id')
         call nc_check(nf90_put_att(ncid, varid, 'long_name', 'particle identifier'), 'att')
         if (n_part_out > 0) call nc_check(nf90_put_var(ncid, varid, p_id), 'put part_id')

         call nc_check(nf90_close(ncid), 'close')
      end if

      call system_clock(clock1)
      wall = real(clock1 - clock0, kind=RKIND) / real(clock_rate, kind=RKIND)
      call mpas_dmpar_sum_int(domain % dminfo, n_cell_out, g_n_cell_out)
      call mpas_dmpar_sum_int(domain % dminfo, n_part_out, g_n_part_out)
      call mpas_dmpar_max_real(domain % dminfo, wall, g_wall)
      call mpas_log_write('PartMC output step $i: wrote $i cells, $i particles, slowest task took $r s', &
                          intArgs=(/itimestep, g_n_cell_out, g_n_part_out/), realArgs=(/g_wall/))

   contains

      subroutine def_put_int(name, dimid, long_name, units, vals)
         character(len=*), intent(in) :: name, long_name, units
         integer, intent(in) :: dimid
         integer, intent(in) :: vals(:)
         integer :: vid
         call nc_check(nf90_def_var(ncid, name, NF90_INT, (/dimid/), vid), 'def ' // name)
         call nc_check(nf90_put_att(ncid, vid, 'long_name', long_name), 'att')
         call nc_check(nf90_put_att(ncid, vid, 'units', units), 'att')
         call nc_check(nf90_put_var(ncid, vid, vals), 'put ' // name)
      end subroutine def_put_int

      subroutine def_put_real(name, dimid, xtype, long_name, units, vals)
         character(len=*), intent(in) :: name, long_name, units
         integer, intent(in) :: dimid, xtype
         real (kind=RKIND), intent(in) :: vals(:)
         integer :: vid
         call nc_check(nf90_def_var(ncid, name, xtype, (/dimid/), vid), 'def ' // name)
         call nc_check(nf90_put_att(ncid, vid, 'long_name', long_name), 'att')
         call nc_check(nf90_put_att(ncid, vid, 'units', units), 'att')
         call nc_check(nf90_put_var(ncid, vid, vals), 'put ' // name)
      end subroutine def_put_real
#endif

   end subroutine atm_partmc_output

end module mpas_atm_partmc
