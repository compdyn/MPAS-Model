! PartMC coupling for MPAS-Atmosphere (development).
!
! A number-mixing-ratio scalar (pmc_num) is advected by MPAS and its total
! time-integrated transport fluxes are captured in atm_advance_scalars_mono_work.
! PartMC particles held in every cell are moved stochastically with those fluxes.
!
! Run-time settings (environment variables):
!   PMC_N_PART        target computational particles per cell (0 = no particles)
!   PMC_COUPLED       1 (default): pmc_num is reset from the particles every step
!   PMC_SCALE_HEIGHT  e-folding height (m) of the initial number mixing ratio
!   PMC_AERO_DATA     PartMC aerosol species file (default aero_data.dat)
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
   use netcdf
   use mpas_atm_partmc_comm
#endif

   implicit none
   private

   public :: atm_partmc_init, atm_partmc_check, atm_partmc_particles_init, atm_partmc_particles_transport, &
             atm_partmc_particles_to_scalar, atm_partmc_output

#ifdef MPAS_PARTMC
   !> Target number of computational particles per cell (0 disables particles). Env: PMC_N_PART.
   integer :: pmc_n_part = 0
   !> Whether pmc_num is reset from the particle counts before every step (Env: PMC_COUPLED, default 1).
   !> If not, particles follow the free-running Eulerian pmc_num (one-way test mode).
   logical :: pmc_coupled = .true.
   !> Aerosol species data.
   type(aero_data_t) :: aero_data
   !> Particle populations, (nVertLevels, nCellsSolve).
   type(aero_state_t), allocatable :: aero_states(:,:)
   !> Particles arriving in each cell during the current transport step.
   type(aero_state_t), allocatable :: aero_incoming(:,:)
   !> MPI tasks that own cells in this block's halo.
   integer, allocatable :: nbr_rank(:)
   !> For each local cell, position of its owner in nbr_rank (0 if owned here).
   integer, allocatable :: nbr_of_cell(:)
   !> Map from global cell ID to local owned-cell index (0 if not owned here).
   integer, allocatable :: global_to_local(:)
   !> Send and receive buffers, one per neighbouring task.
   type(pmc_buf_t), allocatable :: send_buf(:), recv_buf(:)
   !> Number of aerosol species.
   integer :: n_spec = 0
   !> Global mass (kg) of each species at initialization.
   real (kind=RKIND), allocatable :: mass_total_init(:)
   !> Particle output: interval in steps (0 = off), fraction of cells written, directory.
   integer :: pmc_out_interval = 0
   real (kind=RKIND) :: pmc_out_frac = 1.0_RKIND
   character(len=256) :: pmc_out_dir = 'pmc_out'
#endif

contains

   !> Set the initial pmc_num field (kg^-1): uniform background plus a Gaussian
   !> blob in the horizontal, unless the field was read in non-zero.
   subroutine atm_partmc_init(domain)

      type (domain_type), intent(inout) :: domain

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh, state
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      integer, pointer :: index_pmc_num, index_pmc_ref, index_pmc_ref_src2, nCells, nVertLevels
      integer :: iCell, k, stat
      real (kind=RKIND) :: d, cosd, scale_height, z_mid, f_z
      real (kind=RKIND), dimension(:,:), pointer :: zgrid
      character(len=32) :: envval
      real (kind=RKIND), parameter :: q_background = 1.0e8_RKIND
      real (kind=RKIND), parameter :: q_blob = 1.0e9_RKIND
      real (kind=RKIND), parameter :: lat0 = 0.6981317_RKIND   ! 40 N
      real (kind=RKIND), parameter :: lon0 = 1.5707963_RKIND   ! 90 E
      real (kind=RKIND), parameter :: width = 0.25_RKIND       ! radians of arc

      ! optional exponential decay of the mixing ratio with height (m). Env: PMC_SCALE_HEIGHT.
      scale_height = 0.0_RKIND
      call get_environment_variable('PMC_SCALE_HEIGHT', envval, status=stat)
      if (stat == 0) read(envval, *) scale_height

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block % structs, 'state', state)
         nullify(index_pmc_num)
         call mpas_pool_get_dimension(state, 'index_pmc_num', index_pmc_num)
         if (.not. associated(index_pmc_num)) return
         call mpas_pool_get_dimension(state, 'index_pmc_ref', index_pmc_ref)
         call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', index_pmc_ref_src2)
         call mpas_pool_get_dimension(mesh, 'nCells', nCells)
         call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
         call mpas_pool_get_array(mesh, 'latCell', latCell)
         call mpas_pool_get_array(mesh, 'lonCell', lonCell)
         call mpas_pool_get_array(mesh, 'zgrid', zgrid)
         call mpas_pool_get_array(state, 'scalars', scalars, 1)

         if (maxval(scalars(index_pmc_num,:,1:nCells)) <= 0.0_RKIND) then
            do iCell = 1, nCells
               cosd = sin(lat0)*sin(latCell(iCell)) + cos(lat0)*cos(latCell(iCell))*cos(lonCell(iCell)-lon0)
               d = acos(max(-1.0_RKIND, min(1.0_RKIND, cosd)))
               do k = 1, nVertLevels
                  f_z = 1.0_RKIND
                  if (scale_height > 0.0_RKIND) then
                     z_mid = 0.5_RKIND*(zgrid(k,iCell) + zgrid(k+1,iCell))
                     f_z = exp(-z_mid/scale_height)
                  end if
                  ! the blob is the second particle source, the background the first
                  scalars(index_pmc_ref_src2,k,iCell) = f_z * q_blob*exp(-(d/width)**2)
                  scalars(index_pmc_num,k,iCell) = f_z * q_background + scalars(index_pmc_ref_src2,k,iCell)
                  scalars(index_pmc_ref,k,iCell) = scalars(index_pmc_num,k,iCell)
               end do
            end do
            call mpas_log_write('PartMC: initialized pmc_num with background + blob')
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


   !> Physical volume (m^3) of cell (k,iCell).
   pure real (kind=RKIND) function cell_volume(areaCell, rdzw, zz)
      real (kind=RKIND), intent(in) :: areaCell, rdzw, zz
      cell_volume = areaCell / (rdzw * zz)
   end function cell_volume


   !> Create a particle population in every owned cell, consistent with pmc_num.
   subroutine atm_partmc_particles_init(domain)

      type (domain_type), intent(inout) :: domain

#ifdef MPAS_PARTMC
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh, state
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, zz
      integer, pointer :: index_pmc_num, index_src2, nCells, nCellsSolve, nVertLevels
      integer :: iCell, k, i, stat, rank, n2, i_so4, i_nh4, i_oc, i_bc, i_src
      character(len=256) :: envval, aero_data_file
      type(spec_file_t) :: file
      real (kind=dp), allocatable :: vols(:)
      real (kind=dp) :: diam, vol_tot, f2
      type (mpas_pool_type), pointer :: diag
      type (field1DReal), pointer :: owner_field
      real (kind=RKIND), dimension(:), pointer :: owner
      integer, dimension(:), pointer :: indexToCellID
      real (kind=RKIND) :: num_conc
      type(aero_particle_t) :: particle

      call get_environment_variable('PMC_N_PART', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_n_part
      call get_environment_variable('PMC_COUPLED', envval, status=stat)
      if (stat == 0) pmc_coupled = (trim(envval) /= '0')
      if (pmc_n_part <= 0) then
         call mpas_log_write('PartMC: particles disabled (set PMC_N_PART > 0 to enable)')
         return
      end if
      aero_data_file = 'aero_data.dat'
      call get_environment_variable('PMC_AERO_DATA', envval, status=stat)
      if (stat == 0) aero_data_file = envval
      call get_environment_variable('PMC_OUT_INTERVAL', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_out_interval
      call get_environment_variable('PMC_OUT_FRAC', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_out_frac
      call get_environment_variable('PMC_OUT_DIR', envval, status=stat)
      if (stat == 0) pmc_out_dir = envval
      if (pmc_out_interval > 0) call execute_command_line('mkdir -p ' // trim(pmc_out_dir))

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_dimension(state, 'index_pmc_num', index_pmc_num)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)

      call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', index_src2)

      ! species from a PartMC aero_data file; two sources with different composition
      call spec_file_open(trim(aero_data_file), file)
      call spec_file_read_aero_data(file, aero_data)
      call spec_file_close(file)
      call fractal_set_spherical(aero_data%fractal)
      i_src = aero_data_source_by_name(aero_data, 'background')
      i_src = aero_data_source_by_name(aero_data, 'blob')
      i_src = aero_data_weight_class_by_name(aero_data, 'all')
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
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
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

      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            call aero_state_zero(aero_states(k,iCell))
            call aero_state_set_weight(aero_states(k,iCell), aero_data, AERO_STATE_WEIGHT_FLAT)
            call aero_state_set_n_part_ideal(aero_states(k,iCell), real(pmc_n_part, kind=dp))
            call aero_state_zero(aero_incoming(k,iCell))
            call aero_state_set_weight(aero_incoming(k,iCell), aero_data, AERO_STATE_WEIGHT_FLAT)
            ! number concentration (m^-3) = rho * q = rho_zz * zz * q
            num_conc = rho_zz(k,iCell) * zz(k,iCell) * scalars(index_pmc_num,k,iCell)
            aero_states(k,iCell)%awa%weight(1,1)%magnitude = num_conc / real(pmc_n_part, kind=dp)
            ! number of particles from the second source, with randomized rounding
            f2 = real(pmc_n_part, kind=dp) * scalars(index_src2,k,iCell) / scalars(index_pmc_num,k,iCell)
            n2 = int(f2)
            if (pmc_random() < f2 - real(n2, kind=dp)) n2 = n2 + 1
            do i = 1, pmc_n_part
               call aero_particle_zero(particle, aero_data)
               vols = 0.0d0
               if (i <= n2) then
                  ! blob: soot-like, lognormal with geometric mean diameter 50 nm
                  diam = 0.05d-6 * exp(log(1.8d0) * rand_normal(0.0d0, 1.0d0))
                  vol_tot = pii / 6.0d0 * diam**3
                  vols(i_bc) = 0.6d0 * vol_tot
                  vols(i_oc) = 0.4d0 * vol_tot
                  call aero_particle_set_component(particle, 2, 0.0d0)
               else
                  ! background: ammonium sulfate and organics, geometric mean diameter 100 nm
                  diam = 0.1d-6 * exp(log(1.6d0) * rand_normal(0.0d0, 1.0d0))
                  vol_tot = pii / 6.0d0 * diam**3
                  vols(i_so4) = 0.5d0 * vol_tot
                  vols(i_nh4) = 0.2d0 * vol_tot
                  vols(i_oc) = 0.3d0 * vol_tot
                  call aero_particle_set_component(particle, 1, 0.0d0)
               end if
               call aero_particle_set_vols(particle, vols)
               call aero_particle_new_id(particle)
               call aero_particle_set_weight(particle, 1, 1)
               call aero_state_add_particle(aero_states(k,iCell), particle, aero_data)
            end do
         end do
      end do

      call mpas_log_write('PartMC: created $i particles in each of $i cells', &
                          intArgs=(/pmc_n_part, nVertLevels*nCellsSolve/))
      call mpas_log_write('PartMC: $i species from ' // trim(aero_data_file), intArgs=(/n_spec/))
      allocate(mass_total_init(n_spec))
      call species_mass_totals(domain, mass_total_init)
      if (pmc_out_interval > 0) then
         call mpas_log_write('PartMC: particle output every $i steps for cell fraction $r to ' // trim(pmc_out_dir), &
                             intArgs=(/pmc_out_interval/), realArgs=(/pmc_out_frac/))
      end if
      if (pmc_coupled) then
         call mpas_log_write('PartMC: coupled mode, pmc_num is reset from particle counts every step')
      else
         call mpas_log_write('PartMC: one-way mode, particles follow the Eulerian pmc_num')
      end if
#endif

   end subroutine atm_partmc_particles_init


   !> Coupled mode: set pmc_num (time level 1) from the particle counts, so that the
   !> host's scalar transport computes fluxes for the field the particles represent.
   subroutine atm_partmc_particles_to_scalar(domain)

      type (domain_type), intent(inout) :: domain

#ifdef MPAS_PARTMC
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh, state
      type (field3DReal), pointer :: scalars_field
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, zz
      integer, pointer :: index_pmc_num, nCellsSolve, nVertLevels
      integer :: iCell, k

      if (pmc_n_part <= 0 .or. .not. pmc_coupled) return

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_dimension(state, 'index_pmc_num', index_pmc_num)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)

      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            ! number mixing ratio = number concentration / density, density = rho_zz * zz
            scalars(index_pmc_num,k,iCell) = real(aero_state_n_part(aero_states(k,iCell)), kind=RKIND) &
                 * aero_states(k,iCell)%awa%weight(1,1)%magnitude / (rho_zz(k,iCell) * zz(k,iCell))
         end do
      end do

      call mpas_pool_get_field(state, 'scalars', scalars_field, 1)
      call mpas_dmpar_exch_halo_field(scalars_field)
#endif

   end subroutine atm_partmc_particles_to_scalar


   !> Move particles between cells using the captured transport fluxes of pmc_num.
   !>
   !> Each particle in cell c leaves through face f with probability
   !> (outgoing flux through f) / (content of c before transport). A particle
   !> arriving in cell d becomes (W_c / W_d) particles on average, where W is the
   !> number of real particles represented by one computational particle.
   subroutine atm_partmc_particles_transport(domain, itimestep)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep

#ifdef MPAS_PARTMC
      integer, parameter :: MAX_FACES = 12
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh, diag
      real (kind=RKIND), dimension(:,:), pointer :: flux_h, flux_v, c_old, c_new, edgesOnCell_sign, zz
      real (kind=RKIND), dimension(:), pointer :: invAreaCell, areaCell, rdzw
      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, pointer :: nCellsSolve, nVertLevels
      integer :: iCell, k, i, iEdge, n_face, i_face, n_rem, n_move, i_move, i_part
      integer :: dest_cell(MAX_FACES), dest_k(MAX_FACES)
      real (kind=RKIND) :: prob(MAX_FACES), p_rem, w_src
      real (kind=RKIND) :: c_part, err2, norm2, tot_part, tot_euler, vol
      integer :: n_part_total, n_part_min, n_part_max, n_moved
      integer :: g_n_part_total, g_n_part_min, g_n_part_max, g_n_moved, g_n_cells
      real (kind=RKIND) :: g_err2, g_norm2, g_tot_part, g_tot_euler
      integer :: i_nbr, n_bytes, dest_global, comm
      type (mpas_pool_type), pointer :: state
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars_new
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz_new
      integer, pointer :: index_ref, index_src2
      logical :: do_diag
      integer :: n2, i_spec
      real (kind=RKIND) :: c2_part, c2_ref, err2_s, norm2_s, noise2_s, tot2_part, tot2_ref
      real (kind=RKIND) :: g_err2_s, g_norm2_s, g_noise2_s, g_tot2_part, g_tot2_ref, mass_dev
      real (kind=RKIND), allocatable :: mass_total(:)
      real (kind=RKIND) :: p_sum_max, g_p_sum_max, step2, g_step2, c_ref
      real (kind=RKIND) :: p_sum, lost, g_lost
      integer :: n_over, g_n_over
      real (kind=RKIND), allocatable :: err2_k(:), norm2_k(:), g_err2_k(:), g_norm2_k(:)
      real (kind=RKIND) :: ratio_max, ratio_min, g_ratio_max, g_ratio_min, lev_mean, lev_worst
      integer :: n_copy_max, g_n_copy_max, k_worst
      integer, dimension(:), pointer :: indexToCellID
      type(aero_particle_t) :: particle

      if (pmc_n_part <= 0) return

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(diag, 'pmc_flux_h', flux_h)
      call mpas_pool_get_array(diag, 'pmc_flux_v', flux_v)
      call mpas_pool_get_array(diag, 'pmc_content_old', c_old)
      call mpas_pool_get_array(diag, 'pmc_content_new', c_new)
      call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_array(state, 'scalars', scalars_new, 2)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz_new, 2)
      if (pmc_coupled) then
         call mpas_pool_get_dimension(state, 'index_pmc_ref', index_ref)
      else
         call mpas_pool_get_dimension(state, 'index_pmc_num', index_ref)
      end if

      call mpas_pool_get_dimension(state, 'index_pmc_ref_src2', index_src2)
      do_diag = (mod(itimestep, 12) == 0 .or. itimestep <= 2)
      err2_s = 0.0_RKIND
      norm2_s = 0.0_RKIND
      noise2_s = 0.0_RKIND
      tot2_part = 0.0_RKIND
      tot2_ref = 0.0_RKIND

      n_moved = 0
      p_sum_max = 0.0_RKIND
      step2 = 0.0_RKIND
      n_over = 0
      lost = 0.0_RKIND
      ratio_max = 0.0_RKIND
      ratio_min = huge(1.0_RKIND)
      n_copy_max = 0
      allocate(err2_k(nVertLevels), norm2_k(nVertLevels), g_err2_k(nVertLevels), g_norm2_k(nVertLevels))
      err2_k = 0.0_RKIND
      norm2_k = 0.0_RKIND
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            if (c_old(k,iCell) <= 0.0_RKIND) cycle

            ! outgoing probability for each face
            n_face = 0
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
               n_face = n_face + 1
               dest_cell(n_face) = cellsOnCell(i,iCell)
               dest_k(n_face) = k
               prob(n_face) = max(0.0_RKIND, edgesOnCell_sign(i,iCell)*flux_h(k,iEdge)*invAreaCell(iCell)) / c_old(k,iCell)
            end do
            if (k < nVertLevels) then
               n_face = n_face + 1
               dest_cell(n_face) = iCell
               dest_k(n_face) = k + 1
               prob(n_face) = max(0.0_RKIND, rdzw(k)*flux_v(k+1,iCell)) / c_old(k,iCell)
            end if
            if (k > 1) then
               n_face = n_face + 1
               dest_cell(n_face) = iCell
               dest_k(n_face) = k - 1
               prob(n_face) = max(0.0_RKIND, -rdzw(k)*flux_v(k,iCell)) / c_old(k,iCell)
            end if

            ! The limiter bounds the net update, not the gross outflow, so in a noisy field a cell
            ! with much less content than its neighbours can have outflow exceeding its content.
            ! Scale the probabilities to sum to one (as WRF-PartMC does) and record the shortfall.
            p_sum = sum(prob(1:n_face))
            p_sum_max = max(p_sum_max, p_sum)
            if (p_sum > 1.0_RKIND) then
               n_over = n_over + 1
               lost = lost + (p_sum - 1.0_RKIND) * c_old(k,iCell) * zz(k,iCell) &
                             * cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))
               prob(1:n_face) = prob(1:n_face) / p_sum
            end if

            ! multinomial sample of destinations by sequential binomials
            w_src = aero_states(k,iCell)%awa%weight(1,1)%magnitude * cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))
            n_rem = aero_state_n_part(aero_states(k,iCell))
            p_rem = 1.0_RKIND
            do i_face = 1, n_face
               if (n_rem <= 0) exit
               if (prob(i_face) <= 0.0_RKIND) cycle
               n_move = rand_binomial(n_rem, min(1.0_RKIND, prob(i_face)/p_rem))
               p_rem = max(p_rem - prob(i_face), 1.0e-300_RKIND)
               if (n_move == 0) cycle
               i_nbr = nbr_of_cell(dest_cell(i_face))
               do i_move = 1, n_move
                  i_part = pmc_rand_int(n_rem)
                  particle = aero_states(k,iCell)%apa%particle(i_part)
                  call aero_state_remove_particle_no_info(aero_states(k,iCell), i_part)
                  n_rem = n_rem - 1
                  if (i_nbr == 0) then
                     call deliver(particle, w_src, dest_k(i_face), dest_cell(i_face))
                  else
                     ! destination cell is owned by another task: pack for sending
                     n_bytes = 2*pmc_mpi_pack_size_integer(1) + pmc_mpi_pack_size_real(w_src) &
                             + pmc_mpi_pack_size_aero_particle(particle)
                     call pmc_buf_reserve(send_buf(i_nbr), n_bytes)
                     call pmc_mpi_pack_integer(send_buf(i_nbr)%data, send_buf(i_nbr)%pos, indexToCellID(dest_cell(i_face)))
                     call pmc_mpi_pack_integer(send_buf(i_nbr)%data, send_buf(i_nbr)%pos, dest_k(i_face))
                     call pmc_mpi_pack_real(send_buf(i_nbr)%data, send_buf(i_nbr)%pos, w_src)
                     call pmc_mpi_pack_aero_particle(send_buf(i_nbr)%data, send_buf(i_nbr)%pos, particle)
                  end if
               end do
               n_moved = n_moved + n_move
            end do
         end do
      end do

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
               call pmc_mpi_unpack_real(recv_buf(i_nbr)%data, recv_buf(i_nbr)%pos, w_src)
               ! pmc_mpi_unpack_aero_particle allocates component without checking
               if (allocated(particle%component)) deallocate(particle%component)
               call pmc_mpi_unpack_aero_particle(recv_buf(i_nbr)%data, recv_buf(i_nbr)%pos, particle)
               iCell = global_to_local(dest_global)
               if (iCell == 0) call mpas_log_write('PartMC: received particle for cell $i not owned here', &
                                                   intArgs=(/dest_global/), messageType=MPAS_LOG_CRIT)
               call deliver(particle, w_src, k, iCell)
            end do
         end do
      end if

      ! merge arrivals, rebalance, and compare with the Eulerian field
      err2 = 0.0_RKIND
      norm2 = 0.0_RKIND
      tot_part = 0.0_RKIND
      tot_euler = 0.0_RKIND
      n_part_total = 0
      n_part_min = huge(1)
      n_part_max = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            call aero_state_add_particles(aero_states(k,iCell), aero_incoming(k,iCell), aero_data)
            call aero_state_zero(aero_incoming(k,iCell))
            call aero_state_set_weight(aero_incoming(k,iCell), aero_data, AERO_STATE_WEIGHT_FLAT)

            do while (aero_state_n_part(aero_states(k,iCell)) > 2*pmc_n_part)
               call aero_state_halve(aero_states(k,iCell), 1, 1)
            end do
            do while (aero_state_n_part(aero_states(k,iCell)) < pmc_n_part/2 &
                      .and. aero_state_n_part(aero_states(k,iCell)) > 0)
               call aero_state_double(aero_states(k,iCell), aero_data, 1, 1)
            end do

            i = aero_state_n_part(aero_states(k,iCell))
            n_part_total = n_part_total + i
            n_part_min = min(n_part_min, i)
            n_part_max = max(n_part_max, i)
            vol = cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))
            ! particle content in the same units as c_new (number conc / zz)
            c_part = real(i, kind=RKIND) * aero_states(k,iCell)%awa%weight(1,1)%magnitude / zz(k,iCell)
            ! Eulerian reference content, and the one-step expected content c_new
            c_ref = rho_zz_new(k,iCell) * scalars_new(index_ref,k,iCell)
            err2 = err2 + (c_part - c_ref)**2 * vol
            norm2 = norm2 + c_ref**2 * vol
            err2_k(k) = err2_k(k) + (c_part - c_ref)**2 * vol
            norm2_k(k) = norm2_k(k) + c_ref**2 * vol
            step2 = step2 + (c_part - c_new(k,iCell))**2 * vol
            tot_part = tot_part + c_part * zz(k,iCell) * vol
            tot_euler = tot_euler + c_ref * zz(k,iCell) * vol

            ! composition: number from the second source against its Eulerian reference
            if (do_diag) then
               n2 = 0
               do i_part = 1, i
                  if (aero_states(k,iCell)%apa%particle(i_part)%component(1)%source_id == 2) n2 = n2 + 1
               end do
               c2_part = real(n2, kind=RKIND) * aero_states(k,iCell)%awa%weight(1,1)%magnitude / zz(k,iCell)
               c2_ref = rho_zz_new(k,iCell) * scalars_new(index_src2,k,iCell)
               err2_s = err2_s + (c2_part - c2_ref)**2 * vol
               norm2_s = norm2_s + c2_ref**2 * vol
               ! variance of c2_part if the count were Poisson with the reference mean
               noise2_s = noise2_s + c2_ref * aero_states(k,iCell)%awa%weight(1,1)%magnitude / zz(k,iCell) * vol
               tot2_part = tot2_part + c2_part * zz(k,iCell) * vol
               tot2_ref = tot2_ref + c2_ref * zz(k,iCell) * vol
            end if
         end do
      end do

      if (do_diag) then
         call mpas_dmpar_sum_real(domain % dminfo, err2_s, g_err2_s)
         call mpas_dmpar_sum_real(domain % dminfo, norm2_s, g_norm2_s)
         call mpas_dmpar_sum_real(domain % dminfo, noise2_s, g_noise2_s)
         call mpas_dmpar_sum_real(domain % dminfo, tot2_part, g_tot2_part)
         call mpas_dmpar_sum_real(domain % dminfo, tot2_ref, g_tot2_ref)
         call mpas_log_write('PartMC composition step $i: source-2 number L2 rel err = $r, Poisson noise level = $r, ' &
                             // 'total rel diff = $r', intArgs=(/itimestep/), &
                             realArgs=(/sqrt(g_err2_s/g_norm2_s), sqrt(g_noise2_s/g_norm2_s), &
                                        (g_tot2_part - g_tot2_ref)/g_tot2_ref/))
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
         call mpas_dmpar_sum_real(domain % dminfo, err2, g_err2)
         call mpas_dmpar_sum_real(domain % dminfo, norm2, g_norm2)
         call mpas_dmpar_sum_real(domain % dminfo, tot_part, g_tot_part)
         call mpas_dmpar_sum_real(domain % dminfo, tot_euler, g_tot_euler)
         call mpas_dmpar_sum_int(domain % dminfo, n_part_total, g_n_part_total)
         call mpas_dmpar_sum_int(domain % dminfo, n_moved, g_n_moved)
         call mpas_dmpar_sum_int(domain % dminfo, nCellsSolve*nVertLevels, g_n_cells)
         call mpas_dmpar_min_int(domain % dminfo, n_part_min, g_n_part_min)
         call mpas_dmpar_max_int(domain % dminfo, n_part_max, g_n_part_max)
         call mpas_log_write('PartMC particles step $i: L2 rel err = $r, total rel diff = $r, moved frac = $r', &
                             intArgs=(/itimestep/), &
                             realArgs=(/sqrt(g_err2/g_norm2), (g_tot_part - g_tot_euler)/g_tot_euler, &
                                        real(g_n_moved, kind=RKIND)/real(max(1,g_n_part_total), kind=RKIND)/))
         call mpas_log_write('PartMC particles step $i: n_part per cell mean/min/max = $i $i $i', &
                             intArgs=(/itimestep, g_n_part_total/g_n_cells, g_n_part_min, g_n_part_max/))
         call mpas_dmpar_sum_real(domain % dminfo, step2, g_step2)
         call mpas_dmpar_max_real(domain % dminfo, p_sum_max, g_p_sum_max)
         call mpas_log_write('PartMC particles step $i: one-step L2 rel noise = $r, max outflow prob sum = $r', &
                             intArgs=(/itimestep/), realArgs=(/sqrt(g_step2/g_norm2), g_p_sum_max/))
         call mpas_dmpar_sum_real_array(domain % dminfo, nVertLevels, err2_k, g_err2_k)
         call mpas_dmpar_sum_real_array(domain % dminfo, nVertLevels, norm2_k, g_norm2_k)
         call mpas_dmpar_max_real(domain % dminfo, ratio_max, g_ratio_max)
         call mpas_dmpar_min_real(domain % dminfo, ratio_min, g_ratio_min)
         call mpas_dmpar_max_int(domain % dminfo, n_copy_max, g_n_copy_max)
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
         call mpas_log_write('PartMC particles step $i: weight ratio src/dst min = $r, max = $r, max copies of one particle = $i', &
                             intArgs=(/itimestep, g_n_copy_max/), realArgs=(/g_ratio_min, g_ratio_max/))
         call mpas_dmpar_sum_int(domain % dminfo, n_over, g_n_over)
         call mpas_dmpar_sum_real(domain % dminfo, lost, g_lost)
         call mpas_log_write('PartMC particles step $i: cells with outflow > content = $i of $i, unmoved number frac = $r', &
                             intArgs=(/itimestep, g_n_over, g_n_cells/), realArgs=(/g_lost/g_tot_euler/))
      end if

   contains

      !> Put a particle that represented w_from real particles into owned cell (k_to, cell_to),
      !> duplicating or thinning it to match the weight there.
      subroutine deliver(part, w_from, k_to, cell_to)
         type(aero_particle_t), intent(inout) :: part
         real (kind=RKIND), intent(in) :: w_from
         integer, intent(in) :: k_to, cell_to
         real (kind=RKIND) :: r
         integer :: n, j

         r = w_from / (aero_states(k_to,cell_to)%awa%weight(1,1)%magnitude &
                       * cell_volume(areaCell(cell_to), rdzw(k_to), zz(k_to,cell_to)))
         n = int(r)
         if (pmc_random() < r - real(n, kind=RKIND)) n = n + 1
         ratio_max = max(ratio_max, r)
         ratio_min = min(ratio_min, r)
         n_copy_max = max(n_copy_max, n)
         do j = 1, n
            if (j > 1) call aero_particle_new_id(part)
            call aero_state_add_particle(aero_incoming(k_to,cell_to), part, aero_data)
         end do
      end subroutine deliver
#else
   contains
      subroutine deliver_unused()
      end subroutine deliver_unused
#endif

   end subroutine atm_partmc_particles_transport


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
      integer :: iCell, k, i_part

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
               vol_sum = vol_sum + aero_states(k,iCell)%apa%particle(i_part)%vol
            end do
            local = local + vol_sum * aero_data%density * aero_states(k,iCell)%awa%weight(1,1)%magnitude &
                            * cell_volume(areaCell(iCell), rdzw(k), zz(k,iCell))
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


   !> Write the particles of a sample of cells, one NetCDF file per MPI task.
   !>
   !> Layout: per-cell variables on dimension "cell", per-particle variables on
   !> "particle". The particles of cell c are part_start(c) .. part_start(c)+n_part(c)-1
   !> (zero-based). Each computational particle in cell c represents a number
   !> concentration num_conc_per_part(c). num_conc and mass_conc are exact totals
   !> over the whole population of the cell.
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
      integer :: iCell, k, i_part, n_cell_out, n_part_out, ic, ip, n, i_spec
      integer :: g_n_cell_out, g_n_part_out, clock0, clock1, clock_rate
      integer :: ncid, d_cell, d_part, d_spec, varid, stat
      integer, allocatable :: c_id(:), c_lev(:), c_npart(:), c_start(:), p_source(:)
      integer (kind=I8KIND), allocatable :: p_id(:)
      real (kind=RKIND), allocatable :: c_lat(:), c_lon(:), c_z(:), c_vol(:), c_wt(:), c_num(:), c_numref(:), &
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
                  c_lat(n_cell_out), c_lon(n_cell_out), c_z(n_cell_out), c_vol(n_cell_out), c_wt(n_cell_out), &
                  c_num(n_cell_out), c_numref(n_cell_out), c_temp(n_cell_out), c_pres(n_cell_out), &
                  c_qv(n_cell_out), c_rho(n_cell_out), c_mass(n_spec, n_cell_out))
         allocate(p_mass(n_spec, n_part_out), p_source(n_part_out), p_id(n_part_out))

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
               c_wt(ic) = aero_states(k,iCell)%awa%weight(1,1)%magnitude
               c_num(ic) = real(n, kind=RKIND) * c_wt(ic)
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
                  c_mass(:,ic) = c_mass(:,ic) + aero_states(k,iCell)%apa%particle(i_part)%vol * aero_data%density
                  p_mass(:,ip) = real(aero_states(k,iCell)%apa%particle(i_part)%vol * aero_data%density, kind=R4KIND)
                  p_source(ip) = aero_states(k,iCell)%apa%particle(i_part)%component(1)%source_id
                  p_id(ip) = int(aero_states(k,iCell)%apa%particle(i_part)%id, kind=I8KIND)
               end do
               c_mass(:,ic) = c_mass(:,ic) * c_wt(ic)
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
         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'format_version', 1), 'att')
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
         call def_put_real('num_conc_per_part', d_cell, NF90_DOUBLE, &
                           'number concentration represented by one computational particle', 'm-3', c_wt)
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
