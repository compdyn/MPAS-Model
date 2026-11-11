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
!   PMC_LAGR_SUBSTEPS lagrangian mode: trajectory sub-steps per model step (default 1)
!   PMC_TEST_WIND     lagrangian mode: replace the model winds by an analytic flow with a
!                     12-day period: solid (tilted solid-body rotation) or deform
!                     (Nair and Lauritzen 2010 reversing deformational flow)
!   PMC_N_CLASS       flux mode: 1 (default) one weight class; 2 one class per source,
!                     each with its own advected number scalar
!   PMC_COUPLED       1 (default): the class scalars are reset from the particles every step
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
   use netcdf
   use mpas_atm_partmc_comm
#endif

   implicit none
   private

   public :: atm_partmc_init, atm_partmc_check, atm_partmc_particles_init, atm_partmc_particles_transport, &
             atm_partmc_particles_to_scalar, atm_partmc_output

   integer, parameter :: MODE_FLUX = 1, MODE_LAGR = 2
   integer, parameter :: MAX_CLASS = 2
   !> Reals stored per particle in Lagrangian mode: unit vector, zeta, initial unit vector.
   integer, parameter :: NPOS = 7
   integer, parameter :: WIND_MODEL = 0, WIND_SOLID = 1, WIND_DEFORM = 2
   !> Period of the analytic test winds (s).
   real (kind=RKIND), parameter :: test_period = 12.0_RKIND * 86400.0_RKIND

   !> Initial number mixing ratios (kg^-1) and shape of the blob.
   real (kind=RKIND), parameter :: q_background = 1.0e8_RKIND
   real (kind=RKIND), parameter :: q_blob = 1.0e9_RKIND
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
   !> Particle output: interval in steps (0 = off), fraction of cells written, directory.
   integer :: pmc_out_interval = 0
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
   !> Real particles represented by each computational particle (the same for all).
   real (kind=RKIND) :: lagr_w = 0.0_RKIND
   !> Unit vectors of the cell centres (3, nCells), zeta of faces and of level centres.
   real (kind=RKIND), allocatable :: cxyz(:,:), zf(:), zc(:)
   real (kind=RKIND) :: radius = 1.0_RKIND
   !> Source-2 second moments at the first diagnostic step.
   real (kind=RKIND) :: m2_ref0 = 0.0_RKIND, m2_part0 = 0.0_RKIND

contains

   subroutine read_settings()
      character(len=256) :: envval
      integer :: stat

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
      call get_environment_variable('PMC_LAGR_SUBSTEPS', envval, status=stat)
      if (stat == 0) read(envval, *) lagr_substeps
      call get_environment_variable('PMC_N_CLASS', envval, status=stat)
      if (stat == 0) read(envval, *) n_class
      if (pmc_mode == MODE_LAGR) n_class = 1
      n_class = max(1, min(MAX_CLASS, n_class))
      call get_environment_variable('PMC_AERO_DATA', envval, status=stat)
      if (stat == 0) aero_data_file = envval
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
      integer :: iCell, k, i, j, rank, n, n2, n_t, i_so4, i_nh4, i_oc, i_bc, i_src, g_n_cells
      type (field1DReal), pointer :: owner_field
      real (kind=RKIND), dimension(:), pointer :: owner
      type(spec_file_t) :: file
      real (kind=dp), allocatable :: vols(:)
      real (kind=RKIND) :: conc_tot, conc1, conc2, vol, f2, total, g_total, r_max, x(NPOS), e1(3), e2(3), a, b, v(3)
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
      i_src = aero_data_source_by_name(aero_data, 'background')
      i_src = aero_data_source_by_name(aero_data, 'blob')
      if (n_class == 1) then
         i_src = aero_data_weight_class_by_name(aero_data, 'all')
      else
         i_src = aero_data_weight_class_by_name(aero_data, 'background')
         i_src = aero_data_weight_class_by_name(aero_data, 'blob')
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
         if (pmc_mode == MODE_LAGR) then
            ! tangent-plane basis and radius of a disc that covers the cell
            v = (/ 0.0_RKIND, 0.0_RKIND, 1.0_RKIND /)
            if (abs(cxyz(3,iCell)) > 0.9_RKIND) v = (/ 1.0_RKIND, 0.0_RKIND, 0.0_RKIND /)
            e1 = v - dot_product(v, cxyz(:,iCell)) * cxyz(:,iCell)
            e1 = e1 / norm2(e1)
            e2(1) = cxyz(2,iCell)*e1(3) - cxyz(3,iCell)*e1(2)
            e2(2) = cxyz(3,iCell)*e1(1) - cxyz(1,iCell)*e1(3)
            e2(3) = cxyz(1,iCell)*e1(2) - cxyz(2,iCell)*e1(1)
            r_max = 0.0_RKIND
            do j = 1, nEdgesOnCell(iCell)
               i = verticesOnCell(j,iCell)
               v = (/ xVertex(i), yVertex(i), zVertex(i) /)
               v = v / norm2(v)
               r_max = max(r_max, norm2(v - cxyz(:,iCell)))
            end do
            r_max = 1.05_RKIND * r_max
         end if
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
               f2 = conc_tot * vol / lagr_w
               n = int(f2)
               if (pmc_random() < f2 - real(n, kind=RKIND)) n = n + 1
               aero_states(k,iCell)%awa%weight(1,1)%magnitude = lagr_w / vol
               do i = 1, n
                  if (pmc_random() < conc2/conc_tot) then
                     call make_particle(2, 1)
                  else
                     call make_particle(1, 1)
                  end if
                  ! uniform position in the cell
                  do
                     do
                        a = 2.0_RKIND*pmc_random() - 1.0_RKIND
                        b = 2.0_RKIND*pmc_random() - 1.0_RKIND
                        if (a*a + b*b <= 1.0_RKIND) exit
                     end do
                     x(1:3) = cxyz(:,iCell) + r_max*(a*e1 + b*e2)
                     x(1:3) = x(1:3) / norm2(x(1:3))
                     if (find_cell(x(1:3), iCell, nCells, nEdgesOnCell, cellsOnCell) == iCell) exit
                  end do
                  x(4) = zf(k) + pmc_random()*(zf(k+1) - zf(k))
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
                     call make_particle(2, 1)
                  else
                     call make_particle(1, 1)
                  end if
               end do
               cnt(1,k,iCell) = pmc_n_part
            else
               ! one weight class per source, each with half of the particles
               n_t = pmc_n_part / n_class
               aero_states(k,iCell)%awa%weight(1,1)%magnitude = conc1 / real(n_t, kind=dp)
               do i = 1, n_t
                  call make_particle(1, 1)
               end do
               cnt(1,k,iCell) = n_t
               if (conc2 > 1.0e-3_RKIND * conc1) then
                  aero_states(k,iCell)%awa%weight(1,2)%magnitude = conc2 / real(n_t, kind=dp)
                  do i = 1, n_t
                     call make_particle(2, 2)
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

   contains

      !> Add a new particle from source i_source to cell (k,iCell) in weight class i_class.
      subroutine make_particle(i_source, i_class)
         integer, intent(in) :: i_source, i_class
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
         call aero_particle_set_vols(particle, vols)
         call aero_particle_new_id(particle)
         call aero_particle_set_weight(particle, 1, i_class)
         call aero_state_add_particle(aero_states(k,iCell), particle, aero_data)
      end subroutine make_particle
#endif

   end subroutine atm_partmc_particles_init


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
      real (kind=RKIND), dimension(:,:), pointer :: ux, uy, uz, om, ruAvg, wwAvg
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars_new
      real (kind=RKIND), dimension(:), pointer :: invAreaCell, areaCell, rdzw, fzm, fzp
      real (kind=RKIND), pointer :: config_dt
      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, indexToCellID
      integer, pointer :: nCells, nCellsSolve, nEdges, nVertLevels, index_ref, index_src2
      integer :: iCell, k, i, c, i_part, i_nbr, n_bytes, dest_global, comm, i_spec, k_worst
      integer :: n_part_total, n_part_min, n_part_max, n_moved, n_over, n_copy_max, n_clamp, n_t, n_thin, g_n_thin
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
      n_t = pmc_n_part / n_class
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            if (pmc_mode == MODE_LAGR) then
               do i = 1, pos_in(k,iCell) % n
                  call aero_state_add_particle(aero_states(k,iCell), aero_incoming(k,iCell)%apa%particle(i), aero_data)
                  call pos_append(pos(k,iCell), pos_in(k,iCell) % x(:,i))
               end do
               pos_in(k,iCell) % n = 0
            else
               call aero_state_add_particles(aero_states(k,iCell), aero_incoming(k,iCell), aero_data)
            end if
            call aero_state_zero(aero_incoming(k,iCell))
            call aero_state_set_weight(aero_incoming(k,iCell), aero_data, AERO_STATE_WEIGHT_FLAT_SPECIFIED)

            call count_classes(aero_states(k,iCell), counts)
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
               m = aero_states(k,iCell)%awa%weight(1,c)%magnitude / zz(k,iCell)
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
            call mpas_log_write('PartMC lagrangian step $i: particles held at the top or bottom = $i', &
                                intArgs=(/itimestep, g_n_clamp/))
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
         type (field2DReal), pointer :: f_ux, f_uy, f_uz, f_om
         real (kind=RKIND), allocatable :: u_e(:,:), u_zonal(:,:), u_merid(:,:)
         real (kind=RKIND) :: rho_e, rho_f, u_test, v_test
         integer :: iEdge, c1, c2, kk, ic

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
                                   * aero_states(k,iCell)%awa%weight(1,c)%magnitude
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
                  wt = aero_states(k,iCell)%awa%weight(1,c)%magnitude
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
