! PartMC coupling for MPAS-Atmosphere (development).
!
! Stage 0a: a passive number-mixing-ratio scalar (pmc_num) whose total
! time-integrated transport fluxes are captured in atm_advance_scalars_mono_work,
! and a check that those fluxes reproduce the scalar update.

module mpas_atm_partmc

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log, only : mpas_log_write
#ifdef MPAS_PARTMC
   use pmc_aero_data
   use pmc_aero_state
   use pmc_aero_particle
   use pmc_aero_weight
   use pmc_aero_weight_array
   use pmc_fractal
   use pmc_rand
#endif

   implicit none
   private

   public :: atm_partmc_init, atm_partmc_check, atm_partmc_particles_init, atm_partmc_particles_transport

#ifdef MPAS_PARTMC
   !> Target number of computational particles per cell (0 disables particles). Env: PMC_N_PART.
   integer :: pmc_n_part = 0
   !> Aerosol species data.
   type(aero_data_t) :: aero_data
   !> Particle populations, (nVertLevels, nCellsSolve).
   type(aero_state_t), allocatable :: aero_states(:,:)
   !> Particles arriving in each cell during the current transport step.
   type(aero_state_t), allocatable :: aero_incoming(:,:)
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
      integer, pointer :: index_pmc_num, nCells, nVertLevels
      integer :: iCell, k
      real (kind=RKIND) :: d, cosd
      real (kind=RKIND), parameter :: q_background = 1.0e8_RKIND
      real (kind=RKIND), parameter :: q_blob = 1.0e9_RKIND
      real (kind=RKIND), parameter :: lat0 = 0.6981317_RKIND   ! 40 N
      real (kind=RKIND), parameter :: lon0 = 1.5707963_RKIND   ! 90 E
      real (kind=RKIND), parameter :: width = 0.25_RKIND       ! radians of arc

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block % structs, 'state', state)
         nullify(index_pmc_num)
         call mpas_pool_get_dimension(state, 'index_pmc_num', index_pmc_num)
         if (.not. associated(index_pmc_num)) return
         call mpas_pool_get_dimension(mesh, 'nCells', nCells)
         call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
         call mpas_pool_get_array(mesh, 'latCell', latCell)
         call mpas_pool_get_array(mesh, 'lonCell', lonCell)
         call mpas_pool_get_array(state, 'scalars', scalars, 1)

         if (maxval(scalars(index_pmc_num,:,1:nCells)) <= 0.0_RKIND) then
            do iCell = 1, nCells
               cosd = sin(lat0)*sin(latCell(iCell)) + cos(lat0)*cos(latCell(iCell))*cos(lonCell(iCell)-lon0)
               d = acos(max(-1.0_RKIND, min(1.0_RKIND, cosd)))
               do k = 1, nVertLevels
                  scalars(index_pmc_num,k,iCell) = q_background + q_blob*exp(-(d/width)**2)
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

      err_max = 0.0_RKIND
      c_max = 0.0_RKIND
      out_max = 0.0_RKIND
      total_old = 0.0_RKIND
      total_new = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'diag', diag)
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

      call mpas_log_write('PartMC check step $i: recon rel err = $r, max outflow frac = $r, total rel change = $r', &
                          intArgs=(/itimestep/), &
                          realArgs=(/g_err_max/g_c_max, g_out_max, (g_total_new - g_total_old)/g_total_old/))

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
      integer, pointer :: index_pmc_num, nCells, nCellsSolve, nVertLevels
      integer :: iCell, k, i, stat
      character(len=32) :: envval
      real (kind=RKIND) :: num_conc
      type(aero_particle_t) :: particle

      call get_environment_variable('PMC_N_PART', envval, status=stat)
      if (stat == 0) read(envval, *) pmc_n_part
      if (pmc_n_part <= 0) then
         call mpas_log_write('PartMC: particles disabled (set PMC_N_PART > 0 to enable)')
         return
      end if
      if (domain % dminfo % nprocs > 1) then
         call mpas_log_write('PartMC: particle transport is not yet implemented for more than one MPI task', &
                             messageType=MPAS_LOG_CRIT)
      end if

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

      ! minimal one-species aerosol data
      allocate(aero_data%name(1), aero_data%mosaic_index(1), aero_data%density(1), aero_data%num_ions(1), &
               aero_data%molec_weight(1), aero_data%kappa(1), aero_data%abifm_m(1), aero_data%abifm_c(1), &
               aero_data%source_name(1), aero_data%weight_class_name(1), aero_data%wavelengths(n_swbands))
      aero_data%name(1) = 'SO4'
      aero_data%mosaic_index(1) = 0
      aero_data%density(1) = 1800.0d0
      aero_data%num_ions(1) = 0
      aero_data%molec_weight(1) = 0.096d0
      aero_data%kappa(1) = 0.65d0
      aero_data%abifm_m(1) = 0.0d0
      aero_data%abifm_c(1) = 0.0d0
      aero_data%source_name(1) = 'init'
      aero_data%weight_class_name(1) = 'all'
      aero_data%wavelengths = 550.0d-9
      aero_data%i_water = 0
      call fractal_set_spherical(aero_data%fractal)

      call pmc_srand(20261111, domain % dminfo % my_proc_id)

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
            do i = 1, pmc_n_part
               call aero_particle_zero(particle, aero_data)
               call aero_particle_set_vols(particle, (/ 1.0d-21 * (0.5d0 + pmc_random()) /))
               call aero_particle_new_id(particle)
               call aero_particle_set_weight(particle, 1, 1)
               call aero_state_add_particle(aero_states(k,iCell), particle, aero_data)
            end do
         end do
      end do

      call mpas_log_write('PartMC: created $i particles in each of $i cells', &
                          intArgs=(/pmc_n_part, nVertLevels*nCellsSolve/))
#endif

   end subroutine atm_partmc_particles_init


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
      integer :: iCell, k, i, iEdge, n_face, i_face, n_rem, n_move, i_move, i_part, n_copy, i_copy
      integer :: dest_cell(MAX_FACES), dest_k(MAX_FACES)
      real (kind=RKIND) :: prob(MAX_FACES), p_rem, ratio, w_src, w_dst
      real (kind=RKIND) :: c_part, err2, norm2, tot_part, tot_euler, vol
      integer :: n_part_total, n_part_min, n_part_max, n_moved
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

      n_moved = 0
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
               w_dst = aero_states(dest_k(i_face),dest_cell(i_face))%awa%weight(1,1)%magnitude &
                     * cell_volume(areaCell(dest_cell(i_face)), rdzw(dest_k(i_face)), zz(dest_k(i_face),dest_cell(i_face)))
               ratio = w_src / w_dst
               do i_move = 1, n_move
                  i_part = pmc_rand_int(n_rem)
                  particle = aero_states(k,iCell)%apa%particle(i_part)
                  call aero_state_remove_particle_no_info(aero_states(k,iCell), i_part)
                  n_rem = n_rem - 1
                  n_copy = int(ratio)
                  if (pmc_random() < ratio - real(n_copy, kind=RKIND)) n_copy = n_copy + 1
                  do i_copy = 1, n_copy
                     if (i_copy > 1) call aero_particle_new_id(particle)
                     call aero_state_add_particle(aero_incoming(dest_k(i_face),dest_cell(i_face)), particle, aero_data)
                  end do
               end do
               n_moved = n_moved + n_move
            end do
         end do
      end do

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
            err2 = err2 + (c_part - c_new(k,iCell))**2 * vol
            norm2 = norm2 + c_new(k,iCell)**2 * vol
            tot_part = tot_part + c_part * zz(k,iCell) * vol
            tot_euler = tot_euler + c_new(k,iCell) * zz(k,iCell) * vol
         end do
      end do

      if (mod(itimestep, 12) == 0 .or. itimestep <= 2) then
         call mpas_log_write('PartMC particles step $i: L2 rel err = $r, total rel diff = $r, moved frac = $r', &
                             intArgs=(/itimestep/), &
                             realArgs=(/sqrt(err2/norm2), (tot_part - tot_euler)/tot_euler, &
                                        real(n_moved, kind=RKIND)/real(max(1,n_part_total), kind=RKIND)/))
         call mpas_log_write('PartMC particles step $i: n_part per cell mean/min/max = $i $i $i', &
                             intArgs=(/itimestep, n_part_total/(nCellsSolve*nVertLevels), n_part_min, n_part_max/))
      end if
#endif

   end subroutine atm_partmc_particles_transport

end module mpas_atm_partmc
