!This program will calculate the ZND structure of a steady-state detonation using the
!reaction information used by MESA. It is intended to work as Frank Timmes' torch code
!does, where the user specifies the ambient conditions, a reaction network, and a
!detonation velocity. The ZND equations are then integrated from the conditions calculated
!at the Neumann point. The code will also calculate the expected CJ detonation velocity,
!and the user may choose to use that as the detonation velocity.

!Started on 8/28/12 by Kevin Moore (based off of the modified sample_net.f90 in the same
!directory)

!All solvers (CJ, Neumann, ZND) implemented as of 9/10/12

!Next issue (9/11/12): start checking if just using triple-alpha to c12 for our CJ calculation
!and in our reaction network to see if CJ conditions are reached at the end of a ZND integration.
!Also, check if we remain on the Rayleigh line during the ZND integration.

!As of 9/15/15, mass, momentum, and energy are all now being conserved! Working out 
!analytic Jacobian next (for the new composition derivative terms added)

!9/17/12: Analytic Jacobian working (still ~10% differences between the analytic and
!numeric d/dXj() terms, and also large differences between the d_dxdt_dx(:,species) terms
!(these don't seem to affect the overall steps being taken though - not sure why, but
!putting it on the back burner for now).

!9/18/12: Problems found in the analytic Jacobian using larger nets (with prot/neut - 
!specifically approx19.net). Further examination of the net module's derivatives is
!underway in the net_derivs folder.

!9/19/12: Analytic Jacobian now seems correct! (I was missing a dP/dzbar factor in term1 
!of dT/dx) Energy, momentum, and mass fluxes are conserved to ~1e-3

!9/21/12: Revised method calls to work with MESA version 4442 (previous was for 4298)

!10/23/12: Implementation of blowout included (constant mass, momentum, and energy fluxes
!through an expanding slab of material)

!10/24/12: Blowout seems to be consistent (conservation laws satisfied to 1e-3), and is 
!working fine with numerical Jacobian - time for some science! Analytic Jacobian will come
!soon.

!10/26/12: Added functionality to sweep through detonation velocities and record info
!about the sonic point (for investigating the existence of CJ-like solutions in the
!presence of blowout).

!11/7/12: Added and tested analytic Jacobian terms for blowout, everything seems to work
!so far.

!11/13/12: Fixed error with iwork values getting corrupted with large nets (>45 isos).
!Problem was lines containing rpar(11:11+num_vars) = ... where rpar had a pre-set length
!and overflow was presumably corrupting other arrays.

!11/16/12: Added gravity to the radial expansion equation so that the inputs are a WD
!mass and radius (initial scale height is computed from those).

!12/3/12: Vertical velocity and gravitational potential energy treated consistently in
!blowout, with solver conserving total energy used

!12/10/12: Curvature effects and control flag added (assuming steady curvature)

!12/16/12: Added a new (local) prescription for blowout that uses pressure differences
!over a distance determined by the length that a sound wave can travel. This is enabled
!with the do_blowout_local flag (do_blowout used for the standard, scale height dependent
!prescription that didn't seem to match well to FLASH)

!Feb 2013: Added many new running modes for computing sweeps through detonation velocities,
!ZND lengths, and other quantities. Reorganized inlist and updated documentation. Bill W.
!suggests naming the code eZND for EZ-ZND (har-har), although maybe expansive-ZND would
!make sense...

!Apr 2013: Added more code with mapping solutions to WD parameters

!Aug 2014: Fixed some memory leaks with Valgrind so this runs with a 206-isotope network

module znd

	use chem_def
    use chem_lib
    use num_def
    use num_lib
    use net_def
 	use net_lib
 	use mtx_def
 	use mtx_lib
 	use eos_def
 	use eos_lib
    !use alert_lib
    use const_def
    use rates_def
	!use screen_def
	use crlibm_lib
	
	use lane_emden_solver
	use wd_solver
	
	!-------------------------------------------------------------------------------------
	!Specific variables for the underlying implicit solver we'll be calling:
	! args for isolve -- see num_isolve.dek in num/public
      integer :: which_solver ! as defined in num_def.f

      double precision :: x 
         ! input: initial x value
         ! output: x value for which the solution has been computed.
      double precision :: xend ! desired final x value (positive or negative)
      double precision :: h 
         ! input: initial step size guess
         ! output: predicted next step size from the last accepted step
      double precision :: max_step_size
      integer :: max_steps
      
      ! absolute and relative error tolerances
      double precision, pointer :: rtol(:), atol(:) 
      integer :: itol
      double precision :: y_min, y_max

      ! information about the jacobian matrix
      integer :: ijac, nzmax, isparse, mljac, mujac

      ! information about the "mass" matrix
      integer :: imas, mlmas, mumas
      
      ! switch for calling the subroutine solout or nor
      integer :: iout
      
      integer :: lrd, lid
      double precision, pointer :: rpar_decsol(:) ! (lrd)
      integer, pointer :: ipar_decsol(:) ! (lid)

      integer :: caller_id_blk, nvar_blk, nz_blk
      real(dp), dimension(:), pointer :: lblk, dblk, ublk ! =(nvar,nvar,nz)
      real(dp), dimension(:), pointer :: uf_lblk, uf_dblk, uf_ublk ! =(nvar,nvar,nz)
      
      ! work arrays.
      integer :: lwork_isolve, liwork_isolve
      double precision, pointer :: work_isolve(:) ! (lwork_isolve)
      integer, pointer :: iwork_isolve(:) ! (liwork_isolve)
         
      ! parameter arrays.
      integer, parameter :: lrpar = 500, lipar = 1
      double precision, pointer :: rpar(:)
      integer, pointer :: ipar(:)
               
      ! io unit for warnings and errors
      integer :: lout
      
      ! result code
      integer :: idid
      
      ! stiffness parameter for Vdpol equations
      double precision, parameter :: mu = 1d-3
	
	double precision, pointer :: ending_x(:) ! (species)
	integer :: nfcn    ! number of function evaluations
	integer :: njac    ! number of jacobian evaluations
	integer :: nstep   ! number of computed steps
	integer :: naccpt  ! number of accepted steps
	integer :: nrejct  ! number of rejected steps
      
	integer, parameter :: solver_name_len = 10
	character (len=solver_name_len) :: solver_name
	
	!End isolve controls
	!-------------------------------------------------------------------------------------
	
	!Extra variables for Newton solver (for getting CJ conditions)------------------------
	! dimensions
    integer, parameter :: nz = 1 ! number of zones
    integer, parameter :: nvar = 3 ! number of variables per zone
    integer, parameter :: neq = nz*nvar
    integer, parameter :: nvar_neumann = 2 ! number of variables per zone
    integer, parameter :: neq_neumann = nz*nvar_neumann
    integer, parameter :: nvar_hug = 1 ! number of variables per zone
    integer, parameter :: neq_hug = nz*nvar_hug
    
    ! information about the bandwidth of the jacobian matrix
    ! we have a square matrix, so set number zones sub and super both to 0
    integer, parameter :: stencil_zones_subdiagonal = 0
    integer, parameter :: stencil_zones_superdiagonal = 0
    integer, parameter :: m1 = (stencil_zones_subdiagonal+1)*nvar-1 ! number of subdiagonals
    integer, parameter :: m2 = (stencil_zones_superdiagonal+1)*nvar-1 ! number of superdiagonals
    integer, parameter :: m1_neumann = (stencil_zones_subdiagonal+1)*nvar_neumann-1 ! number of subdiagonals
    integer, parameter :: m2_neumann = (stencil_zones_superdiagonal+1)*nvar_neumann-1 ! number of superdiagonals
    integer, parameter :: m1_hug = (stencil_zones_subdiagonal+1)*nvar_hug-1 ! number of subdiagonals
    integer, parameter :: m2_hug = (stencil_zones_superdiagonal+1)*nvar_hug-1 ! number of superdiagonals
    
    ! equation residuals
    double precision, pointer, dimension(:,:) :: equ_newton
    ! equ_newton(i) has the residual for equation i, i.e., the difference between
    ! the left and right hand sides of the equation.  So the goal
    ! of the solver is to find an approximate solution that makes the
    ! magnitude of the residuals acceptably small.
    
    ! the primary variables
	double precision :: x0, x1, x2
    double precision, pointer, dimension(:,:) :: x_newton ! new vector of primaries
    double precision, pointer, dimension(:,:) :: xold_newton ! old vector of primaries
    double precision, pointer, dimension(:,:) :: dx_newton ! increment vector -- on input is initial guess for x - xold
    double precision, pointer, dimension(:,:) :: xscale_newton! typical values
    real(dp), pointer, dimension(:) :: equ1d, x1d, xold1d, dx1d, xscale1d, y1d !(neq)
    
    ! the secondary variables
    integer, parameter :: nsec = 0 ! number of secondaries per zone
    !double precision, pointer :: t1, t2
    integer, parameter :: ldy = nz+nsec ! leading dimension of y, >= nz
    double precision, pointer, dimension(:,:) :: y_newton! the values
	
    logical :: do_numerical_jacobian
    integer :: liwork_newton, lwork_newton, lid_newton, lrd_newton, which_decsol
    integer, dimension(:), pointer :: iwork_newton
    double precision, dimension(:), pointer :: work_newton
    !real(qp), dimension(:), pointer :: qwork_newton
    
    integer, dimension(:), pointer :: ipar_decsol_newton !(lid_newton)
    double precision, dimension(:), pointer :: rpar_decsol_newton !(lrd_newton)
    double precision, dimension(:), pointer :: AF1
    
    logical :: first_step, nonconv, numerical_jacobian, doing_jacobian
    double precision :: tol_correction_norm, tol_max_correction, tol_residual_norm, epsder, &
       dt1_dx0, dt1_dx1, dt1_dx2, dt2_dx0, dt2_dx1, dt2_dx2
    integer :: matrix_type, matrix_type_current
    
    double precision :: q 						!Energy release (erg/g) for finding CJ point
    double precision :: q_hug 					!Energy release (erg/g) while calculating hugoniots
    double precision :: v_det					!detonation velocity
    double precision :: xrho_hug, xt_hug		!Multipliers (eg. xrho = rho/rho0)
    double precision :: xrho, xt, xd			!Multipliers (eg. xrho = rho/rho0)
    double precision :: rho_cj, t_cj, p_cj, e_cj, v_cj!CJ point thermodynamic conditions
    double precision :: rho_neumann, t_neumann, p_neumann, e_neumann, cs_neumann	!Neumann point thermodynamic conditions
    !Vector for the assumed final composition in the CJ calculation
	double precision, pointer, dimension(:) :: xa_end
	!Vector for the  composition in the Hugoniot calculation
	double precision, pointer, dimension(:) :: xa_hug
    
	!End Newton solver variables
	!-------------------------------------------------------------------------------------
	
	!Add all the necessary shared data up here (eg. chem, net, rate info, etc.)
	integer :: ierr, handle, which_rates_choice, num_reactions, species, num_vars
	double precision, pointer, dimension(:) :: vars !(num_vars)
    integer, pointer :: which_rates(:), chem_id(:), net_iso(:)
    character (len=256) :: my_mesa_dir, net_file
    integer, parameter :: file_path_length = 1024
    logical :: dbg	!Whether to turn debugging output on or not
    !type (Net_General_Info), pointer  :: g	!Allocated at end of setup_net
    
    !EOS variables:
    !Common named variables:
	double precision, pointer :: dabar_dx(:), dzbar_dx(:), dmc_dx(:)
	integer :: info
	integer :: eos_handle 	!ID number for calling the EOS module
	character (len=file_path_length) :: eos_file_prefix		!Which EOS tables to use
	logical :: use_cache
	double precision :: zm !mass fraction of metals, Z
	double precision :: Pgas, logPgas
	double precision :: dlnRho_dlnPgas_const_T
    double precision :: dlnRho_dlnT_const_Pgas
	double precision, dimension(num_eos_basic_results) :: res  	!Basic EOS result array
	double precision :: d_dlnRho_const_T(num_eos_basic_results) !Basic EOS derivs array
    double precision :: d_dlnT_const_Rho(num_eos_basic_results) !Basic EOS derivs array
    double precision :: d_dabar_const_TRho(num_eos_basic_results)	!Basic EOS result array
    double precision :: d_dzbar_const_TRho(num_eos_basic_results)	!Basic EOS result array
    double precision, dimension(num_helm_results) :: res_helm   !HELM EOS array
    
	!Specific variables for the burn inlist:
    double precision :: rho0, t0, p0, cs0, g0	!Initial thermodynamic conditions
    double precision :: burn_rho, logRho, burn_t, burn_P, logT, burn_time, burn_u, burn_e
    integer, parameter :: max_num_isos_for_Xinit = 100
    character(len=iso_name_length) :: names_of_isos_for_Xinit(max_num_isos_for_Xinit)
    character(len=iso_name_length) :: names_of_isos_for_Xcj(max_num_isos_for_Xinit)
    character(len=iso_name_length) :: cj_final_iso, cj_final_iso2
    double precision :: values_for_Xinit(max_num_isos_for_Xinit)
    double precision :: values_for_Xcj(max_num_isos_for_Xinit)
    double precision :: h_scale, uy_init
   	integer :: num_isos_for_Xinit, num_isos_for_Xcj
    logical :: use_solar
    logical :: do_constp, do_znd, do_blowout, do_curvature, do_neumann, do_sweep
    logical :: do_blowout_local, do_h_search, do_lznd_vdet, do_reconstruct_lznd_vdet
    logical :: do_vdet_search, use_rc_hs_scaling, use_uy_ux_ratio, use_he_clavin
    logical :: use_const_uy
    logical :: do_cj		!Whether to integrate using the calculated CJ velocity
    logical :: do_rho_sweep
    logical :: do_lznd_wd
    logical :: do_pathological, do_pathological_init
    logical :: do_cjstate_only, do_cjstate_only2	!Only calculate the CJ state - don't do burning
    logical :: wd_env_use_le		!Use Lane-Emden solver for envelope (polytropes only)
    logical :: wd_env_use_mesa		!Use MESA EOS to integrate stellar structure equations
    logical :: use_variable_linearization_ratio
    !double precision :: n_le		!Polytrope index, P = k_le*rho**(1+1/n_le)
    double precision :: env_rho_min, env_rho_max	!Limit on rho0 in WD envelope
    double precision :: m_env_min, m_env_max	!Bounds on envelope mass to use in do_lznd_wd
    double precision :: shock_lead_hp_frac	!Fraction of pressure scale height at envelope base where the leading point on the detonation front is
    double precision :: pathological_linearization_ratio	!Ratio of current x coordinate to use for linearization around pathological point
    double precision :: sonic_limit, sonic_limit_pathological	!When to cut off integration (mach number)
	double precision :: d_hs_scaling			!d = d_hs_scaling*h_scale
	double precision :: rho_sweep_T0			!Initial temerature (for calculating scale height)
	double precision :: rho_sweep_grav			!cm/s^2 (M_sun in R_earth is 3.3d8 cm/s^2)
	double precision :: rho_sweep_log_min		!Lowest density to use
	double precision :: rho_sweep_log_max	!Highest density to use
	integer :: rho_sweep_num_steps		!Number of steps in density to take
    integer :: j_curve !0: plane parallel, 1: cylindrical, 2: spherical
    double precision :: uy_ux_ratio	!Ratio used for uy_init = ux_init*uy_ux_ratio
    double precision :: rc_hs_factor	!R_c = rc_hs_factor*H_s (~3/5, empirically)
    double precision :: l_burn_factor !Fraction of q_tot to use for defining l_burn
	double precision :: burn_constPgas !Actual Pgas to use in a constant pressure burn
	double precision :: h_search_vdet		!Fixed detonation velocity to use in search
	double precision :: h_search_bracket_low, h_search_bracket_high	!search range for H
	integer :: h_search_numsteps	!Max number of steps to try when searching for H
	double precision :: h_search_tol		!H accepted if deltaH/H_prev < h_search_tol
	double precision :: h_search_xmax		!Max location to integrate to in each step
	double precision :: vdet_search_h_scale		!Fixed scale height parameter to use in v_det search
	double precision :: vdet_search_bracket_low, vdet_search_bracket_high	!search range for v_det
	integer :: vdet_search_numsteps	!Max number of steps to try when searching for v_det
	double precision :: vdet_search_tol		!v_det accepted if delta(v_det)/v_det_prev < v_det_search_tol
	double precision :: vdet_search_xmax		!Max location to integrate to in each step
	double precision :: rtol_init, atol_init
    
    !Variables for interacting with the net module:
    !(pointers here will be allocated in the do1_net_burn subroutine)
    integer :: screening_mode
	integer :: caller_id
	logical :: reuse_given_rates
	double precision, pointer, dimension(:) :: xa !(species)
	double precision, pointer, dimension(:) :: Y_mass !(species)
	double precision ::  eta, d_eta_dlnT, d_eta_dlnRho
	double precision ::  t_start, t_end
	double precision, pointer, dimension(:) ::  d_eps_nuc_dx !(species)
	double precision, pointer, dimension(:) ::  rate_factors !(num_reactions)
        double precision :: weak_rate_factor
	double precision ::  xh, xhe
	double precision ::  abar, zbar, z2bar, ye
	double precision ::  xsum
	double precision ::  mass_correction
	double precision ::  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT
	double precision ::  theta_e_for_graboske_et_al
	double precision, pointer, dimension(:) ::  eps_nuc_categories !(num_categories)
	double precision ::  eps_neu_total
	double precision, pointer, dimension(:) ::  dxdt !(species)
	double precision, pointer, dimension(:) ::  d_dxdt_dRho !(species)
	double precision, pointer, dimension(:) ::  d_dxdt_dT !(species)
	double precision, pointer, dimension(:,:) ::  d_dxdt_dx !(species,species)
	double precision ::  approx_abar
	double precision ::  approx_zbar
	integer :: lwork_net
	double precision, pointer :: work_net(:) ! (lwork_net)
	!std_reaction_Qs, etc. are in rates_def
	integer :: num_times ! ending time is times(num_times); starting time is 0
	double precision, pointer :: times(:) ! (num_times) 
	real(dp), dimension(:), pointer :: log10Ts_f1, log10Rhos_f1, etas_f1, log10Ps_f1
	!(4,numtimes) interpolant for log10T(time), log10Rho(time), eta(time), log10Ps(time):
	real(dp), dimension(:,:), pointer :: log10Ts_f, log10Rhos_f, etas_f, log10Ps_f
	double precision, pointer :: dxdt_source_term(:)
	
	!Other utility variables for the burn:
	double precision :: q_init				!For computing q_bind during ZND integration
	double precision :: mass_flux, mom_flux, energy_flux, rayleigh	!Consistency checks
	double precision :: dP_dx, dE_dx, dq_dx, dEtot_dx, dEtot_dx_eqn	!Consistency checks
	double precision :: m_c, r_wd	!mass & radius of WD (for computing g and H_s)
	double precision :: grav_init !Gravitational acceleration used in blowout (cm/s^2)
	double precision :: r_curve !Radius of curvature to use in ZND equations
	double precision :: delta_y_init !Initital length scale for computing dP/dy
	double precision :: sweep_min_mach, sweep_max_mach !Bounds for v_det if do_sweep
	integer :: num_steps !How many x-steps to output the integration
	integer :: sweep_znd_iters !How many steps to take in v_det
	integer :: sys_time_begin, sys_time_end, clock_rate	!Variables for timing the integration
    double precision :: computation_time				!Actual computation time
    logical :: jacobian_out 	!Whether or not to print out the Jacobian during solves
    character (len=file_path_length) :: output_profile_name, output_sweep_name, &
    	reconstruct_lznd_vdet_infile, reconstruct_lznd_vdet_outfile
    double precision :: sonic_loc	!Location where we hit the sonic point (or max ux/cs)
    logical :: gone_sonic			!Whether we hit a sonic point during the integration    
    double precision :: linearization_length, pathological_loc, max_mach, l_burn, q_tot
          
	contains
	
		subroutine read_inlist
			implicit none
			
			character(len=file_path_length) :: filename
			!Insanity - if I have exactly 84 elements (cutoff at cj_final_iso), then I get
			!a segfault during rates_init()...
			namelist /inlist_znd/ t0, rho0, &
				my_mesa_dir, net_file, burn_rho, burn_t, burn_u, burn_time, m_c, r_wd, &
				use_solar, num_isos_for_Xinit, names_of_isos_for_Xinit, values_for_Xinit, &
				num_steps, which_solver, do_constp, do_znd, do_cj, jacobian_out, &
				output_profile_name, ijac, do_blowout, h_scale, v_det, do_neumann, &
				output_sweep_name, do_sweep, sweep_znd_iters, uy_init, grav_init, &
				do_curvature, r_curve, do_blowout_local, delta_y_init, &
				do_h_search, h_search_bracket_low, h_search_bracket_high, h_search_numsteps, &
				h_search_tol, h_search_xmax, h_search_vdet, do_lznd_vdet, sweep_min_mach, &
				sweep_max_mach, do_reconstruct_lznd_vdet, reconstruct_lznd_vdet_infile, &
				reconstruct_lznd_vdet_outfile, &
				do_vdet_search, vdet_search_bracket_low, vdet_search_bracket_high, &
				vdet_search_numsteps, vdet_search_tol, vdet_search_xmax, &
				vdet_search_h_scale, l_burn_factor, use_rc_hs_scaling, rc_hs_factor, &
				j_curve, use_uy_ux_ratio, uy_ux_ratio, use_he_clavin, use_const_uy, &
				do_rho_sweep, d_hs_scaling, rho_sweep_T0, rho_sweep_grav, &
				rho_sweep_log_min, rho_sweep_log_max, rho_sweep_num_steps, do_pathological, &
				sonic_limit, sonic_limit_pathological, pathological_linearization_ratio, &
				do_lznd_wd, shock_lead_hp_frac, m_env_min, m_env_max, env_rho_min, &
				env_rho_max, wd_env_use_le, n_le, use_variable_linearization_ratio, &
				wd_env_use_mesa, cj_final_iso, &
				num_isos_for_Xcj, names_of_isos_for_Xcj, values_for_Xcj, &
				max_steps, rtol_init, atol_init, do_cjstate_only, do_cjstate_only2
				
			filename = 'inlist_znd'
			dbg = .false.
      
      		! set defaults:
      		!-----------------------------------------------------------------------------
      		my_mesa_dir = ''
      		!data_dir = '/Users/Kevin/mesa/data'	!MESA data location
    		net_file = 'basic.net'				!Which nuclear net to use
    		output_profile_name='./diagnosis_runs/rho0_1d5_t0_1d8-he4_to_c12.data'
    		cj_final_iso = 'ni56'
    		
    		rho0 = 1d6
    		t0 = 1d8
    		m_c = msun
    		r_wd = r_earth
    		grav_init = 0d8		!cm/s^2
    		
    		!burn_rho = 1d2		!g/cc
    		!burn_t = 2d7		!K
    		burn_time = 1d7		!sec
    		!burn_u = 3.75d8	!cm/s (1.5d9/4)
    		
    		!l_burn = post-shock distance over which 95% of energy is released
    		l_burn_factor = 0.95
    		
    		!Whether or not to use the analytic Jacobian in the ZND integrator:
    		ijac = 0			!0: finite differences, 1: analytic
    		
    		!Equally spaced in log(x) for sweeping through position to make a profile for
    		!a single detonation:
    		num_steps = 1000
    		
    		do_cjstate_only = .false.
    		!Whether or not to do a const pressure burn (as opposed to one-zone):
    		do_constp = .false.  
    		!Whether or not to use the ZND jacobian and solve for post-shock structure
    		do_znd = .true.
    		!Whether or not to integrate using the calculated CJ velocity
    		do_cj = .true.
    		!Whether to use 1D blowout equations in ZND calculation
    		do_blowout = .false. 
    		h_scale = 1d8		!Scale height for blowout calculation (in cm)
    		uy_init = 0d0		!Initial vertical velocity in blowout calculation (cm/s)
    		!Whether to use a user-specified detonation velocity & integrate from the Neumann pt.
    		do_neumann = .false.
    		v_det = 1.52d9 		!User-specified detonation velocity
    		!Whether to use ZND equations assuming a steady curvature term:
    		do_curvature = .false.
    		r_curve = 1d8		!Constant radius of curvature to use
    		!Whether or not to do blowout (using the local prescription for expansion):
    		do_blowout_local = .false.
    		delta_y_init = 1d0	!cm
    		
    		!Whether to perform a sweep over detonation velocities (outputting sonic pts.)
    		do_sweep = .false.
    		output_sweep_name='./blowout_runs/vdet_sweeps/r1d5_t1d8_he90_c10_hs1d6.data'
    		!Number of steps in v_det when do_sweep = .true.
    		sweep_znd_iters = 10
    		
    		!Whether or not to do a search for H_scale value in blowout, given a v_det
    		do_h_search = .false.
    		h_search_vdet = 1.05d9		!cm/s
    		h_search_bracket_low = 1d5	!cm
    		h_search_bracket_high = 1d9	!cm
    		h_search_numsteps = 100		
    		h_search_tol = 1d-2
    		h_search_xmax = 1d10		!cm
    		
    		wd_env_use_le = .false.
    		n_le = 1.5		!Polytrope index, P = k_le*rho**(1+1/n_le)
    		
    		wd_env_use_mesa = .true.
    		
    		jacobian_out = .false.
    		
    		which_solver = ros3pl_solver
    		
    		!Abundance initialization will take place in initialize()
    		!-----------------------------------------------------------------------------
      		
      		open(unit=10, file=trim(filename), action='read', delim='quote', iostat=ierr)
      		if (ierr /= 0) then
         		write(*, *) 'Failed to open control namelist file ', trim(filename)
         		stop 1
      		else
         		read(10, nml=inlist_znd, iostat=ierr)  
         		close(10)
         		write(*, nml=inlist_znd)
         		if (ierr /= 0) then
            		write(*, *) 'Failed while trying to read control namelist file ',&
            			trim(filename)
            		write(*, '(a)') &
            			'The following runtime error message might help you find the problem'
            		write(*, *) 
            		open(unit=10, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            		read(10, nml=inlist_znd)
            		close(10)
            		stop 1
         		end if  
      		end if
      		
      		do_pathological_init = do_pathological
      		
      		!Check a few flag cases:
      		if(do_znd.and.do_constp) then
      			write(*,*) 'Can only choose either do_znd = .true. or do_constp = .true'
      			stop 1
      		endif
      		     		
		end subroutine read_inlist  
      
		subroutine initialize
			use const_lib, only: const_init
			use chem_lib, only: chem_init
			!use reaclib_lib, only: reaclib_init
			!use weak_lib, only: weak_init
			use rates_lib, only: rates_init
			!use kap_lib, only : kap_init
			use net_lib, only : net_init
			
			!Opacity options
			character (len=256) :: kappa_file_prefix = 'gs98' 	!or eg. 'gn93' or 'gs98' or 'gn93_af94'
			character (len=256) :: kappa_CO_prefix = 'gn93_co'
			character (len=256) :: kappa_lowT_prefix = 'lowT_fa05_gs98' !or eg. 'lowT_Freedman11'
			double precision :: kappa_blend_logT_upper_bdy = 4.1d0
			double precision :: kappa_blend_logT_lower_bdy = 3.93d0
			double precision :: kappa_type2_logT_lower_bdy = 0 ! use default
			logical :: use_cache = .true.
         
          	!my_mesa_dir = '/Users/Kevin/mesa'         
         
         	!EOS options:
         	eos_file_prefix = 'mesa'
         	use_cache = .true.
         
			ierr = 0
			
			call crlibm_init()
			
			call const_init(my_mesa_dir,ierr)
			
			!Initialize my Lane-Emden solver:
			call lane_emden_init()	
			
			!Initialize my WD envelope solver
			call init_mesa_modules()
			call wd_init()
			
			call chem_init('isotopes.data', ierr)
			if (ierr /= 0) then
			   write(*,*) 'chem_init failed'
			   return
			end if
			
			call eos_init(eos_file_prefix, '', '', '', use_cache, info)
			if (info /= 0) then
				write(*,*) 'eos_init failed in Setup_eos'
				stop 1
			end if
		 
			write(*,*) 'loading eos tables...'
		 
			eos_handle = alloc_eos_handle(info)
			if (info /= 0) then
				write(*,*) 'failed trying to allocate eos handle'
				stop 1
			end if
			write(*,*) '1'
			
			!call reaclib_init(ierr)   
			!if (ierr /= 0) then
			!   write(*,*) 'reaclib_init failed'
			!   return
			!end if
			!write(*,*) '2'
			
			!call weak_init('', ierr)   
			!if (ierr /= 0) then
			!   write(*,*) 'weak_init failed'
			!   return
			!end if
			!write(*,*) '3'
			
			 call rates_init('reactions.list', '', .false., '', '', '', ierr)
			if (ierr /= 0) then
			   write(*,*) 'rates_init failed'
			   return
			end if
			write(*,*) '4'
			
			 call net_init(ierr)
			if (ierr /= 0) then
			   write(*,*) 'net_init failed'
			   return
			end if 
			write(*,*) '5'    
			  
			!call kap_init(kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix, &
			!	kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy, &
			!	kappa_type2_logT_lower_bdy, use_cache, '', ierr) 
			!if (ierr /= 0) then
			!   write(*,*) 'kap_init failed'
			!   return
			!end if 
      	
      end subroutine initialize
      
      
      subroutine setup_net
         !use net_lib
         !use rates_def, only: rates_reaction_id_max
         
         ierr = 0
         handle = alloc_net_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_net_handle failed'
            return
         end if
         
         call net_start_def(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_start_def failed'
            return
         end if
         
         write(*,*) 'load ' // trim(net_file)
         ! read_net_file first tries opening the filename in the current directory.
      	 ! if doesn't find that file, then tries the data_dir from the call on net_init.
     	 ! i.e., looks for <data_dir>/net_data/nets/<filename>
      	 ! check net_data/nets/README for info about net files.
         call read_net_file(net_file, handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'read_net_file failed ', trim(net_file)
            return
         end if
         
         call net_finish_def(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_finish_def failed'
            return
         end if
         
         !Turn off the high-T limits on reaction rates to make aprox13 work:
         !call net_turn_off_T_limits(handle, ierr)
   	
      	 allocate(which_rates(rates_reaction_id_max))
         which_rates(:) = which_rates_choice

         call net_set_which_rates(handle, which_rates, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_set_which_rate_f17pg failed'
            return
         end if
         
         call net_setup_tables(handle, 'rate_tables', '', ierr)
         if (ierr /= 0) then
            write(*,*) 'net_setup_tables failed'
            return
         end if
         
         species = net_num_isos(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in net_num_isos'
            return
         end if
         !Here we're adding to the dimension of the system
         !Use species + 3 as the default, and species + 5 to turn on solving the
         !dP/dx and dE/dx equations explicitly for consistency checks
         if(do_blowout) then
         	num_vars = species + 5	!Add scale height & vertical velocity
         else if (do_blowout_local) then
         	num_vars = species + 6  !Add thickness, vertical velocity, and delta_y
         else
         	num_vars = species + 3
         endif
         
         !If we're using my prescription of curvature then add a variable
         if(do_curvature.and.(.not.use_he_clavin)) then
         	num_vars = num_vars + 1
         endif
         
         allocate(rtol(num_vars), atol(num_vars))
         
         call get_chem_id_table_ptr(handle, chem_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_chem_id_table_ptr'
            return
         end if
         
         call get_net_iso_table_ptr(handle, net_iso, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_net_iso_table_ptr'
            return
         end if
         
         num_reactions = net_num_reactions(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in net_num_reactions'
            return
         end if

         lwork_net = net_work_size(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in net_work_size'
            return
         end if
         
         allocate(work_net(lwork_net))
         
         write(*,*) 'species:',species
         write(*,*) 'num_reactions:',num_reactions
         
         !ierr = 0
         !call get_net_ptr(handle, g, ierr)        
      end subroutine setup_net
      
   	!Subroutine to output the Hugoniot curve (in the P-V plane) given a final 
   	!composition (initial composition is taken to be what was specified in the inlist).
   	!Currently, the final composition is 100% whatever isotope final_species is
   	!(eg. 'ni56' or 'c12')
   	!Add a flag, is_initial_state to calculate the q=0 hugoniot (and don't use
   	!final_species) in this case since it will be a mixture in general
   	subroutine output_hugoniot(io, filename, is_initial_state, final_species)
   		implicit none
   		
   		integer, intent(in) :: io
   		character (len=256), intent(in) :: filename, final_species
   		logical, intent(in) :: is_initial_state

   		double precision :: p0, e0			!Initial values
   		double precision :: rho, t, p, e	!Shocked values
   		integer :: i, j
   		
   		open(unit=io,file=trim(filename))
   		
   		!First, set up all the variables needed for the Newton solver
   		!Newton begin-----------------------------------------------------------------
      		
         !Follow mesa/num/test/src/test_newton.f and only allocate the 1d variables (just keep the
         !others as pointers to them - no need to allocate/deallocate them as well)
      	!allocate(equ_newton(nvar_hug,nz), x_newton(nvar_hug,nz), xold_newton(nvar_hug,nz),&
      	!	dx_newton(nvar_hug,nz), xscale_newton(nvar_hug,nz), y_newton(ldy, nsec), stat=ierr)
			!if (ierr /= 0) stop 1
			
			allocate(equ1d(neq_hug), x1d(neq_hug), xold1d(neq_hug), dx1d(neq_hug), &
       			xscale1d(neq_hug), y1d(ldy*nsec), stat=ierr)
    		if (ierr /= 0) stop 1
    		x_newton(1:nvar_hug,1:nz) => x1d(1:neq_hug)
    		xold_newton(1:nvar_hug,1:nz) => xold1d(1:neq_hug)
    		dx_newton(1:nvar_hug,1:nz) => dx1d(1:neq_hug)
    		equ_newton(1:nvar_hug,1:nz) => equ1d(1:neq_hug)
    		xscale_newton(1:nvar_hug,1:nz) => xscale1d(1:neq_hug)
    		y_newton(1:ldy,1:nsec) => y1d(1:ldy*nsec)
      		
			numerical_jacobian = .false.
			which_decsol = lapack
			call lapack_work_sizes(neq_hug, lrd_newton, lid_newton)
			allocate(rpar_decsol_newton(lrd_newton), ipar_decsol_newton(lid_newton), stat=ierr)
			if (ierr /= 0) stop 1
			first_step = .true.
			
			!doing_jacobian = .false.
			matrix_type = square_matrix_type
			
			call newton_work_sizes(m1_hug, m2_hug, nvar_hug, nz, nsec, matrix_type, &
				lwork_newton, liwork_newton, ierr)
			if (ierr /= 0) stop 1
			
			allocate(work_newton(lwork_newton), iwork_newton(liwork_newton), stat=ierr)
			if (ierr /= 0) stop 1
			
			work_newton = 0
			iwork_newton = 0
			!qwork_newton = 0
			
			!iwork_newton(i_try_really_hard) = 1 ! 1: try really hard for first model
			iwork_newton(i_model_number) = 1
			iwork_newton(i_max_tries) = 10000
			iwork_newton(i_itermin) = 2
			!iwork_newton(i_max_iterations_for_jacobian) = 2
			!iwork_newton(i_debug) = 1
			
			! upper limit on magnitude of average scaled correction:
			tol_correction_norm = 1d99 ! lower this to decrease our accuracy
			work_newton(r_tol_max_correction) = 1d-6
			work_newton(r_tol_residual_norm) = 1d-6
			epsder = 1d-12 ! relative variation to compute derivatives
			
			!tol_residual_norm = 1d1
			!tol_max_correction = 1d99
			!Newton end -------------------------------------------------------------------
			
			!Initial composition binding energy:
			q_init = avo*mev_to_ergs*sum(xa*chem_isos% binding_energy(chem_id(1:species))/&
				chem_isos% W(chem_id(1:species)))
			
			if(.not.is_initial_state) then
				!Final composition:
				allocate(xa_hug(species))
				xa_hug = 0d0
				num_isos_for_Xinit = 1
				values_for_Xinit = 0d0
				!names_of_isos_for_Xinit(1:num_isos_for_Xinit) = (/'ni56'/)
				names_of_isos_for_Xinit(1:num_isos_for_Xinit) = (/trim(final_species)/)
				values_for_Xinit(1:num_isos_for_Xinit) = (/1d0/)
				values_for_Xinit = values_for_Xinit/sum(values_for_Xinit)	!Normalize
				q_hug = 0d0
				
				write(*,*) 'Final composition in Hugoniot calculation:'
				do i=1,num_isos_for_Xinit
					j = get_nuclide_index(names_of_isos_for_Xinit(i))
					write(*,*) 'Index in chem_isos: ',j
					write(*,*) names_of_isos_for_Xinit(i), values_for_Xinit(i)
					write(*,*) 'check index (name should match chem_isos result above): ',&
						chem_isos% name(j)
					
					xa_hug(net_iso(j)) = values_for_Xinit(i)
					q_hug = q_hug + avo*mev_to_ergs*xa_hug(net_iso(j))*&
						chem_isos% binding_energy(j)/chem_isos% W(j)
				end do
   			write(*,*)
				
   			q_hug = q_hug - q_init
      	else
      		allocate(xa_hug(species))
      		xa_hug = xa
      		q_hug = 0d0
      	endif
         	
      	write(*,*) 'q_hug:', q_hug
			
			!Loop over the densities to get the Hugoniot curve
			do i=1,1000
			
			xrho_hug = 1d0+(4.8d0*i/1000)
			xt_hug = 3d9/t0

			xold_newton(:,1) = (/ xt_hug /) ! starting model (xt)
			!write(*,*)
			!write(*,'(a15,3ES15.5)') 'guess:', xold_newton(:,1)
			dx_newton(:,1) = 0d0
			!dx_newton(:,1) = (/ rho_hug/10, t_hug/10 /)
			x_newton = xold_newton + dx_newton
			
			AF1 => null()

			call newton(&
				nz, nvar_hug, x1d, xold1d, &
				matrix_type, m1_hug, m2_hug, lapack_decsol, null_decsolblk, &
				lrd_newton, rpar_decsol_newton, lid_newton, ipar_decsol_newton, which_decsol, &
				tol_correction_norm, &
				set_primaries_hug, set_secondaries_hug, default_set_xscale, &
				default_Bdomain, default_xdomain, eval_equations_hug, &
				size_equ, sizeB, default_inspectB, &
				enter_setmatrix_hug, exit_setmatrix, failed_in_setmatrix, default_force_another_iter, &
				xscale1d, equ1d, ldy, nsec, y1d, work_newton, lwork_newton, iwork_newton, &
				liwork_newton, AF1, &
				lrpar, rpar, lipar, ipar, &
				nonconv, ierr)
			if (ierr /= 0) then
				write(*,*) 'ierr = ',ierr
				stop 50
			endif
			if (nonconv) then
				write(*,*) 'hug solver failed to converge'
				write(*,*) 'x_newton: ', x_newton
				stop 2
			end if
			
			rho = xrho_hug*rho0
			t = x_newton(1,1)*t0
			
			logRho = log10(rho)
         	logT = log10(t)
         	call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa_hug, &
        		rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        		d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
			e = exp(res(i_lnE))
			p = exp(res(i_lnPgas)) + 1/3d0*crad*t**4
			
			write(io,*) rho, t, p
			end do
			
			close(io)
			if(associated(AF1)) deallocate(AF1)
			deallocate(equ1d, x1d, xold1d, dx1d, xscale1d, y1d, stat=ierr)	
			deallocate(work_newton, iwork_newton, stat=ierr)
			!deallocate(equ_newton, x_newton, xold_newton, dx_newton, &
			!	xscale_newton, y_newton, stat=ierr)
			deallocate(rpar_decsol_newton, ipar_decsol_newton, stat=ierr)
			deallocate(xa_hug)
      end subroutine output_hugoniot
      
      	!Subroutine to call the newton solver and find the CJ conditions (density, 
      	!temperature, and detonation velocity) given the ambient conditions and a 
      	!composition (this currently assumes a final state of pure 56Ni)
      	!We'll already have created the composition vector xa(:)
		subroutine solve_cj
			implicit none
			
			double precision :: p0, e0, logRho, logT, mu0
			integer :: i,j 								!Loop variables
		
			if(.not.associated(xa)) then
				write(*,*) 'Composition vector xa(:) not initialized before CJ call'
				write(*,*) 'Try again with a valid initial compostion'
				stop 1
			endif
		
			!First, set up all the variables needed for the Newton solver
      		!Newton begin-----------------------------------------------------------------
      		
         !Follow mesa/num/test/src/test_newton.f and only allocate the 1d variables (just keep the
         !others as pointers to them - no need to allocate/deallocate them as well)
      	!allocate(equ_newton(nvar,nz), x_newton(nvar,nz), xold_newton(nvar,nz),&
      	!		dx_newton(nvar,nz), xscale_newton(nvar,nz), y_newton(ldy, nsec), stat=ierr)
			!if (ierr /= 0) stop 1
			
			allocate(equ1d(neq), x1d(neq), xold1d(neq), dx1d(neq), &
       			xscale1d(neq), y1d(ldy*nsec), stat=ierr)
    		if (ierr /= 0) stop 1
    		x_newton(1:nvar,1:nz) => x1d(1:neq)
    		xold_newton(1:nvar,1:nz) => xold1d(1:neq)
    		dx_newton(1:nvar,1:nz) => dx1d(1:neq)
    		equ_newton(1:nvar,1:nz) => equ1d(1:neq)
    		xscale_newton(1:nvar,1:nz) => xscale1d(1:neq)
    		y_newton(1:ldy,1:nsec) => y1d(1:ldy*nsec)
      		
			numerical_jacobian = .false.
			which_decsol = lapack
			call lapack_work_sizes(neq, lrd_newton, lid_newton)
			allocate(rpar_decsol_newton(lrd_newton), ipar_decsol_newton(lid_newton), stat=ierr)
			if (ierr /= 0) stop 1
			first_step = .true.
			
			!doing_jacobian = .false.
			matrix_type = square_matrix_type
			
			call newton_work_sizes(m1, m2, nvar, nz, nsec, matrix_type, &
				lwork_newton, liwork_newton, ierr)
			if (ierr /= 0) stop 1
			
			allocate(work_newton(lwork_newton), iwork_newton(liwork_newton), stat=ierr)
			if (ierr /= 0) stop 1
			
			work_newton = 0
			iwork_newton = 0
			!qwork_newton = 0
			
			!iwork_newton(i_try_really_hard) = 1 ! 1: try really hard for first model
			iwork_newton(i_model_number) = 1
			iwork_newton(i_max_tries) = 500
			iwork_newton(i_itermin) = 2
			!iwork_newton(i_max_iterations_for_jacobian) = 5
			!iwork_newton(i_debug) = 1
			
			! upper limit on magnitude of average scaled correction:
			tol_correction_norm = 1d99 ! lower this to decrease our accuracy
			work_newton(r_tol_max_correction) = 1d-4
			work_newton(r_tol_residual_norm) = 1d-4
			epsder = 1d-12 ! relative variation to compute derivatives
			
			!tol_residual_norm = 1d1
			!tol_max_correction = 1d99
			!Newton end -------------------------------------------------------------------
			
			!Next, set up the final composition and figure out the total q released in
			!a detonation that burns all the way to this composition
			allocate(xa_end(species))
			xa_end = 0d0
			
			!Initial composition binding energy:
			q_init = avo*mev_to_ergs*sum(xa*chem_isos% binding_energy(chem_id(1:species))/&
				chem_isos% W(chem_id(1:species)))
			!write(*,*) 'q_init',q_init
			
			!Final composition:
			!(old - hardcoded composition):
			!num_isos_for_Xcj = 1
			!values_for_Xcj = 0d0
			!names_of_isos_for_Xcj(1:num_isos_for_Xcj) = (/cj_final_iso/)
			!values_for_Xcj(1:num_isos_for_Xcj) = (/1d0/)
						
			values_for_Xcj = values_for_Xcj/sum(values_for_Xcj)	!Normalize
			q = 0d0		
			
			write(*,*) 'Final composition in CJ calculation:'
			do i=1,num_isos_for_Xcj
				j = get_nuclide_index(names_of_isos_for_Xcj(i))
				write(*,*) 'Index in chem_isos: ',j
				write(*,*) names_of_isos_for_Xcj(i), values_for_Xcj(i)
				write(*,*) 'check index (name should match chem_isos result above): ',&
					chem_isos% name(j)
				
         		xa_end(net_iso(j)) = values_for_Xcj(i)
         		
         		!This doesn't work for some reason, so use the next line after the do loop..
         		!q = avo*mev_to_ergs*xa_end(net_iso(j))*chem_isos% binding_energy(j)/&
         		!	chem_isos% W(j)
         	end do
         	q = avo*mev_to_ergs*sum(xa_end*chem_isos% binding_energy(chem_id(1:species))/&
				chem_isos% W(chem_id(1:species)))
         	!write(*,*) 'q',q
         	write(*,*)
         	
         	q = q - q_init
         	!write(*,*) 'q (erg/g): ',q
         	
         	logRho = log10(rho0)
         	logT = log10(T0)
         	
         	call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa, &
        		rho0, logRho, T0, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        		d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
			e0 = exp(res(i_lnE))
			p0 = exp(res(i_lnPgas)) + 1/3d0*crad*t0**4
			mu0 = 1d0/(1/abar + 1/zbar)

			!Output initial values:
			write(*,*) 'Initial values:'
			write(*,'(a7 ES12.4E2 a12)') 'rho0 =',rho0,'g/cm^3'
			write(*,'(a7 ES12.4E2 a12)') 'T0 =',T0,'K'
			write(*,'(a7 ES12.4E2 a12)') 'P0 =',P0,'dyne/cm^2'
			write(*,'(a7 ES12.4E2 a12)') 'mu0 =',mu0
			write(*,'(a7 ES12.4E2 a12)') 'q =',q,'erg/g'
			write(*,*)
         	
         	!Here's the initial conditions for our nondimensional full EOS equations in (xrho, xt, xd):
			xrho = 2d0		!rho/rho0
			xt = 5d9/t0			!T/T0
			xd = 100		!D^2/(P0/rho0)
			xold_newton(:,1) = (/ xrho, xt, xd /) ! starting model (xrho, xt, xd)
			!write(*,'(a15,3ES15.5)') 'guess:', xold_newton(:,1)
			dx_newton(:,1) = 0d0
			!dx(:,1) = (/ -xrho/1d1, xp/1d1,  xd/1d1 /)
			x_newton = xold_newton + dx_newton
			
			AF1 => null()
	
			call newton(&
				nz, nvar, x1d, xold1d, &
				matrix_type, m1, m2, lapack_decsol, null_decsolblk, &
				lrd_newton, rpar_decsol_newton, lid_newton, ipar_decsol_newton, which_decsol, &
				tol_correction_norm, &
				set_primaries, set_secondaries, default_set_xscale, &
				default_Bdomain, default_xdomain, eval_equations, &
				size_equ, sizeB, default_inspectB, &
				enter_setmatrix, exit_setmatrix, failed_in_setmatrix, default_force_another_iter, &
				xscale1d, equ1d, ldy, nsec, y1d, work_newton, lwork_newton, iwork_newton, &
				liwork_newton, AF1, &
				lrpar, rpar, lipar, ipar, &
				nonconv, ierr)
			if (ierr /= 0) then
				write(*,*) 'ierr = ',ierr
				stop 50
			endif
			if (nonconv) then
				write(*, *) 'CJ solver failed to converge'
				stop 2
			end if
			
			rho_cj = rho0*x_newton(1,1)
			t_cj = t0*x_newton(2,1)
			v_cj = sqrt(x_newton(3,1)*P0/rho0)
			
			logRho = log10(rho_cj)
         	logT = log10(T_cj)
         	call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa, &
        		rho_cj, logRho, T_cj, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        		d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
        	e_cj = exp(res(i_lnE))
        	p_cj = exp(res(i_lnPgas)) + 1/3d0*crad*t_cj**4
			
			write(*,'(a15,3ES15.5)') 'CJ solution:', x_newton
			write(*,'(a15 ES12.4E2)') 'rho_cj = ',rho_cj
			write(*,'(a15 ES12.4E2)') 'T_cj = ',t_cj
			write(*,'(a15 ES12.4E2)') 'P_cj = ',p_cj
			write(*,'(a15 ES12.4E2)') 'D_cj = ',v_cj
			write(*,*)
			
			xrho = x_newton(1,1)
			xt = x_newton(2,1)
			xd = x_newton(3,1)
			
			if(associated(AF1)) deallocate(AF1)
			deallocate(xa_end, equ1d, x1d, xold1d, dx1d, xscale1d, y1d, stat=ierr)	
			deallocate(work_newton, iwork_newton, stat=ierr)
			!deallocate(equ_newton, x_newton, xold_newton, dx_newton, &
			!	xscale_newton, y_newton, stat=ierr)
			deallocate(rpar_decsol_newton, ipar_decsol_newton, stat=ierr)
      end subroutine solve_cj
      	
      	!Subroutine to call the newton solver and find the CJ conditions (density, 
      	!temperature, and detonation velocity) given the ambient conditions and a 
      	!composition (this currently assumes a final state of pure 56Ni)
      	!We'll already have created the composition vector xa(:)
		subroutine solve_neumann
			implicit none

			double precision :: p0, e0, logRho, logT, mu0
			integer :: i,j 								!Loop variables
		
			if(.not.associated(xa)) then
				write(*,*) 'Composition vector xa(:) not initialized before CJ call'
				write(*,*) 'Try again with a valid initial compostion'
				stop 1
			endif
		
			!First, set up all the variables needed for the Newton solver
      		!Newton begin-----------------------------------------------------------------
      		
         !Follow mesa/num/test/src/test_newton.f and only allocate the 1d variables (just keep the
         !others as pointers to them - no need to allocate/deallocate them as well)
      	!allocate(equ_newton(nvar_neumann,nz), x_newton(nvar_neumann,nz), xold_newton(nvar_neumann,nz),&
      	!	dx_newton(nvar_neumann,nz), xscale_newton(nvar_neumann,nz), y_newton(ldy, nsec), stat=ierr)
			!if (ierr /= 0) stop 1
			
			allocate(equ1d(neq_neumann), x1d(neq_neumann), xold1d(neq_neumann), dx1d(neq_neumann), &
       			xscale1d(neq_neumann), y1d(ldy*nsec), stat=ierr)
    		if (ierr /= 0) stop 1
    		x_newton(1:nvar_neumann,1:nz) => x1d(1:neq_neumann)
    		xold_newton(1:nvar_neumann,1:nz) => xold1d(1:neq_neumann)
    		dx_newton(1:nvar_neumann,1:nz) => dx1d(1:neq_neumann)
    		equ_newton(1:nvar_neumann,1:nz) => equ1d(1:neq_neumann)
    		xscale_newton(1:nvar_neumann,1:nz) => xscale1d(1:neq_neumann)
    		y_newton(1:ldy,1:nsec) => y1d(1:ldy*nsec)
      		
			numerical_jacobian = .false.
			which_decsol = lapack
			call lapack_work_sizes(neq_neumann, lrd_newton, lid_newton)
			allocate(rpar_decsol_newton(lrd_newton), ipar_decsol_newton(lid_newton), stat=ierr)
			if (ierr /= 0) stop 1
			first_step = .true.
			
			!doing_jacobian = .false.
			matrix_type = square_matrix_type
			
			call newton_work_sizes(m1_neumann, m2_neumann, nvar_neumann, nz, nsec, matrix_type, &
				lwork_newton, liwork_newton, ierr)
			if (ierr /= 0) stop 1
			
			allocate(work_newton(lwork_newton), iwork_newton(liwork_newton), stat=ierr)
			if (ierr /= 0) stop 1
			
			work_newton = 0
			iwork_newton = 0
			!qwork_newton = 0
			
			!iwork_newton(i_try_really_hard) = 1 ! 1: try really hard for first model
			iwork_newton(i_model_number) = 1
			iwork_newton(i_max_tries) = 10000
			iwork_newton(i_itermin) = 2
			!iwork_newton(i_max_iterations_for_jacobian) = 2
			!iwork_newton(i_debug) = 1
			
			! upper limit on magnitude of average scaled correction:
			tol_correction_norm = 1d99 ! lower this to increase our precision
			work_newton(r_tol_max_correction) = 1d-4
			work_newton(r_tol_residual_norm) = 1d-4
			epsder = 1d-12 ! relative variation to compute derivatives
			
			!tol_residual_norm = 1d1
			!tol_max_correction = 1d99
			!Newton end -------------------------------------------------------------------
			
			!Now that the CJ conditions are found, we can also get the Neumann conditions
			!to start the integration at:
			call composition_info(species, chem_id, xa, xh, xhe, zm, abar, zbar, z2bar, &
            	ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
			logRho = log10(rho0)
         	logT = log10(T0)
         	call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa, &
        		rho0, logRho, T0, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        		d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
			e0 = exp(res(i_lnE))
			p0 = exp(res(i_lnPgas)) + 1/3d0*crad*t0**4
			mu0 = 1d0/(1/abar + 1/zbar)
			
			rho_neumann = 4d0
			t_neumann = 3d9/t0
			!If we're calculating the CJ detonation velocity or letting the user
			!specify it. If do_cj is .false. then v_det has already been set in the inlist
			if(do_cj) then
				v_det = sqrt(xd*P0/rho0)		!CJ detonation velocity here
			endif
			
			xold_newton(:,1) = (/ rho_neumann, t_neumann /) ! starting model (xrho, xt)
			!write(*,*)
			!write(*,'(a15,3ES15.5)') 'guess:', xold_newton(:,1)
			dx_newton(:,1) = 0d0
			!dx_newton(:,1) = (/ rho_neumann/10, t_neumann/10 /)
			x_newton = xold_newton + dx_newton
			
			AF1 => null()

			call newton(&
				nz, nvar_neumann, x1d, xold1d, &
				matrix_type, m1_neumann, m2_neumann, lapack_decsol, null_decsolblk, &
				lrd_newton, rpar_decsol_newton, lid_newton, ipar_decsol_newton, which_decsol, &
				tol_correction_norm, &
				set_primaries_neumann, set_secondaries_neumann, default_set_xscale, &
				default_Bdomain, xdomain_neumann, eval_equations_neumann, &
				size_equ, sizeB, default_inspectB, &
				enter_setmatrix_neumann, exit_setmatrix, failed_in_setmatrix, default_force_another_iter, &
				xscale1d, equ1d, ldy, nsec, y1d, work_newton, lwork_newton, iwork_newton, &
				liwork_newton, AF1, &
				lrpar, rpar, lipar, ipar, &
				nonconv, ierr)
			if (ierr /= 0) then
				write(*,*) 'ierr = ',ierr
				stop 50
			endif
			if (nonconv) then
				write(*, *) 'Neumann solver failed to converge'
				write(*, *) x_newton(:,1)
				stop 2
			end if
			
			write(*,'(a20,3ES15.5)') 'Neumann solution:', x_newton
			rho_neumann = x_newton(1,1)*rho0
			t_neumann = x_newton(2,1)*t0
			logRho = log10(rho_neumann)
         logT = log10(T_neumann)
         call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa, &
        		rho_neumann, logRho, T_neumann, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        		d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
        	p_neumann = exp(res(i_lnPgas)) + 1/3d0*crad*t_neumann**4
        	e_neumann = exp(res(i_lnE))
        	cs_neumann = sqrt(res(i_gamma1)*p_neumann/rho_neumann)
			
			if(associated(AF1)) deallocate(AF1)
			deallocate(equ1d, x1d, xold1d, dx1d, xscale1d, y1d, stat=ierr)	
         if (ierr /= 0) then
            write(*,*) 'Deallocation error 1 in solve_neumann()'
            stop 1
         endif
			deallocate(work_newton, iwork_newton, stat=ierr)
         if (ierr /= 0) then
            write(*,*) 'Deallocation error 2 in solve_neumann()'
            stop 1
         endif
			!deallocate(equ_newton, x_newton, xold_newton, dx_newton, &
			!	xscale_newton, y_newton, stat=ierr)
			deallocate(rpar_decsol_newton, ipar_decsol_newton, stat=ierr)	
         if (ierr /= 0) then
            write(*,*) 'Deallocation error 3 in solve_neumann()'
            stop 1
         endif  		
      end subroutine solve_neumann
      	
      	!Begin - Subroutines needed by the newton solver------------------------------------------------
      	!Subroutines for finding the CJ conditions:
    	subroutine set_primaries(nvar, nz, x, lrpar, rpar, lipar, ipar, ierr)
        	integer, intent(in) :: nvar, nz
         	double precision, pointer :: x(:,:)
         	integer, intent(in) :: lrpar, lipar
         	double precision, intent(inout) :: rpar(:)
         	integer, intent(inout) :: ipar(:)
         	integer, intent(out) :: ierr
         	ierr = 0
         	x0 = x(1, 1); x1 = x(2, 1); x2 = x(3, 1) !x0 = rho, x1 = T
         	!write(*, '(a20, 3(a6, 1pe26.16, 3x))') 'primaries', 'x0', x0, 'x1', x1, 'x2', x2
         	
         	if (dbg) write(*, '(a20, 3(a6, 1pe26.16, 3x))') 'primaries', 'x0', x0, 'x1', x1, 'x2', x2
		end subroutine set_primaries
		
		
		subroutine set_secondaries(ivar, lrpar, rpar, lipar, ipar, ierr)
         	integer, intent(in) :: ivar
         	integer, intent(in) :: lrpar, lipar
         	double precision, intent(inout) :: rpar(:)
         	integer, intent(inout) :: ipar(:)
         	integer, intent(out) :: ierr
         	logical, parameter :: skip_partials = .true.
         	call set_sec(ivar, skip_partials, lrpar, rpar, lipar, ipar, ierr)
      	end subroutine set_secondaries
      	
      
      	subroutine set_sec(ivar, skip_partials, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only:pi
         integer, intent(in) :: ivar
         logical, intent(in) :: skip_partials
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(lrpar)
		 integer, intent(inout) :: ipar(lipar)
		 integer, intent(out) :: ierr
		 ierr = 0
		 
		 !No secondaries in this example
      end subroutine set_sec
      
      !For finding the CJ conditions:
      subroutine eval_equations(iter, nvar, nz, x, xscale, equ, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: iter, nvar, nz
         double precision, pointer, dimension(:,:) :: x, xscale, equ
         !double precision, intent(out) :: equ(nvar, nz)
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
      	double precision, dimension(nvar*nz, nvar*nz) :: A ! square matrix for jacobian
			logical, parameter :: skip_partials = .true.			
			call eval_equ(nvar, nz, equ, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)         
      end subroutine eval_equations
      
		!For finding the CJ conditions:
      subroutine eval_equ(nvar, nz, equ, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)
      	implicit none
      	
        integer, intent(in) :: nvar, nz
        double precision, intent(out) :: equ(nvar, nz)
		logical, intent(in) :: skip_partials
      	double precision, dimension(nvar*nz, nvar*nz) :: A 
        integer, intent(in) :: lrpar, lipar
        double precision, intent(inout) :: rpar(:)
        integer, intent(inout) :: ipar(:)
        integer, intent(out) :: ierr  
        integer :: i,j
        double precision :: y0, y1, y2
        !Overwrite the globals here since we're replacing them with the (changing) primaries
        double precision :: e0, e1, p0, p1
        double precision :: rho, t, d, chi_rho, chi_t, dp_dv_hug, dt_drho_hug, xrho, xp, xd, xt
        !First derivatives that we'll need:
        double precision :: dp_drho, de_drho, dp_dT, de_dT
        !Second derivatives that we'll need:
        double precision :: d2p_drhoT, d2p_drho2, d2p_dT2, d2e_drhoT, d2e_drho2, d2e_dT2
        
        ierr = 0
        call set_sec(0, skip_partials, lrpar, rpar, lipar, ipar, ierr); if (ierr /= 0) return
		
		!write(*,*) 'Inside CJ solver'
    	
    	!Override MESA EOS with gamma-law stuff here ----------------------------------------------
    	!Need this part for either choice:
    	!call get_xyz(mass_fracs_init)
		!Y_mass = mass_fracs_init/nuclides% W
		!call nuclides_composition_info(y_mass, nuclides, abar, zbar, z2bar, ye)
    	!mu0 = 1d0/(1d0/abar + 1d0/zbar)
    	!Y_mass = mass_fracs_final/nuclides% W
		!call nuclides_composition_info(y_mass, nuclides, abar, zbar, z2bar, ye)
    	!mu = 1d0/(1d0/abar + 1d0/zbar)
    	
    	!Ideal gas:
    	!gamma1 = 5d0/3d0
    	!p0 = 1d0/V0*kerg*T0/(mu0*mp) 
    	!e0 = kerg*T0/(mu0*mp*(gamma1-1))
    	!p1 = 1d0/V*kerg*T/(mu*mp)
    	!e1 = kerg*T/(mu*mp*(gamma1-1))
    	!dE_dT = e1/T
    	!dE_drho = 0d0
    	!dP_dT = P1/T
    	!dP_drho = P1*V
    	!d2e_drho2 = 0d0
    	!d2e_drhoT = 0d0
    	
    	!Radiation dominated:
    	!gamma1 = 4d0/3d0
    	!p0 = crad*t0**4/3d0
    	!e0 = crad*t0**4 
    	!p1 = crad*t**4/3d0
    	!e1 = crad*t**4 
    	!dE_dT = 4*e1/T
    	!dE_drho = 0d0
    	!dP_dT = 4*P1/T
    	!dP_drho = 0d0
    	!d2e_drho2 = 0d0
    	!d2e_drhoT = 0d0
    	!-----------------------------------------------------------------------------------------
         
        !Begin dimensionless ideal gas equations----------------------------------------------------
        !gamma1 = 5d0/3d0
        !xrho = x0
        !xp = x1
        !xd = x2
        
        !y0 = (xp/xrho-1)/(gamma1-1) - 0.5*(1+xp)*(1-1d0/xrho) - rho0*q/P0
        !y1 = (xp-1)/(1-1d0/xrho) - xd
        !y2 = (1+4*xp)/(1-4d0/xrho) + (xp-1)/(1-1d0/xrho)
        !End dimensionless ideal gas equations------------------------------------------------------
        
        !Dimensionless full EOS equations-----------------------------------------------------------
        xrho = x0
        xt = x1
        xd = x2
        
        !First evaluate the EOS at the initial (rho0, t0) point to get p0:
        logrho = log10(rho0)
        logt = log10(t0)
        call composition_info(species, chem_id, xa, xh, xhe, zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
        call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa, &
        	rho0, logRho, T0, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        	d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
		e0 = exp(res(i_lnE))
		p0 = exp(res(i_lnPgas)) + 1/3d0*crad*t0**4
        
        !Now evaluate the EOS at the current (rho, t) point to get p and necessary derivs:
        rho = xrho*rho0
        t = xt*t0
        logrho = log10(rho)
        logt = log10(t)
        call composition_info(species, chem_id, xa_end, xh, xhe, zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
         call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa_end, &
        	rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        	d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
		e1 = exp(res(i_lnE))
		p1 = exp(res(i_lnPgas)) + 1/3d0*crad*t**4
		xp = p1/p0
		
		de_drho = res(i_dE_dRho)
		de_dt = res(i_Cv)
		chi_rho = res(i_chiRho)
		dp_drho = p1/rho*chi_Rho
		chi_t = res(i_chiT)
		dp_dt = p1/T*chi_T
		
		call eos_get_helm_results( Xh, abar, zbar, Rho, logRho, &
    			T, logT, .true., .false., .false., res_helm, ierr)
    	d2e_drho2 = res_helm(h_dedd)
    	d2e_drhoT = res_helm(h_dedt)
    	d2e_dt2 = res_helm(h_dett)
    	d2p_drho2 = res_helm(h_dpdd)
    	d2p_drhoT = res_helm(h_dpdt)
    	d2p_dt2 = res_helm(h_dptt)
    	
    	!Try alternative form for the second derivatives (use d_dlnRho_const_T, d_dlnT_const_Rho)
    	!d2e_drho2 = d_dlnRho_const_T(i_dE_dRho)/rho
    	!d2e_drhoT = d_dlnT_const_Rho(i_dE_dRho)/t
    	!d2e_dt2 = d_dlnT_const_Rho(i_Cv)/t
    	!d2p_drho2 = p1/rho**2*(chi_rho**2 + d_dlnRho_const_T(i_chiRho) - chi_rho)
    	!d2p_drhoT = p1/(rho*t)*(chi_t*chi_rho + d_dlnT_const_Rho(i_chiRho))
    	!d2p_dt2 = p1/t**2*(chi_t**2 + d_dlnT_const_Rho(i_chiT) - chi_T)
    	
    	!Wow, the numerical jacobian is MORE accurate than the analytic one here (perhaps because
    	!the EOS interpolation is coarser than recomputing the derivatives from finite differences?)
    	
    	y0 = (1+xp)*(1-1/xrho) + 2*rho0*(q+e0-e1)/p0
        y1 = (xp-1)/(1/xrho-1) + xd
        y2 = (xrho-1)*dp_dt*t0/p0*(rho0**2*xrho**2*de_drho - p0*xp) + rho0*xrho*de_dt*t0*&
        	(-(xrho-1)*xrho*dp_drho*rho0/p0 + xp - 1)	
        	
        !End dimensionless full EOS equations-------------------------------------------------------
         
        !write(*, '(a20, 3(a6, 1pe26.16, 3x))') 'primaries', 'x0', x0, 'x1', x1, 'x2', x2
  		!write(*,*) 'allo!',de_drho,de_dt,dp_drho,dp_dt
  	
         equ(:, 1) = (/ y0, y1, y2 /)
         !write(*,*) x0, x1, x2
         !write(*,*) equ(:, 1)
         !write(*,*) d2e_drho2, d2e_drhoT, d2e_dt2, d2p_drho2, d2p_drhoT, d2p_dt2
         !write(*,*)
         
         if (.not. skip_partials) then
            ! A(q, v) = xscale(v) * partial of equation(q) wrt variable(v)
            
            !New (dimensionless) jacobian for ideal gas (checked on Mathematica)--------------------
            !Worked with dimensionless ideal gas equations defined above
            !A(1,1) = -(xp+1)/(2*xrho**2) - xp/(xrho**2*(gamma1-1))
            !A(1,2) = 0.5*(1/xrho-1) + 1/(xrho*(gamma1-1))
            !A(1,3) = 0d0
            
            !A(2,1) = (1-xp)/((xrho-1)**2)
            !A(2,2) = 1/(1-1/xrho)
            !A(2,3) = -1d0

            !A(3,1) = -4*(1+4*xp)/((xrho-4)**2) - (xp-1)/((xrho-1)**2)
            !A(3,2) = 4/(1-4/xrho) + 1/(1-1/xrho)
            !A(3,3) = 0d0
            !End dimensionless jacobian for ideal gas-----------------------------------------------
            
            !Dimensionless jacobian for full EOS----------------------------------------------------
            A(1,1) = (1+xp)/(xrho**2) - 2*rho0**2*de_drho/p0 + (1-1/xrho)*xp*chi_rho/xrho
            A(1,2) = -2*rho0*t0*de_dt/p0 + (1-1/xrho)*xp*chi_t/xt
            A(1,3) = 0d0
            
            A(2,1) = (xp-1-(xrho-1)*xp*chi_rho)/(xrho-1)**2
            A(2,2) = xrho*xp*chi_t/(xt*(1-xrho))
            A(2,3) = 1d0
            
            A(3,1) = dp_dt*t0/p0*(rho0**2*xrho**2*de_drho - p0*xp) + &
            	(xrho-1)*d2p_drhoT*rho0*t0/p0*(rho0*xrho**2*de_drho*rho0 - p0*xp) + &
            	(xrho-1)*dp_dt*t0/p0*(rho0*xrho*(2*de_drho*rho0 + xrho*d2e_drho2*rho0**2) - &
            	p0*dp_drho*rho0/p0) + rho0*de_dt*t0*(-(xrho-1)*xrho*dp_drho*rho0/p0 + &
            	xp - 1) + rho0*xrho*d2e_drhoT*rho0*t0*(-(xrho-1)*xrho*dp_drho*rho0/p0 + &
            	xp - 1) - rho0*(xrho-1)*xrho*de_dt*t0*(2*dp_drho*rho0/p0 + &
            	xrho*d2p_drho2*rho0**2/p0)
            A(3,2) = ((xrho-1)*d2p_dt2*t0**2/p0*(rho0**2*xrho**2*de_drho - p0*xp) + &
            	(xrho-1)*dp_dt*t0/p0*(rho0*xrho**2*d2e_drhoT*rho0*t0 - p0*dp_dt*t0/p0) + &
            	rho0*xrho*d2e_dt2*t0**2*(-(xrho-1)*xrho*dp_drho*rho0/p0 + xp - 1) + &
            	rho0*xrho*de_dt*t0*(dp_dt*t0/p0 - (xrho-1)*xrho*d2p_drhoT*rho0*t0/p0))
            A(3,3) = 0d0
            !End dimensionless jacobian for full EOS------------------------------------------------
            
            A(:, 1) = A(:, 1)*xscale_newton(1,1)
            A(:, 2) = A(:, 2)*xscale_newton(2,1)
            A(:, 3) = A(:, 3)*xscale_newton(3,1)
            
            !write(*,*) A
         end if
      end subroutine eval_equ
      
      subroutine eval_jacobian(ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: ldA ! leading dimension of A
         double precision :: A(ldA, nvar*nz) ! the jacobian matrix
         ! A(idiag+q-v, v) = partial of equation(q) wrt variable(v)
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr         
			logical, parameter :: skip_partials = .false.	
			integer :: i,j
			
         ierr = 0        
		 call eval_equ(nvar, nz, equ_newton, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)
		
		 if(jacobian_out) then	
			 write(*, *)
			 write(*, *) 'result of analytic jacobian'
			 write(*, '(a10,3es20.12)') 'x = ',x_newton
			 write(*, '(a10,3es20.12)') 'xs = ',xscale_newton
			 write(*, '(a10,3es20.12)') 'equ = ',equ_newton
			 do j=1, nvar ! variable
				do i=1, nvar ! equation
				   write(*, '(e20.12, 3x)', advance = 'no') A(i, j)
				end do
				write(*, *)
			 end do
			 write(*, *)
		 endif
			
      end subroutine eval_jacobian
      

      subroutine enter_setmatrix(iter, nvar, nz, neqs, x, xold, xscale, xder, need_solver_to_eval_jacobian, &
            ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: iter, nvar, nz, neqs
         double precision, pointer, dimension(:,:) :: x, xold, xscale, xder
         !double precision, intent(out) :: xder(nvar, nz)
         logical, intent(out) :: need_solver_to_eval_jacobian
         integer, intent(in) :: ldA ! leading dimension of A
         double precision, pointer, dimension(:) :: A !(ldA, neqs) ! the jacobian matrix
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
         integer :: i, j
         if (dbg) write(*, '(/, a)') 'enter_setmatrix'
         
         if (numerical_jacobian) then
            xder=epsder*(xscale+abs(xold))
            need_solver_to_eval_jacobian = .true.
            doing_jacobian = .true.
         else
            call eval_jacobian(ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
            need_solver_to_eval_jacobian = .false.
         end if
      end subroutine enter_setmatrix
      

      subroutine exit_setmatrix(iter, nvar, nz, neqs, dx, ldA, A, idiag, xscale, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: ldA ! leading dimension of A
         integer, intent(in) :: iter, nvar, nz, neqs ! number of equations, 2nd dimension of A
         integer, intent(inout) :: idiag ! row of A with the matrix diagonal entries
         double precision, pointer, dimension(:,:) :: dx
          double precision, pointer, dimension(:) :: A
         !double precision, target :: A(ldA, neqs)
         double precision, pointer, dimension(:,:) :: xscale !(nvar, nz)
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
         integer :: i, j
         if (dbg) write(*, '(a, /)') 'exit_setmatrix'
         ierr = 0
         doing_jacobian = .false.
         
         if(jacobian_out) then
			 write(*, *)
			 write(*, *) 'result of numerical jacobian'
			 write(*, '(a10,3es20.12)') 'x = ',x_newton
			 write(*, '(a10,3es20.12)') 'xs = ',xscale
			 write(*, '(a10,3es20.12)') 'equ = ',equ_newton
			 do i=1, neqs ! equation
				do j=1, nvar ! variable
					write(*, '(e20.12, 3x)', advance = 'no') A((i-1)*nvar + j)
				end do
				write(*, *)
			 end do
			 write(*, *)
		 endif
         
         return
         stop 'exit_setmatrix'
      end subroutine exit_setmatrix
      

      subroutine failed_in_setmatrix(j, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: j
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
         if (dbg) write(*, '(a, /)') 'failed_in_setmatrix'
         ierr = 0
         doing_jacobian = .false.
      end subroutine failed_in_setmatrix
      
      
      subroutine xdomain(nvar, nz, x, dx, xold, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: nvar, nz
         double precision, dimension(nvar, nz), intent(inout) :: x, dx
         double precision, dimension(nvar, nz), intent(in) :: xold
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
         ierr = 0
         ! if the variables must be non-negative, you can do the following:
         x = max(0d0, x)
         dx = x-xold
         x = xold+dx
      end subroutine xdomain
      
      ! the sizeequ routine returns a measure of the magnitude of the residuals.
      subroutine size_equ(iter, nvar, nz, equ, residual_norm, max_residual, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only: dp
         integer, intent(in) :: iter, nvar, nz
         real(dp), pointer :: equ(:,:) ! (nvar, nz)
         real(dp), intent(out) :: residual_norm, max_residual
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
         
         ! for example, you might just use an average magnitude:
          residual_norm = sum(abs(equ(:,:))/dble(neq))
          max_residual = maxval(abs(equ(:,:)))
          
          !write(*,'(99ES10.3)') equ, sum(abs(equ))
          !write(*,'(a20,i10,ES10.3)') 'neq:',neq, dble(neq)
          !write(*,'(a20,ES10.3,a20,ES10.3)') 'residual norm:',residual_norm, 'max residual:',max_residual
          
          ierr = 0
      end subroutine size_equ
      
      ! the "size" of the modification B is an indication of how good the 
      ! solution has become -- small enough B, and we declare victory!
      subroutine sizeB(iter, nvar, nz, x, B, xscale, max_correction, correction_norm, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only: dp
         integer, intent(in) :: iter, nvar, nz
         real(dp), pointer, dimension(:,:) :: x, B, xscale ! (nvar, nz)
         real(dp), intent(out) :: correction_norm ! a measure of the average correction
         real(dp), intent(out) :: max_correction ! magnitude of the max correction
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
         ! for example, 
          max_correction = maxval(abs(b))
          correction_norm = sum(abs(b)/(nvar*nz))
          
          !write(*,'(a20,ES10.3,a20,ES10.3)') 'correction norm:',correction_norm, 'max correction:',max_correction
          !write(*,*)
          
          ierr = 0
      end subroutine sizeB
      
      !Subroutines for finding the Neumann conditions:
       subroutine set_primaries_neumann(nvar, nz, x, lrpar, rpar, lipar, ipar, ierr)
        	integer, intent(in) :: nvar, nz
         	double precision, pointer :: x(:,:)
         	integer, intent(in) :: lrpar, lipar
         	double precision, intent(inout) :: rpar(:)
         	integer, intent(inout) :: ipar(:)
         	integer, intent(out) :: ierr
         	ierr = 0
         	x0 = x(1, 1); x1 = x(2, 1) !x0 = rho, x1 = T
         	!write(*, '(a20, 3(a6, 1pe26.16, 3x))') 'primaries', 'x0', x0, 'x1', x1
         	
         	if (dbg) write(*, '(a20, 3(a6, 1pe26.16, 3x))') 'primaries', 'x0', x0, 'x1', x1, 'x2', x2
		end subroutine set_primaries_neumann
		
		
		subroutine set_secondaries_neumann(ivar, lrpar, rpar, lipar, ipar, ierr)
         	integer, intent(in) :: ivar
         	integer, intent(in) :: lrpar, lipar
         	double precision, intent(inout) :: rpar(:)
         	integer, intent(inout) :: ipar(:)
         	integer, intent(out) :: ierr
         	logical, parameter :: skip_partials = .true.
         	call set_sec_neumann(ivar, skip_partials, lrpar, rpar, lipar, ipar, ierr)
      	end subroutine set_secondaries_neumann
      	
      
      	subroutine set_sec_neumann(ivar, skip_partials, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only:pi
         integer, intent(in) :: ivar
         logical, intent(in) :: skip_partials
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(lrpar)
		 integer, intent(inout) :: ipar(lipar)
		 integer, intent(out) :: ierr
		 ierr = 0
		 
		 !No secondaries in this example
      end subroutine set_sec_neumann
      
      subroutine eval_equations_neumann(iter, nvar, nz, x, xscale, equ, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: iter, nvar, nz
         double precision, pointer, dimension(:,:) :: x, xscale, equ
         !double precision, intent(out) :: equ(nvar, nz)
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
      	double precision, dimension(nvar*nz, nvar*nz) :: A ! square matrix for jacobian
			logical, parameter :: skip_partials = .true.			
			call eval_equ_neumann(nvar, nz, equ, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)         
      end subroutine eval_equations_neumann
      
      subroutine eval_equ_neumann(nvar, nz, equ, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)
      	use eos_def
      	use eos_lib
      	use const_def
      	
      	implicit none
      	
      	include 'formats.dek'
      	
         integer, intent(in) :: nvar, nz
         double precision, pointer, dimension(:,:) :: equ
			logical, intent(in) :: skip_partials
      	double precision, dimension(nvar*nz, nvar*nz) :: A 
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr           
         double precision :: y0, y1, y2
         
         !EOS variables
         double precision :: Z,Xm,Rho,logRho,logT
  		 double precision :: T,V,V0
  		 double precision :: e1, e0, p1, p0, xrho, xp, xd, xt
  		 double precision :: de_drho, de_dt, dp_drho, dp_dt, chi_rho, chi_T
  		 double precision :: q_bind, bind_final, w_final
  		 integer :: info
  		 double precision, dimension(num_eos_basic_results) :: res
		 double precision :: d_dlnRho_const_T(num_eos_basic_results) 
    	 double precision :: d_dlnT_const_Rho(num_eos_basic_results)
  	
         ierr = 0
         call set_sec_neumann(0, skip_partials, lrpar, rpar, lipar, ipar, ierr); if (ierr /= 0) return
         
		!write(*,*) 'Inside Hugoniot solver'
		!write(*,*) rpar
		!write(*,*)
		
		xrho = x0
		xt = x1
		
		!write(*,*) rho0, t0, q
		!write(*,*) x0, x1, v_det
		!write(*,*)
		
		!First evaluate the EOS at the initial (rho0, t0) point to get p0:
		V0 = 1d0/rho0
        logrho = log10(rho0)
        logt = log10(t0)
        call composition_info(species, chem_id, xa, xh, xhe, zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
        call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa, &
        	rho0, logRho, T0, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        	d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
		e0 = exp(res(i_lnE))
		p0 = exp(res(i_lnPgas)) + 1/3d0*crad*t0**4
        
        !Now evaluate the EOS at the current (rho, t) point to get p and necessary derivs:
        rho = xrho*rho0
        !V = 1d0/rho
        t = xt*t0
        logrho = log10(rho)
        logt = log10(t)
        call composition_info(species, chem_id, xa, xh, xhe, zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
         call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa, &
        	rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        	d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
		e1 = exp(res(i_lnE))
		p1 = exp(res(i_lnPgas)) + 1/3d0*crad*t**4
		
		xp = p1/p0
		xd = v_det**2*rho0/p0
		
		de_drho = res(i_dE_dRho)
		de_dt = res(i_Cv)
		chi_rho = res(i_chiRho)
		dp_drho = p1/rho*chi_Rho
		chi_t = res(i_chiT)
		dp_dt = p1/T*chi_T
		
		 !Old (unworking) dimensionfull equations:
		 !y0 = e1-e0-(5d-1)*(V0-V)*(P0+P1)
         !y1 = v_det**2-V0**2*(P1-P0)/(V0-V)
         
         !Dimensionless equations:
         !(We don't need q in here since we want to stay on the initial state Hugoniot)
         y0 = (1+xp)*(1-1/xrho) + 2*rho0*(e0-e1)/p0
         y1 = (xp-1)/(1/xrho-1) + xd
         
         !Rescale eqns to have them behave better:
         !y0 = y0/(V0*p0)
         !y1 = y1/(V0*p0)
  
  		!write(*,*) 'allo!',p1,exp(res(i_lnPgas)),1/3d0*crad*t**4
  		!write(*,*) 'allo!',p1/rho*res(i_chiRho),d_dlnRho_const_T(i_lnPgas)*p1/rho
  		!write(*,*) 'allo!',de_drho,de_dt,dp_drho,dp_dt
  
         equ(:, 1) = (/ y0, y1 /)
         A = 0d0
         if (.not. skip_partials) then
            ! A(q, v) = xscale(v) * partial of equation(q) wrt variable(v)
            
            !Dimensionless Jacobian:
    		A(1,1) = (1+xp)/(xrho**2) - 2*rho0**2*de_drho/p0 + (1-1/xrho)*xp*chi_rho/xrho
            A(1,2) = -2*rho0*t0*de_dt/p0 + (1-1/xrho)*xp*chi_t/xt
            A(2,1) = (xp-1-(xrho-1)*xp*chi_rho)/(xrho-1)**2
            A(2,2) = xrho*xp*chi_t/(xt*(1-xrho))
            
            !Old (unworking) dimensionful Jacobian:
            !A(1, 1) = de_drho - 5d-1*(p0+p1)*V**2 - 5d-1*(V0-V)*dp_drho
            !A(1, 2) = de_dt - 5d-1*(V0-V)*dp_dt
            !A(2, 1) = -V0**2*(dp_drho/(V0-V) - V**2*(p1-p0)/(V0-V)**2)
            !A(2, 2) = -V0**2*(dp_dt/(V0-V))
            
            A(:, 1) = A(:, 1)*xscale_newton(1,1)
            A(:, 2) = A(:, 2)*xscale_newton(2,1)
         end if
      end subroutine eval_equ_neumann
      
      subroutine eval_jacobian_neumann(ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
      	implicit none
      
         integer, intent(in) :: ldA ! leading dimension of A
         double precision :: A(ldA, nvar*nz) ! the jacobian matrix
         ! A(idiag+q-v, v) = partial of equation(q) wrt variable(v)
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr         
			logical, parameter :: skip_partials = .false.	
			integer :: i,j
			
         ierr = 0        
		 call eval_equ_neumann(nvar_neumann, nz, equ_newton, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)
			
		 if(jacobian_out) then
			 write(*, *)
			 write(*, *) 'result of analytic jacobian'
			 write(*, '(a10,3es20.12)') 'x = ',x_newton
			 write(*, '(a10,3es20.12)') 'xs = ',xscale_newton
			 write(*, '(a10,3es20.12)') 'equ = ',equ_newton
			 do j=1, nvar_neumann ! variable
				do i=1, nvar_neumann ! equation
				   write(*, '(e20.12, 3x)', advance = 'no') A(i, j)
				end do
				write(*, *)
			 end do
			 write(*, *)
		 endif
			
      end subroutine eval_jacobian_neumann
      

      subroutine enter_setmatrix_neumann(iter, nvar, nz, neqs, x, xold, xscale, xder, need_solver_to_eval_jacobian, &
            ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
            
            implicit none
            
         integer, intent(in) :: iter, nvar, nz, neqs
         double precision, pointer, dimension(:,:) :: x, xold, xscale, xder
         !double precision, intent(out) :: xder(nvar, nz)
         logical, intent(out) :: need_solver_to_eval_jacobian
         integer, intent(in) :: ldA ! leading dimension of A
         double precision, pointer, dimension(:) :: A !(ldA, neqs) ! the jacobian matrix
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
         integer :: i, j
         if (dbg) write(*, '(/, a)') 'enter_setmatrix'
         
         if (numerical_jacobian) then
            xder=epsder*(xscale+abs(xold))
            need_solver_to_eval_jacobian = .true.
            doing_jacobian = .true.
         else
            call eval_jacobian_neumann(ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
            need_solver_to_eval_jacobian = .false.
         end if
      end subroutine enter_setmatrix_neumann
      
      !Try to limit xrho >= 2 and xt >= 2 so as not to accidentally find the ambient 
      !conditions as a solution:
      subroutine xdomain_neumann(iter, nvar, nz, x, dx, xold, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only: dp
         integer, intent(in) :: iter, nvar, nz
         real(dp), pointer, dimension(:,:) :: x, dx, xold ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: ierr
         ierr = 0
         ! if the variables must be non-negative, you can do the following:
         x = max(0d0, x)
         dx = x-xold
         x = xold+dx
         
         !x(1,1) = max(5d0, x(1,1))  
         !x(2,1) = max(35d0, x(2,1))          
         !dx = x-xold
         !x = xold+dx
      end subroutine xdomain_neumann
      
      !Subroutines for finding the hug conditions:
       subroutine set_primaries_hug(nvar, nz, x, lrpar, rpar, lipar, ipar, ierr)
        	integer, intent(in) :: nvar, nz
         	double precision, pointer :: x(:,:)
         	integer, intent(in) :: lrpar, lipar
         	double precision, intent(inout) :: rpar(:)
         	integer, intent(inout) :: ipar(:)
         	integer, intent(out) :: ierr
         	ierr = 0
         	x0 = x(1, 1); !x1 = x(2, 1) !x0 = rho, x1 = T
         	!write(*, '(a20, 3(a6, 1pe26.16, 3x))') 'primaries', 'x0', x0, 'x1', x1
         	
         	if (dbg) write(*, '(a20, 3(a6, 1pe26.16, 3x))') 'primaries', 'x0', x0, 'x1', x1, 'x2', x2
		end subroutine set_primaries_hug
		
		
		subroutine set_secondaries_hug(ivar, lrpar, rpar, lipar, ipar, ierr)
         	integer, intent(in) :: ivar
         	integer, intent(in) :: lrpar, lipar
         	double precision, intent(inout) :: rpar(:)
         	integer, intent(inout) :: ipar(:)
         	integer, intent(out) :: ierr
         	logical, parameter :: skip_partials = .true.
         	call set_sec_hug(ivar, skip_partials, lrpar, rpar, lipar, ipar, ierr)
      	end subroutine set_secondaries_hug
      	
      
      	subroutine set_sec_hug(ivar, skip_partials, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only:pi
         integer, intent(in) :: ivar
         logical, intent(in) :: skip_partials
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(lrpar)
		 integer, intent(inout) :: ipar(lipar)
		 integer, intent(out) :: ierr
		 ierr = 0
		 
		 !No secondaries in this example
      end subroutine set_sec_hug
      
      subroutine eval_equations_hug(iter, nvar, nz, x, xscale, equ, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: iter, nvar, nz
         double precision, pointer, dimension(:,:) :: x, xscale, equ
         !double precision, intent(out) :: equ(nvar, nz)
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
      	double precision, dimension(nvar*nz, nvar*nz) :: A ! square matrix for jacobian
			logical, parameter :: skip_partials = .true.			
			call eval_equ_hug(nvar, nz, equ, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)         
      end subroutine eval_equations_hug
      
      subroutine eval_equ_hug(nvar, nz, equ, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)
      	use eos_def
      	use eos_lib
      	use const_def
      	
      	implicit none
      	
      	include 'formats.dek'
      	
         integer, intent(in) :: nvar, nz
         double precision, pointer, dimension(:,:) :: equ
			logical, intent(in) :: skip_partials
      	double precision, dimension(nvar*nz, nvar*nz) :: A 
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr           
         double precision :: y0, y1, y2
         
         !EOS variables
         double precision :: Z,Xm,Rho,logRho,logT
  		 double precision :: T,V,V0
  		 double precision e1, e0, p1, p0, xrho, xp, xd, xt
  		 double precision de_drho, de_dt, dp_drho, dp_dt, chi_rho, chi_T
  		 integer :: info
  		 double precision, dimension(num_eos_basic_results) :: res
		 double precision :: d_dlnRho_const_T(num_eos_basic_results) 
    	 double precision :: d_dlnT_const_Rho(num_eos_basic_results)
  	
         ierr = 0
         call set_sec_hug(0, skip_partials, lrpar, rpar, lipar, ipar, ierr); if (ierr /= 0) return
         
		!write(*,*) 'Inside Hugoniot solver'
		!write(*,*) rpar
		!write(*,*)
		
		xt = x0
		
		!First evaluate the EOS at the initial (rho0, t0) point to get p0:
        logrho = log10(rho0)
        logt = log10(t0)
        call composition_info(species, chem_id, xa, xh, xhe, zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
        call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa, &
        	rho0, logRho, T0, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        	d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
		e0 = exp(res(i_lnE))
		p0 = exp(res(i_lnPgas)) + 1/3d0*crad*t0**4
        
        !Now evaluate the EOS at the current (rho, t) point to get p and necessary derivs:
        rho = xrho_hug*rho0
        !V = 1d0/rho
        t = xt*t0
        logrho = log10(rho)
        logt = log10(t)
        call composition_info(species, chem_id, xa_hug, xh, xhe, zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
         call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa_hug, &
        	rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        	d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
		e1 = exp(res(i_lnE))
		p1 = exp(res(i_lnPgas)) + 1/3d0*crad*t**4
		
		xp = p1/p0
		
		de_drho = res(i_dE_dRho)
		de_dt = res(i_Cv)
		chi_rho = res(i_chiRho)
		dp_drho = p1/rho*chi_Rho
		chi_t = res(i_chiT)
		dp_dt = p1/T*chi_T
         
         !Dimensionless equations:
         y0 = (1+xp)*(1-1/xrho_hug) + 2*rho0*(e0-e1+q_hug)/p0
  
  		!write(*,*) 'Inside eval_equ_hug:'
  		!write(*,*) 'y0 = ',y0
  		!write(*,*) 'x0 = ',x0
  		!write(*,*) 'xp = ',xp
  		!write(*,*) 'q_hug = ',q_hug
  		!write(*,*) 'xrho_hug = ',xrho_hug
  		
         equ(:, 1) = (/ y0 /)
         A = 0d0
         if (.not. skip_partials) then
            ! A(q, v) = xscale(v) * partial of equation(q) wrt variable(v)
            
            !Dimensionless Jacobian:
    		A(1,1) = -2*rho0*t0*de_dt/p0 + (1-1/xrho_hug)*xp*chi_t/xt
            
            A(:, 1) = A(:, 1)*xscale_newton(1,1)
         end if
      end subroutine eval_equ_hug
      
      subroutine eval_jacobian_hug(ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
      	implicit none
      
         integer, intent(in) :: ldA ! leading dimension of A
         double precision :: A(ldA, nvar*nz) ! the jacobian matrix
         ! A(idiag+q-v, v) = partial of equation(q) wrt variable(v)
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr         
			logical, parameter :: skip_partials = .false.	
			integer :: i,j
			
         ierr = 0        
		 call eval_equ_hug(nvar_hug, nz, equ_newton, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)
			
		 if(jacobian_out) then
			 write(*, *)
			 write(*, *) 'result of analytic jacobian'
			 write(*, '(a10,3es20.12)') 'x = ',x_newton
			 write(*, '(a10,3es20.12)') 'xs = ',xscale_newton
			 write(*, '(a10,3es20.12)') 'equ = ',equ_newton
			 do j=1, nvar_hug ! variable
				do i=1, nvar_hug ! equation
				   write(*, '(e20.12, 3x)', advance = 'no') A(i, j)
				end do
				write(*, *)
			 end do
			 write(*, *)
		 endif
			
      end subroutine eval_jacobian_hug
      

      subroutine enter_setmatrix_hug(iter, nvar, nz, neqs, x, xold, xscale, xder, need_solver_to_eval_jacobian, &
            ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
            
            implicit none
            
         integer, intent(in) :: iter, nvar, nz, neqs
         double precision, pointer, dimension(:,:) :: x, xold, xscale, xder
         !double precision, intent(out) :: xder(nvar, nz)
         logical, intent(out) :: need_solver_to_eval_jacobian
         integer, intent(in) :: ldA ! leading dimension of A
         double precision, pointer, dimension(:) :: A !(ldA, neqs) ! the jacobian matrix
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(:)
         integer, intent(inout) :: ipar(:)
         integer, intent(out) :: ierr
         integer :: i, j
         if (dbg) write(*, '(/, a)') 'enter_setmatrix'
         
         if (numerical_jacobian) then
            xder=epsder*(xscale+abs(xold))
            need_solver_to_eval_jacobian = .true.
            doing_jacobian = .true.
         else
            call eval_jacobian_hug(ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
            need_solver_to_eval_jacobian = .false.
         end if
      end subroutine enter_setmatrix_hug
      
    !End - Subroutines needed by the newton solver--------------------------------------------------
      
      
      !Begin burn support routines (jacobian, derivs, etc.)
      !-----------------------------------------------------------------------------------
      
      subroutine znd_derivs(n, t, h, x, f, lrpar, rpar, lipar, ipar, ierr)
      	implicit none
      	
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: t, h
         real(dp), intent(inout) :: x(:)
         real(dp), intent(out) :: f(:) ! dxdt
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         integer, parameter :: ld_dfdx = 0
         real(dp) :: dfdx(ld_dfdx,n)
         ierr = 0
         
         !write(*,*) 'entering znd_derivs'
         call znd_jacob(n, t, h, x, f, dfdx, ld_dfdx, lrpar, rpar, lipar, ipar, ierr)
         !write(*,*) 'leaving znd_derivs'
      end subroutine znd_derivs


      subroutine znd_jacob(n, time, ht, x, f, dfdx, ld_dfdx, lrpar, rpar, lipar, ipar, ierr)
      
      	implicit none
         
         integer, intent(in) :: n, ld_dfdx, lrpar, lipar
         real(dp), intent(in) :: time, ht
         real(dp), intent(inout) :: x(:)
         real(dp), intent(out) :: f(:), dfdx(:,:)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         
         double precision :: composition(species)
         double precision :: rho, t, u, ux, uy, H, logT, logRho, P, energy, ye, abar, zbar
         double precision :: bind_final, w_final, q_bind, gamma1, cs
         logical :: just_dxdt
         type (Net_Info), target :: netinfo_target
         type (Net_Info), pointer :: netinfo
         !First derivatives that we'll need:
         double precision :: dp_drho, de_drho, dp_dT, de_dT, dp_da, de_da, &
         	dp_dz, de_dz, dP_de_rho
         !Second derivatives that we'll need:
         double precision :: d2p_drhoT, d2p_drho2, d2p_dT2, d2p_da2, d2p_drhoa, &
         	d2p_dTa, d2e_drhoT, d2e_drho2, d2e_dT2, d2e_da2, d2e_drhoa, &
         	d2e_dTa, d2p_drhoz, d2p_dTz, d2p_daz, d2p_dz2, &
         	d2e_drhoz, d2e_dTz, d2e_daz, d2e_dz2, d_dP_de_rho_drho, d_dP_de_rho_dT
         !Other supplemental variables:
         double precision :: af2, sigma(species), neta, grav, delta_y
         double precision :: f1, f2, f3	!Curvature variables
         double precision :: f1_curvature, f2_curvature, f3_curvature
         double precision :: f1_blowout, f2_blowout, f3_blowout
         double precision :: dneta_drho, dneta_dt, dneta_dux
         double precision, dimension(species) :: df1_curvature_dXi, df2_curvature_dXi, df3_curvature_dXi, df1_dXi, df2_dXi, df3_dXi
         double precision :: df1_curvature_drho, df2_curvature_drho, df3_curvature_drho, df1_drho, df2_drho, df3_drho, &
        	df1_curvature_dT, df2_curvature_dT, df3_curvature_dT, df1_dT, df2_dT, df3_dT, &
        	df1_curvature_dux, df2_curvature_dux, df3_curvature_dux, df1_dux, df2_dux, df3_dux, &
        	df1_curvature_duy, df2_curvature_duy, df3_curvature_duy, df1_duy, df2_duy, df3_duy, &
        	df1_curvature_dH, df2_curvature_dH, df3_curvature_dH, df1_dH, df2_dH, df3_dH
         double precision, dimension(species) :: term1, term2, term3, term1a, term1b, dneta_dx
         double precision, dimension(species) :: dP_dXi, dE_dXi, d_dP_de_rho_dXi, d2P_dTdXi, d2P_drhodXi
         double precision, dimension(species,species) :: d2P_dXidXj, d2E_dXidXj
         double precision, dimension(species) :: d2p_dx2, d2e_dx2	!2nd derivatives components for readability
         double precision :: sum_y, sum_yz
         character :: testchar
         integer :: i
         
         include 'formats.dek'
         
         !write(*,*) 'entering znd_jacob'

         netinfo => netinfo_target
         
         ierr = 0
         f = 0d0
         dfdx = 0d0
         grav = 0d0
         
         !Make sure the abundance pieces of the variable array are >= 0:
         do i=1,species
         	x(i) = max(0d0, min(1d0, x(i)))
         end do
         
         composition = x(1:species) !Array of just the abundances
         
         !Initialize the local thermo variables:
         rho = max(x(species+1), 0d0)
         T = max(x(species+2), 0d0)
         ux = max(x(species+3), 0d0)
         if(do_blowout) then
         	!uy = max(x(species+4), 0d0)	!uy can be negative
         	uy = x(species+4)
         	H = max(x(species+5), 0d0)
         	!grav is acceleration due to gravity at the midpoint of the blowout layer
         	!comment out the following line to turn off gravity evolution:
			!grav = standard_cgrav*m_c/((r_wd + H/2d0)**2)
			grav = grav_init
		 else if (do_blowout_local) then
		 	!uy = max(x(species+4), 0d0)		!uy can be negative
		 	uy = x(species+4)
         	H = max(x(species+5), 0d0)
         	delta_y = max(x(species+6), 0d0)
         endif
         if(do_curvature.and.(.not.use_he_clavin)) then
        	r_curve = max(x(num_vars), 0d0)
         endif
         
         !if(num_vars.gt.species+3) P = x(species+4)
         logT = log10(T)
         logRho = log10(rho)
         
         !write(*,*) 'size(x) - species',size(x) - species
         !write(*,*) 'ux:',ux
         !write(*,*) 'uy:',uy
         !read(*,*) testchar
         
         call composition_info(species, chem_id, composition, xh, xhe, zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
         sum_y = sum(composition/chem_isos% W(chem_id(1:species)))
         sum_yz = sum(composition*chem_isos% Z(chem_id(1:species))/&
         	chem_isos% W(chem_id(1:species)))
            
        logRho = log10(rho)
        logT = log10(T)
        call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, composition, &
        	rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
        	d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
     	if (ierr /= 0) then
            write(*,*) 'failed in eosDT_get', ierr
        endif
        P = exp(res(i_lnPgas)) + (crad*T**4)/3d0			!Total pressure of gas
        energy = exp(res(i_lnE)) 							!Internal erg/g of gas
        gamma1 = res(i_gamma1)
        
        call eos_get_helm_results(Xh, abar, zbar, &
        	Rho, logRho, T, logT, .true., .false., .false., res_helm, ierr)
     	if (ierr /= 0) then
            write(*,*) 'failed in eos_get_helm_results', ierr
        endif
        P = res_helm(h_ptot)			!Total pressure of gas
        energy = res_helm(h_etot)
        !First derivatives:
        dp_drho = res_helm(h_dpd)
        dp_dT = res_helm(h_dpt)
        dp_da = res_helm(h_dpa)
        dp_dz = res_helm(h_dpz)
        de_drho = res_helm(h_ded)
        de_dT = res_helm(h_det)
        de_da = res_helm(h_dea)
        de_dz = res_helm(h_dez)
        dP_de_rho = dP_dT/de_dT
        !Second derivatives:
        d2p_drho2 = res_helm(h_dpdd)
        d2p_dT2 = res_helm(h_dptt)
        d2p_da2 = res_helm(h_dpaa)
        d2p_dz2 = res_helm(h_dpzz)
        d2p_drhoT = res_helm(h_dpdt)
        d2p_drhoa = res_helm(h_dpda)
        d2p_drhoz = res_helm(h_dpdz)
        d2p_dTa = res_helm(h_dpta)
        d2p_dTz = res_helm(h_dptz)
        d2p_daz = res_helm(h_dpaz)
        d2e_drho2 = res_helm(h_dedd)
        d2e_dT2 = res_helm(h_dett)
        d2e_da2 = res_helm(h_deaa)
        d2e_dz2 = res_helm(h_dezz)
        d2e_drhoT = res_helm(h_dedt)
        d2e_drhoa = res_helm(h_deda)
        d2e_drhoz = res_helm(h_dedz)
        d2e_dTa = res_helm(h_deta)
        d2e_dTz = res_helm(h_detz)
        d2e_daz = res_helm(h_deaz)
        d_dP_de_rho_drho = -dP_dT*d2e_drhoT/(de_dT**2) + d2P_drhoT/de_dT
        d_dP_de_rho_dT = -dP_dT*d2e_dT2/(de_dT**2) + d2P_dT2/de_dT
        !Derivatives in terms of X_i & X_j:
        dP_dXi(1:species) = abar/chem_isos% W(chem_id(1:species))*(-abar*dP_da + &
        	dP_dz*(chem_isos% Z(chem_id(1:species)) - zbar))
        dE_dXi(1:species) = abar/chem_isos% W(chem_id(1:species))*(-abar*dE_da + &
        	dE_dz*(chem_isos% Z(chem_id(1:species)) - zbar))
        do j=1, species
        	d2P_dXidXj(1:species,j) = abar/(chem_isos% W(chem_id(1:species))*&
        		chem_isos% W(chem_id(j)))*(abar**2*dP_da + abar**2*d2P_da2 - &
        		abar*((chem_isos% Z(chem_id(1:species)) - zbar) + &
        		(chem_isos% Z(chem_id(j)) - zbar))*dP_dz - &
        		abar**2*((chem_isos% Z(chem_id(1:species)) - zbar) + &
        		(chem_isos% Z(chem_id(j)) - zbar))*d2P_daz + &
        		abar*(chem_isos% Z(chem_id(1:species)) - zbar)*&
        		(chem_isos% Z(chem_id(j)) - zbar)*d2P_dz2)
        	d2E_dXidXj(1:species,j) = abar/(chem_isos% W(chem_id(1:species))*&
        		chem_isos% W(chem_id(j)))*(abar**2*dE_da + abar**2*d2E_da2 - &
        		abar*((chem_isos% Z(chem_id(1:species)) - zbar) + &
        		(chem_isos% Z(chem_id(j)) - zbar))*dE_dz - &
        		abar**2*((chem_isos% Z(chem_id(1:species)) - zbar) + &
        		(chem_isos% Z(chem_id(j)) - zbar))*d2E_daz + &
        		abar*(chem_isos% Z(chem_id(1:species)) - zbar)*&
        		(chem_isos% Z(chem_id(j)) - zbar)*d2E_dz2)
        end do
        !Check if d2P_dXidXj = d2P_dXjdXi (just for diagnosis):
        !do i=1, species
        !	do j=1, species
        !		if(abs((d2P_dXidXj(i,j)-d2P_dXidXj(j,i))/d2P_dXidXj(i,j)).gt.1d-6) then
        !			write(*,*) 'd2P_dXidXj(i,j) != d2P_dXidXj(j,i) in Jacobian',i,j
        !			write(*,*) 'd2P_dXidXj(i,j) = ',d2P_dXidXj(i,j)
        !			write(*,*) 'd2P_dXidXj(j,i) = ',d2P_dXidXj(j,i)
        !			stop 1
        !		endif
        !		if(abs((d2E_dXidXj(i,j)-d2E_dXidXj(j,i))/d2E_dXidXj(i,j)).gt.1d-6) then
        !			write(*,*) 'd2E_dXidXj(i,j) != d2E_dXidXj(j,i) in Jacobian',i,j
        !			write(*,*) 'd2E_dXidXj(i,j) = ',d2E_dXidXj(i,j)
        !			write(*,*) 'd2E_dXidXj(j,i) = ',d2E_dXidXj(j,i)
        !			stop 1
        !		endif
        !	end do
        !end do
        d2P_dTdXi(1:species) = abar/chem_isos% W(chem_id(1:species))*(-abar*d2P_dTa + &
        	(chem_isos% Z(chem_id(1:species)) - zbar)*d2P_dTz)
        d2P_drhodXi(1:species) = abar/chem_isos% W(chem_id(1:species))*(-abar*d2P_drhoa + &
        	(chem_isos% Z(chem_id(1:species)) - zbar)*d2P_drhoz)
        d_dP_de_rho_dXi(1:species) = abar/chem_isos% W(chem_id(1:species))*&
        	((-abar*d2P_dTa + (chem_isos% Z(chem_id(1:species)) - zbar)*d2P_dTz)/de_dT - &
        	dP_dT/de_dT**2*(-abar*d2E_dTa + (chem_isos% Z(chem_id(1:species)) - zbar)*d2E_dTz))
        
        !Consistency checks:
        mass_flux = rho*ux
        mom_flux = P + rho*ux**2
        rayleigh = rho0**2*v_det**2 - (P-P0)/(1/rho0 - 1/rho)
        !Check the binding energy too (1.52d18 for He4 -> Ni56): 
		q_bind = avo*mev_to_ergs*sum(composition*&
			chem_isos% binding_energy(chem_id(1:species))/&
			chem_isos% W(chem_id(1:species))) - q_init
		energy_flux = energy - q_bind + P/rho + ux**2/2
		if(do_blowout.or.do_blowout_local) then
			mass_flux = H*rho*ux
        	mom_flux = H*(P + rho*ux**2)
        	!with gravitational energy:
        	energy_flux = (energy - q_bind + P/rho + (ux**2 + (uy/2)**2)/2 + grav*H/2)
        	!without gravitational energy:
        	!energy_flux = (energy - q_bind + P/rho + ux**2/2)
		endif
		
        just_dxdt = .false.
     	call net_get( &
            handle, just_dxdt, netinfo, species, num_reactions, &
            composition, T, logT, Rho, logRho, & 
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, & 
            std_reaction_Qs, std_reaction_neuQs, .false., .false., &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, & 
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, & 
            screening_mode, theta_e_for_graboske_et_al, &     
            eps_nuc_categories, eps_neu_total, & 
            lwork_net, work_net, ierr)
        if (ierr /= 0) then
            write(*,*) 'failed in net_get', ierr
            write(*,*) sum(composition), lwork_net
            write(*,*) 
            !write(*,*) handle, species, num_reactions, &
            !	x(1), burn_T, logT, burn_Rho, logRho, & 
            !	abar, zbar, z2bar, ye, eta
            return
         end if
         
         !write(*,*) 'net_get successful'
         !write(*,*) dxdt
         !write(*,*)
     
         if (.false. .and. ierr /= 0) then
            write(*,1) 'time', time
            write(*,1) 'T', 10**logT
            write(*,1) 'lgT', logT
            write(*,1) 'rho', 10**logRho
            write(*,1) 'lgRho', logRho
            write(*,1) 'eta', eta
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,1) 'z2bar', z2bar
            write(*,1) 'ye', ye
            write(*,1) 'xh', xh
            write(*,1) 'xhe', xhe
            write(*,1) 'sum(x)', sum(x)
            write(*,1) 'maxval(x)', maxval(x)
            write(*,1) 'minval(x)', minval(x)
            write(*,1) 'eps_nuc', eps_nuc
            write(*,1) 'eps_neu_total', eps_neu_total
            write(*,*)
            if (time > 0) stop
            !write(*,2) 'eval_net ierr', ierr
            !write(*,*) 'failed in znd_jacob'
            !stop '1 zone burn znd_jacob'
         end if
         
        !Set up our secondary variables here:
        !Frozen sound speed squared:
        af2 = dP_drho + (P/rho**2 - de_drho)*dp_dT/de_dT
        !Newer value (with total zbar derivatives):
        sigma(1:species) = ((-abar**2*dp_da + (abar*chem_isos% Z(chem_id(1:species)) - &
        	abar*zbar)*dp_dz) + &
        	(abar**2*de_da - (abar*chem_isos% Z(chem_id(1:species)) - abar*zbar)*de_dz + &
        	avo*chem_isos% binding_energy(chem_id(1:species))*mev_to_ergs)*dp_dt/de_dt)/(rho*af2)
        !New value (with zbar derivatives at constant abar):
        !sigma(1:species) = (abar*chem_isos% Z(chem_id(1:species))*(-abar*dp_da/zbar + dp_dz) + &
        !	(abar*chem_isos% Z(chem_id(1:species))*(abar*de_da/zbar - de_dz) + &
        !	avo*chem_isos% binding_energy(chem_id(1:species))*mev_to_ergs)*dp_dt/de_dt)/(rho*af2)
        !Old value (without zbar derivatives):
        !sigma(1:species) = (-abar**2*dp_da + (abar**2*de_da + &
        !	avo*chem_isos% binding_energy(chem_id(1:species))*mev_to_ergs)*dp_dt/de_dt)/(rho*af2)
        neta = af2 - ux**2
        cs = sqrt(af2)	!local sound speed - or cs = sqrt(gamma1*P/rho)?
        
        !Set up the source terms for the hydro equations here:
        f1 = 0d0	!Mass: 		d(rho*ux)/dx = f1
        f2 = 0d0	!Momentum: 	d(P + rho*ux**2)/dx = f2
        f3 = 0d0	!Energy:	d(E + P/rho + ux**2/2)/dx = f3
        
        f1_blowout = 0d0
        f2_blowout = 0d0
        f3_blowout = 0d0
        f1_curvature = 0d0
        f2_curvature = 0d0
        f3_curvature = 0d0
        
        if(do_blowout) then
        	!duy/dx:
        	!with gravity:
        	if(use_const_uy) then
        		f(species+4) = 0d0
        	else
        		f(species+4) = P/(rho*H*ux) - uy**2/(H*ux) - grav/ux	!new (legit)
        		!f(species+4) = P/(rho*H*ux) - grav/ux					!old
        	endif
        	!without gravity:
        	!f(species+4) = P/(rho*H*ux) - uy**2/(H*ux)
        
        	!dH/dx:
        	f(species+5) = uy/ux
        
        	f1 = f1 - rho*uy/H
        	f2 = f2 - (P+rho*ux**2)*uy/(H*ux)						!newer (legit)
        	!f2 = f2 - rho*ux*uy/H									!new 
        	!f3 = f3 - P*uy/(4*H*rho*ux) - grav*uy/(4*ux) + grav*H*uy/(2*ux*(r_wd + H/2d0))		!old
        	!Need to check if gravity implemented correctly here:
        	f3 = f3 - uy*f(species+4)/4d0 - grav*uy/(2*ux) + grav*H*uy/(2*ux*(r_wd + H/2d0))	!new
        	
        	f1_blowout = -rho*uy/H
        	f2_blowout = -(P+rho*ux**2)*uy/(H*ux)
        	f3_blowout = -uy*f(species+4)/4d0 - grav*uy/(2*ux) + grav*H*uy/(2*ux*(r_wd + H/2d0))	
        endif
        if(do_curvature) then
        	!j_curve = 1		!0: plane parallel, 1: cylindrical, 2: spherical
			f1 = f1 - j_curve*rho*(v_det - ux)/r_curve
			f2 = f2 - j_curve*rho*ux*(v_det - ux)/r_curve
			!(no modification to f3)
			
			f1_curvature = -j_curve*rho*(v_det - ux)/r_curve
			f2_curvature = -j_curve*rho*ux*(v_det - ux)/r_curve
		endif	
		
		!dX_i/dx:
        f(1:species) = dxdt/ux
        
		!Try the equations in a form more amenable to changing the blowout prescription:
		!drho/dx:
        f(species+1) = -rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta) + &  
        	(f2 - 2*ux*f1 - dP_dT/dE_dT*(f3 - f2/rho + ux*f1/rho))/&
        	(-ux**2 + dP_drho + (P/rho**2 - dE_drho)*dP_dT/dE_dT)
        	
        !dT/dx:
        f(species+2) = ((ux**2 - dP_drho)*f(species+1) - &
        	sum((-abar**2*dp_da + (abar*chem_isos% Z(chem_id(1:species)) - &
        	abar*zbar)*dp_dz)*f(1:species)/chem_isos% W(chem_id(1:species))) + &
        	f2 - 2*ux*f1)/dP_dT
        
        !dux/dx:
        f(species+3) = f1/rho - ux*f(species+1)/rho
        
        !For curvature, we also explicitly track R_c(xi)
        if(do_curvature.and.(.not.use_he_clavin)) then
        	f(num_vars) = v_det/ux - 1d0
        endif

!This entire section below was working fine, but is being commented out as part of code 
!streamlining. The equations are now defined more clearly in terms of source terms above.
!To enable, use a regexp find/replace to match "^!" with "" - likewise disable with a
!find/replace on "^" with "!".

!		!Actual ODE equations we're integrating:
!		if(do_blowout) then 
!		!Source terms in the hydro equations
!        f1 = -rho*uy/H
!        f2 = -(P+rho*ux**2)*uy/(H*ux)
!        f3 = -P*uy/(4*H*rho*ux) - grav*uy/(4*ux) + grav*H*uy/(2*ux*(r_wd + H/2d0))
!        
!		!h_scale is approximate scale hight of material (in cm)
!		
!		!dX_i/dx:
!        f(1:species) = dxdt/ux
!        
!		!Try the equations in a form more amenable to changing the blowout prescription:
!		!drho/dx:
!        f(species+1) = -rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta) + &  
!        	(f2 - 2*ux*f1 - dP_dT/dE_dT*(f3 - f2/rho + ux*f1/rho))/&
!        	(-ux**2 + dP_drho + (P/rho**2 - dE_drho)*dP_dT/dE_dT)
!        	
!        !dT/dx:
!        f(species+2) = ((ux**2 - dP_drho)*f(species+1) - &
!        	sum((-abar**2*dp_da + (abar*chem_isos% Z(chem_id(1:species)) - &
!        	abar*zbar)*dp_dz)*f(1:species)/chem_isos% W(chem_id(1:species))) + &
!        	f2 - 2*ux*f1)/dP_dT
!        
!        !dux/dx:
!        f(species+3) = f1/rho - ux*f(species+1)/rho
!         
!        !drho/dx:
!        !With conserved fluxes (energy including H - Pg. 96 in notes, 
!        !f_E(x) is multiplied by dP_dT/dE_dT on 2nd line):
!        !f(species+1) = -rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta) - &  
!        !	((P - rho*ux**2)*uy/(H*ux) + dP_dT/dE_dT*((P - rho*ux**2)*uy/(rho*H*ux) + &
!        !	 ux*uy/H - P*uy/(4*rho*H*ux) - grav*uy/(4*ux)*(R_wd - 3*H/2)/(R_wd + H/2)))/&
!        !	(-ux**2 + dP_drho + (P/rho**2 - dE_drho)*dP_dT/dE_dT)
!        
!        !With conserved fluxes (energy w/o H - Pg. 62 in notes):
!        !f(species+1) = -rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta) - &  
!        !	((P - rho*ux**2)*uy/(H*ux) + dP_dT/dE_dT*(P*uy/(rho*H*ux)))/&
!        !	(-ux**2 + dP_drho + (P/rho**2 - dE_drho)*dP_dT/dE_dT)
!        	
!        !With conserved fluxes, including uy**2/2 into energy flux
!        !f(species+1) = -rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta) + &  
!        !	(dP_dT/dE_dT*(-uy**3/(ux*H)))/&
!        !	(-ux**2 + dP_drho + (P/rho**2 - dE_drho)*dP_dT/dE_dT)
!        
!        !With conserved fluxes, but without including uy**2/2 into energy flux
!        !f(species+1) = -rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta) + &  
!        !	uy/(H*ux)*(dP_dT/dE_dT*(energy - q_bind + ux**2/2) - (P - rho*ux**2))/&
!        !	(-ux**2 + dP_drho + (P/rho**2 - dE_drho)*dP_dT/dE_dT)
!        
!        !dT/dx:
!        !f(species+2) = ((ux**2 - dP_drho)*f(species+1) - &
!        !	sum((-abar**2*dp_da + (abar*chem_isos% Z(chem_id(1:species)) - &
!        !	abar*zbar)*dp_dz)*f(1:species)/chem_isos% W(chem_id(1:species))))/dp_dT - &
!        !	(P - rho*ux**2)*uy/(H*ux*dP_dT)
!        	
!        	!rho**2*ux*uy/(H_scale*dP_dT)
!        
!        !write(*,*) rho**2*ux*uy/(H_scale*dP_dT)
!        !write(*,*) 'E:', energy
!        !write(*,*) 'q_bind:', q_bind
!        !write(*,*)
!        	
!        !dux/dx:
!        !f(species+3) = -ux*f(species+1)/rho - uy/H
!        
!        !duy/dx:
!        !with gravity:
!        f(species+4) = P/(rho*H*ux) - grav/ux
!        !f(species+4) = P/(rho*H*ux) - standard_cgrav*m_c/(ux*(r_wd + H/2d0)**2)
!        !without gravity:
!        !f(species+4) = P/(rho*H*ux) - uy**2/(H*ux)
!        
!        !dH/dx:
!        f(species+5) = uy/ux
!        
!        !Local prescription for blowout (doesn't use scale height to set speed)
!        else if (do_blowout_local) then
!        
!        x(species+4) = cs
!        uy = cs
!        
!        !duy/dx:
!        !f(species+4) = P/(rho*delta_y*ux) - grav/ux
!        f(species+4) = 0d0
!        
!        !dH/dx:
!        f(species+5) = uy/ux
!        
!        !ddelta_y/dx:
!        f(species+6) = cs/ux
!        
!        !Source terms:
!		f1 = -rho*uy/H						!d(rho*ux)/dx = f1
!		f2 = -(P + rho*ux**2)*uy/(H*ux)		!d(P + rho*ux**2)/dx = f2
!		f3 = 0d0
!		!f3 = -uy/4*f(species+4) - grav*f(species+5)*(r_wd - H/2)/(r_wd + H/2)
!		
!		!dX_i/dx:
!        f(1:species) = dxdt/ux
!        
!        !drho/dx:
!        f(species+1) = -rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta) + &  
!        	(f2 - 2*ux*f1 - dP_dT/dE_dT*(f3 - f2/rho + ux*f1/rho))/&
!        	(-ux**2 + dP_drho + (P/rho**2 - dE_drho)*dP_dT/dE_dT)
!        	
!        !dT/dx:
!        f(species+2) = ((ux**2 - dP_drho)*f(species+1) - &
!        	sum((-abar**2*dp_da + (abar*chem_isos% Z(chem_id(1:species)) - &
!        	abar*zbar)*dp_dz)*f(1:species)/chem_isos% W(chem_id(1:species))) + &
!        	f2 - 2*ux*f1)/dP_dT
!        
!        !dux/dx:
!        f(species+3) = f1/rho - ux*f(species+1)/rho
!        
!        
!        !Copied eqns from blowout section (not conserving energy here, although they do above)
!        !dX_i/dx:
!        !f(1:species) = dxdt/ux
!         
!        !drho/dx:
!        !With conserved fluxes (energy including H - Pg. 96 in notes, 
!        !f_E(x) is multiplied by dP_dT/dE_dT on 2nd line):
!        !f(species+1) = -rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta) - &  
!        !	((P - rho*ux**2)*uy/(H*ux) + dP_dT/dE_dT*((P - rho*ux**2)*uy/(rho*H*ux) + &
!        !	 ux*uy/H - P*uy/(4*rho*H*ux) - grav*uy/(4*ux)*(R_wd - 3*H/2)/(R_wd + H/2)))/&
!        !	(-ux**2 + dP_drho + (P/rho**2 - dE_drho)*dP_dT/dE_dT)
!        
!        !dT/dx:
!        !f(species+2) = ((ux**2 - dP_drho)*f(species+1) - &
!        !	sum((-abar**2*dp_da + (abar*chem_isos% Z(chem_id(1:species)) - &
!        !	abar*zbar)*dp_dz)*f(1:species)/chem_isos% W(chem_id(1:species))))/dp_dT - &
!        !	(P - rho*ux**2)*uy/(H*ux*dP_dT)
!        
!        !dux/dx:
!        !f(species+3) = -ux*f(species+1)/rho - uy/H
!		
!		!Use curvature source terms in ZND equations instead of blowout source terms:
!		!See Pg. 106-7 of notes
!		else if(do_curvature) then
!		j_curve = 1		!0: plane parallel, 1: cylindrical, 2: spherical
!		!Curvature source terms in hydrodynamic equations:
!		f1 = -j_curve*rho*ux/r_curve				!d(rho*ux)/dx = f1
!		f2 = -j_curve*rho*ux**2/r_curve					!d(P + rho*ux**2)/dx = f2
!		f3 = 0d0
!		
!		!dX_i/dx:
!        f(1:species) = dxdt/ux
!        
!        !drho/dx:
!        f(species+1) = -rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta) + &  
!        	(f2 - 2*ux*f1 - dP_dT/dE_dT*(f3 - f2/rho + ux*f1/rho))/&
!        	(-ux**2 + dP_drho + (P/rho**2 - dE_drho)*dP_dT/dE_dT)
!        	
!        !dT/dx:
!        f(species+2) = ((ux**2 - dP_drho)*f(species+1) - &
!        	sum((-abar**2*dp_da + (abar*chem_isos% Z(chem_id(1:species)) - &
!        	abar*zbar)*dp_dz)*f(1:species)/chem_isos% W(chem_id(1:species))) + &
!        	f2 - 2*ux*f1)/dP_dT
!        
!        !dux/dx:
!        f(species+3) = f1/rho - ux*f(species+1)/rho
!		
!		else
!		!dX_i/dx:
!        f(1:species) = dxdt/ux
!         
!        !drho/dx:
!        !Equation derived from constant energy flux (rather than via Sharpe 99):
!        !(ended up being the same as the standard 1D ZND, as expected)
!        !f(species+1) = sum(f(1:species)*(abar**2*dp_da - &
!        !	abar*(chem_isos% Z(chem_id(1:species)) - zbar)*dp_dz + &
!        !	dP_dT*(-avo*chem_isos% binding_energy(chem_id(1:species))*mev_to_ergs - &
!        !	abar**2*(dP_da/rho + de_da) + abar*(chem_isos% Z(chem_id(1:species)) - zbar)*&
!        !	(dp_dz/rho + de_dz))/(dP_dT/rho + dE_dT))&
!        !	/chem_isos% W(chem_id(1:species)))/&
!        !	(dP_drho - ux**2 + dP_dT*(P/rho**2 + ux**2/rho - dP_drho/rho - dE_drho)/(dP_dT/rho + dE_dT))
!        !Enhanced 1D ZND (for simulating blowout):
!        !f(species+1) = -1.5*rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta)
!        !Standard 1D ZND:
!        f(species+1) = -rho*af2*sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))/(ux*neta)
!        
!        !dT/dx:
!        !Equation derived from constant energy flux (rather than via Sharpe 99):
!        !f(species+2) = ((P/rho**2 + ux**2/rho - dE_drho - dP_drho/rho)*f(species+1) - &
!        !	sum((-abar**2*(dp_da/rho + de_da) + (abar*chem_isos% Z(chem_id(1:species)) - &
!        !	abar*zbar)*(dp_dz/rho + de_dz) - &
!        !	avo*chem_isos% binding_energy(chem_id(1:species))*mev_to_ergs)*&
!        !	f(1:species)/chem_isos% W(chem_id(1:species))))/(dE_dT + dp_dT/rho)
!        !Newer value (with total zbar derivatives):
!        f(species+2) = ((ux**2 - dP_drho)*f(species+1) - &
!        	sum((-abar**2*dp_da + (abar*chem_isos% Z(chem_id(1:species)) - &
!        	abar*zbar)*dp_dz)*f(1:species)/chem_isos% W(chem_id(1:species))))/dp_dT
!        !New value (with zbar derivatives at constant abar):
!        !f(species+2) = ((ux**2 - dP_drho)*f(species+1) - &
!        !	abar*sum((-abar*dp_da/zbar + dp_dz)*chem_isos% Z(chem_id(1:species))*&
!        !	f(1:species)/chem_isos% W(chem_id(1:species))))/dp_dT
!        !Old value (without zbar derivatives):
!        !f(species+2) = ((ux**2 - dP_drho)*f(species+1) + &
!        !	abar**2*dP_da*sum(f(1:species)/chem_isos% W(chem_id(1:species))))/dp_dT
!        
!        !du/dx:
!        f(species+3) = -ux*f(species+1)/rho
!        
!        !Try overdetermining the system to make sure mass, momentum, and energy are
!        !fully conserved (we're not tracking P and E by default, so let's add them
!        !explicitly here). We won't couple these to the Jacobian so the steps will still
!        !be determined by rho, T, ux, and X_i evolution:
!        !dP/dx:
!        if(num_vars.gt.species+3) f(species+4) = -rho*ux*f(species+3)
!        !dE/dx:
!        if(num_vars.gt.species+4) f(species+5) = P*f(species+1)/rho**2
!        
!    	endif 	!(Blowout, curvature, or standard)
!         
        !write(*,*) f(species+1)
        !write(*,*) size(f), size(dfdx)
		 
        !write(*,*) 'abar = ',abar
        !write(*,*) 'sum_yz:', sum_yz, ' sum_y:',sum_y
        !write(*,*) 'sigma = ',sigma	
        !write(*,*) 'dxdt = ',dxdt
        !write(*,*) 'd_dxdt_dx = ',d_dxdt_dx
        !write(*,*) 'chem_isos% W(chem_id(1:species)) = ',chem_isos% W(chem_id(1:species))
        !write(*,*) 'sum = ',sum(dxdt/chem_isos% W(chem_id(1:species)))
         
        if (ld_dfdx > 0) then
        	!before fix, where dabar/dXi = abar**2/Ai:
        	!dneta_dx(1:species) = abar/chem_isos% W(chem_id(1:species))* &
         	!	(-abar*(d2P_drhoa + (dP_da/rho**2 - d2e_drhoa)*dp_dT/de_dT + &
         	!	(P/rho**2 - de_drho)*(d2P_dta/de_dt - dP_dT*d2e_dta/de_dt**2)) + &
         	!	(chem_isos% Z(chem_id(1:species)) - zbar)*&
         	!	(d2P_drhoz + (dP_dz/rho**2 - d2e_drhoz)*dp_dT/de_dT + &
         	!	(P/rho**2 - de_drho)*(d2P_dtz/de_dt - dP_dT*d2e_dtz/de_dt**2)))
         	!after fix, where dabar/dXi = abar(Ai - abar)/Ai:
         	dneta_dx(1:species) = abar/chem_isos% W(chem_id(1:species))* &
         		((chem_isos% W(chem_id(1:species))-abar)*&
         		(d2P_drhoa + (dP_da/rho**2 - d2e_drhoa)*dp_dT/de_dT + &
         		(P/rho**2 - de_drho)*(d2P_dta/de_dt - dP_dT*d2e_dta/de_dt**2)) + &
         		(chem_isos% Z(chem_id(1:species)) - zbar)*&
         		(d2P_drhoz + (dP_dz/rho**2 - d2e_drhoz)*dp_dT/de_dT + &
         		(P/rho**2 - de_drho)*(d2P_dtz/de_dt - dP_dT*d2e_dtz/de_dt**2)))
        	dneta_drho = d2P_drho2 + (dP_drho/rho**2 - P/rho**3 - d2e_drho2)*dP_dT/de_dT + &
        		(P/rho**2 - de_drho)*(d2P_drhoT/de_dT - dP_dT*d2e_drhot/de_dt**2)
        	dneta_dt = d2P_drhoT + (dP_dT/rho**2 - d2e_drhoT)*dP_dT/de_dT + &
        		(P/rho**2 - de_drho)*(d2P_dT2/de_dT - dP_dT*d2e_dT2/de_dt**2)
        	dneta_dux = -2*ux
        
			!d/dXj (dXi/dx):
         	dfdx(1:species, 1:species) = d_dxdt_dx/ux
         	!d/drho (dXi/dx):
         	dfdx(1:species, species+1) = d_dxdt_dRho/ux
         	!d/dT (dXi/dx):
         	dfdx(1:species, species+2) = d_dxdt_dT/ux
         	!d/dux (dXi/dx):
         	dfdx(1:species, species+3) = -dxdt/ux**2
         	
         	!d/dXj (drho/dx):
         	!First term is -R_i/(ux*neta) * d(rho*af2*sigma_i)/dXj
         	term1a = 0d0
         	term1b = 0d0
         	term1 = 0d0
         	term2 = 0d0
         	term3 = 0d0
         	!The sum over i in my notes is computed in this for loop. The explicit
         	!1:species index is j in my notes.
         	do i=1,species
         		!New terms (with full d/dzbar derivatives), only term1 changes from before
         		!d2P/dXidXj:
         		d2p_dx2 = abar**2*(2*dp_da + abar**2*d2p_da2 - &
         			abar*(chem_isos% Z(chem_id(1:species)) - zbar)*d2p_daz - &
         			(chem_isos% Z(chem_id(1:species)) + chem_isos% Z(chem_id(i)) - 2*zbar)*dp_dz - &
         			abar*(chem_isos% Z(chem_id(i)) - zbar)*d2p_daz + &
         			(chem_isos% Z(chem_id(i)) - zbar)*&
         			(chem_isos% Z(chem_id(1:species)) - zbar)*d2p_dz2)/&
         			(chem_isos% W(chem_id(1:species))*chem_isos% W(chem_id(i)))
         		!d2E/dXidXj:
         		d2e_dx2 = abar**2*(2*de_da + abar**2*d2e_da2 - &
         			abar*(chem_isos% Z(chem_id(1:species)) - zbar)*d2e_daz - &
         			(chem_isos% Z(chem_id(1:species)) + chem_isos% Z(chem_id(i)) - 2*zbar)*de_dz - &
         			abar*(chem_isos% Z(chem_id(i)) - zbar)*d2e_daz + &
         			(chem_isos% Z(chem_id(i)) - zbar)*&
         			(chem_isos% Z(chem_id(1:species)) - zbar)*d2e_dz2)/&
         			(chem_isos% W(chem_id(1:species))*chem_isos% W(chem_id(i)))
         		!d/dXj(rho*af2*sigma(i)):
         		term1(1:species) = term1(1:species) + (-f(i)/(neta*chem_isos% W(chem_id(i))))*&
         			(chem_isos% W(chem_id(1:species))*d2p_dx2 - &
         			chem_isos% W(chem_id(1:species))*d2e_dx2*(dP_dT/dE_dT) - &
         			(-abar**2*de_da + abar*(chem_isos% Z(chem_id(i)) - zbar)*de_dz - &
         			avo*mev_to_ergs*chem_isos% binding_energy(chem_id(i)))*&
         			((-abar**2*d2P_dTa + abar*(chem_isos% Z(chem_id(1:species)) - zbar)*d2P_dtz)/dE_dT - &
         			(-abar**2*d2E_dTa + abar*(chem_isos% Z(chem_id(1:species)) - zbar)*d2E_dtz)*&
         			dP_dT/dE_dT**2)/chem_isos% W(chem_id(1:species)))
         		term2(1:species) = term2(1:species) - rho*af2*sigma(i)*d_dxdt_dx(i,1:species)/&
         			(ux*neta*chem_isos% W(chem_id(i)))
         		term3(1:species) = term3(1:species) - rho*af2*sigma(i)*dxdt(i)/(ux*neta**2)*&
         			abar/chem_isos% W(chem_id(1:species))* &
         			(-abar*(d2P_drhoa + (dP_da/rho**2 - d2e_drhoa)*dp_dT/de_dT + &
         			(P/rho**2 - de_drho)*(d2P_dta/de_dt - dP_dT*d2e_dta/de_dt**2)) + &
         			(chem_isos% Z(chem_id(1:species)) - zbar)*&
         			(d2P_drhoz + (dP_dz/rho**2 - d2e_drhoz)*dp_dT/de_dT + &
         			(P/rho**2 - de_drho)*(d2P_dtz/de_dt - dP_dT*d2e_dtz/de_dt**2)))
         		!Should term3 have the deta_dx as defined above instead? (this looks midding d/dz terms)
         		!term3(1:species) = term3(1:species) - rho*af2*sigma(i)*dxdt(i)/(ux*neta**2)* &
         		!	abar**2/(chem_isos% W(chem_id(1:species))*chem_isos% W(chem_id(i)))* &
         		!	(d2P_drhoa + (dP_da/rho**2 - d2e_drhoa)*dp_dT/de_dT + &
         		!	(P/rho**2 - de_drho)*(d2P_dta/de_dt - dP_dT*d2e_dta/de_dt**2))
         		
         	
         		!Old terms (without full d/dzbar derivatives)
         		!d/dabar(rho*af2**2*sigma(i)):
         		!term1a(1:species) = term1a(1:species) - 2*abar*dP_da - abar**2*d2P_da2 + &
         		!	(chem_isos% Z(chem_id(i)) - zbar)*dP_dz + &
         		!	abar*(chem_isos% Z(chem_id(i)) - zbar)*d2P_daz + (2*abar*de_da + &
         		!	abar**2*d2E_da2 - (chem_isos% Z(chem_id(i)) - zbar)*de_dz - &
         		!	abar*(chem_isos% Z(chem_id(i)) - zbar)*d2E_daz)*dP_dT/de_dT + &
         		!	(abar**2*de_da - abar*(chem_isos% Z(chem_id(i)) - zbar)*de_dz + &
         		!	avo*mev_to_ergs*chem_isos% binding_energy(chem_id(i)))*&
         		!	(d2P_dta/de_dt - dp_dt*d2e_dta/de_dt**2)
         		!d/dzbar(rho*af2**2*sigma(i)):
         		!term1b(1:species) = term1b(1:species) - abar**2*d2P_daz - abar*dP_dz + &
         		!	abar*(chem_isos% Z(chem_id(i)) - zbar)*d2P_dz2 + &
         		!	(abar**2*d2E_daz + abar*dE_dz - abar*(chem_isos% Z(chem_id(i)) - &
         		!	zbar)*d2e_dz2)*dP_dT/de_dT + &
         		!	(abar**2*de_da - abar*(chem_isos% Z(chem_id(i)) - zbar)*de_dz + &
         		!	avo*mev_to_ergs*chem_isos% binding_energy(chem_id(i)))*&
         		!	(d2P_dtz/de_dt - dp_dt*d2e_dtz/de_dt**2)
         		!New term1 (with the d/dzbar derivatives were added):
         		!term1(1:species) = term1(1:species) + (-abar**2*term1a(1:species) + &
         		!	abar*(chem_isos% Z(chem_id(1:species)) - zbar)*term1b(1:species))/&
         		!	chem_isos% W(chem_id(1:species))*&
         		!	(-f(i)/(neta*chem_isos% W(chem_id(i))))
         		!Old term1 (before the d/dzbar derivatives were added):
         		!term1(1:species) = term1(1:species) + abar**2/(ux*neta*chem_isos% W(chem_id(1:species)))* &
         		!	(-2*abar*dP_da - abar**2*d2P_da2 + (2*abar*de_da + abar**2*d2e_da2)*dp_dt/de_dt + &
         		!	(abar**2*de_da + avo*chem_isos% binding_energy(chem_id(i))*mev_to_ergs)* &
         		!	(d2p_dta/de_dT - dp_dT*d2e_dTa/de_dt**2))*dxdt(i)/chem_isos% W(chem_id(i))
         		!term2(1:species) = term2(1:species) - rho*af2*sigma(i)*d_dxdt_dx(i,1:species)/&
         		!	(ux*neta*chem_isos% W(chem_id(i)))
         		!term3(1:species) = term3(1:species) - rho*af2*sigma(i)*dxdt(i)/(ux*neta**2)* &
         		!	abar**2/(chem_isos% W(chem_id(1:species))*chem_isos% W(chem_id(i)))* &
         		!	(d2P_drhoa + (dP_da/rho**2 - d2e_drhoa)*dp_dT/de_dT + &
         		!	(P/rho**2 - de_drho)*(d2P_dta/de_dt - dP_dT*d2e_dta/de_dt**2))
         	end do
         	!write(*,*) chem_isos% binding_energy(chem_id(1))*mev_to_ergs
         	!write(*,*) term1(1), term2(1), term3(1)
         	dfdx(species+1, 1:species) = term1(:) + term2(:) + term3(:)
         
         	!d/dXj (dT/dt):
         	term1a = 0d0
         	term1b = 0d0
         	term1 = 0d0
         	!d/dabar((dP/dY_i)*dY_i/dx):
         	term1a(1:species) = -2*abar*dP_da - abar**2*d2P_da2 + &
         		(chem_isos% Z(chem_id(1:species)) - zbar)*(dP_dz + abar*d2P_daz)
         	!d/dzbar((dP/dY_i)*dY_i/dx):
         	term1b(1:species) = -abar**2*d2P_daz - abar*dP_dz + &
         		abar*(chem_isos% Z(chem_id(1:species)) - zbar)*d2P_dz2
         	!The sum over i in my notes is computed in this for loop. The explicit
         	!1:species index is j in my notes.
         	do i=1,species
         		!Newer term1 (with full d/dabar d/dzbar derivatives):
         		!term1 is the whole sum over i of d/dXj(dP/dXi*dXi/dx)
         		term1(1:species) = term1(1:species) + f(i)*abar**2*(2*dP_da + &
         			abar**2*d2P_da2 - abar*(chem_isos% Z(chem_id(1:species)) - zbar)*d2P_daz - &
         			(chem_isos% Z(chem_id(1:species)) + chem_isos% Z(chem_id(i)) - 2*zbar)*dP_dz - &
         			abar*(chem_isos% Z(chem_id(i)) - zbar)*d2P_daz + d2P_dz2*&
         			(chem_isos% Z(chem_id(i)) - zbar)*(chem_isos% Z(chem_id(1:species)) - zbar))/&
         			(chem_isos% W(chem_id(1:species))*chem_isos% W(chem_id(i))) + &
         			dfdx(i,1:species)*(-abar**2*dP_da + abar*(chem_isos% Z(chem_id(i)) - zbar)*dP_dz)/&
         			chem_isos% W(chem_id(i))
         		!Here's where I figured out the missing dP_dz term:
         		!write(*,*) 'term1(2)', term1(2), f(i)*abar**2*(2*dP_da + &
         		!	abar**2*d2P_da2 - abar*(chem_isos% Z(chem_id(2)) - zbar)*d2P_daz - &
         		!	(chem_isos% Z(chem_id(2)) + chem_isos% Z(chem_id(i)) - 2*zbar)*dP_dz - &
         		!	abar*(chem_isos% Z(chem_id(i)) - zbar)*d2P_daz + d2P_dz2*&
         		!	(chem_isos% Z(chem_id(i)) - zbar)*(chem_isos% Z(chem_id(2)) - zbar))/&
         		!	(chem_isos% W(chem_id(2))*chem_isos% W(chem_id(i))), &
         		!	dfdx(i,2)*(-abar**2*dP_da + abar*(chem_isos% Z(chem_id(i)) - zbar)*dP_dz)/&
         		!	chem_isos% W(chem_id(i))
         	
         		!New term1 (with the d/dzbar terms):
         		!term1(1:species) = term1(1:species) + (-abar**2*term1a(i) + &
         		!	abar*(chem_isos% Z(chem_id(1:species)) - zbar)*term1b(i))*f(i)/&
         		!	(chem_isos% W(chem_id(1:species))*chem_isos% W(chem_id(i))) + &
         		!	(-abar**2*dP_da + abar*(chem_isos% Z(chem_id(i)) - zbar))*&
         		!	dfdx(i,1:species)/chem_isos% W(chem_id(i))
         		!write(*,*) 'term1(2): ',(-abar**2*term1a(i) + &
         		!	abar*(chem_isos% Z(chem_id(2)) - zbar)*term1b(i))*f(i)/&
         		!	(chem_isos% W(chem_id(2))*chem_isos% W(chem_id(i))), &
         		!	(-abar**2*dP_da + abar*(chem_isos% Z(chem_id(i)) - zbar))*&
         		!	dfdx(i,2)/chem_isos% W(chem_id(i))
         		!Old term1 (without the d/dzbar terms added in):
         		!term1(1:species) = term1(1:species) + d_dxdt_dx(i,1:species)/chem_isos% W(chem_id(i))
         	end do
         	
         	!New term (with the d/dzbar terms added in):
         	!write(*,*) 1/(dP_dT)*((abar**2*d2P_drhoa - &
         	!	abar*d2P_drhoz*(chem_isos% Z(chem_id(2)) - zbar))*&
         	!	f(species+1)/chem_isos% W(chem_id(2)) + &
         	!	(ux**2 - dP_drho)*dfdx(species+1,2) - term1(2)), &
         	!	(-abar**2*d2P_dTa + abar*d2P_dTz*(chem_isos% Z(chem_id(2)) - zbar))/&
         	!	(dP_dT**2*chem_isos% W(chem_id(2)))*&
         	!	((ux**2 - dP_drho)*f(species+1) - sum((-abar**2*dp_da + &
         	!	abar*(chem_isos% Z(chem_id(1:species)) - zbar)*dp_dz)*&
         	!	f(1:species)/chem_isos% W(chem_id(1:species))))
         	!write(*,*) 1/(dP_dT)*((abar**2*d2P_drhoa - &
         	!	abar*d2P_drhoz*(chem_isos% Z(chem_id(2)) - zbar))*&
         	!	f(species+1)/chem_isos% W(chem_id(2))), &
         	!	(ux**2 - dP_drho)*dfdx(species+1,2)/dP_dT, & 
         	!	-term1(2)/dP_dT
         	
         	!Newer term (9/19/12) with full d/dabar and d/dzbar derivatives:
         	dfdx(species+2, 1:species) = 1/(dP_dT)*((abar**2*d2P_drhoa - &
         		abar*d2P_drhoz*(chem_isos% Z(chem_id(1:species)) - zbar))*&
         		f(species+1)/chem_isos% W(chem_id(1:species)) + &
         		(ux**2 - dP_drho)*dfdx(species+1,1:species) - term1(1:species)) - &
         		(-abar**2*d2P_dTa + abar*d2P_dTz*(chem_isos% Z(chem_id(1:species)) - zbar))/&
         		(dP_dT**2*chem_isos% W(chem_id(1:species)))*&
         		((ux**2 - dP_drho)*f(species+1) - sum((-abar**2*dp_da + &
         		abar*(chem_isos% Z(chem_id(1:species)) - zbar)*dp_dz)*&
         		f(1:species)/chem_isos% W(chem_id(1:species))))
         		
        	!New term (with the d/dzbar terms added in):	
         	!dfdx(species+2, 1:species) = 1/(dP_dT)*(((abar**2*d2P_drhoa - &
         	!	abar*d2P_drhoz*(chem_isos% Z(chem_id(1:species)) - zbar))*&
         	!	f(species+1)/chem_isos% W(chem_id(1:species))) + &
         	!	(ux**2 - dP_drho)*dfdx(species+1,1:species) - term1(1:species)) + &
         	!	(abar**2*d2P_dta - abar*(chem_isos% Z(chem_id(1:species)) - zbar)*d2P_dtz)/&
         	!	(chem_isos% W(chem_id(1:species))*dP_dT**2)*&
         	!	((ux**2 - dP_drho)*f(species+1) - sum((-abar**2*dp_da + abar*(chem_isos% Z(chem_id(1:species)) - &
        	!	zbar)*dp_dz)*f(1:species)/chem_isos% W(chem_id(1:species))))
         	
         	!Old term (without the d/dzbar terms added in):
         	!dfdx(species+2, 1:species) = 1/(dP_dT)*(abar**2*d2P_drhoa*f(species+1)/chem_isos% W(chem_id(1:species)) + &
         	!	(ux**2 - dP_drho)*dfdx(species+1,1:species) - abar**2/(ux*chem_isos% W(chem_id(1:species)))*&
         	!	(2*abar*dP_da + abar**2*d2P_da2)*sum(dxdt/chem_isos% W(chem_id(1:species))) + &
         	!	abar**2*dP_da*term1(1:species)/ux) + abar**2*d2P_dta/(chem_isos% W(chem_id(1:species))*dP_dT**2)*&
         	!	((ux**2 - dP_drho)*f(species+1) + abar**2*dP_da/ux*sum(dxdt/chem_isos% W(chem_id(1:species))))
         	         	
         	!d_dXj (du/dx):
         	dfdx(species+3, 1:species) = -ux*dfdx(species+1, 1:species)/rho
         	
         	!d_dXj (dP/dx):
         	!dfdx(species+4, 1:species) = -rho*ux*dfdx(species+3, 1:species)
         	
         	!d_drho (dXi/dx):
         	dfdx(1:species, species+1) = d_dxdt_drho/ux
         	
         	!d_drho (drho/dx):
         	!Here we can just create the term arrays as being the pieces of the sum on i
         	!in the equations worked out on Pg. 3 of my notes (will be TeX'd eventually).
         	!This way, we can just use sum(term1) to do the sum over i, rather than using
         	!a do-loop as well.
         	term1 = 0d0
         	term2 = 0d0
         	term3 = 0d0
         	!New term1 (with d/dzbar derivatives):
         	term1(1:species) = term1(1:species) + (-abar**2*d2P_drhoa + abar*&
         		(chem_isos% Z(chem_id(1:species)) - zbar)*d2P_drhoz + (abar**2*d2E_drhoa - &
         		abar*(chem_isos% Z(chem_id(1:species)) - zbar)*d2E_drhoz)*dP_dT/dE_dT + &
         		(abar**2*de_da - abar*(chem_isos% Z(chem_id(1:species)) - zbar)*dE_dz + &
         		avo*mev_to_ergs*chem_isos% binding_energy(chem_id(1:species)))*&
         		(d2P_drhoT/de_dt - dp_dt*d2e_drhot/de_dt**2))*&
         		(-f(1:species)/(neta*chem_isos% W(chem_id(1:species))))
         	!Old term1 (before the d/dzbar derivatives were added):
         	!term1(1:species) = term1(1:species) - (-abar**2*d2P_drhoa + &
         	!	abar**2*d2e_drhoa*dP_dT/de_dT + (abar**2*de_da + &
         	!	avo*chem_isos% binding_energy(chem_id(1:species))*mev_to_ergs)*(d2P_drhot/de_dt - &
         	!	dP_dt*d2e_drhot/de_dt**2))*&
         	!	dxdt(1:species)/(ux*neta*chem_isos% W(chem_id(1:species)))
         	term2(1:species) = term2(1:species) - rho*af2*sigma(1:species)*d_dxdt_drho/&
         		(ux*neta*chem_isos% W(chem_id(1:species)))
         	term3(1:species) = term3(1:species) + rho*af2*sigma(1:species)*dxdt*&
         		(d2P_drho2 + (dP_drho/rho**2 - 2*P/rho**3 - d2e_drho2)*(dP_dT/de_dT) + &
         		(P/rho**2 - de_drho)*(d2P_drhoT/de_dT - dP_dT*d2e_drhoT/de_dT**2))/&
         		(ux*neta**2*chem_isos% W(chem_id(1:species)))
         	dfdx(species+1, species+1) = sum(term1) + sum(term2) + sum(term3)
         	
         	!d_drho (dT/dx):
         	!New term (with d/dzbar derivatives):
         	dfdx(species+2, species+1) = (-d2P_drho2*f(species+1) + &
         		(ux**2 - dP_drho)*dfdx(species+1, species+1) - &
         		abar*sum((-abar*d2P_drhoa + (chem_isos% Z(chem_id(1:species)) - zbar)*&
         		d2P_drhoz)*f(1:species)/chem_isos% W(chem_id(1:species)) + &
         		(-abar*dP_da + (chem_isos% Z(chem_id(1:species)) - zbar)*dP_dz)*&
         		d_dxdt_drho/(ux*chem_isos% W(chem_id(1:species)))))/dP_dT - &
         		d2P_drhoT/dP_dT**2*((ux**2 - dP_drho)*f(species+1) - &
         		sum((-abar**2*dP_da + abar*(chem_isos% Z(chem_id(1:species)) - zbar)*dP_dz)*&
         		f(1:species)/chem_isos% W(chem_id(1:species))))
         	!Old term (without d/dzbar derivatives):
         	!dfdx(species+2, species+1) = (-d2P_drho2*f(species+1) + &
         	!	(ux**2 - dP_drho)*dfdx(species+1, species+1) + &
         	!	abar**2*sum(d2P_drhoa*dxdt/(ux*chem_isos% W(chem_id(1:species))) + &
         	!	dP_da*d_dxdt_drho/(ux*chem_isos% W(chem_id(1:species)))))/dP_dT - &
         	!	d2P_drhoT*((ux**2 - dP_drho)*f(species+1) + abar**2*dP_da/ux*&
         	!	sum(dxdt/chem_isos% W(chem_id(1:species))))/dP_dT**2
         	
         	!d_drho (du/dx):
         	dfdx(species+3, species+1) = ux*f(species+1)/rho**2 - &
         		ux*dfdx(species+1, species+1)/rho
         		
         	!d_drho (dP/dx):
         	!dfdx(species+4, species+1) = -ux*f(species+3) - rho*ux*dfdx(species+3, species+1)
         		
         	!d_dT (dXi/dx):
         	dfdx(1:species, species+2) = d_dxdt_dT/ux
         	
         	!d_dT (drho/dx):
         	!Same as d_drho (drho/dx), but switch rho to T derivs
         	term1 = 0d0
         	term2 = 0d0
         	term3 = 0d0
         	!Here, the 1:species index is i in my notes (will be summed over at the end)
         	!New term (with d/dzbar terms):
         	term1(1:species) = term1(1:species) + (-abar**2*d2P_dTa + abar*&
         		(chem_isos% Z(chem_id(1:species)) - zbar)*d2P_dTz + (abar**2*d2E_dTa - &
         		abar*(chem_isos% Z(chem_id(1:species)) - zbar)*d2E_dTz)*dP_dT/dE_dT + &
         		(abar**2*de_da - abar*(chem_isos% Z(chem_id(1:species)) - zbar)*dE_dz + &
         		avo*mev_to_ergs*chem_isos% binding_energy(chem_id(1:species)))*&
         		(d2P_dT2/de_dt - dp_dt*d2e_dT2/de_dt**2))*&
         		(-f(1:species)/(neta*chem_isos% W(chem_id(1:species))))
         	!Old term (without d/dzbar terms):
         	!term1(1:species) = term1(1:species) - (-abar**2*d2P_dTa + &
         	!	abar**2*d2e_dTa*dP_dT/de_dT + (abar**2*de_da + &
         	!	avo*chem_isos% binding_energy(chem_id(1:species))*mev_to_ergs)*(d2P_dT2/de_dT - &
         	!	dP_dT*d2e_dT2/de_dT**2))*&
         	!	dxdt(1:species)/(ux*neta*chem_isos% W(chem_id(1:species)))
         	term2(1:species) = term2(1:species) - rho*af2*sigma(1:species)*d_dxdt_dT/&
         		(ux*neta*chem_isos% W(chem_id(1:species)))
         	term3(1:species) = term3(1:species) + rho*af2*sigma(1:species)*dxdt*&
         		(d2P_drhoT + (dP_dT/rho**2 - d2e_drhoT)*(dP_dT/de_dT) + &
         		(P/rho**2 - de_drho)*(d2P_dT2/de_dT - dP_dT*d2e_dT2/de_dT**2))/&
         		(ux*neta**2*chem_isos% W(chem_id(1:species)))
         	dfdx(species+1, species+2) = sum(term1) + sum(term2) + sum(term3)
         	
         	!d_dT (dT/dx):
         	!New term (with d/dzbar terms):
         	dfdx(species+2, species+2) = (-d2P_drhoT*f(species+1) + &
         		(ux**2 - dP_drho)*dfdx(species+1, species+2) - &
         		abar*sum((-abar*d2P_dTa + (chem_isos% Z(chem_id(1:species)) - zbar)*&
         		d2P_dTz)*f(1:species)/chem_isos% W(chem_id(1:species)) + &
         		(-abar*dP_da + (chem_isos% Z(chem_id(1:species)) - zbar)*dP_dz)*&
         		d_dxdt_dT/(ux*chem_isos% W(chem_id(1:species)))))/dP_dT - &
         		d2P_dT2/dP_dT**2*((ux**2 - dP_drho)*f(species+1) - &
         		sum((-abar**2*dP_da + abar*(chem_isos% Z(chem_id(1:species)) - zbar)*dP_dz)*&
         		f(1:species)/chem_isos% W(chem_id(1:species))))
         	!dfdx(species+2, species+2) = (-d2P_drhoT*f(species+1) + &
         	!	(ux**2 - dP_drho)*dfdx(species+1, species+2) - &
         	!	abar*sum(abar*d2P_dTa + (chem_isos% Z(chem_id(1:species)) - zbar)*&
         	!	d2P_dTz*f(1:species)/chem_isos% W(chem_id(1:species)) + &
         	!	(abar*dP_da + (chem_isos% Z(chem_id(1:species)) - zbar)*dP_dz)*&
         	!	d_dxdt_dT/(ux*chem_isos% W(chem_id(1:species)))))/dP_dT - &
         	!	d2P_dT2/dP_dT**2*((ux**2 - dP_drho)*f(species+1) - &
         	!	sum((-abar**2*dP_da + abar*(chem_isos% Z(chem_id(1:species)) - zbar)*dP_dz)*&
         	!	f(1:species)/chem_isos% W(chem_id(1:species))))
         	!Old term (without d/dzbar terms):
         	!dfdx(species+2, species+2) = (-d2P_drhoT*f(species+1) + &
         	!	(ux**2 - dP_drho)*dfdx(species+1, species+2) + &
         	!	abar**2*sum(d2P_dTa*dxdt/(ux*chem_isos% W(chem_id(1:species))) + &
         	!	dP_da*d_dxdt_dT/(ux*chem_isos% W(chem_id(1:species)))))/dP_dT - &
         	!	d2P_dT2*((ux**2 - dP_drho)*f(species+1) + abar**2*dP_da/ux*&
         	!	sum(dxdt/chem_isos% W(chem_id(1:species))))/dP_dT**2
         		
         	!d_dT (du/dx):
         	dfdx(species+3, species+2) = -ux*dfdx(species+1, species+2)/rho
         	
         	!d_dT (dP/dx):
         	!dfdx(species+4, species+2) = -rho*ux*dfdx(species+3, species+2)
         		
         	!d/du (dXi_dx):
         	dfdx(1:species, species+3) = -dxdt(1:species)/ux**2
         		
         	!d/du (drho_dx):
         	dfdx(species+1, species+3) = f(species+1)*(2*ux/neta - 1/ux)
         	!dfdx(species+1, species+3) = -2*rho*af2/(neta**2)*&
         	!	sum(sigma*dxdt/chem_isos% W(chem_id(1:species)))
         		
         	!d/du (dT_dx):
         	!New term (with d/dzbar terms):
         	dfdx(species+2, species+3) = (2*ux*f(species+1) + (ux**2 - dP_drho)*&
         		dfdx(species+1, species+3) - abar*sum((-abar*dP_da + &
         		(chem_isos% Z(chem_id(1:species)) - zbar)*dp_dz)/&
         		chem_isos% W(chem_id(1:species))*dfdx(1:species, species+3)))/dP_dT
         	!Old term (without d/dzbar terms):
         	!dfdx(species+2, species+3) = (2*ux*f(species+1) - abar**2*dP_da/ux**2*&
         	!	sum(dxdt/chem_isos% W(chem_id(1:species))) + &
         	!	(ux**2 - dP_drho)*dfdx(species+1, species+3) + &
         	!	abar**2*dP_da/ux*sum(dfdx(1:species, species+3)/&
         	!	chem_isos% W(chem_id(1:species))))/dP_dT
         		
         	!d/du (du_dx):
         	dfdx(species+3, species+3) = -f(species+1)/rho - ux*dfdx(species+1, species+3)/rho
         	
         	!d/du (dP/dx):
         	!dfdx(species+4, species+3) = -rho*f(species+3) - rho*ux*dfdx(species+3, species+3)
         	
         	!If we're using the blowout equations, we'll need to add on the derivatives 
         	!of the extra terms to what we've already computed:
         	if(do_blowout) then
         		f1 = f1_blowout
         		f2 = f2_blowout
         		f3 = f3_blowout
         		
         		!First, just the blowout terms:
         		df1_dXi(1:species) = 0d0
         		df2_dXi(1:species) = -uy*dP_dXi/(H*ux)
         		df3_dXi(1:species) = -uy*dP_dXi/(4*rho*H*ux)
         		
         		df1_drho = -uy/H
         		df2_drho = -ux*uy/H
         		df3_drho = P*uy/(4*rho**2*H*ux)
         		
         		df1_dT = 0d0
         		df2_dT = -uy*dP_dT/(H*ux)
         		df3_dT = -uy*dP_dT/(4*rho*H*ux)
         		
         		df1_dux = 0d0
         		df2_dux = P*uy/(H*ux**2) - rho*uy/H
         		df3_dux = P*uy/(4*rho**H*ux**2) - uy**3/(4*H*ux**2)
         		
         		df1_duy = -rho/H
         		df2_duy = -(P + rho*ux**2)/(H*ux)
         		df3_duy = -P/(4*rho*H*ux) + 3*uy**2/(4*H*ux)
         		
         		df1_dH = rho*uy/H**2
         		df2_dH = (P + rho*ux**2)*uy/(H**2*ux)
         		df3_dH = P*uy/(4*rho*H**2*ux) - uy**3/(4*H**2*ux)
         		
        		if(do_curvature) then
        			!Add on the curvature source terms to the total:
        			f1 = f1 + f1_curvature
        			f2 = f2 + f2_curvature
        			f3 = f3 + f3_curvature
        			
        			!Calculate the derivatives of the curvature source terms:
        			df1_curvature_dXi(1:species) = 0d0
        			df2_curvature_dXi(1:species) = 0d0 
        			df3_curvature_dXi(1:species) = 0d0
        			
        			df1_curvature_drho = -j_curve*(v_det - ux)/r_curve
        			df2_curvature_drho = -j_curve*ux*(v_det - ux)/r_curve
        			df3_curvature_drho = 0d0
        			
        			df1_curvature_dT = 0d0
        			df2_curvature_dT = 0d0
        			df3_curvature_dT = 0d0
        			
        			df1_curvature_dux = j_curve*rho/r_curve
        			df2_curvature_dux = -j_curve*rho*(v_det - 2*ux)/r_curve
        			df3_curvature_dux = 0d0
        			
        			df1_curvature_duy = 0d0
        			df2_curvature_duy = 0d0
        			df3_curvature_duy = 0d0
        			
        			if(use_rc_hs_scaling) then
        				df1_curvature_dH = -j_curve*rho*(v_det - ux)/(r_curve**2/rc_hs_factor)
        				df2_curvature_dH = -j_curve*rho*ux*(v_det - ux)/(r_curve**2/rc_hs_factor)
        				df3_curvature_dH = 0d0
        			else
        				df1_curvature_dH = 0d0
        				df2_curvature_dH = 0d0
        				df3_curvature_dH = 0d0
        			endif
        			
        			!Add on the curvature terms to the previously computed derivatives:
        			df1_dXi(1:species) = df1_dXi + df1_curvature_dXi
        			df2_dXi(1:species) = df2_dXi + df2_curvature_dXi
        			df3_dXi(1:species) = df3_dXi + df3_curvature_dXi
        			
        			df1_drho = df1_drho + df1_curvature_drho
        			df2_drho = df2_drho + df2_curvature_drho
        			df3_drho = df3_drho + df3_curvature_drho
        			
        			df1_dT = df1_dT + df1_curvature_dT
        			df2_dT = df2_dT + df2_curvature_dT
        			df3_dT = df3_dT + df3_curvature_dT
        			
        			df1_dux = df1_dux + df1_curvature_dux
        			df2_dux = df2_dux + df2_curvature_dux
        			df3_dux = df3_dux + df3_curvature_dux
        			
        			df1_duy = df1_duy + df1_curvature_duy
        			df2_duy = df2_duy + df2_curvature_duy
        			df3_duy = df3_duy + df3_curvature_duy
        			
        			df1_dH = df1_dH + df1_curvature_dH
        			df2_dH = df2_dH + df2_curvature_dH
        			df3_dH = df3_dH + df3_curvature_dH
        		endif
         		
         		!New equations (7/10/13):
         		!d_dXj (drho_dx):
         		dfdx(species+1, 1:species) = dfdx(species+1, 1:species) - &
         			(f2 - 2*ux*f1 - dP_de_rho*(f3 - f2/rho + ux*f1/rho))*dneta_dx/neta**2 + &
         			(df2_dXi - 2*ux*df1_dXi - d_dP_de_rho_dXi*(f3 - f2/rho + ux*f1/rho) - &
         			dP_de_rho*(df3_dXi - df2_dXi/rho + ux*df1_dXi/rho))/neta
         		!d_dXj (dT_dx):
         		term1 = 0d0		!Sum along one dimension of multidimensional array products
         		do j=1,species
         			term1 = term1 + d2P_dXidXj(1:species,j)*f(j) + dP_dXi(j)*dfdx(j, 1:species)
         		end do
         		dfdx(species+2, 1:species) = -d2P_dTdXi*f(species+2)/dP_dT + &
         			(dfdx(species+1, 1:species)*(ux**2 - dP_drho) - &
         			f(species+1)*d2P_drhodXi - term1 + df2_dXi - 2*ux*df1_dXi)/dP_dT
         		!d_dXj (dux_dx):
         		dfdx(species+3, 1:species) = -ux*dfdx(species+1, 1:species)/rho + &
         			df1_dXi/rho
         		!d_dXj (duy_dx):
         		dfdx(species+4, 1:species) = dP_dXi/(rho*ux*H)
         		!d_dXj (dH_dx):
         		dfdx(species+5, 1:species) = 0d0
         		
         		!d_drho (drho_dx):
         		dfdx(species+1, species+1) = dfdx(species+1, species+1) + &
         			(df2_drho - 2*ux*df1_drho - d_dP_de_rho_drho*(f3 - f2/rho + ux*f1/rho) - &
         			dP_de_rho*(df3_drho + f2/rho**2 - df2_drho/rho - ux*f1/rho**2 + &
         			ux*df1_drho/rho))/neta - (f2 - 2*ux*f1 - &
         			dP_de_rho*(f3 - f2/rho + ux*f1/rho))*dneta_drho/neta**2
         		!d_drho (dT_dx):
         		dfdx(species+2, species+1) = -d2P_drhoT*f(species+2)/dP_dT + &
         			(dfdx(species+1, species+1)*(ux**2 - dP_drho) - f(species+1)*d2P_drho2 - &
         			(dot_product(d2P_drhodXi,f(1:species)) + &
         			dot_product(dP_dXi,dfdx(1:species,species+1))) + &
         			df2_drho - 2*ux*df1_drho)/dP_dT
         		!d_drho (dux_dx):
         		dfdx(species+3, species+1) = ux*f(species+1)/rho**2 - &
         			ux*dfdx(species+1, species+1)/rho - f1/rho**2 + df1_drho/rho
         		!d_drho (duy_dx):
         		dfdx(species+4, species+1) = dP_drho/(rho*ux*H) - P/(rho**2*ux*H)
         		!d_drho (dH_dx):
         		dfdx(species+5, species+1) = 0d0
         		
         		!d_dT (drho_dx):
         		dfdx(species+1, species+2) = dfdx(species+1, species+2) + &
         			(df2_dT - 2*ux*df1_dT - d_dP_de_rho_dT*(f3 - f2/rho + ux*f1/rho) - &
         			dP_de_rho*(df3_dT - df2_dT/rho + ux*df1_dT/rho))/neta - &
         			(f2 - 2*ux*f1 - dP_de_rho*(f3 - f2/rho + ux*f1/rho))*dneta_dT/neta**2
         		!d_dT (dT_dx):
         		dfdx(species+2, species+2) = -d2P_dT2*f(species+2)/dP_dT + &
         			(dfdx(species+1, species+2)*(ux**2 - dP_drho) - f(species+1)*d2P_drhoT - &
         			(dot_product(d2P_dTdXi,f(1:species)) + &
         			dot_product(dP_dXi,dfdx(1:species, species+2))) + &
         			df2_dT - 2*ux*df1_dT)/dP_dT
         		!d_dT (dux_dx):
         		dfdx(species+3, species+2) = -ux*dfdx(species+1, species+2)/rho + &
         			df1_dT/rho
         		!d_dT (duy_dx):
         		dfdx(species+4, species+2) = dP_dT/(rho*ux*H)
         		!d_dT (dH_dx):
         		dfdx(species+5, species+2) = 0d0
         		
         		!d_dux (drho_dx):
         		dfdx(species+1, species+3) = dfdx(species+1, species+3) + &
         			(df2_dux - 2*f1 - 2*ux*df1_dux - dP_de_rho*(df3_dux - &
         			df2_dux/rho + f1/rho + ux*df1_dux/rho))/neta - &
         			(f2 - 2*ux*f1 - dP_de_rho*(f3 - f2/rho + ux*f1/rho))*dneta_dux/neta**2
         		!d_dux (dT_dx):
         		dfdx(species+2, species+3) = (dfdx(species+1, species+3)*(ux**2 - dP_drho) - &
         			dot_product(dP_dXi,dfdx(1:species, species+3)) + &
         			2*ux*f(species+1) + df2_dux - 2*f1 - 2*ux*df1_dux)/dP_dT
         		!d_dux (dux_dx):
         		dfdx(species+3, species+3) = -ux*dfdx(species+1, species+3)/rho - &
         			f(species+1)/rho + df1_dux/rho
         		!d_dux (duy_dx):
         		dfdx(species+4, species+3) = -P/(rho*ux**2*H) + uy**2/(ux**2*H)
         		!d_dux (dH_dx):
         		dfdx(species+5, species+3) = -uy/ux**2
         		
         		!d_duy (drho_dx):
         		dfdx(species+1, species+4) = (df2_duy - 2*ux*df1_duy - dP_de_rho*(df3_duy - &
         			df2_duy/rho + ux*df1_duy/rho))/neta
         		!d_duy (dT_dx):
         		dfdx(species+2, species+4) = (dfdx(species+1, species+4)*(ux**2 - dP_drho) + &
         			df2_duy - 2*ux*df1_duy)/dP_dT
         		!d_duy (dux_dx):
         		dfdx(species+3, species+4) = -ux*dfdx(species+1, species+4)/rho + df1_duy/rho
         		!d_duy (duy_dx):
         		dfdx(species+4, species+4) = -2*uy/(ux*H)
         		!d_duy (dH_dx):
         		dfdx(species+5, species+4) = 1d0/ux
         		
         		!d_dH (drho_dx):
         		dfdx(species+1, species+5) = (df2_dH - 2*ux*df1_dH - dP_de_rho*(df3_dH - &
         			df2_dH/rho + ux*df1_dH/rho))/neta
         		!d_dH (dT_dx):
         		dfdx(species+2, species+5) = (dfdx(species+1, species+5)*(ux**2 - dP_drho) + &
         			df2_dH - 2*ux*df1_dH)/dP_dT
         		!d_dH (dux_dx):
         		dfdx(species+3, species+5) = -ux*dfdx(species+1, species+5)/rho + df1_dH/rho
         		!d_dH (duy_dx):
         		dfdx(species+4, species+5) = -P/(rho*ux*H**2) + uy**2/(ux*H**2)
         		!d_dH (dH_dx):
         		dfdx(species+5, species+5) = 0d0
         		
         	
         		!Old ish that wasn't working:
         		!d_dXj (drho_dx):
         		!dfdx(species+1, 1:species) = dfdx(species+1, 1:species) + &
!         			uy*(P - rho*ux**2 - P/rho*dP_dT/de_dT)*dneta_dx/(H*ux*neta**2) - &
!         			uy/(H*ux*neta)*abar/chem_isos% W(chem_id(1:species))*&
!         			((-abar)*(dP_da - dP_da*dP_dT/(rho*de_dT) - P/rho*(d2P_dTa/de_dT - &
!         			dP_dT*d2e_dTa/de_dT**2)) + &
!         			(chem_isos% Z(chem_id(1:species)) - zbar)*(dP_dz - &
!         			dP_dz*dP_dT/(rho*de_dT) - P/rho*(d2P_dTz/de_dT - dP_dT*d2e_dTz/de_dT**2)))
!         		!d_drho (drho_dx):
!         		dfdx(species+1, species+1) = dfdx(species+1, species+1) + &
!         			uy*(P - rho*ux**2 - P/rho*dP_dT/de_dT)*dneta_drho/(H*ux*neta**2) - &
!         			uy/(H*ux*neta)*(dP_drho - ux**2 - (dP_drho/rho - P/rho**2)*dP_dT/de_dT - &
!         			P/rho*(d2P_drhoT/de_dT - dP_dT*d2e_drhoT/de_dT**2))
!         		!d_dT (drho_dx):
!         		dfdx(species+1, species+2) = dfdx(species+1, species+2) + &
!         			uy*(P - rho*ux**2 - P/rho*dP_dT/de_dT)*dneta_dT/(H*ux*neta**2) - &
!         			uy/(H*ux*neta)*(dP_dT - dP_dT/rho*dP_dT/de_dT - &
!         			P/rho*(d2P_dT2/de_dT - dP_dT*d2e_dT2/de_dT**2))
!         		!d_dux (drho_dx):
!         		dfdx(species+1, species+3) = dfdx(species+1, species+3) + &
!         			uy*(P - rho*ux**2 - P/rho*dP_dT/de_dT)*(1d0/ux - 2*ux/neta)/(H*ux*neta) + &
!         			2*rho*uy/(H*neta)
!         		!d_duy (drho_dx):
!         		dfdx(species+1, species+4) = -(P - rho*ux**2 - P/rho*dP_dT/de_dT)/(H*ux*neta)
!         		!d_dH (drho_dx):
!         		dfdx(species+1, species+5) = -uy*(P - rho*ux**2 - P/rho*dP_dT/de_dT)/(H**2*ux*neta)
!         		
!         		!d_dXj (dT_dx):
!         		dfdx(species+2, 1:species) = dfdx(species+2, 1:species) + &
!         			(P - rho*ux**2)*uy/(H*ux*dP_dT**2*chem_isos% W(chem_id(1:species)))*&
!         			(abar**2*d2P_dTa + abar*(chem_isos% Z(chem_id(1:species)) - zbar)*d2P_dTz)
!         		!d_drho (dT_dx):
!         		dfdx(species+2, species+1) = dfdx(species+2, species+1) + &
!         			ux*uy/(H*dP_dT) + (P - rho*ux**2)*uy*d2P_drhoT/(H*ux*dP_dT**2)
!         		!d_dT (dT_dx):
!         		dfdx(species+2, species+2) = dfdx(species+2, species+2) + &
!         			(P - rho*ux**2)*uy*d2P_dT2/(H*ux*dP_dT**2)
!         		!d_dux (dT_dx):
!         		dfdx(species+2, species+3) = dfdx(species+2, species+3) + &
!         			(P - rho*ux**2)*uy*d2P_dT2/(H*ux**2*dP_dT) + 2*rho*uy/(H*dP_dT)
!         		!d_duy (dT_dx):
!         		dfdx(species+2, species+4) = -(P - rho*ux**2)/(H*ux*dP_dT)
!         		!d_dH (dT_dx):
!         		dfdx(species+2, species+5) = (P - rho*ux**2)*uy/(H**2*ux*dP_dT)
!         		
!         		!d_duy (dux_dx):
!         		dfdx(species+3, species+4) = 1d0/H
!         		!d_dH (dux_dx):
!         		dfdx(species+3, species+5) = -uy/H**2
!         		
!         		!d_dXj (duy_dx):
!         		dfdx(species+4, 1:species) = abar/(rho*H*ux*chem_isos% W(chem_id(1:species)))*&
!         			(-abar*dP_da + (chem_isos% Z(chem_id(1:species)) - zbar)*dP_dz)
!         		!d_drho (duy_dx):
!         		dfdx(species+4, species+1) = dP_drho/(rho*H*ux) - P/(rho**2*H*ux)
!         		!d_dT (duy_dx):
!         		dfdx(species+4, species+2) = dP_dT/(rho*H*ux)
!         		!d_dux (duy_dx):
!         		!With gravity:
!         		dfdx(species+4, species+3) = -P/(rho*H*ux**2) + standard_cgrav*m_c/&
!         			(ux**2*(r_wd + H/2d0)**2)
!         		!Without gravity:
!         		!dfdx(species+4, species+3) = -P/(rho*H*ux**2) + uy**2/(H*ux**2)
!         		!d_duy (duy_dx):
!         		!With gravity:
!         		dfdx(species+4, species+4) = 0d0
!         		!Without gravity:
!         		!dfdx(species+4, species+4) = -2*uy/(H*ux)
!         		!d_dH (duy_dx):
!         		!With gravity:
!         		dfdx(species+4, species+5) = -P/(rho*H**2*ux) + standard_cgrav*m_c/(ux*(r_wd + H/2d0)**3)
!         		!Without gravity:
!         		!dfdx(species+4, species+5) = -P/(rho*H**2*ux) + uy**2/(H**2*ux)
!         			
!         		!d_dux (dH_dx):
!         		dfdx(species+5, species+3) = -uy/ux**2
!         		!d_duy (dH_dx):
!         		dfdx(species+5, species+4) = 1d0/ux
         	endif !do_blowout
         endif
         
         
         rpar(1) = P 					    !Output total pressure here
         rpar(2) = ux/sqrt(af2) 			!Output the local mach number here
         rpar(3) = ux*f(species+1) + rho*f(species+3)
         rpar(4) = rho*ux*f(species+3) + dp_drho*f(species+1) + dp_dT*f(species+2) - &
         	abar**2*dp_da*sum(f(1:species)/chem_isos% W(chem_id(1:species)))
         !Of course, this is the same as:
         !rpar(4) = rho*ux*f(species+3) + dp_drho*f(species+1) + dp_dT*f(species+2) - &
         !	abar**2*dp_da/ux*sum(dxdt/chem_isos% W(chem_id(1:species)))
         rpar(5) = de_drho*f(species+1) + de_dT*f(species+2) - &
         	abar**2*de_da*sum(f(1:species)/chem_isos% W(chem_id(1:species))) - &
         	avo*sum(f(1:species)*chem_isos% binding_energy(chem_id(1:species))*mev_to_ergs/&
         	chem_isos% W(chem_id(1:species))) - P*f(species+1)/rho**2
         
         !Extra stuff for diagnosing the problems with the 2nd equation (rho*ux*drho/dx + dP/dx = 0):	
         !rpar(6) = rho*ux*f(species+3) + dp_drho*f(species+1) - &
         !	abar**2*dp_da*sum(f(1:species)/chem_isos% W(chem_id(1:species)))
         !rpar(7) = dP_dT*f(species+2)
         
         !Extra stuff for comparing dP/dx (EOS) to dP_dx (eqn II)
         !Old (only d/dabar term) dP/dx term from EOS:
         !rpar(6) = dP_drho*f(species+1) + dP_dT*f(species+2) - &
         !	abar**2*dp_da*sum(f(1:species)/chem_isos% W(chem_id(1:species)))
         !New (d/dabar and d/dzbar with constant other) dP/dx term from EOS:
         !rpar(6) = dP_drho*f(species+1) + dP_dT*f(species+2) + &
         !	abar*sum((-abar*dp_da/zbar + dp_dz)*f(1:species)*&
         !	chem_isos% Z(chem_id(1:species))/chem_isos% W(chem_id(1:species)))
         !New (full d/dabar and d/dzbar terms) dP/dx term from EOS:
         rpar(6) = dP_drho*f(species+1) + dP_dT*f(species+2) + &
         	sum((-abar**2*dp_da + (abar*chem_isos% Z(chem_id(1:species)) - &
         	abar*zbar)*dp_dz)*f(1:species)/chem_isos% W(chem_id(1:species)))
         rpar(7) = -rho*ux*f(species+3)
         
         rpar(8) = energy - q_bind
         rpar(9) = q_bind
         rpar(10) = gamma1
         
         !rpar(11) = f1
         !rpar(12) = f2
         !rpar(13) = f3
         !rpar(14) = 2*ux*f1
         !rpar(15) = f2 - 2*ux*f1
         !rpar(16) = -dP_dT/dE_dT*(f3 - f2/rho + ux*f1/rho)
         
         rpar(11) = f1_blowout
         rpar(12) = f2_blowout
         rpar(13) = f3_blowout
         rpar(14) = f1_curvature
         rpar(15) = f2_curvature
         rpar(16) = f3_curvature

			!Save the derivatives for use in linearization around the pathological point:
			rpar(21:21+num_vars-1) = f(1:num_vars)

         !Don't use these with large nets (>approx19)
         !rpar(11:11+num_vars-1) = f
         !rpar(8+num_vars) = zbar
         
         !More consistency checks:
         dP_dx = dP_drho*f(species+1) + dP_dT*f(species+2) + &
         	sum((-abar**2*dp_da + (abar*chem_isos% Z(chem_id(1:species)) - &
         	abar*zbar)*dp_dz)*f(1:species)/chem_isos% W(chem_id(1:species)))
         dq_dx = avo*sum(chem_isos% binding_energy(chem_id(1:species))*mev_to_ergs*&
        	f(1:species)/chem_isos% W(chem_id(1:species)))
         de_dx = de_drho*f(species+1) + de_dT*f(species+2) + &
         	sum((-abar**2*de_da + (abar*chem_isos% Z(chem_id(1:species)) - &
         	abar*zbar)*de_dz)*f(1:species)/chem_isos% W(chem_id(1:species)))
         detot_dx = de_dx - dq_dx
         dEtot_dx_eqn = P*f(species+1)/rho**2
                           
         !write(*,*) 'f:', f
         !write(*,*) 'dfdx:', dfdx
         !write(*,*)
        
        !write(*,*) 'dfdx(species+3, species+3) = ',dfdx(species+3, species+3)
       	!write(*,*) 'leaving znd_jacob'
         
      end subroutine znd_jacob
      
      subroutine znd_solout(nr, xold, x, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
         ! nr is the step number.
         ! x is the current x value; xold is the previous x value.
         ! y is the current y value.
         ! irtrn negative means terminate integration.
         ! rwork_y and iwork_y hold info for interp_y
         ! note that these are not the same as the rwork and iwork arrays for the solver.
         use const_def, only: dp
         integer, intent(in) :: nr, n, lrpar, lipar
         real(dp), intent(in) :: xold, x
         real(dp), intent(inout) :: y(:)
         ! y can be modified if necessary to keep it in valid range of possible solutions.
         real(dp), intent(inout), target :: rwork_y(*)
         integer, intent(inout), target :: iwork_y(*)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         interface
            include 'num_interp_y.dek'
         end interface
         integer, intent(out) :: irtrn ! < 0 causes solver to return to calling program.
         
         !write(*,*) 'Step: ', nr
         !write(*,*) 'Delta x: ',x-xold
         
      end subroutine znd_solout
      
      !----------End ZND subroutines
      
      subroutine constp_derivs(n, t, h, x, f, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: t, h
         real(dp), intent(inout) :: x(:)
         real(dp), intent(out) :: f(:) ! dxdt
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         integer, parameter :: ld_dfdx = 0
         real(dp) :: dfdx(ld_dfdx,n)
         ierr = 0
         
         !write(*,*) 'entering constp_derivs'
         call constp_jacob(n, t, h, x, f, dfdx, ld_dfdx, lrpar, rpar, lipar, ipar, ierr)
         !write(*,*) 'leaving constp_derivs'
      end subroutine constp_derivs


      subroutine constp_jacob(n, time, h, x, f, dfdx, ld_dfdx, lrpar, rpar, lipar, ipar, ierr)
      
      	implicit none
         
         integer, intent(in) :: n, ld_dfdx, lrpar, lipar
         real(dp), intent(in) :: time, h
         real(dp), intent(inout) :: x(:)
         real(dp), intent(out) :: f(:), dfdx(:,:)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         
         double precision :: composition(species)
         double precision :: rho, cp
         integer :: i

         logical :: just_dxdt
         type (Net_Info), target :: netinfo_target
         type (Net_Info), pointer :: netinfo
         
         include 'formats.dek'

         netinfo => netinfo_target
         
         !write(*,*) 'entering constp_jacob'
         
         ierr = 0
         f = 0
         dfdx = 0
         
         !Normalize all the abundance pieces of the variable array:
         do i=1,species
         	x(i) = max(0d0, min(1d0, x(i)))
         end do
         
         composition = 0d0
         composition = x(1:species) !Array of just the abundances
         
         !write(*,*) composition(1)
         
         call composition_info(species, chem_id, composition, xh, xhe, zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
         
         rho = x(species+1) !local variable
            
        logPgas = log10(burn_constPgas)
        call eosPT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, xa, &
        	burn_constPgas, logPgas, burn_T, logT, rho, logRho, &
        	dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
     		res, d_dlnRho_const_T, d_dlnT_const_Rho, &
     		d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
     	if (ierr /= 0) then
            write(*,*) 'failed in eosPT_get', ierr
        endif
        cp = res(i_Cp)
        
        call eos_get_helm_results(Xh, abar, zbar, &
        	Rho, logRho, burn_T, logT, .true., .false., .false., res_helm, ierr)
     	if (ierr /= 0) then
            write(*,*) 'failed in eosPT_get', ierr
        endif

		!write(*,*) x(1), sum(composition)
		
        just_dxdt = .false.
     	call net_get( &
            handle, just_dxdt, netinfo, species, num_reactions, &
            composition, burn_T, logT, burn_Rho, logRho, & 
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, & 
            std_reaction_Qs, std_reaction_neuQs, .false., .false., &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, & 
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, & 
            screening_mode, theta_e_for_graboske_et_al, &     
            eps_nuc_categories, eps_neu_total, & 
            lwork_net, work_net, ierr)
        if (ierr /= 0) then
            write(*,*) 'failed in net_get', ierr
            write(*,*) x, sum(composition), lwork_net
            write(*,*) 
            !write(*,*) handle, species, num_reactions, &
            !	x(1), burn_T, logT, burn_Rho, logRho, & 
            !	abar, zbar, z2bar, ye, eta
            return
         end if
         
         !write(*,*) 'net_get successful'
         !write(*,*) dxdt
         !write(*,*)
     
         if (.false. .and. ierr /= 0) then
            write(*,1) 'time', time
            write(*,1) 'T', 10**logT
            write(*,1) 'lgT', logT
            write(*,1) 'rho', 10**logRho
            write(*,1) 'lgRho', logRho
            write(*,1) 'eta', eta
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,1) 'z2bar', z2bar
            write(*,1) 'ye', ye
            write(*,1) 'xh', xh
            write(*,1) 'xhe', xhe
            write(*,1) 'sum(x)', sum(x)
            write(*,1) 'maxval(x)', maxval(x)
            write(*,1) 'minval(x)', minval(x)
            write(*,1) 'eps_nuc', eps_nuc
            write(*,1) 'eps_neu_total', eps_neu_total
            write(*,*)
            if (time > 0) stop
            !write(*,2) 'eval_net ierr', ierr
            !write(*,*) 'failed in constp_jacob'
            !stop '1 zone burn constp_jacob'
         end if

         f(1:species) = dxdt
         f(species+1) = eps_nuc/cp !dT/dt
         
         !write(*,*) f(species+1)
         !write(*,*) size(f), size(dfdx)
         
         if (ld_dfdx > 0) then
         	!d/dXi (dXj/dt)
         	dfdx(1:species, 1:species) = d_dxdt_dx
         	
         	!d/dXi (dT/dt):
         	dfdx(species+1, 1:species) = d_eps_nuc_dx/cp + &
         		eps_nuc*abar**2*res_helm(h_dcpda)/(cp**2*chem_isos %W(chem_id(1:species)))
         		
         	!d/dT (dT/dt):
         	dfdx(species+1, species+1) = d_eps_nuc_dT/cp - &
         		eps_nuc*d_dlnT_const_Rho(i_Cp)/(burn_T*cp**2)
         		
         	!d/dT (dXi/dt):
         	dfdx(1:species, species+1) = d_dxdt_dT
         endif
         
       	!write(*,*) 'leaving constp_jacob'
         
      end subroutine constp_jacob
      
      !End Constant P routines
      !------
      
      
      subroutine burn_derivs(n, t, h, x, f, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: t, h
         real(dp), intent(inout) :: x(:)
         real(dp), intent(out) :: f(:) ! dxdt
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         integer, parameter :: ld_dfdx = 0
         real(dp) :: dfdx(ld_dfdx,n)
         ierr = 0
         
         !write(*,*) 'entering burn_derivs'
         call burn_jacob(n, t, h, x, f, dfdx, ld_dfdx, lrpar, rpar, lipar, ipar, ierr)
         !write(*,*) 'leaving burn_derivs'
      end subroutine burn_derivs


      subroutine burn_jacob(n, time, h, x, f, dfdx, ld_dfdx, lrpar, rpar, lipar, ipar, ierr)
      
      	implicit none
         
         integer, intent(in) :: n, ld_dfdx, lrpar, lipar
         real(dp), intent(in) :: time, h
         real(dp), intent(inout) :: x(:)
         real(dp), intent(out) :: f(:), dfdx(:,:)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         
         double precision :: composition(species)
         integer :: i

         logical :: just_dxdt = .false.
         type (Net_Info), target :: netinfo_target
         type (Net_Info), pointer :: netinfo
         
         include 'formats.dek'

         netinfo => netinfo_target
         
         !write(*,*) 'entering burn_jacob'
         
         ierr = 0
         f = 0
         dfdx = 0
         
         !Normalize all the abundance pieces of the variable array:
         do i=1,species
         	x(i) = max(0d0, min(1d0, x(i)))
         end do
         
         composition = 0d0
         composition = x(1:species) !Array of just the abundances
         
         !write(*,*) composition(1)
         
         call composition_info(species, chem_id, composition, xh, xhe, zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)

		!write(*,*) x(1), sum(composition)
		
     	call net_get( &
            handle, just_dxdt, netinfo, species, num_reactions, &
            composition, burn_T, logT, burn_Rho, logRho, & 
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, & 
            std_reaction_Qs, std_reaction_neuQs, .false., .false., &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, & 
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, & 
            screening_mode, theta_e_for_graboske_et_al, &     
            eps_nuc_categories, eps_neu_total, & 
            lwork_net, work_net, ierr)
        if (ierr /= 0) then
            write(*,*) 'failed in net_get', ierr
            write(*,*) x, sum(composition), lwork_net
            write(*,*) 
            !write(*,*) handle, species, num_reactions, &
            !	x(1), burn_T, logT, burn_Rho, logRho, & 
            !	abar, zbar, z2bar, ye, eta
            return
         end if
         
         !write(*,*) 'net_get successful'
         !write(*,*) dxdt
         !write(*,*)
     
         if (.false. .and. ierr /= 0) then
            write(*,1) 'time', time
            write(*,1) 'T', 10**logT
            write(*,1) 'lgT', logT
            write(*,1) 'rho', 10**logRho
            write(*,1) 'lgRho', logRho
            write(*,1) 'eta', eta
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,1) 'z2bar', z2bar
            write(*,1) 'ye', ye
            write(*,1) 'xh', xh
            write(*,1) 'xhe', xhe
            write(*,1) 'sum(x)', sum(x)
            write(*,1) 'maxval(x)', maxval(x)
            write(*,1) 'minval(x)', minval(x)
            write(*,1) 'eps_nuc', eps_nuc
            write(*,1) 'eps_neu_total', eps_neu_total
            write(*,*)
            if (time > 0) stop
            !write(*,2) 'eval_net ierr', ierr
            !write(*,*) 'failed in burn_jacob'
            !stop '1 zone burn burn_jacob'
         end if

         f(1:species) = dxdt
         !f(species+1) = 0d0
         !f(species+2) = 0d0
         f(species+1) = x(species+1)/1d6
         f(species+2) = 2d0*(time/1d7)
         
         !write(*,*) size(f), size(dfdx)
         
         if (ld_dfdx > 0) then
         	dfdx(1:species, 1:species) = d_dxdt_dx
         	dfdx(species+1,:) = 0d0
         	dfdx(:,species+1) = 0d0
         	dfdx(species+1,species+1) = 1/1d6
         	dfdx(species+2,:) = 0d0
         	dfdx(:,species+2) = 0d0
         endif
         
       	!write(*,*) 'leaving burn_jacob'
         
      end subroutine burn_jacob
      
      subroutine set_Aptr(Aptr, dest, n1, n2)
            real(dp), pointer :: Aptr(:, :)
            real(dp), target :: dest(n1, n2) ! reshape work section
            integer, intent(in) :: n1, n2
            Aptr => dest
      end subroutine set_Aptr


      subroutine burn_sjac(n,time,h,y,f,nzmax,ia,ja,values,lrpar,rpar,lipar,ipar,ierr)  
         use mtx_lib, only: dense_to_sparse_with_diag
         use mtx_def
         integer, intent(in) :: n, nzmax, lrpar, lipar
         real(dp), intent(in) :: time, h
         real(dp), intent(inout) :: y(:)
         integer, intent(out) :: ia(:), ja(:)
         real(dp), intent(out) :: f(:), values(:)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means terminate integration
         real(dp), pointer :: dfdv(:,:) ! (n,n)
         integer :: ld_dfdv, nz, i, j, cnt, nnz
      	include 'formats.dek'
      	
      	!write(*,*) 'entering burn_sjac'
      	
      	!write(*,1) 'burn_sjac', x
      	ierr = 0
         ld_dfdv = n
         allocate(dfdv(n,n),stat=ierr)
         if (ierr /= 0) return
         call burn_jacob(n,time,h,y,f,dfdv,ld_dfdv,lrpar,rpar,lipar,ipar,ierr)
         if (ierr /= 0) then
            deallocate(dfdv)
            return
         end if
         ! remove entries with abs(value) < 1d-16
         cnt = 0; nnz = 0
         do i=1,n
            do j=1,n
               if (dfdv(i,j) /= 0) then
                  nnz = nnz + 1
                  if (abs(dfdv(i,j)) < 1d-16) then
                     cnt = cnt+1; dfdv(i,j) = 0
                  end if
               end if
            end do
         end do
         call dense_to_sparse_with_diag(ipar(i_sparse_format),n,n,dfdv,nzmax,nz,ia,ja,values,ierr)
         deallocate(dfdv)
      	!write(*,2) 'done burn_sjac: nz', nz
      	
      	!write(*,*) 'leaving burn_sjac'
      	
      end subroutine burn_sjac
      
      subroutine burn_solout(nr, xold, x, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
         ! nr is the step number.
         ! x is the current x value; xold is the previous x value.
         ! y is the current y value.
         ! irtrn negative means terminate integration.
         ! rwork_y and iwork_y hold info for interp_y
         ! note that these are not the same as the rwork and iwork arrays for the solver.
         use const_def, only: dp
         integer, intent(in) :: nr, n, lrpar, lipar
         real(dp), intent(in) :: xold, x
         real(dp), intent(inout) :: y(n)
         ! y can be modified if necessary to keep it in valid range of possible solutions.
         real(dp), intent(inout), target :: rwork_y(*)
         integer, intent(inout), target :: iwork_y(*)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         interface
            include 'num_interp_y.dek'
         end interface
         integer, intent(out) :: irtrn ! < 0 causes solver to return to calling program.
         
         write(*,*) 'step', nr
         write(*,*) y
      end subroutine burn_solout
      
      !End burn support routines (jacobian, derivs, etc.)
      !-----------------------------------------------------------------------------------

      subroutine burn_finish_substep(nstp, time, y, ierr)
         integer,intent(in) :: nstp
         double precision, intent(in) :: time, y(:)
         integer, intent(out) :: ierr
      end subroutine burn_finish_substep

      
           
      subroutine do1_net_eval
         use num_lib, only:null_solout
         
         implicit none
        
         integer :: iounit, i

         logical :: just_dxdt
         type (Net_Info), target :: netinfo_target
         type (Net_Info), pointer :: netinfo

         double precision :: eps, odescale
         logical :: use_pivoting, trace, burn_dbg

         integer :: burn_lwork, net_lwork
         real(dp), pointer :: burn_work_array(:), net_work_array(:)
         real(dp) :: avg_eps_nuc, eps_neu_total
         
         include "formats.dek"

         netinfo=>netinfo_target
         
         ierr = 0
         
         !Allocate the arrays declared in the beginning of the module:
         allocate (xa(species), y_mass(species), dabar_dx(species), dzbar_dx(species), &
         	d_eps_nuc_dx(species), dmc_dx(species), rate_factors(num_reactions), &
			eps_nuc_categories(num_categories), dxdt(species), &
			d_dxdt_dRho(species), d_dxdt_dT(species), d_dxdt_dx(species,species))
         
         ! set mass fractions -- must add to 1.0
         xa = 0d0
         write(*,*)
         if(use_solar) then
         	write(*,*) 'Solar abundances:'
         	do i=1,species
         		xa(i) = chem_Xsol(chem_isos% name(chem_id(i)))
         		write(*,*) chem_isos% name(chem_id(i)), chem_Xsol(chem_isos% name(chem_id(i)))
         	end do
         else
         	write(*,*) 'User-specified abundances:'
         	!Note, net_iso(ih1) = net_iso(get_nuclide_index('h1')), etc.:
         	!write(*,*) net_iso(ih1), net_iso(ihe4), net_iso(img24)
         	!write(*,*) net_iso(get_nuclide_index('h1')), &
         	!	net_iso(get_nuclide_index('he4')), net_iso(get_nuclide_index('mg24'))
         	do i=1,num_isos_for_Xinit
         		write(*,*) names_of_isos_for_Xinit(i), values_for_Xinit(i)
         		xa(net_iso(get_nuclide_index(names_of_isos_for_Xinit(i)))) = &
         			values_for_Xinit(i)
         	end do
         	!xa(net_iso(ih1)) = 7.5876644280605066d-01
         	!xa(net_iso(ihe4)) = 2.3952230737160904d-01
         	!xa(net_iso(img24)) = 1 - sum(xa(:))
         endif
         xa(:) = xa(:)/sum(xa)	!Normalize the abundances
         y_mass(:) = xa(:)/chem_isos% w(chem_id(1:species))
         
         !write(*,*) net_iso(ih1), net_iso(ihe4)
         !write(*,*) 'Solar abundances:'
         !do i=1,species
         !	write(*,*) chem_isos% name(chem_id(i)), chem_Xsol(chem_isos% name(chem_id(i)))
         !	write(*,*) chem_id(i), net_iso(i), net_iso(chem_id(i))
         !end do
         
         call composition_info(species, chem_id, xa, xh, xhe, zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
            
         logT = log10(burn_t)
         logRho = log10(burn_rho)
         
         eta = 0
         rate_factors(:) = 1
         weak_rate_factor = 1

         theta_e_for_graboske_et_al =  1 ! for nondegenerate
         screening_mode = extended_screening
         reuse_given_rates = .false.

         just_dxdt = .false.
         
         call net_get( &
            handle, just_dxdt, netinfo, species, num_reactions, &
            xa, burn_T, logT, burn_Rho, logRho, & 
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, & 
            std_reaction_Qs, std_reaction_neuQs, .false., .false., &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, & 
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, & 
            screening_mode, theta_e_for_graboske_et_al, &     
            eps_nuc_categories, eps_neu_total, & 
            lwork_net, work_net, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in net_get'
            return
         end if
         
         write(*,*) 'after net_get:', ierr
         write(*,*) xa(1), sum(xa), lwork_net
         write(*,*) 
         
         write(*,1) 'logT', logT
         write(*,1) 'logRho', logRho
         write(*,1) 'eps_nuc', eps_nuc
         write(*,*)
         
         iounit = 6
         call show_net_reactions(handle, iounit, ierr)
         
         num_times = 1
         allocate(ending_x(species), times(num_times), &
        	log10Ts_f1(4*num_times), log10Rhos_f1(4*num_times), &
            etas_f1(4*num_times), log10Ps_f1(4*num_times), stat=ierr)
         if (ierr /= 0) then
            write(*,*) 'allocate failed for Do_One_Test_Burn'
            stop 1
         end if
         
         log10Ts_f(1:4,1:num_times) => log10Ts_f1(1:4*num_times)
         log10Rhos_f(1:4,1:num_times) => log10Rhos_f1(1:4*num_times)
         etas_f(1:4,1:num_times) => etas_f1(1:4*num_times)
         log10Ps_f(1:4,1:num_times) => log10Ps_f1(1:4*num_times)
         
         
         t_start = 0d0 			!Starting time?
         t_end = burn_time			!Ending time (s)?
         !which_solver = 8 		!seulex - defined in num_def
         times = t_end
         log10Ts_f = logT
         log10Rhos_f = logRho
         etas_f = eta
         dxdt_source_term => null()
         !max_steps = 10000
         !rtol(:) = 1d-6
         !atol(:) = 1d-6
         rtol(:) = rtol_init
         atol(:) = atol_init
		 itol = 0
         y_min = -1e100
         y_max = 1e100
         h = 1d-4
         caller_id = 0
         iout = 1
         eps = rtol_init
         odescale = 1e-5
         trace = .false.
         burn_dbg = .false.
         use_pivoting = .true.

         burn_lwork = net_1_zone_burn_work_size(handle,ierr)
         if (ierr /= 0) stop
         net_lwork = net_work_size(handle,ierr)
         if (ierr /= 0) stop
         allocate(net_work_array(net_lwork), burn_work_array(burn_lwork))
         
          ! a 1-zone integrator for nets -- for given temperature and density as functions of time
      	 call net_1_zone_burn(&
              handle, species, num_reactions, t_start, t_end, xa, &
              ! for interpolating T, rho, eta wrt time
              num_times, times, log10Ts_f1, log10Rhos_f1, etas_f1, &
              ! other args for net_get
              dxdt_source_term, rate_factors, weak_rate_factor, std_reaction_Qs, std_reaction_neuQs, &
              screening_mode, theta_e_for_graboske_et_al, &
              ! args to control the solver -- see num/public/num_isolve.dek
              h, max_steps, eps, odescale, &
              .false., use_pivoting, trace, burn_dbg, burn_finish_substep, &
              burn_lwork, burn_work_array, &
              net_lwork, net_work_array, &
              ! results
              ending_x, avg_eps_nuc, eps_neu_total, &
              nfcn, njac, nstep, naccpt, nrejct, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in net_1_zone_burn'
            stop 1
         end if
             
        write(*,*) 
        write(*,*) 'burning complete'
        write(*,*)
        do i=1,species
        	write(*,*) chem_isos% name(chem_id(i)), ending_x(i)
        end do
        
        deallocate(ending_x, times, log10Ts_f1, log10Rhos_f1, etas_f1, log10Ps_f1, stat=ierr)  
        deallocate(net_work_array,burn_work_array)
      end subroutine do1_net_eval
      
      !subroutine lane_emden_init()
      !end subroutine
      
      !subroutine lane_emden_integrate(xi_start, xi_end, num_steps)
      !end subroutine
      
      
    !Subroutine to streamline the ZND integrations below
    subroutine znd_integrate(x_start, x_end, num_steps, output_profile, profile_io, find_lburn)
    	implicit none
    	
    	double precision, intent(in) :: x_start, x_end		!Integration bounds
    	integer, intent(in) :: num_steps	!Number of integration steps
    	logical, intent(in) :: output_profile !Whether to print output during the integration or not
    	integer, intent(in) :: profile_io !IO unit for profile output
    	logical, intent(in) :: find_lburn	!Whether or not to calculate l_burn or not
    	
    	double precision :: burn_time
    	double precision :: t_start, t_end
      !Uh-oh: num_vars isn't known at compile time so we have to allocate memory here. Only allocate if
      !we're actually jumping the pathological point - that's where these are used
    	!double precision :: orig_vars(num_vars), orig_derivs(num_vars)
      double precision, pointer, dimension(:) :: orig_vars, orig_derivs
    	integer :: j, k
      
      if(do_pathological.and.use_variable_linearization_ratio) then
         allocate(orig_vars(num_vars), orig_derivs(num_vars))
         orig_vars = 0d0
         orig_derivs = 0d0
      else
         orig_vars => null()
         orig_derivs => null()
      endif
    
      iwork_isolve = 0
		work_isolve = 0
		t_start = x_start				!Starting time?
		t_end = x_end		!Ending time (s)?
		burn_time = x_end - x_start
		
		vars(1:species) = xa		!Composition
		vars(species+1) = burn_rho	!Density (g/cm^3)
		vars(species+2) = burn_T	!Temperature (K)
		vars(species+3) = burn_u	!Velocity in shock frame (cm/s)
		
		if(do_blowout) then
			vars(species+4) = uy_init	!Initial radial blowout velocity (cm/s)
			if(use_uy_ux_ratio) then
				vars(species+4) = uy_ux_ratio*vars(species+3)
			endif
			!grav = standard_cgrav*m_c/r_wd**2
			!h_scale = P0/(rho0*grav)	!With gravity
			vars(species+5) = h_scale 	!Initial radial scale height (cm)
		else if(do_blowout_local) then
			vars(species+4) = cs0	!Initial radial blowout velocity (cm/s)
			vars(species+5) = h_scale 	!Initial radial scale height (cm)
			vars(species+6) = delta_y_init	!Initial length scale for computing dP/dy (cm)
		endif
		if(do_curvature) then
			if(use_rc_hs_scaling) r_curve = rc_hs_factor*h_scale
			if(.not.use_he_clavin) vars(num_vars) = r_curve
		endif
	
		sonic_loc = 0d0
		max_mach = 0d0
		gone_sonic = .false.
		
		do j=1,num_steps
   		!Evenly spaced timesteps in log
   		!t_start = h_search_xmax*10**(log10(h_search_xmax)*((k-1)*1d0/num_steps-1))
   		!t_end = h_search_xmax*10**(log10(h_search_xmax)*(k*1d0/num_steps-1))
   		!should be equivalent to:
   		t_start = x_start + burn_time**((j-1)*1d0/num_steps)
   		t_end = x_start + burn_time**(j*1d0/num_steps)
		
   		!write(*,*) x_start, x_end

   		call isolve( &
   			which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
   			h, max_step_size, max_steps, &
   			rtol, atol, itol, y_min, y_max, &
   			znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
   			null_mas, imas, mlmas, mumas, &
   			znd_solout, iout, &
   			lapack_decsol, null_decsols, null_decsolblk, &
                        lrd, rpar_decsol, lid, ipar_decsol, &
                        caller_id_blk, nvar_blk, nz_blk, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, &
                        null_fcn_blk_dble, null_jac_blk_dble, &
   			work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
   			lrpar, rpar, lipar, ipar, &
   			lout, idid)
   		if (ierr /= 0) then
   			write(*,*) 'failed in isolve() call'
            if(associated(orig_vars)) deallocate(orig_vars)
            if(associated(orig_derivs)) deallocate(orig_derivs)
   			stop 1
   		end if
   		!Check if we've slowed down too much (good stopping criterion for curvature)
   		if (rpar(2).lt.1d-2) then
   			write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
   			write(*,*)
   			exit
   		endif
		
   		!write(*,*) t_start, t_end, rpar(2)
		
   		if((rpar(2).gt.max_mach).and.(.not.gone_sonic)) then
   			sonic_loc = t_end
   			pathological_loc = t_end
   			max_mach = rpar(2)
			
   			!If we're traversing the pathological point then do that here if we reach the
   			!sonic limit used in pathological detonations:
   			if(do_pathological.and.(rpar(2).ge.sonic_limit_pathological)) then
   				write(*,*) 'pathological point hit - exiting integration...'
   				write(*,*)
   				exit
   			endif
   			if(rpar(2).ge.sonic_limit) then
   				gone_sonic = .true.
   				write(*,*) 'sonic point hit - exiting integration...'
   				write(*,*)
   				exit
   			endif
   		endif
		
   		if(output_profile) then
   			write(profile_io, '(i30,999ES30.20)') j, t_end, vars, rpar(1), &
   				rpar(8), rpar(10), rpar(2), rpar(9), mass_flux, mom_flux, energy_flux, &
   				rpar(11), rpar(12), rpar(13), rpar(14), rpar(15), rpar(16)
   				!dP_dx, dE_dx, dq_dx, detot_dx, dEtot_dx_eqn
   				!rpar(6), rpar(7), rpar(8:8+num_vars-1), rpar(8+num_vars)
   		endif
			
		end do !Integration loop
		
		!Linearize the solution and resume integration:
		if(do_pathological.and.(rpar(2).ge.sonic_limit_pathological)) then
		
			!Try a set of pathological_linearization_ratio values and see whether we
			!successfully jump over the pathological point:
			if(use_variable_linearization_ratio) then
				orig_vars = vars	!Save a copy
				orig_derivs = rpar(21:21+num_vars-1)
            !write(*,*) 'orig_derivs(1)', orig_derivs(1), rpar(21)
            !write(*,*) 'orig_derivs(num_vars)', orig_derivs(num_vars), rpar(21+num_vars), rpar(21+num_vars-1)
				do k=1,20
               if(associated(orig_vars)) write(*,*) 'associated(orig_vars)'
               if(associated(orig_derivs)) write(*,*) 'associated(orig_derivs)'
					pathological_linearization_ratio = 0.1*k
					linearization_length = pathological_loc*pathological_linearization_ratio
					vars(1:num_vars) = orig_vars(1:num_vars) + orig_derivs(1:num_vars)*linearization_length
					
					!Make sure we actually jump over the pathological point:
					if(rpar(2).le.1.0) then
						write(*,*) 'did not make it past the singularity - try again'
						write(*,'(a15,ES25.10)') 'ux/cs = ',rpar(2)
						continue
					endif
					
					!Now try resuming integration:
					do j=1,num_steps
						!Evenly spaced timesteps in log
						t_start = pathological_loc + linearization_length + &
							(burn_time - pathological_loc - linearization_length)**((j-1)*1d0/num_steps)
						t_end = pathological_loc + linearization_length + & 
							(burn_time - pathological_loc - linearization_length)**(j*1d0/num_steps)
						!write(*,*) 'Integrating to ',t_end
						!write(*,*) 'iwork_isolve(2) ',iwork_isolve(2)
						call isolve( &
							which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
							h, max_step_size, max_steps, &
							rtol, atol, itol, y_min, y_max, &
							znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
							null_mas, imas, mlmas, mumas, &
							znd_solout, iout, &
							lapack_decsol, null_decsols, null_decsolblk, &
                                                        lrd, rpar_decsol, lid, ipar_decsol, &
                                                        caller_id_blk, nvar_blk, nz_blk, &
                                                        lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, &
                                                        null_fcn_blk_dble, null_jac_blk_dble, &
							work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
							lrpar, rpar, lipar, ipar, &
							lout, idid)
						if (ierr /= 0) then
							write(*,*) 'failed in isolve() call'
                     if(associated(orig_vars)) deallocate(orig_vars)
                     if(associated(orig_derivs)) deallocate(orig_derivs)
							stop 1
						endif
						!Check if we've hit another singularity. If we're at an odd k,
						!then we should hit it from above (ie. supersonic -> subsonic),
						!while if we're at an even k then we should hit it from below
						!(subsonic -> supersonic)
!						if(mod(k,2).eq.1) then
!							if(rpar(2).le.(2d0-sonic_limit)) then
!								write(*,*) 'hit another singularity (from above) - try again...'
!								write(*,'(a15,ES25.10)') 'ux/cs = ',rpar(2)
!								exit
!							endif
!						else
!							if(rpar(2).ge.sonic_limit) then
!								write(*,*) 'hit another singularity (from below) - try again...'
!								write(*,'(a15,ES25.10)') 'ux/cs = ',rpar(2)
!								exit
!							endif
!						endif
						if(rpar(2).le.(2d0-sonic_limit)) then
							write(*,*) 'hit another singularity (from above) - try again...'
							write(*,'(a15,ES25.10)') 'ux/cs = ',rpar(2)
							exit
						endif
						
						!Check if we've slowed down too much (good stopping criterion for curvature)
						if (rpar(2).lt.1d-2) then
							write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
							write(*,*)
							exit
						endif
						if (rpar(2).gt.max_mach) then
							max_mach = rpar(2)
						endif
			
						if(output_profile) then
							write(profile_io, '(i30,999ES30.20)') j, t_end, vars, rpar(1), &
								rpar(8), rpar(10), rpar(2), rpar(9), mass_flux, mom_flux, energy_flux, &
								rpar(11), rpar(12), rpar(13), rpar(14), rpar(15), rpar(16)
								!dP_dx, dE_dx, dq_dx, detot_dx, dEtot_dx_eqn
								!rpar(6), rpar(7), rpar(8:8+num_vars-1), rpar(8+num_vars)
						endif
					end do !num_steps
					
					if(abs(rpar(2)-1d0).le.1d-1) then
						continue	!Hit the singularity - try again
					else
						exit	!Successful integration
					endif
					
				end do !Pathological linearization attempts
			else
				!Take the derivatives of all the quantities and linearize them over
				!a certain length scale (how to choose this length?) A first choice is to
				!use the time step we would have taken if the integration succeeded.
				!linearization_length = 4d6	!cm
				linearization_length = pathological_loc*pathological_linearization_ratio
				!write(*,*) 'pathological_loc:', pathological_loc
				!write(*,*) 'linearization_length:', linearization_length
				!write(*,*) 'vars pre:', vars
				vars(1:num_vars) = vars(1:num_vars) + rpar(21:21+num_vars-1)*linearization_length
				!write(*,*) 'vars post:', vars
			
				!Now try resuming integration:
				do j=1,num_steps
   				!Evenly spaced timesteps in log
   				t_start = pathological_loc + linearization_length + &
   					(burn_time - pathological_loc - linearization_length)**((j-1)*1d0/num_steps)
   				t_end = pathological_loc + linearization_length + & 
   					(burn_time - pathological_loc - linearization_length)**(j*1d0/num_steps)
   				!write(*,*) 'Integrating to ',t_end
   				!write(*,*) 'iwork_isolve(2) ',iwork_isolve(2)
   				call isolve( &
   					which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
   					h, max_step_size, max_steps, &
   					rtol, atol, itol, y_min, y_max, &
   					znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
   					null_mas, imas, mlmas, mumas, &
   					znd_solout, iout, &
   					lapack_decsol, null_decsols, null_decsolblk, &
                                        lrd, rpar_decsol, lid, ipar_decsol, &
                                        caller_id_blk, nvar_blk, nz_blk, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, &
                                        null_fcn_blk_dble, null_jac_blk_dble, &
   					work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
   					lrpar, rpar, lipar, ipar, &
   					lout, idid)
   				if (ierr /= 0) then
   					write(*,*) 'failed in isolve() call'
                  if(associated(orig_vars)) deallocate(orig_vars)
                  if(associated(orig_derivs)) deallocate(orig_derivs)
   					stop 1
   				endif
   				!Check if we've slowed down too much (good stopping criterion for curvature)
   				if (rpar(2).lt.1d-2) then
   					write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
   					write(*,*)
   					exit
   				endif
   				if (rpar(2).gt.max_mach) then
   					max_mach = rpar(2)
   				endif
			
   				if(output_profile) then
   					write(profile_io, '(i30,999ES30.20)') j, t_end, vars, rpar(1), &
   						rpar(8), rpar(10), rpar(2), rpar(9), mass_flux, mom_flux, energy_flux, &
   						rpar(11), rpar(12), rpar(13), rpar(14), rpar(15), rpar(16)
   						!dP_dx, dE_dx, dq_dx, detot_dx, dEtot_dx_eqn
   						!rpar(6), rpar(7), rpar(8:8+num_vars-1), rpar(8+num_vars)
   				endif
			
				end do !num_steps
			
			endif !use_variable_linearization_ratio
		
		endif !do_pathological
		
		q_tot = rpar(9)
		
		!If we're calculating l_burn, then integrate again until we hit q_tot*l_burn_factor
		if(find_lburn) then
		
			write(*,*) 'Initial integration found q = ',q_tot,' erg/g'
			!write(*,*) x_start, x_end
				
			!Integrate again, but stop when we hit q/2:
			!(we know the Neumann condition, so don't need to find that again)
			iwork_isolve = 0
			work_isolve = 0
			t_start = x_start 	!Starting time?
			t_end = x_end		!Ending time (s)?
			
			vars(1:species) = xa		!Composition
			vars(species+1) = burn_rho	!Density (g/cm^3)
			vars(species+2) = burn_T	!Temperature (K)
			vars(species+3) = burn_u	!Velocity in shock frame (cm/s)
			
			if(do_blowout) then
				vars(species+4) = uy_init	!Initial radial blowout velocity (cm/s)
				if(use_uy_ux_ratio) then
					vars(species+4) = uy_ux_ratio*vars(species+3)
				endif
				!grav = standard_cgrav*m_c/r_wd**2
				!h_scale = P0/(rho0*grav)	!With gravity
				vars(species+5) = h_scale 	!Initial radial scale height (cm)
			else if(do_blowout_local) then
				vars(species+4) = cs0	!Initial radial blowout velocity (cm/s)
				vars(species+5) = h_scale 	!Initial radial scale height (cm)
				vars(species+6) = delta_y_init	!Initial length scale for computing dP/dy (cm)
			endif
			if(do_curvature) then
				if(use_rc_hs_scaling) r_curve = rc_hs_factor*h_scale
				if(.not.use_he_clavin) vars(num_vars) = r_curve
			endif

			!Aha, we don't need to re-check for hitting a sonic point since we're
			!only integrating out halfway:
			!sonic_loc = 0d0
			!max_mach = 0d0
			!gone_sonic = .false.
			
			do j=1,num_steps
   			!Evenly spaced timesteps in log
   			t_start = x_start + burn_time**((j-1)*1d0/num_steps)
   			t_end = x_start + burn_time**(j*1d0/num_steps)
   			call isolve( &
   				which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
   				h, max_step_size, max_steps, &
   				rtol, atol, itol, y_min, y_max, &
   				znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
   				null_mas, imas, mlmas, mumas, &
   				znd_solout, iout, &
   				lapack_decsol, null_decsols, null_decsolblk, &
                                lrd, rpar_decsol, lid, ipar_decsol, &
                                caller_id_blk, nvar_blk, nz_blk, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, &
                                null_fcn_blk_dble, null_jac_blk_dble, &
   				work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
   				lrpar, rpar, lipar, ipar, &
   				lout, idid)
   			if (ierr /= 0) then
   				write(*,*) 'failed in isolve() call'
               if(associated(orig_vars)) deallocate(orig_vars)
               if(associated(orig_derivs)) deallocate(orig_derivs)
   				stop 1
   			end if
			
   			!write(*,*) t_start, t_end, rpar(2)
			
   			!Check if we've hit q_tot*l_burn_factor:
   			if (rpar(9).ge.q_tot*l_burn_factor) then
   				l_burn = t_end
   				write(*,'(a15,ES25.9)') 'l_burn = ',l_burn
   				write(*,*)
   				exit
   			endif
   			!Check if we've slowed down too much (good stopping criterion for curvature)
   			if (rpar(2).lt.1d-2) then
   				write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
   				write(*,*)
   				exit
   			endif
   			!If we're traversing the pathological point then do that here if we reach the
   			!sonic limit used in pathological detonations:
   			if(do_pathological.and.(rpar(2).ge.sonic_limit_pathological)) then
   				pathological_loc = t_end
   				write(*,*) 'pathological point hit - exiting integration...'
   				write(*,*)
   				exit
   			endif
			
			end do !Integration loop
			
			!If we haven't found the ZND length then go through the pathological point:
			if (do_pathological.and.(rpar(9).lt.q_tot*l_burn_factor)) then
				!Use the previously found linearization:
				if(use_variable_linearization_ratio) then
					vars(1:num_vars) = orig_vars(1:num_vars) + orig_derivs*linearization_length
				else
					linearization_length = pathological_loc*pathological_linearization_ratio
					!write(*,*) 'pathological_loc:', pathological_loc
					!write(*,*) 'linearization_length:', linearization_length
					!write(*,*) 'vars pre:', vars
					vars(1:num_vars) = vars(1:num_vars) + rpar(21:21+num_vars-1)*linearization_length
					!write(*,*) 'vars post:', vars
				endif
				
				!Now try resuming integration:
				do j=1,num_steps
				!Evenly spaced timesteps in log
				t_start = pathological_loc + linearization_length + &
					(burn_time - pathological_loc - linearization_length)**((j-1)*1d0/num_steps)
				t_end = pathological_loc + linearization_length + & 
					(burn_time - pathological_loc - linearization_length)**(j*1d0/num_steps)
				!write(*,*) 'Integrating to ',t_end
				!write(*,*) 'iwork_isolve(2) ',iwork_isolve(2)
				call isolve( &
					which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
					h, max_step_size, max_steps, &
					rtol, atol, itol, y_min, y_max, &
					znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
					null_mas, imas, mlmas, mumas, &
					znd_solout, iout, &
					lapack_decsol, null_decsols, null_decsolblk, &
                                        lrd, rpar_decsol, lid, ipar_decsol, &
                                        caller_id_blk, nvar_blk, nz_blk, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, &
                                        null_fcn_blk_dble, null_jac_blk_dble, &
					work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
					lrpar, rpar, lipar, ipar, &
					lout, idid)
				if (ierr /= 0) then
					write(*,*) 'failed in isolve() call'
               if(associated(orig_vars)) deallocate(orig_vars)
               if(associated(orig_derivs)) deallocate(orig_derivs)
					stop 1
				endif
				!Check if we've hit q_tot*l_burn_factor:
				if (rpar(9).ge.q_tot*l_burn_factor) then
					l_burn = t_end
					write(*,'(a15,ES25.9)') 'l_burn = ',l_burn
					write(*,*)
					exit
				endif
				!Check if we've slowed down too much (good stopping criterion for curvature)
				if (rpar(2).lt.1d-2) then
					write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
					write(*,*)
					exit
				endif
				end do !num_steps
            
			endif !do_pathological
		endif !find_lburn
      
      if(associated(orig_vars)) deallocate(orig_vars)
      if(associated(orig_derivs)) deallocate(orig_derivs)
   
   end subroutine znd_integrate
      
      
      !Subroutine to call my local set of burning subroutines rather than net's
      subroutine do_my_burn
      
      	use num_lib, only:null_solout
         
         implicit none
        
         double precision :: time_doing_net, logRho, logT, h_prev, gamma1, grav, rho, t
         double precision :: rho_b, t_b, p_b, h_p, r_wd, DM, m_env
         double precision :: Pgas, logPgas, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
         double precision :: det_mach, v_prev
         double precision :: xi1_low, xi1_high, xi1, dtheta_dxi, theta, xi, m_xi
         double precision :: xi_low, xi_high, m_tot, r_tot, dtheta_dxi1
         double precision :: t_det, t_sound
         double precision :: init_h_search_bracket_high, init_h_search_bracket_low
         double precision :: prev_h_search_bracket_high, prev_h_search_bracket_low
         double precision :: init_vdet_search_bracket_high, init_vdet_search_bracket_low
		 double precision :: prev_vdet_search_bracket_high, prev_vdet_search_bracket_low
		 double precision :: menv_target, mwd_target, rhoc_min, rhoc_max
		 double precision :: rho_c, pb_frac, r_c 
		 logical :: output_profile, find_lburn
         integer :: iounit, i, j, k, io, profile_io, report_io, sweep_io, data_in, ios
         character :: testchar
    	 !character (len=256) :: filestring
         
         include "formats.dek"
         
         ierr = 0
         io = 11 !output file unit for testing the solvers against each other
         profile_io = 12 !output file unit for outputting the x steps
         report_io = 13  !output file for printing the final states as a function of initial conditions
         sweep_io = 14	 !output file for storing sonic point info as a function of v_det
         data_in = 15	 !input file if we need to read in detonation parameters from a file
         
         !Allocate the arrays declared in the beginning of the module:
         allocate (vars(num_vars), xa(species), y_mass(species), dabar_dx(species), dzbar_dx(species), &
         	d_eps_nuc_dx(species), dmc_dx(species), rate_factors(num_reactions), &
			   eps_nuc_categories(num_categories), dxdt(species), &
			   d_dxdt_dRho(species), d_dxdt_dT(species), d_dxdt_dx(species,species))
		   if (ierr /= 0) then
            write(*, *) 'allocate ierr', ierr
            stop 1
         end if
         
         ! set mass fractions -- must add to 1.0
         xa = 0d0
         write(*,*)
         if(use_solar) then
         	write(*,*) 'Initial composition is solar:'
         	do i=1,species
         		xa(i) = chem_Xsol(chem_isos% name(chem_id(i)))
         		write(*,*) chem_isos% name(chem_id(i)), chem_Xsol(chem_isos% name(chem_id(i)))
         	end do
         else
         	write(*,*) 'User-specified initial composition:'
            write(*,*) num_isos_for_Xinit
         	do i=1,num_isos_for_Xinit
         		write(*,*) names_of_isos_for_Xinit(i), values_for_Xinit(i), &
                  net_iso(get_nuclide_index(names_of_isos_for_Xinit(i))), &
                  get_nuclide_index(names_of_isos_for_Xinit(i))
         		xa(net_iso(get_nuclide_index(names_of_isos_for_Xinit(i)))) = values_for_Xinit(i)
         	end do
         	write(*,*)
         endif
         xa(:) = xa(:)/sum(xa)	!Normalize the abundances
         !vars is the total array for all our ODE variables:
         vars = 0d0
         vars(1:species) = xa			!Composition
         vars(species+1) = 1d0
         vars(species+2) = 0d0
         
         write(*,*) 'here?'
         
         call composition_info(species, chem_id, xa, xh, xhe, Zm, abar, zbar, z2bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
         logRho = log10(rho0)
         logT = log10(t0)
         call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, &
         	xa, rho0, logRho, t0, logT, &
         	res, d_dlnRho_const_T, d_dlnT_const_Rho, &
         	d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
         p0 = exp(res(i_lnPgas)) + crad*t0**4/3d0
         cs0 = sqrt(res(i_gamma1)*p0/rho0)
         g0 = standard_cgrav*m_c/(r_wd**2)
         
		 q_init = avo*mev_to_ergs*sum(xa*chem_isos% binding_energy(chem_id(1:species))/&
		 	chem_isos% W(chem_id(1:species)))
		 	
		 write(*,*) 'Ambient conditions:'
		 write(*,'(a20,ES25.9)') 'rho0 (g/cm^3):',rho0
		 write(*,'(a20,ES25.9)') 'T0 (K):',t0
		 write(*,'(a20,ES25.9)') 'P0 (dyne/cm^2):',p0
		 write(*,'(a20,ES25.9)') 'free_e:',exp(res(i_lnfree_e))
		 write(*,'(a20,ES25.9)') 'cs0 (cm/s):',cs0
		 !write(*,'(a20,ES25.9)') 'g0 (cm/s^2):',g0
		 write(*,'(a20,ES25.9)') 'grav_init (cm/s^2):',grav_init
		 !write(*,'(a20,ES25.9)') 'H_s0 (cm):',P0/(rho0*standard_cgrav*m_c/r_wd**2)
		 write(*,*)
		 write(*,'(a20,i25)') 'num_vars:',num_vars
		 write(*,*)
		 
		 !write(*,*) 'Testing some sums used in the ZND solver:'
		 !write(*,*) 'species = ',species
		 !write(*,*) chem_isos% W(chem_id(1:species))
		 !write(*,*) chem_isos% binding_energy(chem_id(1:species))
         
         !Here's where we'd call the solve_cj subroutine to find the initial conditions
         !for the ZND integrator:
         if(do_cj) then
         	call solve_cj
         	call solve_neumann
         	write(*,'(a20,ES25.9)') 'rho_neumann:',rho_neumann
         	write(*,'(a20,ES25.9)') 't_neumann:',t_neumann
         	write(*,'(a20,ES25.9)') 'p_neumann:',p_neumann
         	write(*,'(a20,ES25.9)') 'e_neumann:',e_neumann
         	write(*,'(a20,ES25.9)') 'u_neumann:',(rho0/rho_neumann)*v_det
         	write(*,'(a20,ES25.9)') 'cs_neumann:',cs_neumann
         	write(*,*)
         	burn_rho = rho_neumann
         	burn_t = t_neumann
         	burn_u = (rho0/rho_neumann)*v_det
         	burn_P = p_neumann
         	burn_e = e_neumann
         	!burn_u = 1.258*(rho0/rho_neumann)*v_det
         	!burn_u = 1.5*(rho0/rho_neumann)*v_det
         endif
         
         !If we just specify the detonation velocity and initial conditions, then we
         !just need to find the neumann point:
         if(do_neumann) then
         	call solve_neumann
         	write(*,'(a20,ES25.9)') 'rho_neumann:',rho_neumann
         	write(*,'(a20,ES25.9)') 't_neumann:',t_neumann
         	write(*,'(a20,ES25.9)') 'p_neumann:',p_neumann
         	write(*,'(a20,ES25.9)') 'e_neumann:',e_neumann
         	write(*,'(a20,ES25.9)') 'u_neumann:',(rho0/rho_neumann)*v_det
         	write(*,'(a20,ES25.9)') 'cs_neumann:',cs_neumann
         	write(*,*)
         	burn_rho = rho_neumann
         	burn_t = t_neumann
         	burn_u = (rho0/rho_neumann)*v_det
         	burn_P = p_neumann
         	burn_e = e_neumann
         	!burn_u = 1.258*(rho0/rho_neumann)*v_det
         	!burn_u = 1.5*(rho0/rho_neumann)*v_det
         endif
            
         logT = log10(burn_t)
         logRho = log10(burn_rho)
         
         eta = 0
         rate_factors(:) = 1
         weak_rate_factor = 1
         theta_e_for_graboske_et_al =  1 ! for nondegenerate
         screening_mode = extended_screening
         reuse_given_rates = .false.
         dxdt_source_term => null()
		 itol = 0
         iout = 1
         
         t_start = 0d0 			!Starting time?
         t_end = burn_time		!Ending time (s)?
         h = 1d-10

         max_step_size = 0 
         !max_steps = 10000
         
         !itol = 1			!O: rtol & atol are scalars, 1: they're vectors
         !rtol(1:species) = 1d-4
         !atol(1:species) = 1d-4
         !rtol(species+1:num_vars) = 1d-4
         !atol(species+1:num_vars) = 1d4
         !write(*,*) 'rtol:', rtol

         itol = 1			!O: rtol & atol are scalars, 1: they're vectors
         rtol(:) = rtol_init
         atol(:) = atol_init
         y_min = -1e100
         y_max = 1e100
         
         !ijac = 1			!0: finite differences, 1: analytic (in inlist now)
         nzmax = 0
         isparse = 0
         mljac = num_vars ! square matrix
         mujac = num_vars

         imas = 0
         mlmas = 0
         mumas = 0        

         caller_id_blk = 0
         nvar_blk = 0
         nz_blk = 0
         nullify( lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk )
         
         iout = 1
         
         lid = 0
         lrd = 0

		 allocate(rpar(lrpar), ipar(lipar))
		 if (ierr /= 0) then
            write(*, *) 'allocate ierr', ierr
            stop 1
         end if
         ipar = 0
         rpar = 0         

         lout = 6
         
         write(*,*)
         !write(*,*) species, lrd, lid

         call lapack_work_sizes(num_vars, lrd, lid)

         call isolve_work_sizes(num_vars, nzmax, imas, mljac, mujac, &
         	mlmas, mumas, liwork_isolve, lwork_isolve)
         
         allocate(iwork_isolve(liwork_isolve), work_isolve(lwork_isolve), &
         	ipar_decsol(lid), rpar_decsol(lrd), stat=ierr)
         if (ierr /= 0) then
            write(*, *) 'allocate ierr', ierr
            stop 1
         end if
      
         iwork_isolve = 0
         work_isolve = 0d0
         write(*,'(a20,i25)') 'liwork_isolve',liwork_isolve
         write(*,'(a20,i25)') 'lwork_isolve',lwork_isolve
         write(*,'(a20,ES25.9)') 't_start',t_start
         write(*,'(a20,ES25.9)') 't_end',t_end
         write(*,'(a20,i25)') 'max_steps',max_steps
         
         if(do_constp) then
         	call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, &
         		xa, burn_Rho, logRho, burn_T, logT, &
         		res, d_dlnRho_const_T, d_dlnT_const_Rho, &
         		d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
         	burn_constPgas = exp(res(i_lnPgas))! + crad*T**4/3d0
         endif
         
         !open(unit=io, file='znd_solvers.data')
         !write(io,'(999a25)') 'burn length (cm)', 'rtol', 'atol', 'solver', &
         !	chem_isos% name(chem_id(1:species)), 'rho (g/cm^3)', 'T (K)', 'u (cm/s)', 'u/c_s'
         
         !If we're doing a single integration, then set up the headers appropriately:
         if(do_constp.or.do_znd.or.do_vdet_search.or.do_h_search) then
			 open(unit=profile_io, file=output_profile_name)
			 write(profile_io,'(999a25)') 'rho0 (g/cm^3)', 'T0 (K)', 'P0 (dyne/cm^2)', &
         		'cs0 (cm/s)', chem_isos% name(chem_id(1:species)), &
         		'H_scale (cm)', 'v_det (cm/s)', 'x_end (cm)'
         	 write(profile_io,'(999ES25.9)') rho0, t0, p0, cs0, xa(1:species), &
         	 	h_scale, v_det, burn_time
         	 write(profile_io,*)
			 
			 if(do_blowout) then
				write(profile_io,'(999a30)',advance='no') 'Step', 'Position (cm)', &
				chem_isos% name(chem_id(1:species)), 'rho (g/cm^3)', 'T (K)', 'u (cm/s)',&
				'u_y (cm/s)', 'H_scale (cm)'
				
				if(do_curvature.and.(.not.use_he_clavin)) then
					write(profile_io,'(999a30)',advance='no') 'R_curve (cm)'
				endif
				
				write(profile_io,'(999a30)') &
				'P (dyne/cm^2)', 'Etot (erg/g)', 'gamma1', 'u/c_s', 'q (erg/g)', &
				!'p_eqn (dyne/cm^2)', 'p_eos (dyne/cm^2)', 'u/c_s', &
				!'eqn1', 'eqn2', 'eqn3', 'Etot_EOS (erg/g)', &
				 'mass flux', 'mom. flux', 'energy flux'
				!'dP_dx (EOS)', 'dE_dx (EOS)', 'dq_dx (EOS)', 'dEtot_dx (eqn)', 'dEtot_dx (eqn)'
				!'dP_dx (EOS)', 'dP_dx (eqn)', &
				!'f(1)', 'f(2)', 'f(3)', 'f(4)', 'f(5)', 'f(6)'
			 else if(do_blowout_local) then
				write(profile_io,'(999a30)',advance='no') 'Step', 'Position (cm)', &
				chem_isos% name(chem_id(1:species)), 'rho (g/cm^3)', 'T (K)', 'u (cm/s)',&
				'u_y (cm/s)', 'H_scale (cm)', 'delta_y (cm)'
				
				if(do_curvature.and.(.not.use_he_clavin)) then
					write(profile_io,'(999a30)',advance='no') 'R_curve (cm)'
				endif
				
				write(profile_io,'(999a30)') &
				'P (dyne/cm^2)', 'Etot (erg/g)', 'gamma1', 'u/c_s', 'q (erg/g)', &
				 'mass flux', 'mom. flux', 'energy flux'
			 else
				write(profile_io,'(999a30)') 'Step', 'Position_cm', &
				chem_isos% name(chem_id(1:species)), 'rho_g_p_cm_3', 'T_K', 'u_cm_p_s',&
				'P_dyne_p_cm_2', 'Etot_erg_p_g', 'gamma1', 'u_over_c_s', 'q_erg_p_g', &
				 'mass_flux', 'mom._flux', 'energy_flux', &
				 'f1_blowout', 'f2_blowout', 'f3_blowout', &
				 'f1_curvature', 'f2_curvature', 'f3_curvature'
			 endif
		 endif
         	
         !open(unit=report_io, file='znd_report.data', access='APPEND')
         !write(report_io,'(99a25)') 'rho0 (g/cm^3)', 'T0 (K)', 'rho_cj (g/cm^3)', &
         !	'T_cj (K)', 'v_cj (cm/s)', 'rho_neumann (g/cm^3)', 'T_neumann (K)', &
         !	'u_burn (cm/s)', 'l_burn (cm)', &
         !	chem_isos% name(chem_id(1:species)), &
         !	'rho_final (g/cm^3)', 'T_final (K)', 'u_final (cm/s)', '(u/c_s)_final'
         
         
        if(do_cjstate_only) then 
        	write(*,*) 'Sweeping through rho - finding naive CJ state only' 
        	        
        	!t0 = rho_sweep_t0
			logT = log10(t0)
			!grav_init = rho_sweep_grav
       		
       		open(unit=sweep_io, file=output_sweep_name)
			write(sweep_io,'(999a25)') 'T0 (K)', &
				'grav0 (cm/s^2)', 'cs0 (cm/s)', chem_isos% name(chem_id(1:species))
         	write(sweep_io,'(999ES25.9)') t0, grav_init, &
         		cs0, xa(1:species)
         	write(sweep_io,*)
			write(sweep_io,'(999a25)') 'rho0 (g/cm^3)', 'T0 (K)', 'P0 (dyne/cm^2)', &
				'q_tot (erg/g)', 'v_CJ (cm/s)', 'v_CJ/cs0', &
				'rho_neumann (g/cm^3)', 'T_neumann (K)', 'P_neumann (dyne/cm^2)', &
				'e_neumann (erg/g)', 'burn_u (cm/s)', &
				'rho_CJ (g/cm^3)', 'T_CJ (K)', 'P_CJ (dyne/cm^2)', 'e_CJ (erg/g)'
				
			call system_clock(sys_time_begin,clock_rate)
        
        	!Sweep through densities:
        	do i=1,rho_sweep_num_steps
        		rho0 = 10**(rho_sweep_log_min + &
        			(rho_sweep_log_max-rho_sweep_log_min)*(i-1d0)/(rho_sweep_num_steps-1d0))
				logRho = log10(rho0)
         		call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, &
         			xa, rho0, logRho, t0, logT, &
         			res, d_dlnRho_const_T, d_dlnT_const_Rho, &
         			d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
         		p0 = exp(res(i_lnPgas)) + crad*t0**4/3d0
         		cs0 = sqrt(res(i_gamma1)*p0/rho0)
				h_scale = p0/(rho0*grav_init)
				
				
				call solve_cj
				q_tot = q
				
				!Find the Neumann point:
				call solve_neumann
				burn_rho = rho_neumann
				burn_t = t_neumann
				burn_u = (rho0/rho_neumann)*v_det
				burn_P = p_neumann
				burn_e = e_neumann
				
				write(sweep_io,'(999ES25.9)') rho0, t0, p0, q_tot, v_cj, v_cj/cs0, &
					rho_neumann, t_neumann, p_neumann, e_neumann, burn_u, rho_cj, &
					t_cj, p_cj, e_cj
			end do
        else if(do_constp) then
			do i=1,8
         		t_start = 0d0 			!Starting time?
         		t_end = burn_time		!Ending time (s)?
         		vars(1:species) = xa
         		vars(species+1) = burn_T
         		vars(species+2) = 0d0
         
         		which_solver = i		!Defined in num_def.f
         		call system_clock(sys_time_begin,clock_rate)
      	 		call isolve( &
      	 			which_solver, num_vars, constp_derivs, t_start, vars, t_end, &
            		h, max_step_size, max_steps, &
            		rtol, atol, itol, y_min, y_max, &
            		constp_jacob, ijac, burn_sjac, nzmax, isparse, mljac, mujac, &
            		null_mas, imas, mlmas, mumas, &
            		null_solout, iout, &
            		lapack_decsol, null_decsols, null_decsolblk, &
                        lrd, rpar_decsol, lid, ipar_decsol, &
                        caller_id_blk, nvar_blk, nz_blk, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, &
                        null_fcn_blk_dble, null_jac_blk_dble, &
            		work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
            		lrpar, rpar, lipar, ipar, &
            		lout, idid)
         		if (ierr /= 0) then
            		write(*,*) 'failed in isolve() call'
            		stop 1
         		end if
         		call system_clock(sys_time_end,clock_rate)
             
        		write(*,*) 
        		write(*,*) 'isolve complete with local subroutines', idid
        		write(*,*) 'Computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate
        		write(*,*)
        		do j=1,species
        			write(*,*) chem_isos% name(chem_id(j)), vars(j)
        		end do
        		write(*,*) 'vars(species+1)', vars(species+1)
        		write(*,*) 'vars(species+2)', vars(species+2)
        
      			write(io,'(3ES25.9,i25,99ES25.9)') t_end, rtol(1), atol(1), &
      				which_solver, vars(1:species)
        
        	end do
        else if(do_znd) then
			t_start = 0d0 			!Starting location
         	t_end = burn_time		!Ending location
         	output_profile = .true.
         	find_lburn = .true.
			call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
			call system_clock(sys_time_end,clock_rate)			
		 
			write(*,*) 
			write(*,*) 'isolve complete with local subroutines', idid
			write(*,'(a20,ES15.4,a2)') 'Computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate, 's'
			write(*,*)
			do j=1,species
				write(*,'(a15,999ES25.9)') chem_isos% name(chem_id(j)), vars(j)
			end do
			write(*,'(a15,ES25.9)') 'rho (g/cm^3)', vars(species+1)
			write(*,'(a15,ES25.9)') 'T (K)', vars(species+2)
			write(*,'(a15,ES25.9)') 'u (cm/s)', vars(species+3)
			write(*,'(a15,ES25.9)') 'u/c_s', rpar(2)
			write(*,'(a15,ES25.9)') 'max(u/cs)', max_mach
			write(*,'(a15,ES25.9)') 'maxloc(u/cs)', sonic_loc
			write(*,'(a15,ES25.9)') 'q_final', rpar(9)
			write(*,'(a15,ES25.9)') 'gamma1', rpar(10)
			write(*,'(a15,ES25.9)') 'l_burn', l_burn
	
			!write(io,'(3ES25.9,i25,999ES25.9)') t_end, rtol(1), atol(1), &
			!	which_solver, vars, rpar(1)
        	
        	!write(report_io,'(99ES25.9)') rho0, t0, rho_cj, t_cj, v_cj, rho_neumann, &
        	!	t_neumann, burn_u, burn_time, vars, rpar(1)
        	
        !This section will automate finding a H_scale value such that the blowout causes
        !the detonation velocity specified to be the 'generalized CJ' velocity. This
        !means that the sonic point is reached as the burning terminates.
        else if(do_h_search) then
        	write(*,*) 'Searching for H_scale...'
        	write(*,*) 'v_det = ',h_search_vdet
        	h_prev = 0d0
        	v_det = h_search_vdet
        	call system_clock(sys_time_begin,clock_rate)
        	
        	do i=1,h_search_numsteps
        		h_scale = (h_search_bracket_high + h_search_bracket_low)/2d0
         		
         		!Find the Neumann point:
         		call solve_neumann
				burn_rho = rho_neumann
				burn_t = t_neumann
				burn_u = (rho0/rho_neumann)*v_det
				burn_P = p_neumann
				burn_e = e_neumann
         		
         		t_start = 0d0 				!Starting location
         		t_end = h_search_xmax		!Ending location
         		output_profile = .false.
         		find_lburn = .false.
         		do_pathological = .false.
				call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
         		
        		write(*,'(a15,ES25.9)') 'u_x (cm/s)', vars(species+3)
        		write(*,'(a15,ES25.9)') 'u_x/c_s', rpar(2)
        		write(*,'(a15,ES25.9)') 'max(u/cs)', max_mach
        		write(*,'(a15,ES25.9)') 'maxloc(u/cs)', sonic_loc
        		write(*,'(a15,ES25.9)') 'q_final', rpar(9)
        		write(*,'(a15,ES25.9)') 'h_scale', h_scale
        		write(*,'(a15,ES25.9)') 'delta_h', abs(h_scale - h_prev)/h_scale
        		
        		!Check if we've found a solution:
        		if(abs(h_scale - h_prev)/h_scale.le.h_search_tol) then
        			write(*,*) 'Solution found:'
        			if(do_pathological) then !take the case that will hit the sonic point
        				h_scale = h_search_bracket_high
        			else !take the case that won't hit the sonic point
        				h_scale = h_search_bracket_low
        			endif
        			write(*,'(a15,ES25.9)') 'h_scale (cm)', h_scale
        			exit
        		endif
        		
        		!Update the search paramaters:
        		!For finding a solution where the burning stops at the sonic point
        		!If we hit a sonic point then decrease the scale height, otherwise increase it
        		if(max_mach.gt.sonic_limit) then
        			h_search_bracket_high = h_scale
        		else
        			h_search_bracket_low = h_scale
        		endif
        		
        		write(*,*) 'Brackets now:'
        		write(*,'(a15,ES25.9)') 'high (cm)', h_search_bracket_high
        		write(*,'(a15,ES25.9)') 'low (cm)', h_search_bracket_low
        		
        		h_prev = h_scale
        		
        	end do !sweep through v_det 
        	
        	!Integrate again to create a profile at our given h_scale
			call solve_neumann
			burn_rho = rho_neumann
			burn_t = t_neumann
			burn_u = (rho0/rho_neumann)*v_det
			burn_P = p_neumann
			burn_e = e_neumann
			
			t_start = 0d0 				!Starting location
			t_end = h_search_xmax		!Ending location
			output_profile = .true.
			find_lburn = .true.
			do_pathological = do_pathological_init
			call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
				
			write(*,'(a15,ES25.9)') 'u_x (cm/s)', vars(species+3)
			write(*,'(a15,ES25.9)') 'u_x/c_s', rpar(2)
			write(*,'(a15,ES25.9)') 'max(u/cs)', max_mach
			write(*,'(a15,ES25.9)') 'maxloc(u/cs)', sonic_loc
			write(*,'(a15,ES25.9)') 'q_final', q_tot
			write(*,'(a15,ES25.9)') 'l_burn', l_burn
			write(*,'(a15,ES25.9)') 'v_det', v_det
			write(*,'(a15,ES25.9)') 'h_scale', h_scale
			write(*,'(a15,ES25.9)') 'delta_h', abs(h_scale - h_prev)/h_scale
        	
        	call system_clock(sys_time_end,clock_rate)
        	write(*,'(a20,ES15.4,a2)') 'Computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate, 's'
        
        !This section will automate finding a v_det value such that the blowout causes
        !the detonation velocity specified to be the 'generalized CJ' velocity. This
        !means that the sonic point is reached as the burning terminates. 
        !*** Be careful about choosing your vdet_search_bracket_low since this search will
        !assume you're on the upper (stable) branch of v_det vs. h_scale. The lower
        !bracket must be a velocity such that the integration will hit a sonic point.
        !If your bracket is too low, then this may accidentally find the lower 
        !(unstable) velocity. If you're unsure where this lies for your initial 
        !conditions, you should run a full search using do_lznd_vdet *** 
        else if(do_vdet_search) then
        	write(*,*) 'Searching for v_det...'
        	write(*,*) 'H_scale = ',vdet_search_h_scale
        	v_prev = 0d0
        	h_scale = vdet_search_h_scale
        	call system_clock(sys_time_begin,clock_rate)
        	
        	do i=1,vdet_search_numsteps
        		v_det = (vdet_search_bracket_high + vdet_search_bracket_low)/2d0

         		!Find the Neumann point:
         		call solve_neumann
				burn_rho = rho_neumann
				burn_t = t_neumann
				burn_u = (rho0/rho_neumann)*v_det
				burn_P = p_neumann
				burn_e = e_neumann
				
				t_start = 0d0 				!Starting location
				t_end = vdet_search_xmax	!Ending location
				output_profile = .false.
				find_lburn = .false.
				do_pathological = .false.
				write(*,*) 'integrating with vdet =', v_det
				call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)

        		write(*,'(a15,ES25.9)') 'u_x (cm/s)', vars(species+3)
        		write(*,'(a15,ES25.9)') 'u_x/c_s', rpar(2)
        		write(*,'(a15,ES25.9)') 'max(u/cs)', max_mach
        		write(*,'(a15,ES25.9)') 'maxloc(u/cs)', sonic_loc
        		write(*,'(a15,ES25.9)') 'q_final', rpar(9)
        		write(*,'(a15,ES25.9)') 'v_det', v_det
        		write(*,'(a15,ES25.9)') 'delta_vdet', abs(v_det - v_prev)/v_det
        		
        		!Check if we've found a solution:
        		if(abs(v_det - v_prev)/v_det.le.vdet_search_tol) then
        			write(*,*) 'Generalized CJ solution found:'
        			if(do_pathological_init) then !take the case that will hit the sonic point
        				v_det = vdet_search_bracket_low
        			else !take the case that won't hit the sonic point
        				v_det = vdet_search_bracket_high
        			endif
        			write(*,'(a15,ES25.9)') 'v_det (cm/s)', v_det
        			exit
        		endif
        		
        		!Update the search paramaters:
        		!For finding a solution where the burning stops at the sonic point
        		!If we hit a sonic point then increase v_det, otherwise decrease it
        		if(max_mach.ge.sonic_limit) then
        			vdet_search_bracket_low = v_det
        		else
        			vdet_search_bracket_high = v_det
        		endif
        		
        		write(*,*) 'Brackets now:'
        		write(*,'(a15,ES25.9)') 'high (cm/s)', vdet_search_bracket_high
        		write(*,'(a15,ES25.9)') 'low (cm/s)', vdet_search_bracket_low
        		
        		v_prev = v_det
        		
        	end do !sweep through v_det 
        	
        	!Do an actual integration at the found v_det to save time:
			call solve_neumann
			burn_rho = rho_neumann
			burn_t = t_neumann
			burn_u = (rho0/rho_neumann)*v_det
			burn_P = p_neumann
			burn_e = e_neumann
			
			t_start = 0d0 				!Starting location
			t_end = vdet_search_xmax	!Ending location
			output_profile = .true.
			find_lburn = .true.
			do_pathological = do_pathological_init
			call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
			
			!Check the free electron ratio at the end to see if there were electron captures
			call composition_info(species, chem_id, vars(1:species), xh, xhe, Zm, abar, zbar, z2bar, &
            	ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
			Rho = vars(species+1)
			logRho = log10(Rho)
			T = vars(species+2)
			logT = log10(T)
			call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, &
         		vars(1:species), Rho, logRho, T, logT, &
         		res, d_dlnRho_const_T, d_dlnT_const_Rho, &
         		d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
			
			write(*,'(a15,ES25.9)') 'Y_e', ye
			write(*,'(a15,ES25.9)') 'free_e', exp(res(i_lnfree_e))
			write(*,'(a15,ES25.9)') 'u_x (cm/s)', vars(species+3)
			write(*,'(a15,ES25.9)') 'u_x/c_s', rpar(2)
			write(*,'(a15,ES25.9)') 'max(u/cs)', max_mach
			write(*,'(a15,ES25.9)') 'maxloc(u/cs)', sonic_loc
			write(*,'(a15,ES25.9)') 'q_final', q_tot
			write(*,'(a15,ES25.9)') 'v_det', v_det
			
        	call system_clock(sys_time_end,clock_rate)
        	write(*,'(a20,ES15.4,a2)') 'Computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate, 's'

        !This will automate a series of ZND calculations. It assumes initial conditions
        !(rho0, T0, composition, and h_scale if necessary) are fixed, and sweeps through
        !detonation velocities, calculating maxloc(ux/cs) and max(ux/cs) for each
        !value of v_det
        else if(do_sweep) then
        	open(unit=sweep_io, file=output_sweep_name)
        	write(sweep_io,'(99a25)') 'rho0 (g/cm^3)', 'T0 (K)', 'P0 (dyne/cm^2)', &
         		'cs0 (cm/s)', chem_isos% name(chem_id(1:species)), 'x_end (cm)'
         	write(sweep_io,'(99ES25.9)') rho0, t0, p0, cs0, xa(1:species), burn_time
         	write(sweep_io,*)
        	write(sweep_io,'(99a25)') 'v_det (cm/s)', 'v_det/cs0', &
        		'rho_neumann (g/cm^3)', 'T_neumann (K)', 'P_neumann (dyne/cm^2)', &
        		'e_neumann (erg/g)', 'burn_u (cm/s)', &
         		'maxloc(ux/cs) (cm)', 'max(ux/cs)'
        	call system_clock(sys_time_begin,clock_rate)
			do i=1,sweep_znd_iters
         		v_det = (sweep_min_mach + sweep_max_mach*(i-1)/(sweep_znd_iters-1))*cs0
         		
         		!Find the Neumann point:
         		call solve_neumann
				burn_rho = rho_neumann
				burn_t = t_neumann
				burn_u = (rho0/rho_neumann)*v_det
				burn_P = p_neumann
				burn_e = e_neumann
				
				t_start = 0d0 				!Starting location
				t_end = burn_time			!Ending location
				output_profile = .false.
				find_lburn = .true.
				call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
        
        		write(*,'(99ES25.9)') v_det, v_det/cs0, sonic_loc, max_mach
      			write(sweep_io,'(99ES25.9)') v_det, v_det/cs0, &
      				rho_neumann, t_neumann, p_neumann, e_neumann, burn_u, &
      				sonic_loc, max_mach
        	end do !sweep through v_det 
        	call system_clock(sys_time_end,clock_rate)
        	write(*,'(a20,ES15.4,a2)') 'Computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate, 's'
        	close(sweep_io)
        	
        !This section will sweep through detonation velocities, for each one finding the
        !minimum H_scale where a sonic point is not hit and integrate twice - once to
        !find the total energy release, abundances, etc. and again to find the ZND length
        else if(do_lznd_vdet) then
       		write(*,*) 'Finding l_burn as a function of v_det...'
       		
       		open(unit=sweep_io, file=output_sweep_name)
			write(sweep_io,'(999a25)') 'rho0 (g/cm^3)', 'T0 (K)', 'P0 (dyne/cm^2)', &
         		'cs0 (cm/s)', chem_isos% name(chem_id(1:species)), 'x_end (cm)'
         	write(sweep_io,'(999ES25.9)') rho0, t0, p0, cs0, xa(1:species), burn_time
         	write(sweep_io,*)
			write(sweep_io,'(999a25)') 'v_det (cm/s)', 'v_det/cs0', 'h_scale (cm)', &
				'rho_neumann (g/cm^3)', 'T_neumann (K)', 'P_neumann (dyne/cm^2)', &
				'e_neumann (erg/g)', 'burn_u (cm/s)', 'maxloc(ux/cs) (cm)', &
				'max(ux/cs)', 'q_tot (erg/g)', 'l_burn (cm)', &
				chem_isos% name(chem_id(1:species))
        
        	call system_clock(sys_time_begin,clock_rate)
        	init_h_search_bracket_high = h_search_bracket_high
			init_h_search_bracket_low = h_search_bracket_low
			prev_h_search_bracket_high = h_search_bracket_high
			prev_h_search_bracket_low = h_search_bracket_low
			
        	do i=1,sweep_znd_iters
				write(*,*) 'Searching for H_scale...'
				v_det = (sweep_min_mach + sweep_max_mach*(i-1)/(sweep_znd_iters-1))*cs0
				write(*,*) 'v_det = ',v_det
				
				!Find the Neumann point (just once per v_det):
				call solve_neumann
				burn_rho = rho_neumann
				burn_t = t_neumann
				burn_u = (rho0/rho_neumann)*v_det
				burn_P = p_neumann
				burn_e = e_neumann
				
				!Reset all the values used in finding h_scale:
				h_prev = 0d0
				!Use what we know about the h_scale for the previous v_det to save some
				!effort (we know the h_scale for a neighboring v_det will be close)
				h_search_bracket_high = min(prev_h_search_bracket_high*1d1, init_h_search_bracket_high)
				h_search_bracket_low = max(prev_h_search_bracket_low/1d1, init_h_search_bracket_low)
			
				do j=1,h_search_numsteps
					h_scale = (h_search_bracket_high + h_search_bracket_low)/2d0
					
					t_start = 0d0 				!Starting location
					t_end = h_search_xmax		!Ending location
					output_profile = .false.
					find_lburn = .false.
					do_pathological = .false.
					call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
					
					write(*,'(a15,ES25.9)') 'u_x (cm/s)', vars(species+3)
					write(*,'(a15,ES25.9)') 'u_x/c_s', rpar(2)
					write(*,'(a15,ES25.9)') 'max(u/cs)', max_mach
					write(*,'(a15,ES25.9)') 'maxloc(u/cs)', sonic_loc
					write(*,'(a15,ES25.9)') 'q_final', q_tot
					write(*,'(a15,ES25.9)') 'h_scale', h_scale
					write(*,'(a15,ES25.9)') 'delta_h', abs(h_scale - h_prev)/h_scale
			
					!write(*,'(99ES25.9)') v_det, v_det/cs0, sonic_loc, max_mach
					
					!Check if we've found a solution:
					if(abs(h_scale - h_prev)/h_scale.le.h_search_tol) then
						write(*,*) 'Solution found:'
						if(do_pathological_init) then !take the case that will hit the sonic point
							h_scale = h_search_bracket_high
						else !take the case that won't hit the sonic point
							h_scale = h_search_bracket_low
						endif
						write(*,'(a15,ES25.9)') 'h_scale (cm)', h_scale
						exit
					endif
					
					!Update the search paramaters:
					!For finding a solution where the burning stops at the sonic point
					!If we hit a sonic point then decrease the scale height, otherwise increase it
					if(gone_sonic) then
						h_search_bracket_high = h_scale
					else
						h_search_bracket_low = h_scale
					endif
					
					write(*,*) 'Brackets now:'
					write(*,'(a15,ES25.9)') 'high (cm)', h_search_bracket_high
					write(*,'(a15,ES25.9)') 'low (cm)', h_search_bracket_low
					
					h_prev = h_scale
					
				end do !search for h_scale
				
				!h_scale = h_search_bracket_low !To ensure we don't hit the sonic point, gives frozen subsonic solution
				!h_scale = h_search_bracket_high !To ensure we can go through the pathological point, gives frozen supersonic solution
				prev_h_search_bracket_high = h_search_bracket_high
				prev_h_search_bracket_low = h_search_bracket_low
				write(*,*)
				write(*,*) 'H_scale found: ',h_scale,' cm'
				
				!Finally, check if we've hit the upper limit of H_scale during our sweep.
      			!If so, then we can stop our sweep, as we'll only need larger H_scale
      			!values as we increase v_det from here
      			if(h_search_bracket_high.eq.init_h_search_bracket_high) then
      				write(*,*) 'Maximum scale height of ',init_h_search_bracket_high,' cm hit - stopping sweep'
      				exit
      			endif
				
				!Now that we've found a valid h_scale, we need to find the ZND length
				write(*,*) 'Calculating ZND length...'
				
				t_start = 0d0 				!Starting location
				t_end = h_search_xmax		!Ending location
				output_profile = .false.
				find_lburn = .true.
				do_pathological = do_pathological_init
				call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
				
				write(sweep_io,'(999ES25.9)') v_det, v_det/cs0, h_scale, &
      				rho_neumann, t_neumann, p_neumann, e_neumann, burn_u, &
      				!sonic_loc, max_mach, q_tot, t_end, vars(1:species)
      				pathological_loc, max_mach, q_tot, l_burn, vars(1:species)
								
        	end do !sweep over v_det
        	
        	call system_clock(sys_time_end,clock_rate)
			write(*,'(a20,ES15.4,a2)') 'Computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate, 's'

       		!write(*,*) 'Finding l_burn as a function of v_det...'
!       		
!       		open(unit=sweep_io, file=output_sweep_name)
!			write(sweep_io,'(999a25)') 'rho0 (g/cm^3)', 'T0 (K)', 'P0 (dyne/cm^2)', &
!         		'cs0 (cm/s)', chem_isos% name(chem_id(1:species)), 'x_end (cm)'
!         	write(sweep_io,'(999ES25.9)') rho0, t0, p0, cs0, xa(1:species), burn_time
!         	write(sweep_io,*)
!			write(sweep_io,'(999a25)') 'v_det (cm/s)', 'v_det/cs0', 'h_scale (cm)', &
!				'rho_neumann (g/cm^3)', 'T_neumann (K)', 'P_neumann (dyne/cm^2)', &
!				'e_neumann (erg/g)', 'burn_u (cm/s)', 'maxloc(ux/cs) (cm)', &
!				'max(ux/cs)', 'q_tot (erg/g)', 'l_burn (cm)', &
!				chem_isos% name(chem_id(1:species))
!        
!        	call system_clock(sys_time_begin,clock_rate)
!        	init_h_search_bracket_high = h_search_bracket_high
!			init_h_search_bracket_low = h_search_bracket_low
!			prev_h_search_bracket_high = h_search_bracket_high
!			prev_h_search_bracket_low = h_search_bracket_low
!			
!        	do i=1,sweep_znd_iters
!				write(*,*) 'Searching for H_scale...'
!				v_det = (sweep_min_mach + sweep_max_mach*(i-1)/(sweep_znd_iters-1))*cs0
!				write(*,*) 'v_det = ',v_det
!				
!				!Find the Neumann point (just once per v_det):
!				call solve_neumann
!				burn_rho = rho_neumann
!				burn_t = t_neumann
!				burn_u = (rho0/rho_neumann)*v_det
!				burn_P = p_neumann
!				burn_e = e_neumann
!				
!				!Reset all the values used in finding h_scale:
!				h_prev = 0d0
!				!Use what we know about the h_scale for the previous v_det to save some
!				!effort (we know the h_scale for a neighboring v_det will be close)
!				h_search_bracket_high = min(prev_h_search_bracket_high*1d1, init_h_search_bracket_high)
!				h_search_bracket_low = max(prev_h_search_bracket_low/1d1, init_h_search_bracket_low)
!			
!				do j=1,h_search_numsteps
!					h_scale = (h_search_bracket_high + h_search_bracket_low)/2d0
!					iwork_isolve = 0
!					work_isolve = 0
!					t_start = 0d0 				!Starting time?
!					t_end = h_search_xmax		!Ending time (s)?
!					
!					vars(1:species) = xa		!Composition
!					vars(species+1) = burn_rho	!Density (g/cm^3)
!					vars(species+2) = burn_T	!Temperature (K)
!					vars(species+3) = burn_u	!Velocity in shock frame (cm/s)
!					
!					if(do_blowout) then
!						vars(species+4) = uy_init	!Initial radial blowout velocity (cm/s)
!						if(use_uy_ux_ratio) then
!         					vars(species+4) = uy_ux_ratio*vars(species+3)
!         				endif
!						!grav = standard_cgrav*m_c/r_wd**2
!						!h_scale = P0/(rho0*grav)	!With gravity
!						vars(species+5) = h_scale 	!Initial radial scale height (cm)
!					else if(do_blowout_local) then
!						vars(species+4) = cs0	!Initial radial blowout velocity (cm/s)
!						vars(species+5) = h_scale 	!Initial radial scale height (cm)
!						vars(species+6) = delta_y_init	!Initial length scale for computing dP/dy (cm)
!					endif
!					if(do_curvature) then
!						if(use_rc_hs_scaling) r_curve = rc_hs_factor*h_scale
!						if(.not.use_he_clavin) vars(num_vars) = r_curve
!					endif
!	
!					sonic_loc = 0d0
!					max_mach = 0d0
!					gone_sonic = .false.
!					
!					do k=1,num_steps
!					!Evenly spaced timesteps in log
!					!t_start = burn_time*10**(log10(h_search_xmax)*((k-1)*1d0/num_steps-1))
!					!t_end = burn_time*10**(log10(h_search_xmax)*(k*1d0/num_steps-1))
!					!should be equivalent to:
!					t_start = h_search_xmax**((k-1)*1d0/num_steps)
!					t_end = h_search_xmax**(k*1d0/num_steps)
!					call isolve( &
!						which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
!						h, max_step_size, max_steps, &
!						rtol, atol, itol, &
!						znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
!						null_mas, imas, mlmas, mumas, &
!						znd_solout, iout, &
!						lapack_decsol, null_decsols, lrd, rpar_decsol, lid, ipar_decsol, &
!						work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
!						lrpar, rpar, lipar, ipar, &
!						lout, idid)
!					if (ierr /= 0) then
!						write(*,*) 'failed in isolve() call'
!						stop 1
!					end if
!					!Check if we've slowed down too much (good stopping criterion for curvature)
!					if (rpar(2).lt.1d-2) then
!						write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
!						write(*,*)
!						exit
!					endif
!					if((rpar(2).gt.max_mach).and.(.not.gone_sonic)) then
!						sonic_loc = t_end
!						max_mach = rpar(2)
!						if((rpar(2).ge.sonic_limit).and.(.not.gone_sonic)) then
!							gone_sonic = .true.
!							write(*,*) 'sonic point hit - exiting integration...'
!							write(*,*)
!							exit
!						endif
!					endif
!					end do !Integration loop
!					
!					write(*,'(a15,ES25.9)') 'u_x (cm/s)', vars(species+3)
!					write(*,'(a15,ES25.9)') 'u_x/c_s', rpar(2)
!					write(*,'(a15,ES25.9)') 'max(u/cs)', max_mach
!					write(*,'(a15,ES25.9)') 'maxloc(u/cs)', sonic_loc
!					write(*,'(a15,ES25.9)') 'q_final', rpar(9)
!					write(*,'(a15,ES25.9)') 'h_scale', h_scale
!					write(*,'(a15,ES25.9)') 'delta_h', abs(h_scale - h_prev)/h_scale
!			
!					!write(*,'(99ES25.9)') v_det, v_det/cs0, sonic_loc, max_mach
!					
!					!Check if we've found a solution:
!					if(abs(h_scale - h_prev)/h_scale.le.h_search_tol) then
!						write(*,*) 'Solution found:'
!						write(*,'(a15,ES25.9)') 'h_scale (cm)', h_scale
!						exit
!					endif
!					
!					!Update the search paramaters:
!					!For finding a solution where the burning stops at the sonic point
!					!If we hit a sonic point then decrease the scale height, otherwise increase it
!					if(gone_sonic) then
!						h_search_bracket_high = h_scale
!					else
!						h_search_bracket_low = h_scale
!					endif
!					
!					!For finding a solution where the burning stops @ 1d7cm (doesn't work...)
!					!if((max_mach.gt.0.99).and.(sonic_loc.lt.1d7)) then
!					!	h_search_bracket_high = h_scale
!					!else
!					!	h_search_bracket_low = h_scale
!					!endif
!					
!					write(*,*) 'Brackets now:'
!					write(*,'(a15,ES25.9)') 'high (cm)', h_search_bracket_high
!					write(*,'(a15,ES25.9)') 'low (cm)', h_search_bracket_low
!					
!					h_prev = h_scale
!					
!				end do !search for h_scale
!				
!				!h_scale = h_search_bracket_low !To ensure we don't hit the sonic point, gives frozen subsonic solution
!				h_scale = h_search_bracket_high !To ensure we can go through the pathological point, gives frozen supersonic solution
!				prev_h_search_bracket_high = h_search_bracket_high
!				prev_h_search_bracket_low = h_search_bracket_low
!				write(*,*)
!				write(*,*) 'H_scale found: ',h_scale,' cm'
!				
!				!Finally, check if we've hit the upper limit of H_scale during our sweep.
!      			!If so, then we can stop our sweep, as we'll only need larger H_scale
!      			!values as we increase v_det from here
!      			if(h_search_bracket_high.eq.init_h_search_bracket_high) then
!      				write(*,*) 'Maximum scale height of ',init_h_search_bracket_high,' cm hit - stopping sweep'
!      				exit
!      			endif
!				
!				!Now that we've found a valid h_scale, we need to find the ZND length
!				write(*,*) 'Calculating ZND length...'
!				
!				iwork_isolve = 0
!				work_isolve = 0
!				t_start = 0d0 				!Starting time?
!				t_end = h_search_xmax		!Ending time (s)?
!				
!				vars(1:species) = xa		!Composition
!				vars(species+1) = burn_rho	!Density (g/cm^3)
!				vars(species+2) = burn_T	!Temperature (K)
!				vars(species+3) = burn_u	!Velocity in shock frame (cm/s)
!				
!				if(do_blowout) then
!					vars(species+4) = uy_init	!Initial radial blowout velocity (cm/s)
!					if(use_uy_ux_ratio) then
!         				vars(species+4) = uy_ux_ratio*vars(species+3)
!         			endif
!					!grav = standard_cgrav*m_c/r_wd**2
!					!h_scale = P0/(rho0*grav)	!With gravity
!					vars(species+5) = h_scale 	!Initial radial scale height (cm)
!				else if(do_blowout_local) then
!					vars(species+4) = cs0	!Initial radial blowout velocity (cm/s)
!					vars(species+5) = h_scale 	!Initial radial scale height (cm)
!					vars(species+6) = delta_y_init	!Initial length scale for computing dP/dy (cm)
!				endif
!				if(do_curvature) then
!					if(use_rc_hs_scaling) r_curve = rc_hs_factor*h_scale
!					if(.not.use_he_clavin) vars(num_vars) = r_curve
!				endif
!
!				sonic_loc = 0d0
!				max_mach = 0d0
!				gone_sonic = .false.
!				
!				do j=1,num_steps
!				!Evenly spaced timesteps in log
!				!t_start = h_search_xmax*10**(log10(h_search_xmax)*((k-1)*1d0/num_steps-1))
!				!t_end = h_search_xmax*10**(log10(h_search_xmax)*(k*1d0/num_steps-1))
!				!should be equivalent to:
!				t_start = h_search_xmax**((j-1)*1d0/num_steps)
!				t_end = h_search_xmax**(j*1d0/num_steps)
!				
!				!Try linear spaced steps to get better resolution on l_burn since we know it's
!				!going to be near the end of integration for our fiducial of 95% energy generation:
!				!t_start = (j-1)*burn_time/num_steps
!				!t_end = j*burn_time/num_steps
!				call isolve( &
!					which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
!					h, max_step_size, max_steps, &
!					rtol, atol, itol, &
!					znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
!					null_mas, imas, mlmas, mumas, &
!					znd_solout, iout, &
!					lapack_decsol, null_decsols, lrd, rpar_decsol, lid, ipar_decsol, &
!					work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
!					lrpar, rpar, lipar, ipar, &
!					lout, idid)
!				if (ierr /= 0) then
!					write(*,*) 'failed in isolve() call'
!					stop 1
!				end if
!				!Check if we've slowed down too much (good stopping criterion for curvature)
!				if (rpar(2).lt.1d-2) then
!					write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
!					write(*,*)
!					exit
!				endif
!				if((rpar(2).gt.max_mach).and.(.not.gone_sonic)) then
!					sonic_loc = t_end
!					pathological_loc = t_end
!					max_mach = rpar(2)
!					if((rpar(2).ge.sonic_limit_pathological).and.(.not.gone_sonic)) then
!						gone_sonic = .true.
!						write(*,*) 'sonic point hit - exiting integration...'
!						write(*,*)
!						exit
!					endif
!				endif
!				
!				end do !Integration loop
!				
!				!Take the derivatives of all the quantities and linearize them over
!				!a certain length scale (how to choose this length?) A first choice is to
!				!use the time step we would have taken if the integration succeeded.
!				!linearization_length = 4d6	!cm
!				linearization_length = pathological_loc*pathological_linearization_ratio
!				write(*,*) 'vars pre:', vars
!				vars(1:num_vars) = vars(1:num_vars) + rpar(21:21+num_vars)*linearization_length
!				write(*,*) 'vars post:', vars
!				
!				!Now try resuming integration:
!				sonic_loc = 0d0
!				max_mach = 0d0
!				gone_sonic = .false.
!				do j=1,num_steps
!				!Evenly spaced timesteps in log
!				t_start = pathological_loc + linearization_length + &
!					(burn_time - pathological_loc + linearization_length)**((j-1)*1d0/num_steps)
!				t_end = pathological_loc + linearization_length + & 
!					(burn_time - pathological_loc + linearization_length)**(j*1d0/num_steps)
!				!write(*,*) 'Integrating to ',t_end
!				!write(*,*) 'iwork_isolve(2) ',iwork_isolve(2)
!				call isolve( &
!					which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
!					h, max_step_size, max_steps, &
!					rtol, atol, itol, &
!					znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
!					null_mas, imas, mlmas, mumas, &
!					znd_solout, iout, &
!					lapack_decsol, null_decsols, lrd, rpar_decsol, lid, ipar_decsol, &
!					work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
!					lrpar, rpar, lipar, ipar, &
!					lout, idid)
!				if (ierr /= 0) then
!					write(*,*) 'failed in isolve() call'
!					stop 1
!				endif
!				!Check if we've slowed down too much (good stopping criterion for curvature)
!				if (rpar(2).lt.1d-2) then
!					write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
!					write(*,*)
!					exit
!				endif
!				
!				!write(profile_io, '(i30,999ES30.20)') k, t_end, vars, rpar(1), &
!				!	rpar(8), rpar(10), rpar(2), rpar(9), mass_flux, mom_flux, energy_flux, &
!				!	rpar(11), rpar(12), rpar(13), rpar(14), rpar(15), rpar(16)
!					
!					!dP_dx, dE_dx, dq_dx, detot_dx, dEtot_dx_eqn
!					!rpar(6), rpar(7), rpar(8:8+num_vars-1), rpar(8+num_vars)
!				end do !num_steps
!				
!				q_tot = rpar(9)
!				write(*,*) 'Initial integration found q = ',q_tot,' erg/g'
!				
!				!Integrate again, but stop when we hit q/2:
!				!(we know the Neumann condition, so don't need to find that again)
!				iwork_isolve = 0
!				work_isolve = 0
!				t_start = 0d0 				!Starting time?
!				t_end = h_search_xmax		!Ending time (s)?
!				
!				vars(1:species) = xa		!Composition
!				vars(species+1) = burn_rho	!Density (g/cm^3)
!				vars(species+2) = burn_T	!Temperature (K)
!				vars(species+3) = burn_u	!Velocity in shock frame (cm/s)
!				
!				if(do_blowout) then
!					vars(species+4) = uy_init	!Initial radial blowout velocity (cm/s)
!					if(use_uy_ux_ratio) then
!         				vars(species+4) = uy_ux_ratio*vars(species+3)
!         			endif
!					!grav = standard_cgrav*m_c/r_wd**2
!					!h_scale = P0/(rho0*grav)	!With gravity
!					vars(species+5) = h_scale 	!Initial radial scale height (cm)
!				else if(do_blowout_local) then
!					vars(species+4) = cs0	!Initial radial blowout velocity (cm/s)
!					vars(species+5) = h_scale 	!Initial radial scale height (cm)
!					vars(species+6) = delta_y_init	!Initial length scale for computing dP/dy (cm)
!				endif
!				if(do_curvature) then
!					if(use_rc_hs_scaling) r_curve = rc_hs_factor*h_scale
!					if(.not.use_he_clavin) vars(num_vars) = r_curve
!				endif
!
!				!Aha, we don't need to re-check for hitting a sonic point since we're
!				!only integrating out halfway:
!				!sonic_loc = 0d0
!				!max_mach = 0d0
!				!gone_sonic = .false.
!				
!				do j=1,num_steps
!				!Evenly spaced timesteps in log
!				t_start = burn_time**((j-1)*1d0/num_steps)
!				t_end = burn_time**(j*1d0/num_steps)
!				call isolve( &
!					which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
!					h, max_step_size, max_steps, &
!					rtol, atol, itol, &
!					znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
!					null_mas, imas, mlmas, mumas, &
!					znd_solout, iout, &
!					lapack_decsol, null_decsols, lrd, rpar_decsol, lid, ipar_decsol, &
!					work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
!					lrpar, rpar, lipar, ipar, &
!					lout, idid)
!				if (ierr /= 0) then
!					write(*,*) 'failed in isolve() call'
!					stop 1
!				end if
!				!Check if we've hit q_tot*l_burn_factor:
!				if (rpar(9).ge.q_tot*l_burn_factor) then
!					write(*,'(a15,ES25.9)') 'l_burn = ',t_end
!					write(*,*)
!					exit
!				endif
!				!Check if we've slowed down too much (good stopping criterion for curvature)
!				if (rpar(2).lt.1d-2) then
!					write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
!					write(*,*)
!					exit
!				endif
!				if((rpar(2).gt.max_mach).and.(.not.gone_sonic)) then
!					sonic_loc = t_end
!					pathological_loc = t_end
!					max_mach = rpar(2)
!					if((rpar(2).ge.sonic_limit_pathological).and.(.not.gone_sonic)) then
!						gone_sonic = .true.
!						write(*,*) 'sonic point hit - exiting integration...'
!						write(*,*)
!						exit
!					endif
!				endif
!				
!				end do !Integration loop
!				
!				!If we haven't found the ZND length then go through the pathological point:
!				if (rpar(9).lt.q_tot*l_burn_factor) then
!					linearization_length = pathological_loc*pathological_linearization_ratio
!					write(*,*) 'vars pre:', vars
!					vars(1:num_vars) = vars(1:num_vars) + rpar(21:21+num_vars)*linearization_length
!					write(*,*) 'vars post:', vars
!					
!					!Now try resuming integration:
!					sonic_loc = 0d0
!					max_mach = 0d0
!					gone_sonic = .false.
!					do j=1,num_steps
!					!Evenly spaced timesteps in log
!					t_start = pathological_loc + linearization_length + &
!						(burn_time - pathological_loc + linearization_length)**((j-1)*1d0/num_steps)
!					t_end = pathological_loc + linearization_length + & 
!						(burn_time - pathological_loc + linearization_length)**(j*1d0/num_steps)
!					!write(*,*) 'Integrating to ',t_end
!					!write(*,*) 'iwork_isolve(2) ',iwork_isolve(2)
!					call isolve( &
!						which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
!						h, max_step_size, max_steps, &
!						rtol, atol, itol, &
!						znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
!						null_mas, imas, mlmas, mumas, &
!						znd_solout, iout, &
!						lapack_decsol, null_decsols, lrd, rpar_decsol, lid, ipar_decsol, &
!						work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
!						lrpar, rpar, lipar, ipar, &
!						lout, idid)
!					if (ierr /= 0) then
!						write(*,*) 'failed in isolve() call'
!						stop 1
!					endif
!					!Check if we've hit q_tot*l_burn_factor:
!					if (rpar(9).ge.q_tot*l_burn_factor) then
!						write(*,'(a15,ES25.9)') 'l_burn = ',t_end
!						write(*,*)
!						exit
!					endif
!					!Check if we've slowed down too much (good stopping criterion for curvature)
!					if (rpar(2).lt.1d-2) then
!						write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
!						write(*,*)
!						exit
!					endif
!					end do !num_steps
!				endif
!				
!				write(sweep_io,'(999ES25.9)') v_det, v_det/cs0, h_scale, &
!      				rho_neumann, t_neumann, p_neumann, e_neumann, burn_u, &
!      				!sonic_loc, max_mach, q_tot, t_end, vars(1:species)
!      				pathological_loc, rpar(2), q_tot, t_end, vars(1:species)
!								
!        	end do !sweep over v_det
!        	
!        	call system_clock(sys_time_end,clock_rate)
!			write(*,'(a20,ES15.4,a2)') 'Computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate, 's'

		!This section will sweep through M_env values (mapped to h_scale), for each one 
		!finding the minimum H_scale where a sonic point is not hit and integrate 
		!twice - once to find the total energy release, abundances, etc. and again to 
		!find the ZND length
        else if(do_lznd_wd) then
       		write(*,*) 'Finding l_burn as a function of M_env...'
       		
       		t_b = t0
       		call composition_info(species, chem_id, xa, xh, xhe, Zm, abar, zbar, z2bar, &
            	ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
            !First calculate r_wd and grav given m_c:
			r_wd = 7.62d-3*rsun*m_c**(-1.616)		!WD mass-radius relation
			grav = (standard_cgrav*msun/rsun**2)*(1.72d4*m_c**4.23)	!surface gravity (cm/s^2)
       		
       		!write(*,*) 'm_c = ',m_c
       		!write(*,*) 'r_wd = ',r_wd
			!write(*,*) 'grav = ',grav
			
			if(wd_env_use_le)  then
				xi1_low = 1.0
				xi1_high = 5.0
				call find_xi1(n_le, xi1_low, xi1_high, xi1, dtheta_dxi1)
				write(*,*) 'Lane-Emden solution:'
				write(*,'(a10,ES20.10)') 'xi1:', xi1
				write(*,'(a10,ES20.10)') 'dtheta_dxi1:', dtheta_dxi1
				write(*,*)
				call get_aux_vars(m_c, r_wd, xi1, dtheta_dxi1)
			endif
       		
       		open(unit=sweep_io, file=output_sweep_name)
			write(sweep_io,'(999a25)') 'm_c (M_sun)', 'R_wd (M_sun)', 'T_b (K)', &
         		'g_b (cm/s^2)', chem_isos% name(chem_id(1:species)), 'x_end (cm)'
         	write(sweep_io,'(999ES25.9)') m_c, r_wd, t0, grav, xa(1:species), burn_time
         	write(sweep_io,*)
			write(sweep_io,'(999a25)') 'rho_b (g/cm^3)', 'rho0 (g/cm^3)', &
				'P_b (dyne/cm^2)', 'P0 (dyne/cm^2)', &
				'M_env (Msun)', 't_sound (s)', 'v_min (cm/s)', &
				'v_det (cm/s)', 'v_det/cs0', 'h_scale (cm)', &
				'rho_neumann (g/cm^3)', 'T_neumann (K)', 'P_neumann (dyne/cm^2)', &
				'e_neumann (erg/g)', 'burn_u (cm/s)', 'maxloc(ux/cs) (cm)', &
				'max(ux/cs)', 'q_tot (erg/g)', 'l_burn (cm)', &
				chem_isos% name(chem_id(1:species))
				
        	call system_clock(sys_time_begin,clock_rate)
        	init_vdet_search_bracket_high = vdet_search_bracket_high
			init_vdet_search_bracket_low = vdet_search_bracket_low
			prev_vdet_search_bracket_high = vdet_search_bracket_high
			prev_vdet_search_bracket_low = vdet_search_bracket_low
			
			!If we're using the MESA EOS to build our envelopes, then initialize our
			!guesses for rho_c and pb_frac here:
			if(wd_env_use_mesa) then
				menv_target = m_env_min
				mwd_target = m_c + menv_target
				rhoc_min = 2d6
				rhoc_max = 6d10
				call find_rhoc(mwd_target, rhoc_min, rhoc_max, rho_c)
				pb_frac = 1d-4
			endif
			
        	do i=1,sweep_znd_iters
				write(*,*) 'Searching for v_det...'
				
				!Envelope mass to use (in msun)
				m_env = m_env_min*(m_env_max/m_env_min)**((i-1)*1d0/(sweep_znd_iters-1))
				
				if(wd_env_use_le)  then
					!In the polytrope model, we construct a WD with total mass m_c, then
					!integrate the lane-emden equations to find where the xi coordinate
					!where the enclosed mass is equal to the core mass (WD - env)
					m_tot = m_c + m_env
					r_tot = 7.62d-3*rsun*m_tot**(-1.616)
					!m_c = m_c - m_env
		
					call get_aux_vars(m_tot, r_tot, xi1, dtheta_dxi1)
					xi_low = 1d-4
					xi_high = xi1
					call find_xi_mc(m_c, xi_low, xi_high, xi, dtheta_dxi)
		
					theta = theta_le(xi,dtheta_dxi,lrpar_le,rpar_le,lipar_le,ipar_le,ierr_le)
					m_xi = -4*pi*r_n**3*rho_c*xi**2*dtheta_dxi/msun
					p_b = k_le*rho_c**(1+1d0/n_le)*theta**(n_le+1)
					Pgas = p_b - crad*T_b**4/3d0
					grav = standard_cgrav*m_xi*msun/(r_n*xi)**2
					rho_b = rho_c*theta**n_le
					
					gamma1 = 1 + 1d0/n_le
				elseif(wd_env_use_mesa) then
					write(*,*) 'Searching for envelope solution for'
					write(*,'(a10,ES20.10)') 'm_c', m_c
					write(*,'(a10,ES20.10)') 'm_env', m_env
					write(*,*)
					call find_envelope_newton(m_c, m_env, rho_c, pb_frac, &
						rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
						
					write(*,*) 'Envelope solution found'
					Pgas = p_b - crad*T_b**4/3d0
					gamma1 = 5d0/3d0
				else
					p_b = m_env*msun*grav/(4*pi*r_wd**2)
					Pgas = p_b - crad*T_b**4/3d0
					logPgas = log10(Pgas)
					logT = log10(T_b)
					call eosPT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, &
						xa, Pgas, logPgas, T_b, logT, &
						rho_b, logRho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
						res, d_dlnRho_const_T, d_dlnT_const_Rho, &
						d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
				
					gamma1 = res(i_gamma1)	
				endif
				
				write(*,*) 'm_env = ',m_env
				write(*,*) 'p_b = ',p_b
				write(*,*) 'Pgas = ',Pgas
				
				rho0 = rho_b*(1-(gamma1-1)*shock_lead_hp_frac/gamma1)**(1d0/(gamma1-1))
				p0 = p_b*(1 + (1d0/gamma1 - 1)*shock_lead_hp_frac)**(gamma1/(gamma1-1))
				h_scale = p_b/(2*rho_b*grav)
				cs0 = sqrt(gamma1*p0/rho0)
				
				write(*,*) 'rho0 = ',rho0
				write(*,*) 'h_scale = ',h_scale
				
				!If our h_scale is too large (typically > 1d9cm) then it's no longer a
				!thin shell so we can't do our calculation - stop
				if(h_scale.ge.h_search_bracket_high) then
					write(*,*) 'exiting - h_scale became larger than ',h_search_bracket_high
					exit	!Going to higher M_env only makes things worse so exit
				endif
				
				!If M_tot = m_c + M_env > 1.4 then stop since we're above the Chandrasekhar mass:
				if((m_c + m_env).ge.1.5) then
					write(*,*) 'exiting - total mass above M_ch: ',(m_c + m_env)
					exit	!Going to higher M_env only makes things worse so exit
				endif
				
				!If our rho0 is too small (typically < 1d5cm) then it's no longer a
				!thin shell so we can't do our calculation - stop
				if(rho0.le.env_rho_min) then
					write(*,*) 'exiting - rho0 less than env_rho_min'
					cycle	!Going to higher M_env helps so skip to next loop iteration
				endif
				
				!If our rho0 is too large (typically > 5d6cm) then the integration
				!takes way too long - stop
				if(rho0.ge.env_rho_max) then
					write(*,*) 'exiting - rho0 greater than env_rho_max'
					cycle	!Going to higher M_env helps so skip to next loop iteration
				endif
				
				!Reset all the values used in finding h_scale:
				v_prev = 0d0
				!Use what we know about the h_scale for the previous v_det to save some
				!effort (we know the h_scale for a neighboring v_det will be close)
				!h_search_bracket_high = min(prev_h_search_bracket_high*1d1, init_h_search_bracket_high)
				!h_search_bracket_low = max(prev_h_search_bracket_low/1d1, init_h_search_bracket_low)
				
				vdet_search_bracket_low = init_vdet_search_bracket_low
				vdet_search_bracket_high = init_vdet_search_bracket_high
				do j=1,vdet_search_numsteps
					v_det = (vdet_search_bracket_high + vdet_search_bracket_low)/2d0

					!Find the Neumann point:
					call solve_neumann
					burn_rho = rho_neumann
					burn_t = t_neumann
					burn_u = (rho0/rho_neumann)*v_det
					burn_P = p_neumann
					burn_e = e_neumann
				
					t_start = 0d0 				!Starting location
					t_end = vdet_search_xmax	!Ending location
					output_profile = .false.
					find_lburn = .false.
					do_pathological = .false.
					call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)

					write(*,'(a15,ES25.9)') 'u_x (cm/s)', vars(species+3)
					write(*,'(a15,ES25.9)') 'u_x/c_s', rpar(2)
					write(*,'(a15,ES25.9)') 'max(u/cs)', max_mach
					write(*,'(a15,ES25.9)') 'maxloc(u/cs)', sonic_loc
					write(*,'(a15,ES25.9)') 'q_final', rpar(9)
					write(*,'(a15,ES25.9)') 'v_det', v_det
					write(*,'(a15,ES25.9)') 'delta_vdet', abs(v_det - v_prev)/v_det
				
					!Check if we've found a solution:
					if(abs(v_det - v_prev)/v_det.le.vdet_search_tol) then
						write(*,*) 'Generalized CJ solution found:'
						if(do_pathological_init) then !take the case that will hit the sonic point
							v_det = vdet_search_bracket_low
						else !take the case that won't hit the sonic point
							v_det = vdet_search_bracket_high
						endif
						write(*,'(a15,ES25.9)') 'v_det (cm/s)', v_det
						exit
					endif
				
					!Update the search paramaters:
					!For finding a solution where the burning stops at the sonic point
					!If we hit a sonic point then increase v_det, otherwise decrease it
					if(max_mach.ge.sonic_limit) then
						vdet_search_bracket_low = v_det
					else
						vdet_search_bracket_high = v_det
					endif
				
					write(*,*) 'Brackets now:'
					write(*,'(a15,ES25.9)') 'high (cm/s)', vdet_search_bracket_high
					write(*,'(a15,ES25.9)') 'low (cm/s)', vdet_search_bracket_low
				
					v_prev = v_det
				
				end do !sweep through v_det
				
				!Check if we have a valid solution. If we never hit a sonic point during
				!our search, then there's not a propagating solution
				if(v_det.eq.init_vdet_search_bracket_low) then
					write(*,*) 'exiting - v_det hit our lower bracket - no propagating solution'
					cycle 	!since going to higher m_env will help
				endif
				
				!h_scale = h_search_bracket_low !To ensure we don't hit the sonic point, gives frozen subsonic solution
				!h_scale = h_search_bracket_high !To ensure we can go through the pathological point, gives frozen supersonic solution
				!prev_h_search_bracket_high = h_search_bracket_high
				!prev_h_search_bracket_low = h_search_bracket_low
				write(*,*)
				write(*,*) 'v_det found: ',v_det,' cm/s'
				
				!Now that we've found a valid h_scale, we need to find the ZND length
				write(*,*) 'Calculating ZND length...'
				
				t_start = 0d0 					!Starting location
				t_end = burn_time				!Ending location
				output_profile = .false.
				find_lburn = .true.
				do_pathological = do_pathological_init
				call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
				
				!Final step is calculating some supplementary information:
				
				!Time it takes the detonation to sweep across the WD to the antipode
				t_det = pi*r_wd/v_det
				
				!Time it takes for a sound wave to cross the WD:
				!t_sound = 2*
				
				write(sweep_io,'(999ES25.9)') rho_b, rho0, p_b, p0, &
					m_env, t_sound, r_c*pi/t_sound, &
					v_det, v_det/cs0, h_scale, &
      				rho_neumann, t_neumann, p_neumann, e_neumann, burn_u, &
      				!sonic_loc, max_mach, q_tot, t_end, vars(1:species)
      				!pathological_loc, rpar(2), q_tot, l_burn, vars(1:species)
      				pathological_loc, max_mach, q_tot, l_burn, vars(1:species)
								
        	end do !sweep over v_det
        	
        	call system_clock(sys_time_end,clock_rate)
			write(*,'(a20,ES15.4,a2)') 'Computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate, 's'
		
		else if(do_reconstruct_lznd_vdet) then
       		write(*,*) 'Reconstructing lznd_vdet data file: ',reconstruct_lznd_vdet_infile
       		
       		!Input file:
       		open(unit=data_in, file=reconstruct_lznd_vdet_infile)
       		!Skip the header line:
       		read(data_in,*) 
       		read(data_in,*) rho0, t0, p0, cs0, xa(1:species), burn_time
       		read(data_in,*)
       		read(data_in,*)
       		
       		!Output file:
       		open(unit=sweep_io, file=reconstruct_lznd_vdet_outfile)
			write(sweep_io,'(999a25)') 'rho0 (g/cm^3)', 'T0 (K)', 'P0 (dyne/cm^2)', &
         		'cs0 (cm/s)', chem_isos% name(chem_id(1:species)), 'x_end (cm)'
         	write(sweep_io,'(999ES25.9)') rho0, t0, p0, cs0, xa(1:species), burn_time
         	write(sweep_io,*)
         	!First two lines are for adding white dwarf (WD) values:
			write(sweep_io,'(999a25)') 'rho_b (g/cm^3)', 'P_b (dyne/cm^2)', 'g (cm/s^2)', &
				'm_c (Msun)', 'R_wd (cm)', 'DM (Msun)', &
				'v_det (cm/s)', 'v_det/cs0', 'h_scale (cm)', &
				'rho_neumann (g/cm^3)', 'T_neumann (K)', 'P_neumann (dyne/cm^2)', &
				'e_neumann (erg/g)', 'burn_u (cm/s)', 'maxloc(ux/cs) (cm)', &
				'max(ux/cs)', 'q_tot (erg/g)', 'l_burn (cm)', &
				chem_isos% name(chem_id(1:species))
				
			ios = 0
			!Loop through input file, reading off v_det and h_scale values to use:
			do while(ios.eq.0)    
				!Read off the relevant parameters (just the first three)
				read(data_in,*,iostat=ios) v_det, det_mach, h_scale, burn_rho, &
					burn_T, p_neumann, e_neumann, burn_u, &
					sonic_loc, max_mach, q_tot, l_burn, vars(1:species)
				xa(1:species) = vars(1:species)
				
				if(ios.ne.0) exit
					
				!All of this code is for adding WD variables (mapping (rho_0,H) to (M,DM))
				!to a data file-----------------------------------------------------------
				write(*,*) 'Mapping (rho0,H) to (M,DM) using'
				write(*,'(a15,ES25.9)') 'rho0 (g/cm^3):', rho0
				write(*,'(a15,ES25.9)') 'H (cm):', h_scale
				
				rho_b = rho0/0.69	!Base of envelope relation from gamma=5/3
				T_b = T0
				
				call composition_info(species, chem_id, xa, xh, xhe, Zm, abar, zbar, z2bar, &
					ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
				logRho = log10(rho_b)
				logT = log10(t_b)
				call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, &
					xa, rho_b, logRho, t_b, logT, &
					res, d_dlnRho_const_T, d_dlnT_const_Rho, &
					d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
				p_b = exp(res(i_lnPgas)) + crad*t_b**4/3d0
				
				H_p = 2*h_scale			!Pressure scale height (H_p) from layer thickness (h_scale)
				g0 = P_b/(rho_b*H_p)	!constant gravity in envelope
				
				!From Bill W.'s plot, it looks like
				m_c = ((g0/(standard_cgrav*msun/(rsun**2)))*5.80d-5)**(1d0/4.23)	!(msun)
				r_wd = 7.62d-3*(m_c**(-1.616))*rsun	!(cm)
				DM = 4*pi*r_wd**2*p_b/(g0*msun)		!(msun)
				
				write(sweep_io,'(999ES25.9)') rho_b, P_b, g0, m_c, r_wd, DM, &
					v_det, det_mach, h_scale, &
      				burn_rho, burn_t, p_neumann, e_neumann, burn_u, &
      				sonic_loc, max_mach, q_tot, l_burn, vars(1:species)
				
				!-------------------------------------------------------------------------
					
				!All of this code is for adding/changing ZND length to a data file--------
				!write(*,*) 'Calculating ZND length using'
				!write(*,'(a15,ES25.9)') 'v_det (cm/s):', v_det
				!write(*,'(a15,ES25.9)') 'h_scale (cm):', h_scale
				
				!t_start = 0d0 				!Starting location
				!t_end = burn_time			!Ending location
				!output_profile = .false.
				!find_lburn = .true.
				!call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
				
				!Variables are a little weird here since some are read from a file:
				!write(sweep_io,'(999ES25.9)') v_det, det_mach, h_scale, &
      			!	burn_rho, burn_t, p_neumann, e_neumann, burn_u, &
      			!	sonic_loc, max_mach, q_tot, l_burn, vars(1:species)
				!-------------------------------------------------------------------------
				
      		end do !Reading from infile
      		
      	else if(do_rho_sweep) then
       		write(*,*) 'Sweeping through densities...'
       		
       		!do_rho_sweep = .false.
			!d_hs_scaling = 5.0				!d = d_hs_scaling*h_scale
			!rho_sweep_T0 = 1e8				!Initial temerature (for calculating scale height)
			!rho_sweep_grav = 3d8			!cm/s^2 (M_sun in R_earth is 3.3d8 cm/s^2)
			!rho_sweep_log_min = 5.0		!Lowest density to use
			!rho_sweep_log_max = 7.0		!Highest density to use
			!rho_sweep_num_steps = 100		!Number of steps in density to take
			
			!t0 = rho_sweep_t0
			logT = log10(t0)
			!grav_init = rho_sweep_grav
       		
       		open(unit=sweep_io, file=output_sweep_name)
			write(sweep_io,'(999a25)') 'T0 (K)', &
				'grav0 (cm/s^2)', 'cs0 (cm/s)', chem_isos% name(chem_id(1:species)), &
				'x_end (cm)'
         	write(sweep_io,'(999ES25.9)') t0, grav_init, &
         		cs0, xa(1:species), burn_time
         	write(sweep_io,*)
			write(sweep_io,'(999a25)') 'rho0 (g/cm^3)', 'P0', &
				'v_det (cm/s)', 'v_det/cs0', 'h_scale (cm)', &
				'rho_neumann (g/cm^3)', 'T_neumann (K)', 'P_neumann (dyne/cm^2)', &
				'e_neumann (erg/g)', 'burn_u (cm/s)', 'maxloc(ux/cs) (cm)', &
				'max(ux/cs)', 'q_tot (erg/g)', 'l_burn (cm)', &
				chem_isos% name(chem_id(1:species))
				
			init_vdet_search_bracket_high = vdet_search_bracket_high
			init_vdet_search_bracket_low = vdet_search_bracket_low
			prev_vdet_search_bracket_high = vdet_search_bracket_high
			prev_vdet_search_bracket_low = vdet_search_bracket_low
			
			!Sweep through densities:
        	do i=1,rho_sweep_num_steps
        		rho0 = 10**(rho_sweep_log_min + &
        			(rho_sweep_log_max-rho_sweep_log_min)*(i-1d0)/(rho_sweep_num_steps-1d0))
				logRho = log10(rho0)
         		call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species, chem_id, net_iso, &
         			xa, rho0, logRho, t0, logT, &
         			res, d_dlnRho_const_T, d_dlnT_const_Rho, &
         			d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
         		p0 = exp(res(i_lnPgas)) + crad*t0**4/3d0
         		cs0 = sqrt(res(i_gamma1)*p0/rho0)
				h_scale = p0/(rho0*grav_init)
				write(*,*) 'Thickness d = ',h_scale
				h_scale = h_scale/d_hs_scaling
				write(*,*) 'H_scale = ',h_scale
				
				v_prev = 0d0
				call system_clock(sys_time_begin,clock_rate)
				
				if(do_cj) then
					call solve_cj
					q_tot = q
					
					!Find the Neumann point:
					call solve_neumann
					burn_rho = rho_neumann
					burn_t = t_neumann
					burn_u = (rho0/rho_neumann)*v_det
					burn_P = p_neumann
					burn_e = e_neumann
				else
					write(*,*) 'Searching for v_det...'
					
					vdet_search_bracket_high = init_vdet_search_bracket_high
					vdet_search_bracket_low = init_vdet_search_bracket_low
					
					do j=1,vdet_search_numsteps
						v_det = (vdet_search_bracket_high + vdet_search_bracket_low)/2d0
						
						!Find the Neumann point:
						call solve_neumann
						burn_rho = rho_neumann
						burn_t = t_neumann
						burn_u = (rho0/rho_neumann)*v_det
						burn_P = p_neumann
						burn_e = e_neumann
										
						t_start = 0d0 				!Starting location
						t_end = vdet_search_xmax	!Ending location
						output_profile = .false.
						find_lburn = .false.
						do_pathological = .false.
						call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
						
						write(*,'(a15,ES25.9)') 'u_x (cm/s)', vars(species+3)
						write(*,'(a15,ES25.9)') 'u_x/c_s', rpar(2)
						write(*,'(a15,ES25.9)') 'max(u/cs)', max_mach
						write(*,'(a15,ES25.9)') 'maxloc(u/cs)', sonic_loc
						write(*,'(a15,ES25.9)') 'q_final', rpar(9)
						write(*,'(a15,ES25.9)') 'v_det', v_det
						write(*,'(a15,ES25.9)') 'delta_vdet', abs(v_det - v_prev)/v_det
						
						!Check if we've found a solution:
						if(abs(v_det - v_prev)/v_det.le.vdet_search_tol) then
							write(*,*) 'Generalized CJ solution found:'
							if(do_pathological_init) then !take the case that will hit the sonic point
								v_det = vdet_search_bracket_low
							else !take the case that won't hit the sonic point
								v_det = vdet_search_bracket_high
							endif
							write(*,'(a15,ES25.9)') 'v_det (cm/s)', v_det
							exit
						endif
						
						!Update the search paramaters:
						!For finding a solution where the burning stops at the sonic point
						!If we hit a sonic point then increase v_det, otherwise decrease it
						if(max_mach.gt.sonic_limit) then
							vdet_search_bracket_low = v_det
						else
							vdet_search_bracket_high = v_det
						endif
						
						write(*,*) 'Brackets now:'
						write(*,'(a15,ES25.9)') 'high (cm/s)', vdet_search_bracket_high
						write(*,'(a15,ES25.9)') 'low (cm/s)', vdet_search_bracket_low
						
						v_prev = v_det
						
					end do !sweep through v_det
				endif
					
				!Now that we've found a valid v_det, we need to find the ZND length
				write(*,*) 'Calculating ZND length...'
				
				!Find the Neumann point:
				call solve_neumann
				burn_rho = rho_neumann
				burn_t = t_neumann
				burn_u = (rho0/rho_neumann)*v_det
				burn_P = p_neumann
				burn_e = e_neumann
				
				t_start = 0d0 				!Starting location
				t_end = burn_time			!Ending location
				output_profile = .false.
				find_lburn = .true.
				do_pathological = do_pathological_init
				call znd_integrate(t_start, t_end, num_steps, output_profile, profile_io, find_lburn)
						
				write(sweep_io,'(999ES25.9)') rho0, p0, &
					v_det, v_det/cs0, h_scale, &
					rho_neumann, t_neumann, p_neumann, e_neumann, burn_u, &
					sonic_loc, max_mach, q_tot, l_burn, vars(1:species)
								
        	end do !sweep over rho0
        	
        	call system_clock(sys_time_end,clock_rate)
			write(*,'(a20,ES15.4,a2)') 'Computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate, 's'
			
		!Here we will try to linearize the ZND equations around the pathological point and 
		!try to move to the frozen supersonic solution that FLASH reveals
		!else if(do_pathological) then
!        	work_isolve = 0d0
!        	iwork_isolve = 0
!			t_start = 0d0 			!Starting time?
!			t_end = burn_time		!Ending time (s)?
!			
!			vars(1:species) = xa		!Composition
!			vars(species+1) = burn_rho	!Density (g/cm^3)
!			vars(species+2) = burn_T	!Temperature (K)
!			vars(species+3) = burn_u	!Velocity in shock frame (cm/s)
!			if(do_blowout) then
!				vars(species+4) = uy_init		!Initial radial blowout velocity (cm/s)
!				if(use_uy_ux_ratio) then
!         			vars(species+4) = uy_ux_ratio*vars(species+3)
!         		endif
!				!vars(species+4) = cs_neumann 	!Initial radial velocity = sound speed
!				!grav = standard_cgrav*m_c/r_wd**2
!				!h_scale = P0/(rho0*grav)	!With gravity
!				vars(species+5) = h_scale 	!Initial radial scale height (cm)
!			else if(do_blowout_local) then
!				vars(species+4) = cs0	!Initial radial blowout velocity (cm/s)
!				vars(species+5) = h_scale 	!Initial radial scale height (cm)
!				vars(species+6) = delta_y_init	!Initial length scale for computing dP/dy (cm)
!			else
!				!Pressure in shock frame (dyne/cm^2):
!				if(num_vars.gt.species+3) vars(species+4) = burn_P
!				!Energy in shock frame (erg/g):
!				if(num_vars.gt.species+4) vars(species+5) = burn_e
!			endif
!			if(do_curvature) then
!				if(use_rc_hs_scaling) r_curve = rc_hs_factor*h_scale
!				if(.not.use_he_clavin) vars(num_vars) = r_curve
!			endif
!	 
!			!write(*,*) 'initial vars:', num_vars-species
!			!write(*,*) vars
!			!read(*,*) testchar
!			
!			sonic_loc = 0d0
!			max_mach = 0d0
!			gone_sonic = .false.
!			write(*,*) 'Entering isolve...'
!			call system_clock(sys_time_begin,clock_rate)
!			
!			do j=1,num_steps
!			!Evenly spaced timesteps in log
!			t_start = burn_time**((j-1)*1d0/num_steps)
!			t_end = burn_time**(j*1d0/num_steps)
!			!t_start = burn_time*10**(12*((j-1)*1d0/num_steps - 1))
!			!t_end = burn_time*10**(12*(j*1d0/num_steps - 1))
!			!write(*,*) 'Integrating to ',t_end
!			!write(*,*) 'iwork_isolve(2) ',iwork_isolve(2)
!			call isolve( &
!				which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
!				h, max_step_size, max_steps, &
!				rtol, atol, itol, &
!				znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
!				null_mas, imas, mlmas, mumas, &
!				znd_solout, iout, &
!				lapack_decsol, null_decsols, lrd, rpar_decsol, lid, ipar_decsol, &
!				work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
!				lrpar, rpar, lipar, ipar, &
!				lout, idid)
!			if (ierr /= 0) then
!				write(*,*) 'failed in isolve() call'
!				stop 1
!			endif
!			!Check if we've slowed down too much (good stopping criterion for curvature)
!			if (rpar(2).lt.1d-2) then
!				write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
!				write(*,*)
!				exit
!			endif
!			
!			!We also want to record the location of the maximum mach number, sonic_loc
!			!as well as the maximum mach number itself, max_mach:
!			if((rpar(2).gt.max_mach).and.(.not.gone_sonic)) then
!				sonic_loc = t_end
!				pathological_loc = t_end
!				max_mach = rpar(2)
!				if((rpar(2).ge.sonic_limit_pathological).and.(.not.gone_sonic)) then
!					write(*,*) 'ux/cs = ',rpar(2)
!					write(*,*) 'pathological point hit - linearizing ZND equations...'
!					write(*,*)
!					gone_sonic = .true.
!					k = j	!Save the time step so we can pick up again
!					exit
!				endif
!			endif
!			
!			write(profile_io, '(i30,999ES30.20)') j, t_end, vars, rpar(1), &
!				rpar(8), rpar(10), rpar(2), rpar(9), mass_flux, mom_flux, energy_flux, &
!				rpar(11), rpar(12), rpar(13), rpar(14), rpar(15), rpar(16)
!				
!				!dP_dx, dE_dx, dq_dx, detot_dx, dEtot_dx_eqn
!				!rpar(6), rpar(7), rpar(8:8+num_vars-1), rpar(8+num_vars)
!			end do !num_steps
!			
!			!Take the derivatives of all the quantities and linearize them over
!			!a certain length scale (how to choose this length?) A first choice is to
!			!use the time step we would have taken if the integration succeeded.
!			!linearization_length = 4d6	!cm
!			linearization_length = pathological_loc*pathological_linearization_ratio
!			write(*,*) 'vars pre:', vars
!			vars(1:num_vars) = vars(1:num_vars) + rpar(21:21+num_vars)*linearization_length
!			write(*,*) 'vars post:', vars
!			
!			!Now try resuming integration:
!			sonic_loc = 0d0
!			max_mach = 0d0
!			gone_sonic = .false.
!			do j=1,num_steps
!			!Evenly spaced timesteps in log
!			t_start = pathological_loc + linearization_length + &
!				(burn_time - pathological_loc + linearization_length)**((j-1)*1d0/num_steps)
!			t_end = pathological_loc + linearization_length + & 
!				(burn_time - pathological_loc + linearization_length)**(j*1d0/num_steps)
!			!write(*,*) 'Integrating to ',t_end
!			!write(*,*) 'iwork_isolve(2) ',iwork_isolve(2)
!			call isolve( &
!				which_solver, num_vars, znd_derivs, t_start, vars, t_end, &
!				h, max_step_size, max_steps, &
!				rtol, atol, itol, &
!				znd_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
!				null_mas, imas, mlmas, mumas, &
!				znd_solout, iout, &
!				lapack_decsol, null_decsols, lrd, rpar_decsol, lid, ipar_decsol, &
!				work_isolve, lwork_isolve, iwork_isolve, liwork_isolve, &
!				lrpar, rpar, lipar, ipar, &
!				lout, idid)
!			if (ierr /= 0) then
!				write(*,*) 'failed in isolve() call'
!				stop 1
!			endif
!			!Check if we've slowed down too much (good stopping criterion for curvature)
!			if (rpar(2).lt.1d-2) then
!				write(*,*) 'ux/cs dropped below 0.01 - exiting integration...'
!				write(*,*)
!				exit
!			endif
!			
!			write(profile_io, '(i30,999ES30.20)') j, t_end, vars, rpar(1), &
!				rpar(8), rpar(10), rpar(2), rpar(9), mass_flux, mom_flux, energy_flux, &
!				rpar(11), rpar(12), rpar(13), rpar(14), rpar(15), rpar(16)
!				
!				!dP_dx, dE_dx, dq_dx, detot_dx, dEtot_dx_eqn
!				!rpar(6), rpar(7), rpar(8:8+num_vars-1), rpar(8+num_vars)
!			end do !num_steps
!			
!			
!			call system_clock(sys_time_end,clock_rate)
!			
!		 
!			write(*,*) 
!			write(*,*) 'isolve complete with local subroutines', idid
!			write(*,'(a20,ES15.4,a2)') 'Computation time:', (sys_time_end-sys_time_begin)*1d0/clock_rate, 's'
!			write(*,*)
!			do j=1,species
!				write(*,'(a15,999ES25.9)') chem_isos% name(chem_id(j)), vars(j)
!			end do
!			write(*,'(a15,ES25.9)') 'rho (g/cm^3)', vars(species+1)
!			write(*,'(a15,ES25.9)') 'T (K)', vars(species+2)
!			write(*,'(a15,ES25.9)') 'u (cm/s)', vars(species+3)
!			write(*,'(a15,ES25.9)') 'u/c_s', rpar(2)
!			write(*,'(a15,ES25.9)') 'max(u/cs)', max_mach
!			write(*,'(a15,ES25.9)') 'maxloc(u/cs)', sonic_loc
!			write(*,'(a15,ES25.9)') 'q_final', rpar(9)
!			write(*,'(a15,ES25.9)') 'gamma1', rpar(10)
!	
!			write(io,'(3ES25.9,i25,999ES25.9)') t_end, rtol(1), atol(1), &
!				which_solver, vars, rpar(1)
!				
!			q_tot = rpar(9)
!        	
!        	!write(report_io,'(99ES25.9)') rho0, t0, rho_cj, t_cj, v_cj, rho_neumann, &
!        	!	t_neumann, burn_u, burn_time, vars, rpar(1)
			
        endif !Main running mode selector
        
        close(io)
        close(profile_io)
        close(sweep_io)
        close(data_in)
        !close(report_io)
        
      end subroutine do_my_burn
      
      subroutine cleanup
      	!use weak_lib, only: weak_shutdown
      	
      	implicit none
      
      	ierr = 0
      	deallocate(xa, vars, y_mass, dabar_dx, dzbar_dx, d_eps_nuc_dx, dmc_dx, rate_factors, &
			eps_nuc_categories, dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, stat=ierr)
		if (ierr /= 0) STOP "*** Deallocation error 1 ***"
        deallocate(rpar, ipar, work_isolve, iwork_isolve, work_net, &
        	ipar_decsol, rpar_decsol, stat=ierr)
        if (ierr /= 0) STOP "*** Deallocation error 2 ***"
      	deallocate(which_rates, rtol, atol, chem_id, net_iso, stat=ierr)
      	if (ierr /= 0) STOP "*** Deallocation error 3 ***"
      	
      	call lane_emden_cleanup()
      	call wd_cleanup()
      	
      	!Newton variables:
      	!deallocate(work_newton, iwork_newton, qwork_newton, stat=ierr)
      	!deallocate(rpar_decsol_newton, ipar_decsol_newton, stat=ierr)
      	
      	call free_eos_handle(eos_handle)
	 	call eos_shutdown
      	!call weak_shutdown
      end subroutine cleanup
      
end module znd


program main
	
	use znd
	
	implicit none
      	
    !If true, then use the burn subroutines defined in this module. If false, then
    !use the mesa subroutine net_1_zone_burn
    logical :: do_my_solver, intial_state
	character (len=256) :: filename, final_iso
      
    do_my_solver = .true.  	
    which_rates_choice = rates_NACRE_if_available

    ierr = 0
    
    call read_inlist
	if (ierr /= 0) then
		write(*,*) 'Error during read_inlist' 
		stop 1
	endif
    
	call initialize
	if (ierr /= 0) then
		write(*,*) 'Error during initialize' 
		stop 1
	endif
	 
	call setup_net
	if (ierr /= 0) then
		write(*,*) 'Error during setup_net' 
		stop 1
	endif
	
	!Subroutine just to test out using isolve:
	!call solve_vdpol
	!if (ierr /= 0) then
	!	write(*,*) 'Error during solve_vdpol' 
	!	stop 1
	!endif 
	
	if(do_my_solver) then
		call do_my_burn
		if (ierr /= 0) then
			write(*,*) 'Error during do_my_burn' 
			stop 1
		endif
	else 
		call do1_net_eval
		if (ierr /= 0) then
			write(*,*) 'Error during do1_net_eval' 
			stop 1
		endif
	endif 
	
	!Hugoniot constructors:
	!filename = 'initial_hugoniot-approx19.data'
	!final_iso = 'he4'
	!intial_state = .true.
	!call output_hugoniot(8,filename,intial_state,final_iso)
	!filename = 'final_hugoniot-approx19.data'
	!intial_state = .false.
	!final_iso = 'ni56'
	!call output_hugoniot(8,filename,intial_state,final_iso)
	
	call cleanup
	if (ierr /= 0) then
		write(*,*) 'Error during cleanup' 
		stop 1
	endif

end program main
