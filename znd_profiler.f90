!Module for profiling various aspects of the znd code.
!Written by: Kevin Moore - 8/6/14

module znd_test

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
	use const_def
	use rates_def
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

      ! information about the jacobian matrix
      integer :: ijac, nzmax, isparse, mljac, mujac

      ! information about the "mass" matrix
      integer :: imas, mlmas, mumas
      
      ! switch for calling the subroutine solout or nor
      integer :: iout
      
      integer :: lrd, lid
      double precision, pointer :: rpar_decsol(:) ! (lrd)
      integer, pointer :: ipar_decsol(:) ! (lid)
      
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
    character (len=100) :: my_mesa_dir, net_file
    integer, parameter :: file_path_length = 256
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
	logical :: reuse_given_rates, clip
	double precision, pointer, dimension(:) :: xa !(species)
	double precision, pointer, dimension(:) :: Y_mass !(species)
	double precision ::  eta, d_eta_dlnT, d_eta_dlnRho
	double precision ::  t_start, t_end
	double precision, pointer, dimension(:) ::  d_eps_nuc_dx !(species)
	double precision, pointer, dimension(:) ::  rate_factors !(num_reactions)
	double precision, pointer, dimension(:) ::  category_factors !(num_categories)
	double precision, pointer, dimension(:,:) ::  rate_screened !(num_rvs, num_reactions)
	double precision, pointer, dimension(:,:) ::  reaction_eps_nuc !(num_rvs, num_reactions)
	double precision, pointer , dimension(:,:)::  rate_raw !(num_rvs, num_reactions)
	double precision ::  xh, xhe
	double precision ::  abar, zbar, z2bar, ye
	double precision ::  xsum
	double precision ::  mass_correction
	double precision ::  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT
	double precision ::  theta_e_for_graboske_et_al
	double precision, pointer, dimension(:,:) ::  eps_nuc_categories !(num_rvs, num_categories)
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
			
			character(len=256) :: filename
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
      		my_mesa_dir = '/Users/Kevin/mesa'
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
			
			call eos_init(eos_file_prefix, '', '', use_cache, info)
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
			
			 call rates_init('reactions.list', '', ierr)
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
		
      subroutine cleanup
      	!use weak_lib, only: weak_shutdown
      	
      	implicit none
      
      	ierr = 0
      	!deallocate(xa, vars, y_mass, dabar_dx, dzbar_dx, d_eps_nuc_dx, dmc_dx, rate_factors, &
         !	category_factors, rate_screened, reaction_eps_nuc, rate_raw, &
			!   eps_nuc_categories, dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, stat=ierr)
		   !if (ierr /= 0) STOP "*** Deallocation error 1 ***"
         !deallocate(rpar, ipar, work_isolve, iwork_isolve, &
        	!ipar_decsol, rpar_decsol, stat=ierr)
         !if (ierr /= 0) STOP "*** Deallocation error 2 ***"
      	!deallocate(which_rates, rtol, atol, work_net, chem_id, net_iso, stat=ierr)
      	!if (ierr /= 0) STOP "*** Deallocation error 3 ***"
      	
      	call lane_emden_cleanup()
      	call wd_cleanup()
      	
      	!Newton variables:
      	!deallocate(work_newton, iwork_newton, qwork_newton, stat=ierr)
      	!deallocate(rpar_decsol_newton, ipar_decsol_newton, stat=ierr)
      	
      	call free_eos_handle(eos_handle)
	 	   call eos_shutdown
      	!call weak_shutdown
      end subroutine cleanup
		
end module znd_test


program main
	
	use znd_test
	
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
	 
	!call setup_net
	!if (ierr /= 0) then
	!	write(*,*) 'Error during setup_net' 
	!	stop 1
   !endif
	
	!Subroutine just to test out using isolve:
	!call solve_vdpol
	!if (ierr /= 0) then
	!	write(*,*) 'Error during solve_vdpol' 
	!	stop 1
	!endif 
	
   !Don't call the solvers yet (still profiling):
	!if(do_my_solver) then
	!	call do_my_burn
	!	if (ierr /= 0) then
	!		write(*,*) 'Error during do_my_burn' 
	!		stop 1
	!	endif
	!else 
	!	call do1_net_eval
	!	if (ierr /= 0) then
	!		write(*,*) 'Error during do1_net_eval' 
	!		stop 1
	!	endif
	!endif 
	
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