!Module to integrate the stellar structure equations using the MESA EOS and assuming 
!an isothermal WD with the possibility of an adiabatic envelope of differing composition.
!Originally designed to work as a module in the ZND code to provide more realistic
!mappings from (M_c, M_env) to (rho_b, H_p) when calculating detonation structure.

!Written by Kevin Moore on 5/7/13
!Newton solver implemented 5/10/13, improved by 5/13/13

module wd_solver

	use const_def
	use const_lib
	use chem_def
	use chem_lib
	use eos_def
	use eos_lib
	use num_def
	use num_lib
	use mtx_def
	use mtx_lib
	
	!I want this stuff to be completely self-contained, so all the variables should
	!be inaccessible to an outside program. This lets the outside program have public
	!variables with the same name, cutting down on confusion (eg. having to name this
	!work array "work_wd" to avoid conflicting with the "work" array in my ZND code. 
	PRIVATE	
	
	integer, parameter :: nv_wd = 4  ! the number of variables in the white dwarf system of ODEs
	real(dp) :: rtol_wd(1) ! relative error tolerance(s)
	real(dp) :: atol_wd(1) ! absolute error tolerance(s)
	real(dp) :: x_wd ! starting value for the interval of integration
	real(dp) :: xend_wd ! ending value for the interval of integration
	real(dp) :: expect(nv_wd), yprime(nv_wd)
	character (len=64) :: str
	character (len=256) :: dir, fname
	integer, parameter :: lrpar_wd = 3, lipar_wd = 1
	real(dp) :: max_abs_yp2, h_wd, max_step_size_wd
	integer :: io_unit, i, lout_wd, iout_wd, idid_wd, itol_wd, j
	integer :: liwork_wd, lwork_wd, max_steps_wd, ierr_wd
	real(dp), pointer :: work_wd(:)
	integer, pointer :: iwork_wd(:)
	real(dp), target :: y_ary(nv_wd)
	real(dp), pointer :: y_wd(:)
	real(dp), target :: rpar_ary(lrpar_wd)
	integer, target :: ipar_ary(lipar_wd)
	integer, pointer :: ipar_wd(:) ! (lipar)
	real(dp), pointer :: rpar_wd(:) ! (lrpar)
	
	integer, parameter :: file_path_length = 256
	
	!Composition variables:
	integer, pointer :: chem_id(:), net_iso(:)
	double precision, pointer :: xa_wd(:), xa_env(:)
	integer, parameter :: max_num_isos_for_Xinit = 100
	integer :: num_isos_for_Xinit, species_wd
	character(len=iso_name_length) :: names_of_isos_for_Xinit(max_num_isos_for_Xinit)
    double precision :: values_for_Xinit(max_num_isos_for_Xinit)
	double precision ::  xh, xhe, zm, abar, zbar, z2bar, ye, xsum, mass_correction
	 
	!EOS variables:
    !Common named variables:
	double precision, pointer :: dabar_dx(:), dzbar_dx(:), dmc_dx(:)
	integer :: info
	integer :: eos_handle 	!ID number for calling the EOS module
	character (len=file_path_length) :: eos_file_prefix		!Which EOS tables to use
	logical :: use_cache
	double precision :: Pgas, logPgas, logRho, logT, t0, p0, e0, p, e, rho, t, t_env
	double precision :: dlnRho_dlnPgas_const_T
    double precision :: dlnRho_dlnT_const_Pgas
	double precision, dimension(num_eos_basic_results) :: res  	!Basic EOS result array
	double precision :: d_dlnRho_const_T(num_eos_basic_results) !Basic EOS derivs array
    double precision :: d_dlnT_const_Rho(num_eos_basic_results) !Basic EOS derivs array
    double precision :: d_dabar_const_TRho(num_eos_basic_results)	!Basic EOS result array
    double precision :: d_dzbar_const_TRho(num_eos_basic_results)	!Basic EOS result array
    
    !Newton solver variables:
    real(dp), parameter :: one=1
	integer, parameter :: nz = 1, nvar = 2 !use odd number of zones for problem symmetry
	integer, parameter :: nsec = 0 ! number of secondaries per zone
	integer, parameter :: ldy = nz 

	integer, parameter :: i_conc=1, i_flux=2, equ_conc=1, equ_flux=2

	logical :: do_numerical_jacobian
	integer :: neq, which_decsol
	integer :: matrix_type
	real(dp) :: xp0, xp1 	!Primaries
	real(dp), pointer, dimension(:) :: equ1, x1, xold1, dx1, xscale1, y1
	real(dp), pointer, dimension(:,:) :: equ, x, xold, dx, xscale, y
	real(dp), pointer, dimension(:,:,:) :: ublk, dblk, lblk
	
	!Turn on various debugging output
	logical, parameter :: dbg = .true.
	!if true, P is our independent variable for integration (otherwise it's r):
	logical, parameter :: indep_var_P = .true.
    
    !Specify which subroutines will be callable by outside prorams:
    public :: init_mesa_modules, wd_init, wd_integrate, wd_envelope, find_envelope, &
    	set_initial_conditions, find_rhoc, find_envelope_newton, wd_cleanup, &
    	find_wd_params

	contains
	
		!Initializes all the variables used in the newton solver
		subroutine newton_init()
			implicit none
			
			character (len=64) :: decsol_option_name
			integer :: ierr
			
			do_numerical_jacobian = .true.   
      	which_decsol = lapack
      		
			call decsol_option_str(which_decsol, decsol_option_name, ierr)
			if (ierr /= 0) return
         	
         write(*,*) 'Newton solver using ' // trim(decsol_option_name)
         write(*,*)
         	
         neq = nvar*nz
			allocate(equ1(neq), x1(neq), xold1(neq), dx1(neq), &
				xscale1(neq), y1(ldy*nsec), stat=ierr)
			if (ierr /= 0) stop 1
         
			x(1:nvar,1:nz) => x1(1:neq)
      	xold(1:nvar,1:nz) => xold1(1:neq)
      	dx(1:nvar,1:nz) => dx1(1:neq)
      	equ(1:nvar,1:nz) => equ1(1:neq)
      	xscale(1:nvar,1:nz) => xscale1(1:neq)
      	y(1:ldy,1:nsec) => y1(1:ldy*nsec)
		end subroutine newton_init
		
		subroutine init_mesa_modules()
			implicit none
			
			write(*,*) 'Initializing WD builder'
			write(*,*)
		
			!mesa_dir = '/Users/Kevin/mesa'
			!mesa_dir = '/Users/Kevin/mesa_5271'
			!mesa_dir = '/Users/Kevin/mesa_5819'
			mesa_dir = ''
			call const_init(mesa_dir, ierr_wd)

			!EOS options:
			eos_file_prefix = 'mesa'
			use_cache = .true.
         species_wd = 3

			call chem_init('isotopes.data', ierr_wd)
			if (ierr_wd /= 0) then
			   write(*,*) 'chem_init failed'
			   return
			end if
         
			write(*,*) 'num_chem_isos', num_chem_isos
			write(*,*) 'species_wd', species_wd

			allocate(net_iso(num_chem_isos), chem_id(species_wd), dabar_dx(species_wd), &
				dzbar_dx(species_wd), dmc_dx(species_wd), stat=ierr_wd)
				if (ierr_wd /= 0) stop 'allocate failed'

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

			call newton_init()
		end subroutine init_mesa_modules
      
      subroutine wd_init()
			implicit none  	 

			ipar_wd => ipar_ary
			rpar_wd => rpar_ary

			y_wd => y_ary

			t0 = 1d7		!Constant core temperature (K)
			t_env = 1d8		!Envelope temperature typically higher than T_c

			!Next, set up the final composition and figure out the total q released in
			!a detonation that burns all the way to this composition
			allocate(xa_wd(species_wd),xa_env(species_wd))
			xa_wd = 0d0
			xa_env = 0d0

			!Set up possible isotopes
			chem_id(:) = 0; net_iso(:) = 0
			chem_id(1) = ihe4; net_iso(ihe4) = 1
			chem_id(2) = ic12; net_iso(ic12) = 2
			chem_id(3) = io16; net_iso(io16) = 3

			!Core composition:
			num_isos_for_Xinit = 2
			values_for_Xinit = 0d0
			names_of_isos_for_Xinit(1:num_isos_for_Xinit) = (/'c12', 'o16'/)
			values_for_Xinit(1:num_isos_for_Xinit) = (/0.5, 0.5/)
			values_for_Xinit = values_for_Xinit/sum(values_for_Xinit)	!Normalize
			write(*,*) 'WD core composition'
			do i=1,num_isos_for_Xinit
				j = get_nuclide_index(names_of_isos_for_Xinit(i))
				write(*,*) 'Index in chem_isos: ',j
				write(*,*) names_of_isos_for_Xinit(i), values_for_Xinit(i)
				write(*,*) 'check index (name should match chem_isos result above): ',&
					chem_isos% name(j)

				xa_wd(net_iso(j)) = values_for_Xinit(i)
			end do
			write(*,*)
			
			!Envelope composition:
			num_isos_for_Xinit = 2
			values_for_Xinit = 0d0
			names_of_isos_for_Xinit(1:num_isos_for_Xinit) = (/'he4', 'c12'/)
			values_for_Xinit(1:num_isos_for_Xinit) = (/1.0, 0.0/)
			values_for_Xinit = values_for_Xinit/sum(values_for_Xinit)	!Normalize
			write(*,*) 'WD envelope composition'
			do i=1,num_isos_for_Xinit
				j = get_nuclide_index(names_of_isos_for_Xinit(i))
				write(*,*) 'Index in chem_isos: ',j
				write(*,*) names_of_isos_for_Xinit(i), values_for_Xinit(i)
				write(*,*) 'check index (name should match chem_isos result above): ',&
					chem_isos% name(j)

				xa_env(net_iso(j)) = values_for_Xinit(i)
			end do
			write(*,*)

			call composition_info(species_wd, chem_id, xa_wd, xh, xhe, zm, abar, zbar, z2bar, &
				ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
			!write(*,*) abar, zbar, zm

			lout_wd = 6
			max_steps_wd = 10000
			max_step_size_wd = 0

			itol_wd = 0 ! scalar tolerances
			iout_wd = 0 ! no intermediate output

			rtol_wd(1) = 1d-6
			atol_wd(1) = 0d0
			h_wd = 1d0

			!vdp:
			!rpar_wd(1) = eps

			!WD structure equations:
			rpar_wd = 0
			ipar_wd = 0

			call cash_karp_work_sizes(nv_wd,liwork_wd,lwork_wd)
			allocate(work_wd(lwork_wd), iwork_wd(liwork_wd))

			iwork_wd = 0
			work_wd = 0
		end subroutine wd_init
		
		subroutine set_initial_conditions(rho_c)
			implicit none
			
			double precision, intent(in) :: rho_c
			double precision :: x_wd
			
			x_wd = 1d0
			y_wd(1) = rho_c							!rho(0) in g/cm^3
			y_wd(2) = 4*pi*y_wd(1)*x_wd**3/3d0		!M(0) in g
			y_wd(3) = 0d0							!t_sound in s
			y_wd(4) = x_wd							!r(0), if using P as indep variable
			
			call composition_info(species_wd, chem_id, xa_wd, xh, xhe, zm, abar, zbar, z2bar, &
				ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
			!write(*,*) abar, zbar, zm

			rho = y_wd(1)
			logRho = log10(rho)
			T = T0
			logT = log10(T)
			call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species_wd, chem_id, net_iso, xa_wd, &
				rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
				d_dabar_const_TRho, d_dzbar_const_TRho, ierr_wd)
			e0 = exp(res(i_lnE))
			p0 = exp(res(i_lnPgas)) + 1/3d0*crad*t**4
			
			!rpar_wd = 0
			!ipar_wd = 0
			
			!iwork_wd = 0
			!work_wd = 0
		end subroutine set_initial_conditions
      
		subroutine wd_derivs(n,x,h,y,f,lrpar,rpar,lipar,ipar,ierr)
			integer, intent(in) :: n, lrpar, lipar
			real(dp), intent(in) :: x, h
			real(dp), intent(inout) :: y(:)
			real(dp), intent(out) :: f(:)
			integer, intent(inout), pointer :: ipar(:) ! (lipar)
			real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
			integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.

			double precision :: rho, T, e, p, dP_drho, gamma1, grav

			ierr = 0

			!if(y(1).lt.0.0) then
			!	y(1) = 0.0
			!endif

			rho = y(1)
			logRho = log10(rho)
			T = T0
			logT = log10(T)
			call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species_wd, chem_id, net_iso, xa_wd, &
				rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
				d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
			e = exp(res(i_lnE))
			gamma1 = res(i_gamma1)

			!write(*,*) x, y, dP_drho, p, exp(res(i_lnPgas))

			if(indep_var_P) then
				grav = standard_cgrav*y(2)/y(4)**2
				p = x     
				dP_drho = res(i_chiRho)*p/rho
				
				f(1) = 1d0/(dP_drho)							!f(1) = drho/dP
				f(2) = -4*pi*y(4)**2/grav						!f(2) = dM/dP   
				f(3) = -2*(gamma1*p/rho)**(-0.5)/(rho*grav)		!f(3) = d(t_sound)/dP
				f(4) = -1d0/(rho*grav)							!f(4) = dr/dP
			else
				grav = standard_cgrav*y(2)/x**2
				p = exp(res(i_lnPgas)) + 1/3d0*crad*t**4     
				dP_drho = res(i_chiRho)*p/rho
				
				f(1) = -rho*grav/dP_drho						!f(1) = drho/dr
				f(2) = 4*pi*x**2*rho							!f(2) = dM/dr   
				f(3) = 2*(gamma1*p/rho)**(-0.5)					!f(3) = d(t_sound)/dr
				f(4) = 0d0										!unused here
			endif

			!write(*,'(99ES20.10)') y(1), y(2), y(3), y(4)

		end subroutine wd_derivs
		
		subroutine envelope_derivs(n,x,h,y,f,lrpar,rpar,lipar,ipar,ierr)
			implicit none
		
			integer, intent(in) :: n, lrpar, lipar
			real(dp), intent(in) :: x, h
			real(dp), intent(inout) :: y(:)
			real(dp), intent(out) :: f(:)
			integer, intent(inout), pointer :: ipar(:) ! (lipar)
			real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
			integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.

			double precision :: rho, T, e, p, dP_drho, rho_b, gamma1, grav

			ierr = 0

			!if(y(1).lt.0.0) then
			!	y(1) = 0.0
			!endif
			
			rho_b = rpar(1)
			gamma1 = 5d0/3d0
						
			rho = y(1)
			logRho = log10(rho)
			T = T_env*(rho/rho_b)**(gamma1-1)
			!T = T0
			logT = log10(T)
			call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species_wd, chem_id, net_iso, xa_env, &
				rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
				d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
			e = exp(res(i_lnE))
			p = exp(res(i_lnPgas)) + 1/3d0*crad*t**4     
			!dP_drho = res(i_chiRho)*p/rho

			!write(*,*) x, y, dP_drho, p, exp(res(i_lnPgas))

			if(indep_var_P) then
				grav = standard_cgrav*y(2)/y(4)**2
				p = x
				dP_drho = res(i_chiRho)*p/rho
				
				f(1) = 1d0/(dP_drho + &							!f(1) = drho/dP
					(gamma1-1)*T_env/rho_b*(rho/rho_b)**(gamma1-2))	
				!f(1) = 1d0/dP_drho	
				f(2) = -4*pi*y(4)**2/grav						!f(2) = dM/dP   
				f(3) = 0d0										!f(3) = d(t_sound)/dP
				f(4) = -1d0/(rho*grav)							!f(4) = dr/dP
			else
				grav = standard_cgrav*y(2)/x**2
				p = exp(res(i_lnPgas)) + 1/3d0*crad*t**4    
				dP_drho = res(i_chiRho)*p/rho
			
				f(1) = -rho*grav/&									!f(1) = drho/dr
					(dP_drho + (gamma1-1)*T_env/rho_b*(rho/rho_b)**(gamma1-2))		
				f(2) = 4*pi*x**2*rho								!f(2) = dM/dr 
				f(3) = 0d0			  								!f(3) = d(t_sound)/dr
				f(4) = 0d0											!unused here
			endif

			!write(*,'(99ES20.10)') f(1), f(2), f(3), f(4)

		end subroutine envelope_derivs

		subroutine wd_integrate(m_wd)
			implicit none
			
			double precision, intent(inout) :: m_wd

			integer :: i, num_steps

			ierr_wd = 0
			num_steps = 100

			if(dbg) then
				write(*,*) 'Starting integration'
				write(*,*)
			endif
			
			do i=1,num_steps
				if(indep_var_P) then
					x_wd = p0*10**(-8d0*(i-1)/num_steps)
					xend_wd = p0*10**(-8d0*i/num_steps)
				else
					x_wd = 1d0 + 10**(5d0 + 5d0*(i-1)/num_steps)
					xend_wd = 1d0 + 10**(5d0 + 5d0*i/num_steps)
				endif
				
				!write(*,'(99ES20.10)') x_wd, xend_wd
				!write(*,'(99ES20.10)') y_wd(1), y_wd(2), y_wd(3), y_wd(4)
				call cash_karp( &
				   nv_wd,wd_derivs,x_wd,y_wd,xend_wd, &
				   h_wd,max_step_size_wd,max_steps_wd, &
				   rtol_wd,atol_wd,itol_wd, &
				   null_solout,iout_wd,work_wd,lwork_wd,iwork_wd,liwork_wd, &
				   lrpar_wd,rpar_wd,lipar_wd,ipar_wd,lout_wd,idid_wd)
				!write(*,'(99ES20.10)') y_wd(1), y_wd(2), y_wd(3), y_wd(4)

				if (idid_wd /= 1) then ! trouble
					write(*,*) 'idid', idid_wd
					stop 1
				end if

				rho = y_wd(1)
				logRho = log10(rho)
				T = T0
				logT = log10(T)
				call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species_wd, chem_id, net_iso, xa_wd, &
					rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
					d_dabar_const_TRho, d_dzbar_const_TRho, ierr_wd)
				e = exp(res(i_lnE))
				p = exp(res(i_lnPgas)) + 1/3d0*crad*t**4   

				if((p/p0.le.1d-8).or.((i.eq.num_steps).and.(indep_var_P))) then
					if(dbg) then
						write(*,*) 'Pressure has dropped to 10^(-6)*P_c, stop integration'
						if(indep_var_P) then
							write(*,'(a15,ES25.10)') 'R_wd (rsun)', y_wd(4)/rsun
							write(*,'(a15,ES25.10)') 'R_wd (cm)', y_wd(4)
							write(*,'(a15,ES25.10)') 'rho (g/cc)', y_wd(1)
							write(*,'(a15,ES25.10)') 'M_wd (msun)', y_wd(2)/msun
							write(*,'(a15,ES25.10)') 't_sound', y_wd(3)
						else
							write(*,'(a15,ES25.10)') 'R_wd (rsun)', x_wd/rsun
							write(*,'(a15,ES25.10)') 'rho (g/cc)', y_wd(1)
							write(*,'(a15,ES25.10)') 'M_wd (msun)', y_wd(2)/msun
							write(*,'(a15,ES25.10)') 't_sound', y_wd(3)
						endif
					endif
					
					m_wd = y_wd(2)/msun
					exit
				endif

				!write(*,*) x_wd, p, p0
			end do

      	end subroutine wd_integrate
      	
      	subroutine wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
			implicit none
			
			double precision, intent(in) :: rho_c, pb_frac
			double precision, intent(inout) :: m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound

			integer :: i, num_steps

			ierr_wd = 0
			num_steps = 100
			
			if(dbg) then
				write(*,*) 'Starting core integration'
				write(*,*) rho_c, pb_frac
			endif
			
			do i=1,num_steps
				if(indep_var_P) then
					x_wd = p0*10**(log10(pb_frac)*(i-1)/num_steps)
					xend_wd = p0*10**(log10(pb_frac)*i/num_steps)
				else
					x_wd = 1d0 + 10**(5d0 + 5d0*(i-1)/num_steps)
					xend_wd = 1d0 + 10**(5d0 + 5d0*i/num_steps)
				endif
				call cash_karp( &
				   nv_wd,wd_derivs,x_wd,y_wd,xend_wd, &
				   h_wd,max_step_size_wd,max_steps_wd, &
				   rtol_wd,atol_wd,itol_wd, &
				   null_solout,iout_wd,work_wd,lwork_wd,iwork_wd,liwork_wd, &
				   lrpar_wd,rpar_wd,lipar_wd,ipar_wd,lout_wd,idid_wd)

				if (idid_wd /= 1) then ! trouble
					write(*,*) 'idid', idid_wd
					stop 1
				end if

				rho = y_wd(1)
				logRho = log10(rho)
				T = T0
				logT = log10(T)
				call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species_wd, chem_id, net_iso, xa_wd, &
					rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
					d_dabar_const_TRho, d_dzbar_const_TRho, ierr_wd)
				e = exp(res(i_lnE))
				Pgas = exp(res(i_lnPgas))
				p = exp(res(i_lnPgas)) + 1/3d0*crad*t**4   

				if((p/p0.le.pb_frac).or.((i.eq.num_steps).and.(indep_var_P))) then
					if(dbg) then
						write(*,*) 'Pressure has dropped to 10^(-2)*P_c, stop core integration'
						if(indep_var_P) then
							write(*,'(a15,ES25.10)') 'R_wd (rsun)', y_wd(4)/rsun
							write(*,'(a15,ES25.10)') 'rho (g/cc)', y_wd(1)
							write(*,'(a15,ES25.10)') 'M_wd (msun)', y_wd(2)/msun
							write(*,'(a15,ES25.10)') 't_sound', y_wd(3)
							write(*,'(a15,ES25.10)') 'P', x_wd
							write(*,'(a15,ES25.10)') 'P_eos', p
						else
							write(*,'(a15,ES25.10)') 'R_wd (rsun)', x_wd/rsun
							write(*,'(a15,ES25.10)') 'rho (g/cc)', y_wd(1)
							write(*,'(a15,ES25.10)') 'M_wd (msun)', y_wd(2)/msun
							write(*,'(a15,ES25.10)') 't_sound', y_wd(3)
						endif
					endif
					t_sound = y_wd(3)
					m_c = y_wd(2)/msun
					if(indep_var_P) then
						r_c = y_wd(4)
						p_b = x_wd
					else
						r_c = x_wd
						p_b = p
					endif
					exit
				endif
			end do
			
			!Find base density in envelope (rho_b):
			call composition_info(species_wd, chem_id, xa_env, xh, xhe, zm, abar, zbar, z2bar, &
				ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
			logPgas = log10(Pgas)
			T = t_env	!Envelope temperature higher than T_c typically
			logT = log10(T)
			call eosPT_get(eos_handle, Zm, Xh, abar, zbar, species_wd, chem_id, net_iso, &
				xa_env, Pgas, logPgas, T, logT, rho_b, logRho, dlnRho_dlnPgas_const_T, &
				dlnRho_dlnT_const_Pgas, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
				d_dabar_const_TRho, d_dzbar_const_TRho, ierr_wd)
			if(dbg) then
				write(*,*) 'Base density in core:', y_wd(1)
				write(*,*) 'Base density in envelope:', rho_b
				write(*,*)
			endif
			
			y_wd(1) = rho_b
			rpar_wd(1) = rho_b
			
			!Now resume integration with a isentropic envelope:
			if(dbg) then
				write(*,*) 'Starting envelope integration'
				write(*,*) y_wd
			endif
			
			do i=1,num_steps
				if(indep_var_P) then
					x_wd = p0*10**(log10(pb_frac) - (log10(pb_frac) + 8d0)*(i-1)/num_steps)
					xend_wd = p0*10**(log10(pb_frac) - (log10(pb_frac) + 8d0)*i/num_steps)
				else
					x_wd = r_c + 10**(5d0 + 5.0*(i-1)/num_steps) - 1d0
					xend_wd = r_c + 10**(5d0 + 5.0*i/num_steps) - 1d0
				endif
				
				call cash_karp( &
				   nv_wd,envelope_derivs,x_wd,y_wd,xend_wd, &
				   h_wd,max_step_size_wd,max_steps_wd, &
				   rtol_wd,atol_wd,itol_wd, &
				   null_solout,iout_wd,work_wd,lwork_wd,iwork_wd,liwork_wd, &
				   lrpar_wd,rpar_wd,lipar_wd,ipar_wd,lout_wd,idid_wd)

				if (idid_wd /= 1) then ! trouble
					write(*,*) 'idid', idid_wd
					stop 1
				end if

				rho = y_wd(1)
				logRho = log10(rho)
				T = T_env*(rho/rho_b)**(2d0/3)
				logT = log10(T)
				call eosDT_get(eos_handle, Zm, Xh, abar, zbar, species_wd, chem_id, net_iso, xa_env, &
					rho, logRho, T, logT, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
					d_dabar_const_TRho, d_dzbar_const_TRho, ierr_wd)
				e = exp(res(i_lnE))
				p = exp(res(i_lnPgas)) + 1/3d0*crad*t**4   

				if((p/p0.le.1d-8).or.((i.eq.num_steps).and.(indep_var_P))) then
					if(dbg) then
						write(*,*) 'Pressure has dropped to 10^(-6)*P_c, stop integration'
						write(*,'(a15,ES25.10)') 'R_wd (rsun)', r_wd
						write(*,'(a15,ES25.10)') 'rho (g/cc)', y_wd(1)
						write(*,'(a15,ES25.10)') 'M_c (msun)', m_c
						write(*,'(a15,ES25.10)') 'M_env (msun)', (y_wd(2)/msun - m_c)
						write(*,'(a15,ES25.10)') 'M_wd (msun)', y_wd(2)/msun
						write(*,'(a15,ES25.10)') 't_sound', y_wd(3)
						write(*,'(a15,ES25.10)') 'min v_det', pi*r_c/y_wd(3)
						write(*,*)
					endif
					exit
				endif
				
				!write(*,*) x_wd, p, p0
			end do
			
			m_env = y_wd(2)/msun - m_c
			if(indep_var_P) then
				r_wd = y_wd(4)
			else
				r_wd = x_wd
			endif
			
      	end subroutine wd_envelope
      	
		double precision function find_rhoc_func(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
			implicit none 
			! returns with ierr = 0 if was able to evaluate f and df/dx at x
			! if df/dx not available, it is okay to set it to 0
			integer, intent(in) :: lrpar, lipar
			real(dp), intent(in) :: x
			real(dp), intent(out) :: dfdx
			integer, intent(inout), pointer :: ipar(:) ! (lipar)
			real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
			integer, intent(out) :: ierr
			
			double precision :: m_wd, target_mwd
			
			ierr = 0
			target_mwd = rpar(1)
			
			call set_initial_conditions(x)
			call wd_integrate(m_wd)
			
			!write(*,*) 'm_wd, target_mwd', m_wd, target_mwd
			
			find_rhoc_func = m_wd - target_mwd
			dfdx = 0
			
			!write(*,*) 'find_rhoc_func=',find_rhoc_func
		end function find_rhoc_func
      	
      	!Solves for the central density (rho_c) given a total WD mass (m_wd).
      	subroutine find_rhoc(target_mwd, rhoc_min, rhoc_max, rho_c)
      		implicit none
      		
      		double precision, intent(in) :: target_mwd, rhoc_min, rhoc_max
      		double precision, intent(inout) :: rho_c
      		
      		integer, parameter :: lrpar = 1, lipar = 1
      		integer, pointer :: ipar(:)
      		double precision, pointer :: rpar(:)
      		double precision :: x1, x3, y1, y3, dfdx, epsx, epsy
      		integer :: imax, ierr
      		
      		imax = 100
      		epsx = 1d3
      		epsy = 1d-4
      		x1 = rhoc_min
      		x3 = rhoc_max
      		
      		allocate(rpar(lrpar), ipar(lipar))
      		
      		!Target WD mass:
      		rpar(1) = target_mwd
      		
      		y1 = find_rhoc_func(x1, dfdx, lrpar, rpar, lipar, ipar, ierr)
      		y3 = find_rhoc_func(x3, dfdx, lrpar, rpar, lipar, ipar, ierr)
      		
      		!write(*,*) 'x1, x3', x1, x3
      		!write(*,*) 'y1, y3', y1, y3
      		
      		if(y1.gt.0d0) then
      			rho_c = x1
      		elseif(y3.lt.0d0) then
      			rho_c = x3
      		else
      			rho_c = safe_root(find_rhoc_func, x1, x3, y1, y3, imax, &
      				epsx, epsy, lrpar, rpar, lipar, ipar, ierr)
      		endif
      			
      		write(*,*) 'WD solution, rho_c = ', rho_c
      		
      		deallocate(rpar,ipar)
      		
      	end subroutine find_rhoc
      	
      	!Given a core mass (M_c) and an envelope mass (M_env), find the integration
      	!parameters rho_c and pb_frac that create a WD with that structure.
      	subroutine find_envelope(mc_target, menv_target)
      		implicit none
      		
      		double precision, intent(in) :: mc_target, menv_target
      		double precision :: rho_c_min, rho_c_max, pb_frac_min, pb_frac_max
      		double precision :: rho_c_upper, rho_c_lower, pb_frac_upper, pb_frac_lower
      		double precision :: mc_tol, menv_tol, pb_frac_tol, rho_c_tol
      		double precision :: rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound
      		integer :: i,j
      		
      		rho_c_min = 1d5
      		rho_c_max = 1d10
      		pb_frac_min = 1d-5
      		pb_frac_max = 1d-1
      		
      		!Tolerances on solutions:
      		mc_tol = 1d-1
      		menv_tol = 1d-1
      		
      		!Tolerances on searches:
      		pb_frac_tol = 1d-1
      		rho_c_tol = 1d-1
      		
      		rho_c_upper = rho_c_max
      		rho_c_lower = rho_c_min
      		pb_frac_upper = pb_frac_max
      		pb_frac_lower = pb_frac_min
      		
      		do i=1,100
      		
      			!First check whether the solution even lies in this pb_frac bracketing:
      			x_wd = 1d0
				rho_c = (rho_c_upper + rho_c_lower)/2d0
				y_wd(1) = rho_c							!rho(0) in g/cm^3
				y_wd(2) = 4*pi*y_wd(1)*x_wd**3/3d0		!M(0) in g
				y_wd(3) = 0d0							!t_sound(0) in s
				t0 = 1d7	!Constant temperature (K)
				pb_frac = pb_frac_upper
			
				call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
				write(*,'(a15,ES25.10)') 'm_c (msun)', m_c
				write(*,'(a15,ES25.10)') 'm_env (msun)', m_env
				write(*,*)
				
				!If pb_frac_upper is a true bracket, then the core mass should be too
				!small, and the envelope mass should be too big
				if ((m_c.gt.mc_target).or.(m_env.lt.menv_target)) then
					rho_c_upper = rho_c
					continue
				endif
				
				x_wd = 1d0
				rho_c = (rho_c_upper + rho_c_lower)/2d0
				y_wd(1) = rho_c							!rho(0) in g/cm^3
				y_wd(2) = 4*pi*y_wd(1)*x_wd**3/3d0		!M(0) in g
				y_wd(3) = 0d0							!t_sound(0) in s
				t0 = 1d7	!Constant temperature (K)
				pb_frac = pb_frac_lower
			
				call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
				write(*,'(a15,ES25.10)') 'm_c (msun)', m_c
				write(*,'(a15,ES25.10)') 'm_env (msun)', m_env
				write(*,*)
				
				!If pb_frac_lower is a true bracket, then the core mass should be too
				!big, and the envelope mass should be too small
				if ((m_c.lt.mc_target).or.(m_env.gt.menv_target)) then
					rho_c_lower = rho_c
					continue
				endif
      			
      			do j=1,100
      			
      				write(*,*) 'Brackets are:'
      				write(*,'(a20,ES25.10)') 'rho_c_lower (g/cc)', rho_c_lower
      				write(*,'(a20,ES25.10)') 'rho_c_upper (g/cc)', rho_c_upper
					write(*,'(a20,ES25.10)') 'pb_frac_lower', pb_frac_lower
					write(*,'(a20,ES25.10)') 'pb_frac_upper', pb_frac_upper
					write(*,*)
      			
      				x_wd = 1d0
      				rho_c = (rho_c_upper + rho_c_lower)/2d0
					y_wd(1) = rho_c							!rho(0) in g/cm^3
					y_wd(2) = 4*pi*y_wd(1)*x_wd**3/3d0		!M(0) in g
					y_wd(3) = 0d0							!t_sound(0) in s
					t0 = 1d7	!Constant temperature (K)
      				pb_frac = (pb_frac_upper + pb_frac_lower)/2d0
      				
      				!Try our integral solution
      				call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
      				write(*,'(a15,ES25.10)') 'm_c (msun)', m_c
					write(*,'(a15,ES25.10)') 'm_env (msun)', m_env
					write(*,*)
      				
      				!Check if we've hit the solution:
      				if((abs(m_c - mc_target)/mc_target.le.mc_tol).and.&
      					(abs(m_env - menv_target)/menv_target.le.menv_tol)) then
      					write(*,*) 'Solution found:'
						write(*,'(a15,ES25.10)') 'rho_c (g/cc)', rho_c
						write(*,'(a15,ES25.10)') 'pb_frac', pb_frac
						write(*,'(a15,ES25.10)') 'rho_b', rho_b
						write(*,'(a15,ES25.10)') 'P_b', p_b
						write(*,'(a15,ES25.10)') 'R_c (cm)', r_c
						write(*,'(a15,ES25.10)') 'g_b', standard_cgrav*m_c*msun/r_c**2
						write(*,'(a15,ES25.10)') 't_sound', y_wd(3)
						write(*,*)
						write(*,*) 'Comparisons to thin-shell limit:'
						write(*,'(a15,ES25.10)') 'Pb_thin', &
							m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4)
						write(*,'(a15,ES25.10)') 'Delta P_b', (p_b - &
							m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4))/p_b
						exit
					endif
      				
      				!We now have four cases, adjust the brackets accordingly:
      				if((m_c.lt.mc_target).and.(m_env.lt.menv_target)) then
      					rho_c_lower = rho_c
      					exit
      				elseif ((m_c.lt.mc_target).and.(m_env.gt.menv_target)) then
      					pb_frac_upper = pb_frac
      				elseif ((m_c.gt.mc_target).and.(m_env.lt.menv_target)) then
      					pb_frac_lower = pb_frac
      				elseif ((m_c.gt.mc_target).and.(m_env.gt.menv_target)) then
      					rho_c_upper = rho_c
      					exit
      				endif
      				
      				!Check if our brackets are too close:
      				if(abs(pb_frac_lower - pb_frac_upper)/pb_frac.le.pb_frac_tol) then
      					if((m_c.lt.mc_target).and.(m_env.gt.menv_target)) then
      						rho_c_lower = rho_c
      					elseif((m_c.gt.mc_target).and.(m_env.lt.menv_target)) then
      						rho_c_upper = rho_c
      					endif
      					
      					exit
      				endif
      				
      			end do !inner loop (j)
      			
      			!Check if we've hit the solution (need to check in outer loop too):
				if((abs(m_c - mc_target)/mc_target.le.mc_tol).and.&
					(abs(m_env - menv_target)/menv_target.le.menv_tol)) then
					exit
				endif
				
				!Reset bounds on pb_frac (since rho_c has changed)
				pb_frac_lower = pb_frac_min
      			pb_frac_upper = pb_frac_max
	
      		end do !outer loop (i)
      		
      	end subroutine find_envelope
      	
      	!Given a core mass (M_c) and an envelope mass (M_env), find the integration
      	!parameters rho_c and pb_frac that create a WD with that structure.
      	subroutine find_envelope_newton(mc_target, menv_target, rhoc_guess, pbfrac_guess, &
      		rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
      		implicit none
      		
      		double precision, intent(in) :: mc_target, menv_target, rhoc_guess, pbfrac_guess
      		double precision, intent(inout) :: rho_c, pb_frac, m_c, m_env, rho_b, &
      			p_b, r_c, r_wd, t_sound
      		
   			!Begin Newton solver variables -----------------------------------------------
   			integer :: liwork, lwork!, lqwork
   			integer, dimension(:), pointer :: iwork
   			real(dp), dimension(:), pointer :: work
   			!real(qp), dimension(:), pointer :: qwork

   			integer, parameter :: lrpar = 3, lipar = 1
   			integer, target :: ipar_target(lipar)
   			real(dp), target :: rpar_target(lrpar)
   			integer, pointer :: ipar(:)
   			real(dp), pointer :: rpar(:)

   			integer :: lrd, lid
   			integer, pointer :: ipar_decsol(:) ! (lid)
   			real(dp), pointer :: rpar_decsol(:) ! (lrd)
   			real(dp), pointer :: AF1(:)

   			integer :: mljac, mujac
   			real(dp) :: tol_correction_norm, tol_max_correction, tol_residual_norm
   			logical :: nonconv
   			integer :: ierr
   			!End Newton solver variables -------------------------------------------------
			
   			double precision :: mwd_target, rhoc_min, rhoc_max
			
   			ierr = 0

   			ipar => ipar_target
   			rpar => rpar_target         

   			if (do_numerical_jacobian) then
   				ipar(1) = 1
   			else
   				ipar(1) = 0
   			end if
			
   			xold(:,1) = (/ rhoc_guess, pbfrac_guess /)
   			rpar(2) = mc_target
   			rpar(3) = menv_target
   			dx = 0d0
   			x = xold

   			mljac = 2*nvar-1 ! number of subdiagonals
   			mujac = mljac ! number of superdiagonals
   			if (which_decsol == lapack) then
   				call lapack_work_sizes(nz*nvar, lrd, lid)
   				if (do_numerical_jacobian) then
   				   matrix_type = square_matrix_type
   				   mljac = nz*nvar-1
   				   mujac = mljac
   				   if (nz*nvar > 51) then
   					  write(*,*) 'numerical jac is very slow for large nz*nvar'
   					  stop 1
   				   end if
   				else
   			   		matrix_type = banded_matrix_type
   				end if
   			else if (which_decsol == block_thomas_dble) then
   				call block_thomas_dble_work_sizes(nvar, nz, lrd, lid)
   				matrix_type = block_tridiag_dble_matrix_type
   			!else if (which_decsol == block_thomas_quad) then
   			!	call block_thomas_quad_work_sizes(nvar, nz, lrd, lid)
   			!		matrix_type = block_tridiag_quad_matrix_type
   			else if (which_decsol == bcyclic_dble) then
   				call bcyclic_dble_work_sizes(nvar, nz, lrd, lid)
   				matrix_type = block_tridiag_dble_matrix_type
   			!else if (which_decsol == block_dc_mt_dble) then
   			!	call block_dc_mt_dble_work_sizes(nvar, nz, lrd, lid)
   			!		matrix_type = block_tridiag_dble_matrix_type
   			!else if (which_decsol == block_dc_mt_quad) then
   			!	call block_dc_mt_quad_work_sizes(nvar, nz, lrd, lid)
   			!		matrix_type = block_tridiag_quad_matrix_type
   			else
   				write(*,*) 'bad value for which_decsol'
   				stop 1
   			end if
			
   			allocate(rpar_decsol(lrd), ipar_decsol(lid), stat=ierr)
            if (ierr /= 0) stop 1
         
         	call newton_work_sizes( &
            	mljac, mujac, nvar, nz, nsec, matrix_type, lwork, liwork, ierr)
         	if (ierr /= 0) stop 1
      
         	allocate(work(lwork), iwork(liwork), stat=ierr)
         	if (ierr /= 0) stop 1
         
           	work = 0
            iwork = 0
         	
            tol_correction_norm = 1d-1 ! upper limit on magnitude of average scaled correction
   			work(r_tol_residual_norm) = 1d-3
   			work(r_tol_max_residual) = 1d-2
         	
            AF1 => null()
            call newton(&
               nz, nvar, x1, xold1, matrix_type, mljac, mujac, &
               lapack_decsol, null_decsolblk, &
               !decsol, decsolblk, decsolblk_quad, &
               lrd, rpar_decsol, lid, ipar_decsol, which_decsol, tol_correction_norm, &
               set_primaries, default_set_secondaries, set_xscale, &
               Bdomain, default_xdomain, eval_equations, &
               size_equ, default_sizeB, default_inspectB, &
               enter_setmatrix, exit_setmatrix, failed_in_setmatrix, &
               default_force_another_iter, xscale1, equ1, ldy, nsec, y1, &
               work, lwork, iwork, liwork, &
               AF1, lrpar, rpar, lipar, ipar, nonconv, ierr)
            if (ierr /= 0) stop 1
         
            write(*,*) ierr, nonconv
         
            rho_c = x(1,1)
            pb_frac = x(2,1)
         
            call set_initial_conditions(rho_c)
         
            !x_wd = 1d0
   			!xend_wd = rsun
   			!y_wd(1) = rho_c							!rho(0) in g/cm^3
   			!y_wd(2) = 0d0							!M(0) in g
   			!y_wd(3) = 0d0							!t_sound(0) in s

   			call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
            
            !write(*,*) 'Solution found:'
            !write(*,'(a15,ES25.10)') 'm_c (msun)', m_c
   			!write(*,'(a15,ES25.10)') 'm_env (msun)', m_env
            !write(*,'(a15,ES25.10)') 'rho_c (g/cc)', rho_c
   			!write(*,'(a15,ES25.10)') 'pb_frac', pb_frac			
   			!write(*,'(a15,ES25.10)') 'rho_b', rho_b
   			!write(*,'(a15,ES25.10)') 'P_b', p_b
   			!write(*,'(a15,ES25.10)') 'R_c (cm)', r_c
   			!write(*,'(a15,ES25.10)') 'g_b', standard_cgrav*m_c*msun/r_c**2
   			!write(*,'(a15,ES25.10)') 't_sound', t_sound
   			!write(*,*)
   			!write(*,*) 'Comparisons to thin-shell limit:'
   			!write(*,'(a15,ES25.10)') 'Pb_thin', &
   			!	m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4)
   			!write(*,'(a15,ES25.10)') 'Delta P_b', (p_b - &
   			!	m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4))/p_b
            
            if (associated(AF1)) deallocate(AF1)
            deallocate(rpar_decsol, ipar_decsol, stat=ierr)
         	if (ierr /= 0) stop 1
         	deallocate(work, iwork, stat=ierr)
         	if (ierr /= 0) stop 1
      	
         end subroutine find_envelope_newton
      	
      	subroutine set_primaries(nvar, nz, x, lrpar, rpar, lipar, ipar, ierr)
        	integer, intent(in) :: nvar, nz
         	double precision, pointer :: x(:,:)
         	integer, intent(in) :: lrpar, lipar
         	double precision, intent(inout) :: rpar(:)
         	integer, intent(inout) :: ipar(:)
         	integer, intent(out) :: ierr
         	ierr = 0
         	xp0 = x(1, 1); xp1 = x(2, 1)
         	!write(*, '(a20, 3(a6, 1pe26.16, 3x))') 'primaries', 'x0', x0, 'x1', x1, 'x2', x2
         	
         	if (dbg) write(*, '(a20, 2(a6, 1pe26.16, 3x))') 'primaries', 'xp0', xp0, 'xp1', xp1
		end subroutine set_primaries
      	
		subroutine eval_equations(iter, nvar, nz, x, xscale, equ, lrpar, rpar, lipar, ipar, ierr)
			implicit none
			
			integer, intent(in) :: iter, nvar, nz
			real(dp), pointer, dimension(:,:) :: x, xscale, equ ! (nvar, nz)
			integer, intent(in) :: lrpar, lipar
			real(dp), intent(inout) :: rpar(:) ! (lrpar)
			integer, intent(inout) :: ipar(:) ! (lipar)
			integer, intent(out) :: ierr
			
			double precision :: rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound
			double precision :: mc_target, menv_target
			double precision :: x_wd, xend_wd

			ierr = 0
			equ = 0d0
			
			rho_c = xp0
			pb_frac = xp1
			mc_target = rpar(2)
			menv_target = rpar(3)
			
			call set_initial_conditions(rho_c)
			
			!WD:
			!x_wd = 1d0
			!xend_wd = rsun
			!y_wd(1) = rho_c						!rho(0) in g/cm^3
			!y_wd(2) = 4*pi*y_wd(1)*x_wd**3/3d0		!M(0) in g
			!y_wd(3) = 0d0							!t_sound(0) in s

			call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
			
			!equ(i,j) is equ(var,zone) and we're using one zone only
			equ(1,1) = m_c - mc_target
			equ(2,1) = m_env - menv_target
			
			!write(*,*) equ(:,1)

		end subroutine eval_equations

		subroutine eval_jacobian(ldA, A1, idiag, lrpar, rpar, lipar, ipar, ierr)
			integer, intent(in) :: ldA ! leading dimension of A
			real(dp), pointer :: A1(:) ! (ldA, nvar*nz) ! the jacobian matrix
			! A(idiag+q-v, v) = partial of equation(q) wrt variable(v)
			integer, intent(inout) :: idiag 
			integer, intent(in) :: lrpar, lipar
			real(dp), intent(inout) :: rpar(:) ! (lrpar)
			integer, intent(inout) :: ipar(:) ! (lipar)
			integer, intent(out) :: ierr         

		end subroutine eval_jacobian
		
		subroutine size_equ(iter, nvar, nz, equ, residual_norm, residual_max, lrpar, rpar, lipar, ipar, ierr)
			implicit none
			
			integer, intent(in) :: iter, nvar, nz
			real(dp), pointer :: equ(:,:) ! (nvar, nz)
			real(dp), intent(out) :: residual_norm, residual_max
			integer, intent(in) :: lrpar, lipar
			real(dp), intent(inout) :: rpar(:) ! (lrpar)
			integer, intent(inout) :: ipar(:) ! (lipar)
			integer, intent(out) :: ierr
			ierr = 0
			
			residual_norm = sum(abs(equ))/(neq)
			residual_max = maxval(abs(equ))
			
			!write(*,*) 'Residual norm:',residual_norm
			!write(*,*) 'Residual max:',residual_max
      	end subroutine size_equ
      	
      	! you might want to use a different value of xscale_min for this
		subroutine set_xscale(nvar, nz, xold, xscale, lrpar, rpar, lipar, ipar, ierr)
			implicit none
			
			integer, intent(in) :: nvar, nz
			real(dp), pointer :: xold(:,:) ! (nvar, nz)
			real(dp), pointer :: xscale(:,:) ! (nvar, nz)
			integer, intent(in) :: lrpar, lipar
			real(dp), intent(inout) :: rpar(:) ! (lrpar)
			integer, intent(inout) :: ipar(:) ! (lipar)
			integer, intent(out) :: ierr
			real(dp), parameter :: xscale_min = 1d-10
			
			xscale = max(xscale_min, abs(xold))
			ierr = 0
		end subroutine set_xscale
		
		! the proposed change to x is B*xscale*correction_factor
		! edit correction_factor and/or B as necessary so that the new x will be valid.
		! set ierr nonzero if things are beyond repair.
		subroutine Bdomain(iter, nvar, nz, B, x, xscale, correction_factor, lrpar, rpar, lipar, ipar, ierr)
			implicit none
			
			integer, intent(in) :: iter, nvar, nz
			real(dp), pointer, dimension(:,:) :: x, xscale, B ! (nvar, nz)
			real(dp), intent(inout) :: correction_factor
			integer, intent(in) :: lrpar, lipar
			real(dp), intent(inout) :: rpar(:) ! (lrpar)
			integer, intent(inout) :: ipar(:) ! (lipar)
			integer, intent(out) :: ierr
			
			double precision :: rho_c, pb_frac, delta_rhoc, delta_pbfrac
			
			rho_c = x(1,1)
			delta_rhoc = B(1,1)*xscale(1,1)*correction_factor
			pb_frac = x(2,1)
			delta_pbfrac = B(2,1)*xscale(2,1)*correction_factor
			
			!write(*,*) 'Proposed changes:'
			!write(*,*) delta_rhoc, delta_pbfrac
			!write(*,*)
			
			!Make sure our variables can't become negative:
			
			if((rho_c + delta_rhoc).le.0d0) then
				B(1,1) = rho_c/(2d0*xscale(1,1)*correction_factor)
			endif
			
			if((pb_frac + delta_pbfrac).le.0d0) then
				B(2,1) = pb_frac/(2d0*xscale(2,1)*correction_factor)
			endif
			
		end subroutine Bdomain

		subroutine enter_setmatrix( &
		iter, nvar, nz, neqs, x, xold, xscale, xder, need_solver_to_eval_jacobian, &
		ldA, A1, idiag, lrpar, rpar, lipar, ipar, ierr)
			integer, intent(in) :: iter, nvar, nz, neqs
			real(dp), pointer, dimension(:,:) :: x, xold, xscale, xder ! (nvar, nz)
			logical, intent(out) :: need_solver_to_eval_jacobian
			integer, intent(in) :: ldA ! leading dimension of A
			real(dp), pointer, dimension(:) :: A1 ! =(ldA, neqs)
			integer, intent(inout) :: idiag 
			integer, intent(in) :: lrpar, lipar
			real(dp), intent(inout) :: rpar(:) ! (lrpar)
			integer, intent(inout) :: ipar(:) ! (lipar)
			integer, intent(out) :: ierr
			real(dp) :: epsder
			include 'formats.dek'
			if (dbg) write(*, '(/, a)') 'enter_setmatrix'
			if (ipar(1) /= 0) then ! do numerical jacobian
				epsder = 1d-6 ! relative variation to compute numerical derivatives
				!xder = epsder*(xscale+abs(xold))
				xder = epsder*abs(xold)
				need_solver_to_eval_jacobian = .true.
			else
				call eval_jacobian(ldA, A1, idiag, lrpar, rpar, lipar, ipar, ierr)
					need_solver_to_eval_jacobian = .false.
			end if
		end subroutine enter_setmatrix

		subroutine exit_setmatrix( &
			iter, nvar, nz, neqs, dx, ldA, A1, idiag, xscale, lrpar, rpar, lipar, ipar, ierr)
			integer, intent(in) :: ldA ! leading dimension of A
			integer, intent(in) :: iter, nvar, nz, neqs ! number of equations, 2nd dimension of A
			integer, intent(inout) :: idiag ! row of A with the matrix diagonal entries
			real(dp), pointer, dimension(:,:) :: dx
			real(dp), pointer, dimension(:) :: A1
			real(dp), pointer, dimension(:,:) :: xscale ! (nvar, nz)
			integer, intent(in) :: lrpar, lipar
			real(dp), intent(inout) :: rpar(:) ! (lrpar)
			integer, intent(inout) :: ipar(:) ! (lipar)
			integer, intent(out) :: ierr
			integer :: i, j
			if (dbg) write(*, '(a, /)') 'exit_setmatrix'
			ierr = 0
		end subroutine exit_setmatrix

		subroutine failed_in_setmatrix(j, lrpar, rpar, lipar, ipar, ierr)
			integer, intent(in) :: j
			integer, intent(in) :: lrpar, lipar
			real(dp), intent(inout) :: rpar(:) ! (lrpar)
			integer, intent(inout) :: ipar(:) ! (lipar)
			integer, intent(out) :: ierr
			if (dbg) write(*, '(a, /)') 'failed_in_setmatrix'
			ierr = 0
		end subroutine failed_in_setmatrix
		
		!Given an envelope base density (rho_b) and base gravity (g_b), find the integration
      	!parameters rho_c and pb_frac that create a WD with that structure.
      	subroutine find_wd_params(rhob_target, gb_target, rhoc_guess, pbfrac_guess, &
      		rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
      		implicit none
      		
      		double precision, intent(in) :: rhob_target, gb_target, rhoc_guess, pbfrac_guess
      		double precision, intent(inout) :: rho_c, pb_frac, m_c, m_env, rho_b, &
      			p_b, r_c, r_wd, t_sound
      		
			!Begin Newton solver variables -----------------------------------------------
			integer :: liwork, lwork!, lqwork
			integer, dimension(:), pointer :: iwork
			real(dp), dimension(:), pointer :: work
			!real(qp), dimension(:), pointer :: qwork

			integer, parameter :: lrpar = 3, lipar = 1
			integer, target :: ipar_target(lipar)
			real(dp), target :: rpar_target(lrpar)
			integer, pointer :: ipar(:)
			real(dp), pointer :: rpar(:)

			integer :: lrd, lid
			integer, pointer :: ipar_decsol(:) ! (lid)
			real(dp), pointer :: rpar_decsol(:) ! (lrd)
			real(dp), pointer :: AF1(:)

			integer :: mljac, mujac
			real(dp) :: tol_correction_norm, tol_max_correction, tol_residual_norm
			logical :: nonconv
			integer :: ierr
			!End Newton solver variables -------------------------------------------------
			
			double precision :: mwd_target, rhoc_min, rhoc_max
			
			ierr = 0

			ipar => ipar_target
			rpar => rpar_target         

			if (do_numerical_jacobian) then
				ipar(1) = 1
			else
				ipar(1) = 0
			end if
			
			xold(:,1) = (/ rhoc_guess, pbfrac_guess /)
			rpar(2) = rhob_target
			rpar(3) = gb_target
			dx = 0d0
			x = xold

			mljac = 2*nvar-1 ! number of subdiagonals
			mujac = mljac ! number of superdiagonals
			if (which_decsol == lapack) then
				call lapack_work_sizes(nz*nvar, lrd, lid)
				if (do_numerical_jacobian) then
				   matrix_type = square_matrix_type
				   mljac = nz*nvar-1
				   mujac = mljac
				   if (nz*nvar > 51) then
					  write(*,*) 'numerical jac is very slow for large nz*nvar'
					  stop 1
				   end if
				else
			   		matrix_type = banded_matrix_type
				end if
			else if (which_decsol == block_thomas_dble) then
				call block_thomas_dble_work_sizes(nvar, nz, lrd, lid)
					matrix_type = block_tridiag_dble_matrix_type
			!else if (which_decsol == block_thomas_quad) then
			!	call block_thomas_quad_work_sizes(nvar, nz, lrd, lid)
			!		matrix_type = block_tridiag_quad_matrix_type
			else if (which_decsol == bcyclic_dble) then
				call bcyclic_dble_work_sizes(nvar, nz, lrd, lid)
					matrix_type = block_tridiag_dble_matrix_type
			!else if (which_decsol == block_dc_mt_dble) then
			!	call block_dc_mt_dble_work_sizes(nvar, nz, lrd, lid)
			!		matrix_type = block_tridiag_dble_matrix_type
			!else if (which_decsol == block_dc_mt_quad) then
			!	call block_dc_mt_quad_work_sizes(nvar, nz, lrd, lid)
			!		matrix_type = block_tridiag_quad_matrix_type
			else
				write(*,*) 'bad value for which_decsol'
				stop 1
			end if
			
			allocate(rpar_decsol(lrd), ipar_decsol(lid), stat=ierr)
         		if (ierr /= 0) stop 1
         
         	call newton_work_sizes( &
            	mljac, mujac, nvar, nz, nsec, matrix_type, lwork, liwork, ierr)
         	if (ierr /= 0) stop 1
         
         	allocate(work(lwork), iwork(liwork), stat=ierr)
         	if (ierr /= 0) stop 1
         
        	work = 0
         	iwork = 0
         	
         	tol_correction_norm = 1d-2 ! upper limit on magnitude of average scaled correction
			work(r_tol_residual_norm) = 1d-4
			work(r_tol_max_residual) = 1d-2
         	
            AF1 => null()
            call newton(&
               nz, nvar, x1, xold1, matrix_type, mljac, mujac, &
               lapack_decsol, null_decsolblk, &
               !decsol, decsolblk, decsolblk_quad, &
               lrd, rpar_decsol, lid, ipar_decsol, which_decsol, tol_correction_norm, &
               set_primaries, default_set_secondaries, set_xscale, &
               Bdomain, default_xdomain, eval_wd_params_equations, &
               size_equ, default_sizeB, default_inspectB, &
               enter_setmatrix, exit_setmatrix, failed_in_setmatrix, &
               default_force_another_iter, xscale1, equ1, ldy, nsec, y1, &
               work, lwork, iwork, liwork, &
               AF1, lrpar, rpar, lipar, ipar, nonconv, ierr)
            if (ierr /= 0) stop 1
            
            write(*,*) ierr, nonconv
            
            rho_c = x(1,1)
            pb_frac = x(2,1)
            
            call set_initial_conditions(rho_c)
            
            !x_wd = 1d0
			!xend_wd = rsun
			!y_wd(1) = rho_c						!rho(0) in g/cm^3
			!y_wd(2) = 0d0							!M(0) in g
			!y_wd(3) = 0d0							!t_sound(0) in s

			call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
            
            !write(*,*) 'Solution found:'
            !write(*,'(a15,ES25.10)') 'm_c (msun)', m_c
			!write(*,'(a15,ES25.10)') 'm_env (msun)', m_env
            !write(*,'(a15,ES25.10)') 'rho_c (g/cc)', rho_c
			!write(*,'(a15,ES25.10)') 'pb_frac', pb_frac			
			!write(*,'(a15,ES25.10)') 'rho_b', rho_b
			!write(*,'(a15,ES25.10)') 'P_b', p_b
			!write(*,'(a15,ES25.10)') 'R_c (cm)', r_c
			!write(*,'(a15,ES25.10)') 'g_b', standard_cgrav*m_c*msun/r_c**2
			!write(*,'(a15,ES25.10)') 't_sound', t_sound
			!write(*,*)
			!write(*,*) 'Comparisons to thin-shell limit:'
			!write(*,'(a15,ES25.10)') 'Pb_thin', &
			!	m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4)
			!write(*,'(a15,ES25.10)') 'Delta P_b', (p_b - &
			!	m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4))/p_b
            
            if (associated(AF1)) deallocate(AF1)
            deallocate(rpar_decsol, ipar_decsol, stat=ierr)
         		if (ierr /= 0) stop 1
         	deallocate(work, iwork, stat=ierr)
         		if (ierr /= 0) stop 1
      	
      	end subroutine find_wd_params
      	
      	subroutine eval_wd_params_equations(iter, nvar, nz, x, xscale, equ, lrpar, rpar, lipar, ipar, ierr)
			implicit none
			
			integer, intent(in) :: iter, nvar, nz
			real(dp), pointer, dimension(:,:) :: x, xscale, equ ! (nvar, nz)
			integer, intent(in) :: lrpar, lipar
			real(dp), intent(inout) :: rpar(:) ! (lrpar)
			integer, intent(inout) :: ipar(:) ! (lipar)
			integer, intent(out) :: ierr
			
			double precision :: rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound
			double precision :: rhob_target, gb_target
			double precision :: x_wd, xend_wd

			ierr = 0
			equ = 0d0
			
			rho_c = xp0
			pb_frac = xp1
			rhob_target = rpar(2)
			gb_target = rpar(3)
			
			call set_initial_conditions(rho_c)
			
			!WD:
			!x_wd = 1d0
			!xend_wd = rsun
			!y_wd(1) = rho_c						!rho(0) in g/cm^3
			!y_wd(2) = 4*pi*y_wd(1)*x_wd**3/3d0		!M(0) in g
			!y_wd(3) = 0d0							!t_sound(0) in s

			call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
			
			!equ(i,j) is equ(var,zone) and we're using one zone only
			equ(1,1) = rho_b - rhob_target
			equ(2,1) = standard_cgrav*m_c*msun/r_c**2 - gb_target
			
			!write(*,*) equ(:,1)
			!write(*,*)

		end subroutine eval_wd_params_equations
    
		subroutine wd_cleanup()
			call free_eos_handle(eos_handle)
			call eos_shutdown
			
			deallocate(work_wd, iwork_wd, xa_wd, xa_env)
			deallocate(net_iso, chem_id, dabar_dx, dzbar_dx, dmc_dx)
			deallocate(equ1, x1, xold1, dx1, xscale1, y1)
		end subroutine wd_cleanup
      
end module wd_solver
