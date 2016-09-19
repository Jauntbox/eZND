program main
	use const_def
	use const_lib
	use wd_solver
	
	implicit none
	
	integer, parameter :: nv_wd = 60
	double precision, parameter :: m_ch = 1.38	!Chandrasekhar mass
	integer :: i, num_steps
	double precision :: rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, t_sound, m_wd, r_wd
	double precision :: mc_target, menv_target, mwd_target, rhoc_min, rhoc_max
	double precision :: rhoc_guess, pbfrac_guess, rhob_target, gb_target
    
    !(should output 60, not 2 since nv_wd is private in the wd_solver module)
    !write(*,*) 'In wd_driver nv_wd:', nv_wd
	
	call init_mesa_modules()
	
	rho_c =  3.1838d7
	pb_frac =  7.9117d-4
	call wd_init()
	
	!call set_initial_conditions(rho_c)
	!call wd_integrate(m_wd)
	
	!1D Solver for rho_c given a total WD mass, mwd_target (in msun units)
	!mwd_target = 1.0
	!rhoc_min = 1d6
	!rhoc_max = 1d10
	!call find_rhoc(mwd_target, rhoc_min, rhoc_max, rho_c)
	!call set_initial_conditions(rho_c)
	!call wd_integrate(m_wd)
	
	!rho_c = 1d8
	!pb_frac = 1d-2
	!call set_initial_conditions(rho_c)
	!call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
	!call set_initial_conditions(rho_c)
	!call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
	!call set_initial_conditions(rho_c)
	!call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
	!write(*,*) m_c, m_env, rho_b, p_b, r_c
	
	mc_target = 0.6
	menv_target = 0.01
	rhoc_guess = 3d6
	pbfrac_guess = 1d-4
	call find_envelope_newton(mc_target, menv_target, rhoc_guess, pbfrac_guess, &
			rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
	call set_initial_conditions(rho_c)
	call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
		
	!rhob_target = 2d5	!g/cm^3
	!gb_target = 2d8		!cm/s^2
	!rhoc_guess = 1d8
	!pbfrac_guess = 1d-3
	!call find_wd_params(rhob_target, gb_target, rhoc_guess, pbfrac_guess, &
	!	rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
	!call set_initial_conditions(rho_c)
	!call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
	
	write(*,*) 'Solution found:'
	write(*,'(a15,ES25.10)') 'm_c (msun)', m_c
	write(*,'(a15,ES25.10)') 'm_env (msun)', m_env
	write(*,'(a15,ES25.10)') 'R_c (cm)', r_c
	write(*,'(a15,ES25.10)') 'R_wd (cm)', r_wd
	write(*,'(a15,ES25.10)') 'rho_c (g/cc)', rho_c
	write(*,'(a15,ES25.10)') 'pb_frac', pb_frac			
	write(*,'(a15,ES25.10)') 'rho_b', rho_b
	write(*,'(a15,ES25.10)') 'P_b', p_b
	write(*,'(a15,ES25.10)') 'g_b', standard_cgrav*m_c*msun/r_c**2
	write(*,'(a15,ES25.10)') 't_sound', t_sound
	write(*,'(a15,ES25.10)') 'min v_det', pi*r_c/t_sound
	write(*,'(a15,ES25.10)') 'rho_0', &
		rho_b*(1-(5d0/3-1)*0.5/(5d0/3))**(1d0/(5d0/3-1))
	write(*,*)
	write(*,*) 'Comparisons to thin-shell limit:'
	write(*,'(a15,ES25.10)') 'Pb_thin', &
		m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4)
	write(*,'(a15,ES25.10)') 'Delta P_b', (p_b - &
		m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4))/p_b
	
!	open(unit=7, file='wd_models/mc100_tenv2d8.data')
!	write(7,'(99a25)') 'm_c (msun)', 'm_env (msun)','R_c (cm)', 'R_wd (cm)', &
!		'rho_c (g/cc)', 'pb_frac', 'rho_b', &
!		'P_b', 'g_b', 't_sound', 'min v_det', 'rho_0', 'Pb_thin', 'Delta P_b'
!	write(7,*)
	
	num_steps = 100
	mc_target = 1.0
	menv_target = (m_ch - mc_target)*1*1d0/(num_steps)
	!rho_c = 1d8
			
	!We'll start off our solution by finding the rho_c that will support the
	!total WD mass:
	!mwd_target = mc_target + menv_target
	!rhoc_min = 1d6
	!rhoc_max = 6d10
	!call find_rhoc(mwd_target, rhoc_min, rhoc_max, rho_c)
	!pb_frac = 1d-4
	
	do i=1,0!num_steps
		!Code for solving for a specific core/envelope mass:
		menv_target = (m_ch - mc_target)*i*1d0/(num_steps)
		write(*,*) 'menv_target', menv_target
		
		call find_envelope_newton(mc_target, menv_target, rho_c, pb_frac, &
			rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, r_wd, t_sound)
		
		!Code for just sweeping through pb_frac, given a fixed rho_c
	!	pb_frac = 10**(-5d0 + 4d0*i/num_steps)
	!	write(*,*) 'Using pb_frac = ',pb_frac
	!	call set_initial_conditions(rho_c)
	!	call wd_envelope(rho_c, pb_frac, m_c, m_env, rho_b, p_b, r_c, t_sound)
		
		write(*,*) 'Solution found:'
		write(*,'(a15,ES25.10)') 'm_c (msun)', m_c
		write(*,'(a15,ES25.10)') 'm_env (msun)', m_env
		write(*,'(a15,ES25.10)') 'R_c (cm)', r_c
		write(*,'(a15,ES25.10)') 'R_wd (cm)', r_wd
		write(*,'(a15,ES25.10)') 'rho_c (g/cc)', rho_c
		write(*,'(a15,ES25.10)') 'pb_frac', pb_frac			
		write(*,'(a15,ES25.10)') 'rho_b', rho_b
		write(*,'(a15,ES25.10)') 'P_b', p_b
		write(*,'(a15,ES25.10)') 'g_b', standard_cgrav*m_c*msun/r_c**2
		write(*,'(a15,ES25.10)') 't_sound', t_sound
		write(*,'(a15,ES25.10)') 'min v_det', pi*r_c/t_sound
		write(*,'(a15,ES25.10)') 'rho_0', &
			rho_b*(1-(5d0/3-1)*0.5/(5d0/3))**(1d0/(5d0/3-1))
		write(*,*)
		write(*,*) 'Comparisons to thin-shell limit:'
		write(*,'(a15,ES25.10)') 'Pb_thin', &
			m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4)
		write(*,'(a15,ES25.10)') 'Delta P_b', (p_b - &
			m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4))/p_b
			
		write(7,'(99ES25.10)') m_c, m_env, r_c, r_wd, rho_c, pb_frac, rho_b, P_b, &
			standard_cgrav*m_c*msun/r_c**2, t_sound, pi*r_c/t_sound, &
			rho_b*(1-(5d0/3-1)*0.5/(5d0/3))**(1d0/(5d0/3-1)), &
			m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4), &
			(p_b - m_env*msun*standard_cgrav*m_c*msun/(4*pi*r_c**4))/p_b
			
		!if(rho_b*(1-(5d0/3-1)*0.5/(5d0/3))**(1d0/(5d0/3-1)).gt.5d6) then
		!	exit
		!endif
	end do
	
	close(7)
	call wd_cleanup()
end program