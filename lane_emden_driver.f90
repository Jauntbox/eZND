program main
	use lane_emden_solver
	
	implicit none

	logical :: show_all
	double precision :: xi, theta, xi1, dtheta_dxi, xi1_low, xi1_high, index
	double precision :: xi_low, xi_high
	double precision :: m_xi, h_p, grav, rho_b, p_b, p_thin, dfdx
	double precision :: m_wd, r_wd, m_env, m_c
	integer :: nloop
	
	mesa_dir = '/Users/Kevin/mesa'
	call const_init(mesa_dir, ierr_le)
	
	call lane_emden_init()
	
	show_all = .true.
	call test_cash_karp(show_all)
	
	index = 1.5
	xi1_low = 1.0
	xi1_high = 5.0
	call find_xi1(index, xi1_low, xi1_high, xi1, dtheta_dxi)
	write(*,*) 'Lane-Emden solution:'
	write(*,'(a10,ES20.10)') 'xi1:', xi1
	write(*,'(a10,ES20.10)') 'dtheta_dxi:', dtheta_dxi
	write(*,*)
	
	m_wd = 1.0	!solar mass
	r_wd = 7.62d-3*rsun*m_wd**(-1.616)
	
	!p_c = (m_wd*msun)**2/(4*pi*r_wd**4*dtheta_dxi**2)
	!rho_c = sqrt((nle+1)*p_c*xi1**2/(4*pi*standard_cgrav*r_wd**2))
	call get_aux_vars(m_wd, r_wd, xi1, dtheta_dxi)
	write(*,'(a10,ES25.10)') 'rho_c', rho_c
	write(*,'(a10,ES25.10)') 'P_c', p_c
	write(*,'(a10,ES25.10)') 'r_n', r_n
	write(*,'(a10,ES25.10)') 'xi1*rn', xi1*r_n
	write(*,'(a10,ES25.10)') 'r_wd', r_wd
	write(*,*)
	
	nloop = 100
	write(*,'(3a18)') 'M_c (msun)', 'M_env (msun)', 'Delta P_b' !'P_b (dyne/cm^2)', 'P_thin (dyne/cm^2)'
	do i=1,nloop
		xi = xi1*10**(0.25*((i-1)*1d0/nloop-1.0))
		m_env = 10**(3.0*(-1d0 + (i-1)*1d0/nloop))
		m_c = m_wd - m_env
		
		xi_low = 1d-4
		xi_high = xi1
		call find_xi_mc(m_c, xi_low, xi_high, xi, dtheta_dxi)
		
		theta = theta_le(xi,dtheta_dxi,lrpar_le,rpar_le,lipar_le,ipar_le,ierr_le)
		m_xi = -4*pi*r_n**3*rho_c*xi**2*dtheta_dxi/msun
		p_b = k_le*rho_c**(1+1d0/nle)*theta**(nle+1)
		grav = standard_cgrav*m_xi*msun/(r_n*xi)**2
		rho_b = rho_c*theta**nle
		m_env = m_wd - m_xi
		h_p = p_b/(rho_b*grav)
		p_thin = m_env*msun*grav/(4*pi*r_wd**2)
		write(*,'(3ES18.10)') m_c, m_env, abs(p_b-p_thin)/p_thin
		!write(*,'(3ES18.10)') m_xi, m_env, abs(p_b-p_thin)/p_thin
	end do
	
	call lane_emden_cleanup()
end program