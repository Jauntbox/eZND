module lane_emden_solver

	use const_def
	use const_lib
	use num_def
	use num_lib
	
	integer, parameter :: nv_le = 2  ! the number of variables in the van der Pol system of ODEs
	real(dp) :: n_le	!Polytropic index gamma = 1 + 1/n_le
	 real(dp), parameter :: eps = 1d-3 ! stiffness parameter for van der Pol
	 real(dp) :: rtol_le(1) ! relative error tolerance(s)
	 real(dp) :: atol_le(1) ! absolute error tolerance(s)
	 real(dp) :: x_le ! starting value for the interval of integration
	 real(dp) :: xend_le ! ending value for the interval of integration
	 real(dp) :: expect(nv_le), yprime(nv_le)
	 character (len=64) :: str
	 character (len=256) :: dir, fname
	 integer, parameter :: lrpar_le = 2, lipar_le = 1
	 real(dp) :: max_abs_yp2, h_le, max_step_size_le
	 integer :: io_unit, i, lout_le, iout_le, idid_le, itol_le, j
	 integer :: liwork_le, lwork_le, max_steps_le, ierr_le
	 real(dp), pointer :: work_le(:)
	 integer, pointer :: iwork_le(:)
	 real(dp), target :: y_ary(nv_le)
	 real(dp), pointer :: y_le(:)
	 real(dp), target :: rpar_ary(lrpar_le)
	 integer, target :: ipar_ary(lipar_le)
	 integer, pointer :: ipar_le(:) ! (lipar)
	 real(dp), pointer :: rpar_le(:) ! (lrpar)
	 
	 real(dp) :: r_n, rho_c, p_c, k_le

	contains
	
	subroutine show_results(nv,y,expect,show_all)
         integer, intent(in) :: nv
         real(dp), dimension(nv), intent(in) :: y, expect
         logical, intent(in) :: show_all
         integer :: i
         if (show_all) then
            write(*,'(/,a5,99a20)') 'i', 'calculated    ', 'reference    ', 'lg(abs rel diff)'
            do i=1,nv
               write(*,'(i5,2e20.10,f20.10)') i, y(i), expect(i), log10(abs(y(i)-expect(i))/max(1d-299,abs(expect(i))))
            end do
         else
            write(*,'(/,a5,99a20)') 'i', 'calculated    ', 'reference    '
            do i=1,nv
               write(*,'(i5,2e20.10,f20.10)') i, y(i), expect(i)
            end do
         end if
         write(*,*)
      end subroutine show_results
	
	subroutine van_der_Pol_derivs(n,x,y,f,lrpar,rpar,lipar,ipar,ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x
         real(dp), intent(inout) :: y(n)
         real(dp), intent(out) :: f(n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
         ierr = 0
         f(1) = y(2)
         f(2) = ((1 - y(1)**2) * y(2) - y(1))/rpar(1)
         ! the derivatives do not depend on x         
      end subroutine van_der_Pol_derivs
      
      subroutine lane_emden_derivs(n,x,h,y,f,lrpar,rpar,lipar,ipar,ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:)
         real(dp), intent(out) :: f(:)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
         ierr = 0
         
         if(y(1).lt.0.0) then
         	y(1) = 0.0
         endif
         
         f(1) = y(2)
         f(2) = -(y(1)**n_le) - 2*f(1)/x
         
         !rpar(1) = f(2)       

      end subroutine lane_emden_derivs
      
      subroutine lane_emden_init()
      	implicit none
      
      	ipar_le => ipar_ary
         rpar_le => rpar_ary

         write(*,*)        
         write(*,*) 'vdpol'
         write(*,*) 'cash_karp'
         
         y_le => y_ary

         !vdp:
         !x = 0
         !xend = 2.0
         !y(1) = 2
         !y(2) = 0
         
         !le:
         x_le = 1d-6
         xend_le = 3.65375
         !xend = 3.653
         n_le = 1.5
         y_le(1) = 1.0
         y_le(2) = 0.0
         
         lout_le = 6
         max_steps_le = 10000
         max_step_size_le = 9

         itol_le = 0 ! scalar tolerances
         iout_le = 0 ! no intermediate output
         
         rtol_le(1) = 1d-6
         atol_le(1) = 1d-6
         h_le = 1d-10
         
         !vdp:
        !rpar_le(1) = eps
         
         !lane-emden
         rpar_le = 0
         ipar_le = 0

         call cash_karp_work_sizes(nv_le,liwork_le,lwork_le)
         allocate(work_le(lwork_le), iwork_le(liwork_le))
         
         iwork_le = 0
         work_le = 0
      end subroutine lane_emden_init

	subroutine test_cash_karp(show_all)
		implicit none
		
         logical, intent(in) :: show_all
         
         ierr_le = 0
         call cash_karp( &
               !nv_le,van_der_Pol_derivs,x,y,xend, &
               nv_le,lane_emden_derivs,x_le,y_le,xend_le, &
               h_le,max_step_size_le,max_steps_le, &
               rtol_le,atol_le,itol_le, &
               null_solout,iout_le,work_le,lwork_le,iwork_le,liwork_le, &
               lrpar_le,rpar_le,lipar_le,ipar_le,lout_le,idid_le)

         if (idid_le /= 1) then ! trouble
            write(*,*) 'idid', idid_le
            stop 1
         end if
		
		!vdp
         !expect(1:2) = (/ 1.7632345401889102d+00, -8.3568868191466206d-01 /)
         
        !lane-emden
         expect(1:2) = (/ 0.0, -0.2033 /)
         
         call show_results(nv_le,y_le,expect,show_all)
         
         if (.not. show_all) return
         
         write (6,91) (iwork_le(j),j=1,4)
 91      format(' fcn=',i8,'     step=',i6,'     accpt=',i6,'     rejct=',i5)
 
         write(*,*)
      
      end subroutine test_cash_karp
      
      	!Lane-Emden function: theta(xi)
      	real(dp) function theta_le(x,dfdx,lrpar,rpar,lipar,ipar,ierr)
      		implicit none
      		
         integer, intent( in ) :: lrpar, lipar
         real(dp), intent( in ) :: x
         real(dp), intent( out ) :: dfdx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent( out ) :: ierr
         
         real(dp) :: xstart
         
         ierr = 0
         xstart = 1d-8	!xend = x, here
         y_le(1) = 1d0
         y_le(2) = 0d0
         
         call cash_karp( &
               !nv_le,van_der_Pol_derivs,x,y,xend, &
               nv_le,lane_emden_derivs,xstart,y_le,x, &
               h_le,max_step_size_le,max_steps_le, &
               rtol_le,atol_le,itol_le, &
               null_solout,iout_le,work_le,lwork_le,iwork_le,liwork_le, &
               lrpar_le,rpar_le,lipar_le,ipar_le,lout_le,idid_le)

         if (idid_le /= 1) then ! trouble
            write(*,*) 'idid', idid_le
            stop 1
         end if
         
         theta_le = y_le(1)
         dfdx = y_le(2)
      end function theta_le
      
      !Lane-Emden function: theta(xi)
      	real(dp) function m_xi_func(x,dfdx,lrpar,rpar,lipar,ipar,ierr)
      		implicit none
      		
         integer, intent( in ) :: lrpar, lipar
         real(dp), intent( in ) :: x
         real(dp), intent( out ) :: dfdx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent( out ) :: ierr
         
         real(dp) :: xstart
         
         ierr = 0
         xstart = 1d-8	!xend = x, here
         y_le(1) = 1d0
         y_le(2) = 0d0
         
         call cash_karp( &
               !nv_le,van_der_Pol_derivs,x,y,xend, &
               nv_le,lane_emden_derivs,xstart,y_le,x, &
               h_le,max_step_size_le,max_steps_le, &
               rtol_le,atol_le,itol_le, &
               null_solout,iout_le,work_le,lwork_le,iwork_le,liwork_le, &
               lrpar,rpar,lipar,ipar,lout_le,idid_le)

         if (idid_le /= 1) then ! trouble
            write(*,*) 'idid', idid_le
            stop 1
         end if
         
         m_xi_func = -4*pi*r_n**3*rho_c*x**2*y_le(2) - rpar(1)*msun
         dfdx = -y_le(1)**n_le - 2*y_le(2)/x
         !write(*,*) 'm_xi_func', m_xi_func
      end function m_xi_func
      
      	!Searches for the first zero of the Lane-Emden function, theta(xi1) = 0:
		subroutine find_xi1(n, xi1_low, xi1_high, xi1, dtheta_dxi)
      		implicit none
      		real(dp), intent(in) :: n	!Index of Lane-Emden eqn (eg. 3/2 for NR degeneracy)
      		real(dp), intent(in) :: xi1_low, xi1_high
      		real(dp), intent(inout) :: xi1, dtheta_dxi
      		
      		real(dp) :: y_low, y_high, epsx, epsy, dfdx, theta
      		!integer, parameter :: lrpar = 0, lipar = 0
      		!real(dp), pointer :: rpar(lrpar)
      		!integer, pointer :: ipar(lipar)
      		integer :: imax, ierr
      		
      		n_le = n
      		y_low = theta_le(xi1_low,dfdx,lrpar_le,rpar_le,lipar_le,ipar_le,ierr)
      		y_high = theta_le(xi1_high,dfdx,lrpar_le,rpar_le,lipar_le,ipar_le,ierr)
      		epsx = 1d-6
      		epsy = 1d-6
      		
      		ierr = 0
      		xi1 = safe_root(theta_le, xi1_low, xi1_high, y_low, y_high, imax, epsx, epsy, &
      			lrpar_le, rpar_le, lipar_le, ipar_le, ierr)
      		!find dtheta_dxi:
      		theta = theta_le(xi1,dtheta_dxi,lrpar_le,rpar_le,lipar_le,ipar_le,ierr)
      		
      		if (ierr /= 0) then ! trouble
            	write(*,*) 'ierr', ierr
            	stop 1
         	end if
         	
      	end subroutine find_xi1
      	
      	!Given a WD mass and radius, calculate the auxiliary variables:
      	subroutine get_aux_vars(m_wd, r_wd, xi1, dtheta_dxi)
      		real(dp), intent(in) :: m_wd, r_wd, xi1, dtheta_dxi
      	
      		rho_c = -1/3d0 * m_wd*msun/(4*pi*r_wd**3/3d0) * (xi1/dtheta_dxi)
			p_c = standard_cgrav*(m_wd*msun)**2/(r_wd**4)*1d0/(4*pi*(n_le+1)*dtheta_dxi**2)
			k_le = p_c/(rho_c**(1+1d0/n_le))
			r_n = sqrt((n_le+1)*p_c/(4*pi*standard_cgrav*rho_c**2))
      	end subroutine get_aux_vars
      	
      	!Given an M_c, find the xi such that M(xi) = M_c
		subroutine find_xi_mc(mc, xi_low, xi_high, xi, dtheta_dxi)
      		implicit none
      		!real(dp), intent(in) :: n	!Index of Lane-Emden eqn (eg. 3/2 for NR degeneracy)
      		real(dp), intent(in) :: xi_low, xi_high, mc
      		real(dp), intent(inout) :: xi, dtheta_dxi
      		
      		real(dp) :: y_low, y_high, epsx, epsy, dfdx, theta
      		!integer, parameter :: lrpar = 0, lipar = 0
      		!real(dp), pointer :: rpar(lrpar)
      		!integer, pointer :: ipar(lipar)
      		integer :: imax, ierr
      		
      		rpar_le(1) = mc
      		y_low = m_xi_func(xi_low,dfdx,lrpar_le,rpar_le,lipar_le,ipar_le,ierr)
      		y_high = m_xi_func(xi_high,dfdx,lrpar_le,rpar_le,lipar_le,ipar_le,ierr)
      		epsx = 1d-6
      		epsy = 1d-6
      		
      		ierr = 0
      		xi = safe_root(m_xi_func, xi_low, xi_high, y_low, y_high, imax, epsx, epsy, &
      			lrpar_le, rpar_le, lipar_le, ipar_le, ierr)
      		!find dtheta_dxi:
      		theta = theta_le(xi,dtheta_dxi,lrpar_le,rpar_le,lipar_le,ipar_le,ierr)
      		
      		!write(*,*) xi, dtheta_dxi
      		
      		if (ierr /= 0) then ! trouble
            	write(*,*) 'ierr', ierr
            	stop 1
         	end if
         	
      	end subroutine find_xi_mc
      
	  subroutine lane_emden_cleanup
		deallocate(work_le, iwork_le)
	  end subroutine lane_emden_cleanup
      
end module lane_emden_solver
