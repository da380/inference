program S2_point_bound

  use nrtype
  use module_S2_point
  use module_maps
  use ran_state
  use module_gsht
  use module_interp
  use module_opt
  use module_crust
  implicit none
  

  character(len = 256) :: efile,pref
  
  integer(i4b) :: np,lmax,io,nparm,iparm,jparm,info,lwork,nth,ith,lt,i
  integer(i4b), dimension(:), allocatable :: la,ma
  

  real(dp) :: s,mu,rtol,U2_min,U2_true,th,dth,fac1,fac2,U2
  
  real(dp), dimension(:,:), allocatable :: point_loc
  real(dp), dimension(:,:), allocatable :: point_dat

  real(dp), dimension(:), allocatable :: u
  real(dp), dimension(:), allocatable :: ut
  real(dp), dimension(:), allocatable :: v
  real(dp), dimension(:), allocatable :: wt
  real(dp), dimension(:), allocatable :: w0
  real(dp), dimension(:), allocatable :: w
  real(dp), dimension(:), allocatable :: wh
  real(dp), dimension(:), allocatable :: work
  real(dp), dimension(:), allocatable :: eval
  real(dp), dimension(:), allocatable :: eval0
  

  real(dp), dimension(:,:), allocatable :: bbs
  real(dp), dimension(:,:), allocatable :: bbs0

  

  !====================================================!
  !          set the model and synthetic data          !
  !====================================================!


  ! set up the GSHT routines
  lmax = 128
  call setup_gsht(lmax = lmax,nmax = 0)
  
  ! set up random numbers
!  call setup_random()

  ! set the sobolev parameters
  s_S2_point =  2.0_dp
  mu_S2_point = 0.25_dp

  ! set up the crustal model
  call get_crust

  np = 400
  call random_locations_topo(np,point_loc,sign = -1)
    
  ! set the point locations
  m_S2_point = np
  allocate(x_S2_point(np,2))
  x_S2_point(:,1) = (90.0_dp-point_loc(:,1))*deg2rad
  x_S2_point(:,2) = point_loc(:,2)*deg2rad

  ! allocate the data vector
  allocate(v_S2_point(m_S2_point,2))
  
  ! plot the coast lines
  call write_coast('hammer','coast.out')
    
  ! set dimension of model vector 
  n_S2_point = nth_gsht*nph_gsht

  ! allocate the model vector 
  allocate(u(n_S2_point))
  allocate(ut(n_S2_point))

!========================================================++!  
!  ! make random model
!  allocate(qq_S2_point(lmax+1))
!  call covariance_given_norm(lmin = 1,lmax = lmax,t = 3.0_dp, &
  !       un = 1.0_dp,qq = qq_S2_point,mu = 0.25_dp)
!==========================================================!

  open(newunit = io,file='norm.out')
  do i = 1,10000
     print *, i
     call random_model(lmin = 1,lmax = lmax,qq = qq_S2_point,u = u)
     call sobolev_product_S2(u,u,U2_true)
     write(io,*)  U2_true
  end do
  close(io)
  
  stop
  
  call plot_function_S2('fun.in',360,720,u)




  
  
  ! make synthetic data
  allocate(v(m_S2_point))
  call point_operator_S2(u,v_S2_point(:,1))
  call write_point_data('hoggard.dat',np,point_loc,v_S2_point)  


  ! reset the sobolev parameters
  s_S2_point =  2.00_dp
  mu_S2_point = 0.5_dp

!  ! reset up the GSHT routines
!  lmax = 128
!  call setup_gsht(lmax = lmax,nmax = 0)
  
  ! get the true norm of the model
  call sobolev_product_S2(u,u,U2_true)
  
  !=============================================!
  !          set up the property space          !
  !=============================================!

  ! set the number of parameters
  lt = 2
  nparm_S2_point = 2*lt+1

  ! allocate spherical harmonic index arrays
  allocate(la_S2_point(nparm_S2_point))
  allocate(ma_S2_point(nparm_S2_point))

  ! set index values
  do iparm = 1,2*lt+1
     la_S2_point(iparm) =  lt
     ma_S2_point(iparm) =  -(lt+1)+iparm
    end do


  
  ! allocate some property vectors
  allocate(w(nparm_S2_point))
  allocate(w0(nparm_S2_point))
  allocate(wt(nparm_S2_point))
  allocate(wh(nparm_S2_point))

  ! get property vector for true model
  call spharm_operator_S2(u,w0)


  
  !================================================================!
  !                 set up the optimization codes                  !
  !================================================================!

  ! set dimension of model vector 
  n_S2_point = nth_gsht*nph_gsht

  ! set the function to optimise
  fun1 => point_misfit_S2
  
  ! set the line search method
  call setup_cvsrch(n_S2_point,1.0e-6_dp,1.e-6_dp,10.0_dp)
  lins => cvsrch

  ! set the inner product
  iprod => sobolev_product_S2

  ! set the preconditioner
  precon => preid

  ! set the tolerance for optimisations
  rtol = 1.e-4_dp

  
  !=============================================!
  !        find the minimum norm solution       !
  !=============================================!

  ! setup the mask for the constrained coefficients
  allocate(mask_S2_point(ncoef_gsht))
  mask_S2_point(:) = .true.
  mask_S2_point(ilm(0,0)) = .false.
  
  ! solve the optimisation problem
  u = 0.0_dp
  call lbfgs(u,5,verb = .true.,rtol = rtol)
  call plot_function_S2('fun.min',360,720,u)

  ! get the minimum norm property vector
  call spharm_operator_S2(u,wt)

  ! get the squared norm of the model
  call sobolev_product_S2(u,u,U2_min)


  
  !=========================================!
  !   compute  BB^{*} and B_{0}B_{0}^{*}    !
  !=========================================!

  allocate(bbs(nparm_S2_point,nparm_S2_point))
  allocate(bbs0(nparm_S2_point,nparm_S2_point))

  
  do iparm = 1,nparm_S2_point

     ! compute action of B^{*} on a basis vector
     w = 0.0_dp
     w(iparm) = 1.0_dp
     call adjoint_spharm_operator_S2(w,ut)


     if(nparm_S2_point == 1) then
        call plot_function_S2('fun.real',360,720,ut)
     end if

     ! compute the action of B on the result
     call spharm_operator_S2(ut,bbs0(:,iparm))
     
     ! compute the associated data
     call point_operator_S2(ut,v_S2_point(:,1))

     ! solve the minimum norm problem
     u = 0.0_dp
     call lbfgs(u,5,verb = .true.,rtol = rtol)
     
     ! compute the projection
     ut = ut-u

     if(nparm_S2_point == 1) then
        call plot_function_S2('fun.proj',360,720,ut)
     end if

     ! compute the action of B on the result
     call spharm_operator_S2(ut,bbs(:,iparm))
     
  end do


  ! print results to file
  open(newunit = io,file='S2_point_bound.dat')
  write(io,*) la_S2_point
  write(io,*) ma_S2_point
  write(io,*) w0
  write(io,*) wt
  close(io)

  open(newunit = io,file='S2_point_bound.mat')
  write(io,*) U2_min,w0(2:)*0.0_dp
  write(io,*) w0
  do iparm = 1,nparm_S2_point
     write(io,*) bbs0(iparm,:)
  end do
  write(io,*) wt
  do iparm = 1,nparm_S2_point
     write(io,*) bbs(iparm,:)
  end do     
  close(io)

  !-------------------------------------!
  !       diagonalise the matrices      !
  !-------------------------------------!


  
  lwork = 3*nparm_S2_point-1
  allocate(work(lwork))
  allocate(eval(nparm_S2_point))
  allocate(eval0(nparm_S2_point))
  call dsyev('V', 'U', nparm_S2_point, bbs0, nparm_S2_point, &
             eval0, work, lwork, info)
  call dsyev('V', 'U', nparm_S2_point, bbs, nparm_S2_point, &
              eval, work, lwork, info)  


!  do iparm = 1,nparm_S2_point
!     eval0(iparm) = bbs0(iparm,iparm)
!     bbs0(iparm,:) = 0.0_dp
!     bbs0(iparm,iparm) = 1.0_dp
!  end do


  
  

  !---------------------------------------------!
  !            plot bounding ellipse            !
  !---------------------------------------------!

  ! set norm bound
  U2 = 1.1_dp*sqrt(U2_min)
  U2 = U2**2
  
  nth = 360
  dth = twopi_d/(nth-1)
  
  do iparm = 1,nparm_S2_point
     do jparm = iparm+1,nparm_S2_point

        call string_cat_int('S2_point_bound.out.',iparm,pref)
        call string_cat_int(trim(pref)//'.',jparm,efile)

        open(newunit = io,file=trim(efile))

        ! initialise property vector
        w = 0.0_dp
        
        do ith = 1,nth

           th = (ith-1)*dth
           w(iparm) = cos(th)
           w(jparm) = sin(th)

           fac1 = dot_product(bbs(:,iparm),w)**2/eval(iparm) + & 
                  dot_product(bbs(:,jparm),w)**2/eval(jparm)

           fac2 = dot_product(bbs0(:,iparm),w)**2/eval0(iparm) + & 
                  dot_product(bbs0(:,jparm),w)**2/eval0(jparm) 

           fac1 = sqrt((U2-U2_min)/fac1)
           fac2 = sqrt(U2/fac2)


           write(io,*) wt(iparm) + fac1*w(iparm),wt(jparm) + fac1*w(jparm), &
                       fac2*w(iparm),fac2*w(jparm)
           
        end do
        
        close(io)


        
     end do
  end do

end program S2_point_bound
