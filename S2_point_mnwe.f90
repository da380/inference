program S2_point_mnwe

  use nrtype
  use module_S2_point
  use module_gsht
  use module_opt
  use ran_state
  implicit none

  logical(lgt) :: go
  character(len=256) :: string
  integer(i4b) :: lmax,lmin,l,narg,nsamp,isamp,i,maxit,it
  real(dp) :: rtol,alpha,Jc,eta1,eta2
  real(dp), dimension(:,:), allocatable :: point_loc
  real(dp), dimension(:), allocatable :: u,samp
  

  narg = command_argument_count()
  if(narg == 0 .or. narg /= 6) then
     print *, '[lmin] [lmax] [s] [mu] [rtol] [alpha] '
     stop
  end if
  call get_command_argument(1,string)
  read(string,*) lmin
  call get_command_argument(2,string)
  read(string,*) lmax
  call get_command_argument(3,string)
  read(string,*) s_S2_point
  call get_command_argument(4,string)
  read(string,*) mu_S2_point
  call get_command_argument(5,string)
  read(string,*) rtol
  call get_command_argument(6,string)
  read(string,*) alpha

  ! seed the random numbers
  call setup_random()
  
  ! set up the GSHT routines
  call setup_gsht(lmax = lmax,nmax = 0)
  
  ! read in the data
  call read_point_data('syn.dat',m_S2_point,point_loc,v_S2_point)
  allocate(x_S2_point(m_S2_point,2))
  x_S2_point(:,1) = (90.0_dp-point_loc(:,1))*deg2rad
  x_S2_point(:,2) = point_loc(:,2)*deg2rad

  ! set dimension of model vector 
  n_S2_point = nth_gsht*nph_gsht
  allocate(u(n_S2_point))

  ! setup the mask for the constrained coefficients
  allocate(mask_S2_point(ncoef_gsht))
  mask_S2_point(:) = .true.
  do l = 0,lmin-1
     mask_S2_point(ilm(l,-l):ilm(l,l)) = .false.
  end do

  ! get the squared radius of the confidence set
  call gaussian_confidence_set(m_S2_point,alpha,ss_S2_point)

  print *, ss_S2_point

  
  ! set the function to optimise
  fun1 => point_misfit_we_S2
  
  ! set the line search method
  call setup_cvsrch(n_S2_point,1.0e-6_dp,1.e-6_dp,10.0_dp)
  lins => cvsrch

  ! set the inner product
  iprod => sobolev_product_S2

  ! set the preconditioner
  precon => preid


  ! check to see if the zero-model is compatible with the data
  u = 0.0_dp
  call point_misfit_we_con_S2(u,Jc)

  if(Jc <= 0.0_dp) then
     print *, ' zero-model sufficient'
     call plot_function_S2('min.mod',360,720,u)
     stop
  end if

  ! set an initial guess for the Lagrange multiplier
  eta1 = 0.01_dp
  eta2 = 0.0_dp
  
  ! find a lower bound
  print *, ' finding lower bound for multipler'
  do 
     eta_S2_point = eta1
     u = 0.0_dp
     call lbfgs(u,5,verb = .false.,rtol = rtol)
     call point_misfit_we_con_S2(u,Jc)
     print *, eta1,Jc
     if(Jc > 0.0_dp) exit
     eta2 = eta1
     eta1 = 0.5_dp*eta1
  end do


  ! find an upper bound
  print *, ' finding upper bound for multipler'
  if(eta2 == 0.0_dp) then
     eta2 = 2.0_dp*eta1
     do
        eta_S2_point = eta2
        u = 0.0_dp
        call lbfgs(u,5,verb = .false.,rtol = rtol)
        call point_misfit_we_con_S2(u,Jc)
        print *, eta2,Jc
        if(Jc < 0.0_dp) exit
        eta1 = eta2
        eta2 = 2.0_dp*eta2
     end do
  end if

  
  ! bisect to find the root
  print *, ' bisecting to find multipler'
  do
     eta_S2_point = 0.5_dp*(eta1+eta2)
     u = 0.0_dp
     call lbfgs(u,5,verb = .false.,rtol = rtol)
     call point_misfit_we_con_S2(u,Jc)
     print *, eta_S2_point,Jc/ss_S2_point
     if(abs(Jc)/ss_S2_point < rtol) exit
     if(Jc > 0.0_dp) eta1 = eta_S2_point
     if(Jc < 0.0_dp) eta2 = eta_S2_point
  end do

  ! plot the minimum norm solution
  call plot_function_S2('min.mod',360,720,u)
  




  
end program S2_point_mnwe
