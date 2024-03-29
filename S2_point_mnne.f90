program S2_point_mnne

  use nrtype
  use module_S2_point
  use module_gsht
  use module_opt
  implicit none

  character(len=256) :: string
  integer(i4b) :: lmax,lmin,l,narg
  real(dp) :: rtol
  real(dp), dimension(:,:), allocatable :: point_loc
  real(dp), dimension(:), allocatable :: u
  

  narg = command_argument_count()
  if(narg == 0 .or. narg /= 5) then
     print *, '[lmin] [lmax] [s] [mu] [rtol] '
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
  
  ! set the function to optimise
  fun1 => point_misfit_ne_S2
  
  ! set the line search method
  call setup_cvsrch(n_S2_point,1.0e-6_dp,1.e-6_dp,10.0_dp)
  lins => cvsrch

  ! set the inner product
  iprod => sobolev_product_S2

  ! set the preconditioner
  precon => preid

  ! solve the optimisation problem
  u = 0.0_dp
  call lbfgs(u,5,verb = .true.,rtol = rtol)
  call plot_function_S2('min.mod',360,720,u)


  
end program S2_point_mnne
