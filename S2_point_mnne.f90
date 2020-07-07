program S2_point_mnne

  use nrtype
  use module_S2_point
  use module_maps
  use module_gsht
  use module_interp
  use module_opt
  implicit none


  integer(i4b) :: lmax,lmin,l
  real(dp) :: rtol
  real(dp), dimension(:,:), allocatable :: point_loc
  real(dp), dimension(:), allocatable :: u
  

  ! set up the GSHT routines
  lmin =  0
  lmax = 128
  call setup_gsht(lmax = lmax,nmax = 0)


  ! set the sobolev parameters
  s_S2_point =  2.0_dp
  mu_S2_point = 0.25_dp


  ! read in the data
  call read_point_data('syn.dat',m_S2_point,point_loc,v_S2_point)
  allocate(x_S2_point(m_S2_point,2))
  x_S2_point(:,1) = (90.0_dp-point_loc(:,1))*deg2rad
  x_S2_point(:,2) = point_loc(:,2)*deg2rad


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


  ! setup the mask for the constrained coefficients
  allocate(mask_S2_point(ncoef_gsht))
  mask_S2_point(:) = .true.
  do l = 0,lmin-1
     mask_S2_point(ilm(l,-l):ilm(l,l)) = .false.
  end do
  
  ! solve the optimisation problem
  allocate(u(n_S2_point))
  u = 0.0_dp
  call lbfgs(u,5,verb = .true.,rtol = rtol)
  call plot_function_S2('mnne.mod',360,720,u)


  
end program S2_point_mnne
