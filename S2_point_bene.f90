program S2_point_bene


  use nrtype
  use module_S2_point
  use module_gsht
  use module_opt
  implicit none


  integer(i4b) :: lmax,lmin,l,lt,iparm,io
  real(dp) :: rtol,utn
  real(dp), dimension(:,:), allocatable :: point_loc
  real(dp), dimension(:), allocatable :: u,ut,w,wt
  real(dp), dimension(:,:), allocatable :: bb,bb0
  

  !==========================================!
  !     set up the model and data spaces     !
  !==========================================!

  
  ! set up the GSHT routines
  lmin =  0
  lmax = 64
  call setup_gsht(lmax = lmax,nmax = 0)


  ! set the sobolev parameters
  s_S2_point =  2.0_dp
  mu_S2_point = 0.2_dp


  ! read in the data
  call read_point_data('syn.dat',m_S2_point,point_loc,v_S2_point)
  allocate(x_S2_point(m_S2_point,2))
  x_S2_point(:,1) = (90.0_dp-point_loc(:,1))*deg2rad
  x_S2_point(:,2) = point_loc(:,2)*deg2rad

  ! set dimension of model vector 
  n_S2_point = nth_gsht*nph_gsht
  allocate(u(n_S2_point))
  allocate(ut(n_S2_point))


  ! setup the mask for the constrained coefficients
  allocate(mask_S2_point(ncoef_gsht))
  mask_S2_point(:) = .true.
  do l = 0,lmin-1
     mask_S2_point(ilm(l,-l):ilm(l,l)) = .false.
  end do

  !========================================!
  !        set up the property space       !
  !========================================!

  
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
  allocate(wt(nparm_S2_point))


  !================================================================!
  !                 set up the optimization codes                  !
  !================================================================!

  ! set dimension of model vector 
  n_S2_point = nth_gsht*nph_gsht

  ! set the function to optimise
  fun1 => point_misfit_ne_S2
  
  ! set the line search method
  call setup_cvsrch(n_S2_point,1.0e-6_dp,1.e-6_dp,10.0_dp)
  lins => cvsrch

  ! set the inner product
  iprod => sobolev_product_S2

  ! set the preconditioner
  precon => preid

  ! set the tolerance for optimisations
  rtol = 1.e-4_dp
  

  !=========================================!
  !      solve the minimum norm problem     !
  !=========================================!

  ! solve the optimisation problem
  ut = 0.0_dp
  call lbfgs(ut,5,verb = .true.,rtol = rtol)

  ! get the minimum norm value
  call sobolev_product_S2(ut,ut,utn)

  ! get the minimum norm property vector
  call spharm_operator_S2(ut,wt)


  !=========================================!
  !   compute  BB^{*} and B_{0}B_{0}^{*}    !
  !=========================================!

  allocate(bb(nparm_S2_point,nparm_S2_point))
  allocate(bb0(nparm_S2_point,nparm_S2_point))

  
  do iparm = 1,nparm_S2_point

     ! compute action of B^{*} on a basis vector
     w = 0.0_dp
     w(iparm) = 1.0_dp
     call adjoint_spharm_operator_S2(w,u)

     ! compute the action of B on the result
     call spharm_operator_S2(u,bb0(:,iparm))
     
     ! compute the associated data
     call point_operator_S2(u,v_S2_point(:,1))

     ! solve the minimum norm problem
     ut = 0.0_dp
     call lbfgs(ut,5,verb = .true.,rtol = rtol)
     
     ! compute the projection
     u = u-ut

     ! compute the action of B on the result
     call spharm_operator_S2(u,bb(:,iparm))
     
  end do


  ! write out the results for plotting 
  open(newunit = io,file='bene.dat')
  write(io,*) utn,wt(2:)*0.0_dp
  do iparm = 1,nparm_S2_point
     write(io,*) bb0(iparm,:)
  end do
  write(io,*) wt
  do iparm = 1,nparm_S2_point
     write(io,*) bb(iparm,:)
  end do
  close(io)

end program S2_point_bene
