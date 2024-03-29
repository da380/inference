program S2_point_data

  use nrtype
  use module_S2_point
  use module_maps
  use ran_state
  use module_gsht
  use module_interp
  use module_crust
  implicit none
  

  character(len = 256) :: efile,pref,string
  
  integer(i4b) :: np,lmax,io,l,m,i,nth,lmin,nsamp,narg,aopt

  real(dp) :: sigma1,sigma2,ran,pow,t,lam,un,up,alpha, &
              th,dth,cor,aval
  
  real(dp), dimension(:,:), allocatable :: point_loc
  real(dp), dimension(:,:), allocatable :: point_dat

  real(dp), dimension(:), allocatable :: u
  real(dp), dimension(:), allocatable :: v
  real(dp), dimension(:), allocatable :: ulm



  !----------------------------------------------------!
  !                    get arguments                   !
  !----------------------------------------------------!

  narg = command_argument_count()
  if(narg == 0 .or. narg /= 12) then
     print *, '[lmin] [lmax] [s] [mu] [t] [lam] [aopt] [aval] [nsamp] [npoint] [sigma1] [sigma2] '
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
  read(string,*) t
  call get_command_argument(6,string)
  read(string,*) lam
  call get_command_argument(7,string)
  read(string,*) aopt
  call get_command_argument(8,string)
  read(string,*) aval
  call get_command_argument(9,string)
  read(string,*) nsamp
  call get_command_argument(10,string)
  read(string,*) m_S2_point
  call get_command_argument(11,string)
  read(string,*) sigma1
  call get_command_argument(12,string)
  read(string,*) sigma2

  

  !====================================================!
  !          set the model and synthetic data          !
  !====================================================!


  ! set up the GSHT routines
  call setup_gsht(lmax = lmax,nmax = 0)
  
  ! set up random numbers
  call setup_random()

  ! set up the crustal model
  call get_crust

  ! set the point locations
  call random_locations_topo(m_S2_point,point_loc,sign = -1)
  allocate(x_S2_point(m_S2_point,2))
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


  ! make random model
  if(aopt == 1) then     
     alpha = aval
  else if(aopt == 2) then
     un = aval
  else if(aopt == 3) then
     up = aval
  else
     stop 'S2_point_data: bad input for aopt'
  end if
  allocate(qq_S2_point(lmax+1))  
  call covariance_parm(lmin,lmax,alpha,t,lam,qq_S2_point,up = up)
  call random_model(lmin,lmax,qq_S2_point,u)
  call  sobolev_product_S2(u,u,un)

  
  ! write out the two-point correlation function
  nth = 180
  dth = pi_d/(nth-1)
  open(newunit = io,file='tpc.out')
  do i = 1,nth
     th = (i-1)*dth
     call  point_correlation_S2(lmin,lmax,qq_S2_point,th,cor)
     write(io,*) th,cor
  end do
  close(io)
  

  ! write out model for plotting
  call plot_function_S2('syn.mod',360,720,u)

  
  ! get the spectral coefficients and write out to file
  allocate(ulm(ncoef_gsht))
  call real_coefs_from_model(u,ulm)

  open(newunit = io,file='syn.spec')
  do l = 0,lmax
     do m = -l,l
        write(io,*) l,m,ulm(ilm(l,m))
     end do
  end do
  close(io)

  ! write out the power spectrum to file
  open(newunit = io,file='syn.pow')
  do l = 0,lmax
     pow = 0.0_dp
     do m = -l,l
        pow = pow + ulm(ilm(l,m))**2
     end do
     write(io,*) l,pow/(2*l+1)
  end do
  close(io)

  ! make synthetic data
  allocate(v(m_S2_point))
  call point_operator_S2(u,v_S2_point(:,1))
  
  ! set the standard deviations  
  do i = 1,m_S2_point
     call random_number(ran)     
     v_S2_point(i,2) =  ran*sigma1+(1.0_dp-ran)*sigma2
  end do

  ! add in random errors to synthetic data
  do i = 1,m_S2_point
     call gasdev(ran)
     v_S2_point(i,1) =  v_S2_point(i,1) + ran*v_S2_point(i,2)     
  end do

  ! write out the data
  call write_point_data('syn.dat',m_S2_point,point_loc,v_S2_point)

  ! write out the data using map projection
  call write_point_data('syn.dat.proj',m_S2_point,point_loc,v_S2_point,proj = .true.)  

  
  ! sample from the model distribution
  if(nsamp > 0) then
     open(newunit = io,file='syn.mod.samp')
     do i = 1,nsamp
        print *, 'doing sample ',i
        call random_model(lmin,lmax,qq_S2_point,u)
        call  sobolev_product_S2(u,u,un)
        call point_operator_S2(u,v_S2_point(:,1))
        write(io,*) un,v_S2_point(:,1)
     end do
     close(io)
  end if
     
  



  
end program S2_point_data
