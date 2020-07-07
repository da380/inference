program delta_test

  use nrtype
  use module_util
  use module_function

  implicit none

  integer(i4b), parameter :: np = 10
  integer(i4b) :: l,nth,ith,lmax,io,ip
  real(dp) :: th1,th2,dth,th,s,mu,mu1,mu2,s1,s2,ds,dmu,ft,fac
  real(dp), dimension(2) :: x,xc,xp
  real(dp), dimension(np) :: f

  ! set plot parameters
  lmax = 2000
  nth = 1001
  th1 = -0.01_dp
  th2 =  0.01_dp
  dth =(th2-th1)/(nth-1)

  ! set Sobolev parameters
  s1  = 1.1_dp
  s2  = 2.1_dp
  mu1 = 0.1_dp
  mu2 = 0.1_dp

  ds  = (s2-s1)/(np-1)
  dmu = (mu2-mu1)/(np-1)


  open(newunit = io,file='delta.out')
  
  ! loop over angle
  do ith = 1,nth

     ! set th value
     th = th1+(ith-1)*dth

     ! initialise the sum
     f = 0.0_dp

     ! loop over l
     do l = 0,lmax

        ! compute Legendre polynomials
        call legendre(th,l,1,x,xp,xc)
        
        do ip = 1,np

           mu = mu1+(ip-1)*dmu
           s  = s1+(ip-1)*ds

           
           ! compute pre-factor
           fac = 1.0+mu**2*l*(l+1)
           fac = sqrt((2*l+1)/fourpi_d)*fac**(-s)

           ! sum up the series
           f(ip) = f(ip) + fac*x(1)

        end do
           
     end do

     write(io,*) th,f
     
  end do

  close(io)

end program delta_test
