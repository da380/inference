module module_S2_point

  use nrtype
  implicit none

  integer(i4b), save :: n_S2_point
  integer(i4b), save :: m_S2_point

  logical(lgt), dimension(:), allocatable, save :: mask_S2_point
  
  real(dp), save :: s_S2_point
  real(dp), save :: mu_S2_point
  real(dp), save :: ss_S2_point
  real(dp), save :: eta_S2_point

  real(dp), dimension(:,:), allocatable, save :: x_S2_point
  real(dp), dimension(:,:), allocatable, save :: v_S2_point
  real(dp), dimension(:), allocatable, save :: up_S2_point
  real(dp), dimension(:), allocatable, save :: qq_S2_point

  integer(i4b), save :: nparm_S2_point
  integer(i4b), dimension(:), allocatable, save :: la_S2_point
  integer(i4b), dimension(:), allocatable, save :: ma_S2_point
  
contains

  !==============================================================!
  !             spherical harmonic operator routines             !
  !==============================================================!

  subroutine spharm_operator_S2(u,w)
    use nrtype
    use module_gsht
    use module_interp
    implicit none

    real(dp), dimension(:), intent(in) :: u
    real(dp), dimension(:), intent(out) :: w

    integer(i4b) :: iparm
    real(dp), dimension(ncoef_gsht) :: ulm

    call real_coefs_from_model(u,ulm)
    do iparm = 1,nparm_S2_point
       w(iparm) = ulm(ilm(la_S2_point(iparm),ma_S2_point(iparm)))
    end do
    
  end subroutine spharm_operator_S2

  subroutine adjoint_spharm_operator_S2(w,u)
    use nrtype
    use module_gsht
    use module_function
    implicit none

    real(dp), dimension(:), intent(in) :: w
    real(dp), dimension(:), intent(out) :: u

    integer(i4b) :: iparm,l,m,k,ith,iph
    real(dp) :: fac
    real(dp), dimension(ncoef_gsht) :: ulmr
    
    complex(dpc), dimension(ncoef_gsht) :: ulm
    complex(dpc), dimension(ncoef_gsht) :: ulmt
    complex(dpc), dimension(nth_gsht,nph_gsht) :: uc
    
    
    ! initialise the coefficients
    ulmr = 0.0_dp
    do iparm = 1,nparm_S2_point

       ! set the degree
       l = la_S2_point(iparm)
       
       ! set the Sobolev factor
       fac = (1.0_dp + mu_S2_point**2*l*(l+1))**s_S2_point       
       fac = 1.0_dp/fac

       ! set the coefficient
       ulmr(ilm(la_S2_point(iparm),ma_S2_point(iparm))) = w(iparm)*fac

       ! project orthogonal to constraints
       do l = 0,lmax_gsht
          do m = -l,l
             if(.not.mask_S2_point(ilm(l,m))) then
                ulmr(ilm(l,m)) = 0.0_dp
             end if
          end do
       end do

       ! convert to complex coefficients
       call real_to_complex_SH(ulmr,ulm)
       
       ! perform the inverse spherical harmonic transform
       call igsht(0,ulm,uc)
       
       ! package the function as a vector
       k = 0
       do ith = 1,nth_gsht
          do iph = 1,nph_gsht
             k = k+1
             u(k) = real(uc(ith,iph))
          end do
       end do
       
    end do


  end subroutine adjoint_spharm_operator_S2
    
  !==============================================================!
  !                    point operator routines                   !
  !==============================================================!



  
  subroutine point_operator_S2(u,v)
    use nrtype
    use module_gsht
    use module_interp
    implicit none

    real(dp), dimension(:), intent(in) :: u
    real(dp), dimension(:), intent(out) :: v

    integer(i4b) :: ith,iph,k,ip
    real(dp) :: th,ph,ue
    real(dp), dimension(nth_gsht,nph_gsht) :: ua

    ! unpack the function vector
    k = 0
    do ith = 1,nth_gsht
       do iph = 1,nph_gsht
          k = k+1
          ua(ith,iph) = u(k)
       end do
    end do

    ! loop over the points evaluating the function
    do ip = 1,m_S2_point
       th = x_S2_point(ip,1)
       ph = x_S2_point(ip,2)
       call polin2_sub(th_gsht,ph_gsht,ua,th,ph,5,5,v(ip),ue)
    end do
    
    
    return
  end subroutine point_operator_S2



  subroutine adjoint_point_operator_S2(v,u)
    use nrtype
    use module_gsht
    use module_function
    implicit none

    real(dp), dimension(:), intent(in) :: v
    real(dp), dimension(:), intent(out) :: u

    integer(i4b) :: ip,l,k,ith,iph,m
    real(dp) :: th,ph,fac
    real(dp), dimension(ncoef_gsht) :: ulmr
    
    complex(dpc), dimension(ncoef_gsht) :: ulm
    complex(dpc), dimension(ncoef_gsht) :: ulmt
    complex(dpc), dimension(nth_gsht,nph_gsht) :: uc
    

    ! add up the delta function coefficients
    ulm = 0.0_dp
    do ip = 1,m_S2_point

       th = x_S2_point(ip,1)
       ph = x_S2_point(ip,2)

       call ylm_array(lmax_gsht,th,ph,ulmt)

       ! add contribution to the kernel
       ulm = ulm + v(ip)*conjg(ulmt)
       
    end do
    
    ! do the sobolev processing
    do l = 0,lmax_gsht
       fac = (1.0_dp + mu_S2_point**2*l*(l+1))**s_S2_point
       fac = 1.0_dp/fac       
       ulm(ilm(l,-l):ilm(l,l)) = fac*ulm(ilm(l,-l):ilm(l,l))
    end do
    
    ! project orthogonal to constraints
    call complex_to_real_SH(ulm,ulmr)
    do l = 0,lmax_gsht
       do m = -l,l
          if(.not.mask_S2_point(ilm(l,m))) then
             ulmr(ilm(l,m)) = 0.0_dp
          end if
       end do
    end do
    call real_to_complex_SH(ulmr,ulm)
    
    ! perform the inverse spherical harmonic transform
    call igsht(0,ulm,uc)

    ! package the function as a vector
    k = 0
    do ith = 1,nth_gsht
       do iph = 1,nph_gsht
          k = k+1
          u(k) = real(uc(ith,iph))
       end do
    end do
    
    return
  end subroutine adjoint_point_operator_S2


  subroutine point_misfit_ne_S2(u,J,dJ)
    use nrtype
    implicit none
    real(dp), dimension(:), intent(in) :: u
    real(dp), intent(out) :: J
    real(dp), dimension(:), intent(out), optional :: dJ

    real(dp), dimension(m_S2_point) :: v


    call point_operator_S2(u,v)
    v = v-v_S2_point(:,1)
    J = dot_product(v,v)

    if(present(dJ)) then
       call adjoint_point_operator_S2(2.0_dp*v,dJ)       
    end if
    

    return
  end subroutine point_misfit_ne_S2


  subroutine point_misfit_we_S2(u,J,dJ)
    use nrtype
    implicit none
    real(dp), dimension(:), intent(in) :: u
    real(dp), intent(out) :: J
    real(dp), dimension(:), intent(out), optional :: dJ
    
    real(dp) :: un2
    real(dp), dimension(m_S2_point) :: v,vt

    call sobolev_product_S2(u,u,un2)
    call point_operator_S2(u,v)
    v = v_S2_point(:,1) - v
    vt = v/v_S2_point(:,2)**2
    J = 0.5_dp*dot_product(vt,v)    
    J = 0.5_dp*un2 + eta_S2_point*(J-ss_S2_point)
    
    
    if(present(dJ)) then
       call adjoint_point_operator_S2(vt,dJ)
       dJ = u - eta_S2_point*dJ
    end if
    

    return
  end subroutine point_misfit_we_S2


  subroutine point_misfit_we_con_S2(u,Jc)
    use nrtype
    implicit none
    real(dp), dimension(:), intent(in) :: u
    real(dp), intent(out) :: Jc
   
    real(dp), dimension(m_S2_point) :: v,vt


    call point_operator_S2(u,v)
    v = v_S2_point(:,1) - v
    vt = v/v_S2_point(:,2)**2
    Jc = 0.5_dp*dot_product(vt,v)-ss_S2_point

    return
  end subroutine point_misfit_we_con_S2


  subroutine covariance_operator_S2(u,uq)
    use nrtype
    use module_gsht
    implicit none
    real(dp), dimension(:), intent(in) :: u
    real(dp), dimension(:), intent(out) :: uq

    integer(i4b) :: ith,iph,k,l
    real(dp) :: th,ph,ue
    complex(dpc), dimension(nth_gsht,nph_gsht) :: ua
    complex(dpc), dimension(ncoef_gsht) :: ulm

    ! unpack the function vector
    k = 0
    do ith = 1,nth_gsht
       do iph = 1,nph_gsht
          k = k+1
          ua(ith,iph) = u(k)
       end do
    end do

    ! get spherical harmonic coefficients
    call gsht(0,ua,ulm)

    ! act the covariance operator
    do l = 0,lmax_gsht
       ulm(ilm(l,-l):ilm(l,l)) = qq_S2_point(l+1)*ulm(ilm(l,-l):ilm(l,l)) 
    end do


    ! pass back into the spatial domain
    call igsht(0,ulm,ua)

    ! repack the function
    k = 0
    do ith = 1,nth_gsht
       do iph = 1,nph_gsht
          k = k+1
          uq(k) = real(ua(ith,iph))
       end do
    end do

    return
  end subroutine covariance_operator_S2
  
  
  !==============================================================!
  !                 routines for model generation                !
  !==============================================================!

  subroutine model_from_real_coefs(ulm,u)
    use nrtype
    use module_gsht
    implicit none
    real(dp), dimension(:), intent(in) :: ulm
    real(dp), dimension(:), intent(out) :: u

    integer(i4b) :: k,ith,iph,l,m    
    complex(dpc), dimension(nth_gsht,nph_gsht) :: uc
    complex(dpc), dimension(ncoef_gsht) :: ulm_tmp

    ! convert to complex coefficients
    call real_to_complex_SH(ulm,ulm_tmp)

    ! set the model
    call model_from_coefs(ulm_tmp,u)
    
    return
  end subroutine model_from_real_coefs
  
  subroutine model_from_coefs(ulm,u)
    use nrtype
    use module_gsht
    implicit none
    complex(dpc), dimension(:), intent(in) :: ulm
    real(dp), dimension(:), intent(out) :: u

    integer(i4b) :: k,ith,iph,l,m
    complex(dpc), dimension(nth_gsht,nph_gsht) :: uc

    call igsht(0,ulm,uc)
    k = 0
    do ith = 1,nth_gsht
       do iph = 1,nph_gsht
          k = k+1
          u(k) = real(uc(ith,iph))
       end do
    end do
    
    return
  end subroutine model_from_coefs


  subroutine random_model(lmin,lmax,qq,u)
    use nrtype
    use ran_state
    use module_gsht
    implicit none
    integer(i4b), intent(in) :: lmin
    integer(i4b), intent(in) :: lmax
    real(dp), dimension(lmax+1), intent(in) :: qq
    real(dp), dimension(:), intent(out) :: u

    complex(dpc), dimension(ncoef_gsht) :: ulm
    complex(dpc), dimension(nth_gsht,nph_gsht) :: uc
    
    integer(i4b) :: l,m,k,ith,iph
    real(dp) :: rana,ranp,fac1,fac2,norm,ustd


    
    ulm(:) = 0.0_dp
    do l = lmin,lmax

       ! get the standard deviation for the lth coefficients
       fac1 = 1.0_dp/((1.0_dp + (mu_S2_point**2)*l*(l+1))**s_S2_point)
       fac2 = qq(l+1)
       ustd = sqrt(fac1*fac2)
       
       
       ! do m = 0 coefficient
       call gasdev(rana)
       ulm(ilm(l,0)) = rana*ustd

       ! do m /= 0 coefficients
       do m = 1,l

          call gasdev(rana)
          call ran1(ranp)
          
          ! positive m
          ulm(ilm(l,m)) = rana*ustd*exp(ii*twopi_d*ranp)

          ! negative m
          ulm(ilm(l,-m)) = (-1)**m*conjg(ulm(ilm(l,m)))
          
       end do

    end do

    ! perform the inverse spherical harmonic transform
    call igsht(0,ulm,uc)

    k = 0
    do ith = 1,nth_gsht
       do iph = 1,nph_gsht
          k = k+1
          u(k) = real(uc(ith,iph))
       end do
    end do
    
    
    return
  end subroutine random_model


  subroutine covariance_parm(lmin,lmax,alpha,t,lam,qq,un,up)
    ! produces a covariance operator of the form
    !
    ! Q = \sum_{l\in \mathbb{N}} q_{l} \mathmm{P}_{l}
    !
    ! where the q_{l} have the parametric form
    !
    ! q_{l} = \alpha \mult{l}_{\lambda}^{-t}
    !
    ! where t > 2 in order for the trace to be finite
    !
    ! optional argument un given if \alpha is to be determined   
    ! through the requirement that:
    !
    ! \tr Q = un**2
    !
    ! optional argument up given if \alpha is to be determined 
    ! through the requirement that:
    !
    ! Var[u(x)] = up**2
    !
    ! for each x \in \mathbb{S}^{2}
    !
    use nrtype
    use ran_state
    implicit none
    integer(i4b), intent(in) :: lmin
    integer(i4b), intent(in) :: lmax
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: t
    real(dp), intent(in) :: lam
    real(dp), dimension(lmax+1)  :: qq
    real(dp), intent(in), optional  :: un
    real(dp), intent(in), optional  :: up

    integer(i4b) :: l
    real(dp) :: fac,sum


    ! check value of t
    if(t <= 2.0_dp) stop ' covariance_given_norm: t is too small!'
    

    ! check optional arguments do not conflict
    if(present(un) .and. present(up)) stop 'covariance_given_norm: only one optional argument!'

    ! initialise the covariance
    qq = 0.0_dp

    ! set the unnormalised version of the coefficients
    do l = lmin,lmax
       fac = (1.0_dp+lam**2*l*(l+1))**(-0.5_dp*t)
       qq(l+1) = fac
    end do

    
    
    ! set the normalisation
    if(present(un)) then

       ! compute the trace
       sum = 0.0_dp
       do l = lmin,lmax
          sum = sum + (2*l+1)*qq(l+1)
       end do

       ! normalise the terms
       qq = un**2*qq/sum

    else if(present(up)) then
       
       ! compute the variance at a point
       sum = 0.0_dp
       do l = lmin,lmax
          fac = (1.0_dp+mu_S2_point**2*l*(l+1))**(-s_S2_point)
          sum = sum  + fac*(2*l+1)/fourpi_d*qq(l+1)
       end do

       qq = up**2*qq/sum
       
    else
       qq = alpha*qq
    end if
    
    
    
    return
  end subroutine covariance_parm




  
  !==============================================================!
  !            routines for data reading or generation           !
  !==============================================================!
  
  subroutine hoggard_locations(np,points,data)
    use nrtype
    implicit none    
    integer(i4b), intent(out) :: np
    real(dp), dimension(:,:), intent(inout), allocatable :: points
    real(dp), dimension(:,:), intent(inout), allocatable, optional :: data
    
    integer(i4b) :: ios,ip,io
    real(dp) :: lon,lat,dat,err

    ! read in the data file
    open(newunit = io, file='/home/da380/raid1/codes/work/hoggard_2017_data/nst.llze')
    np = 0
    do 
       read(io,*,iostat=ios) lat,lon,dat,err
       if(ios /= 0) exit
       np = np+1
    end do
    
    rewind(io)
    allocate(points(np,2))
    if(present(data)) allocate(data(np,2))
    do ip = 1,np
       read(io,*,iostat=ios) lat,lon,dat,err
       points(ip,1) = lat
       if(lon < 0.0_dp) lon = lon+360.0_dp
       points(ip,2) = lon
       if(present(data)) then
          data(ip,1) = dat
          data(ip,2) = err          
       end if
    end do
    close(io)
    
    return
  end subroutine hoggard_locations



  subroutine random_locations(np,points)
    use nrtype
    use ran_state
    implicit none    
    integer(i4b), intent(in) :: np
    real(dp), dimension(:,:), intent(inout), allocatable :: points

    integer(i4b) :: ios,ip,io
    real(dp) :: lon,lat,ran

    allocate(points(np,2))
    ip = 0
    do 

       call ran1(ran)
       lat = 2.0_dp*(ran-0.5_dp)*85.0_dp
       call ran1(ran)
       lon =  ran*358.0_dp
       
       ip = ip+1
       points(ip,1) = lat 
       points(ip,2) = lon

       if(ip == np) exit
       
    end do
    
    
    return
  end subroutine random_locations


  subroutine random_locations_topo(np,points,sign,amp)
    use nrtype    
    use ran_state
    use module_crust
    use module_maps
    implicit none    
    integer(i4b), intent(in) :: np
    real(dp), dimension(:,:), intent(inout), allocatable :: points
    integer(i4b), intent(in) :: sign
    real(dp), intent(in), optional :: amp

    integer(i4b) :: ios,ip,io
    real(dp) :: lon,lat,ran,lont,topo,r,th,ph
    real(dp), dimension(3) :: x
    real(dp) :: amp_loc

    if(present(amp)) then
       amp_loc = amp
    else
       amp_loc = 0.0_dp
    end if


    allocate(points(np,2))
    ip = 0
    do 

       call random_sample_SN(2,x)
       call cart_2_sph(x(1),x(2),x(3),r,th,ph)

       lat = (0.5_dp*pi_d-th)*rad2deg
       lon = ph*rad2deg
       if(lon < 0.0_dp) lon = 360.0_dp+lon

       call get_topo(lat,lon,topo)
       if(sign /= 0) then
          if(sign*topo < amp_loc) cycle
       end if
       
       ip = ip+1
       points(ip,1) = lat 
       points(ip,2) = lon
       

       if(ip == np) exit
       
    end do
    
    
    return
  end subroutine random_locations_topo



  !==============================================================!
  !                     routines for plotting                    !
  !==============================================================!

  subroutine plot_function_S2(ufile,nth,nph,u)
    use nrtype
    use module_gsht
    use module_maps
    use module_interp
    implicit none
    character(len=*), intent(in) :: ufile
    integer(i4b), intent(in) :: nth
    integer(i4b), intent(in) :: nph
    real(dp), dimension(:), intent(in) :: u


    integer(i4b) :: io,ith,iph,k
    real(dp) :: lat,lon,x,y,dth,dph,th,ph,ui,uie,th1,th2
    real(dp), dimension(nth_gsht,nph_gsht) :: ua

    k = 0
    do ith = 1,nth_gsht
       do iph = 1,nph_gsht
          k = k+1
          ua(ith,iph) = u(k)
       end do
    end do

    th1 = th_gsht(1)
    th2 = th_gsht(nth_gsht)
    dth = (th2-th1)/(nth-1)
    dph = twopi_d/(nph-1)
    
    open(newunit = io,file = trim(ufile))
    write(io,*) nth,nph,0
    do ith = 1,nth
       th = th1+(ith-1)*dth
       lat = 0.5_dp*pi_d-th
       do iph = 1,nph
          ph = (iph-1)*dph
          lon = ph
          call pr_hammer_aitoff(lat,lon,x,y)
          call polin2_sub(th_gsht,ph_gsht,ua,th,ph,5,5,ui,uie)
          write(io,*) x,y,ui
       end do
    end do
    close(io)
    
    return
  end subroutine plot_function_S2

  
  subroutine write_points(file_name,n,x)
    use nrtype
    use module_maps
    implicit none
    character(len=*), intent(in) :: file_name
    integer(i4b), intent(in) :: n
    real(dp), dimension(:,:), intent(in) :: x

    integer(i4b) :: i,io
    real(dp) :: xp,yp

    open(newunit = io, file=trim(file_name))
    do i = 1,n
       call pr_hammer_aitoff(deg2rad*x(i,1), & 
                             deg2rad*x(i,2),xp,yp)
       write(io,*) xp,yp
    end do
    close(io)    

    return
  end subroutine write_points


  subroutine write_point_data(file_name,n,x,data,proj)
    use nrtype
    use module_maps
    implicit none
    character(len=*), intent(in) :: file_name
    integer(i4b), intent(in) :: n
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:,:), intent(in) :: data
    logical(lgt), intent(in), optional :: proj

    logical(lgt) :: proj_loc
    integer(i4b) :: i,io
    real(dp) :: xp,yp

    if(present(proj)) then
       proj_loc = proj
    else
       proj_loc = .false.
    end if

    open(newunit = io, file=trim(file_name))
    do i = 1,n
       if(proj_loc) then
          call pr_hammer_aitoff(deg2rad*x(i,1), & 
                                deg2rad*x(i,2),xp,yp)
          write(io,*) xp,yp,data(i,1),data(i,2)
       else
          write(io,*) x(i,1),x(i,2),data(i,1),data(i,2)
       end if
    end do
    close(io)    

    return
  end subroutine write_point_data

  
  subroutine read_point_data(file_name,n,x,data)
    use nrtype
    use module_maps
    implicit none
    character(len=*), intent(in) :: file_name
    integer(i4b), intent(out) :: n
    real(dp), dimension(:,:), allocatable, intent(inout) :: x
    real(dp), dimension(:,:), allocatable, intent(inout) :: data

    integer(i4b) :: i,io,ios
    real(dp) :: xp,yp

    open(newunit = io, file=trim(file_name))
    n = 0
    do 
       read(io,*,iostat = ios) xp
       if(ios /= 0) exit
       n = n+1
    end do
    rewind(io)
    allocate(x(n,2),data(n,2))
    do i = 1,n
       read(io,*) x(i,1),x(i,2),data(i,1),data(i,2)
    end do
    close(io)    

    return
  end subroutine read_point_data


  !==============================================================!
  !                        utility routines                      !
  !==============================================================!

  

  

  subroutine sobolev_product_S2(u1,u2,prod)
    use nrtype
    use module_gsht
    real(dp), dimension(:), intent(in) :: u1
    real(dp), dimension(:), intent(in) :: u2
    real(dp), intent(out) :: prod

    integer(i4b) :: ith,iph,k,l
    real(dp) :: fac
    complex(dpc), dimension(nth_gsht,nph_gsht) :: u1c,u2c
    complex(dpc), dimension(ncoef_gsht) :: ulm1,ulm2

    k = 0
    do ith = 1,nth_gsht
       do iph = 1,nph_gsht
          k = k+1
          u1c(ith,iph) = u1(k)
          u2c(ith,iph) = u2(k)
       end do
    end do

    ! perform the spherical harmonic transforms
    call gsht(0,u1c,ulm1)
    call gsht(0,u2c,ulm2)

    prod = 0.0_dp
    do l = 0,lmax_gsht

       ! set the sobolev factor
       fac = (1.0_dp + mu_S2_point**2*l*(l+1))**s_S2_point

       ! add contribution to the inner product
       prod = prod + fac*real(dot_product(ulm1(ilm(l,-l):ilm(l,l)),ulm2(ilm(l,-l):ilm(l,l))))
       
    end do

    
    return
  end subroutine sobolev_product_S2

  subroutine real_coefs_from_model(u,ulm)
    use nrtype
    use module_gsht
    implicit none
    real(dp), dimension(:), intent(in) :: u
    real(dp), dimension(:), intent(out) :: ulm

    complex(dpc), dimension(ncoef_gsht) :: ulmc

    call coefs_from_model(u,ulmc)
    call complex_to_real_SH(ulmc,ulm)

  end subroutine real_coefs_from_model
    
  subroutine coefs_from_model(u,ulm)
    use nrtype
    use module_gsht
    implicit none
    real(dp), dimension(:), intent(in) :: u
    complex(dpc), dimension(:), intent(out) :: ulm

    integer(i4b) :: k,ith,iph,l,m
    complex(dpc), dimension(nth_gsht,nph_gsht) :: uc

    ! unpack function vector
    k = 0
    do ith = 1,nth_gsht
       do iph = 1,nph_gsht
          k = k+1
          uc(ith,iph) = u(k)
       end do
    end do
    call gsht(0,uc,ulm)
    
    return
  end subroutine coefs_from_model


  subroutine quad_coefs_1D(w0,f0,w,f,a)
    use nrtype
    implicit none
    real(dp), intent(in) :: w0
    real(dp), intent(in) :: f0
    real(dp), intent(in) :: w
    real(dp), intent(in) :: f
    real(dp), intent(out) :: a
    a = (f-f0)/(w-w0)**2    
    return
  end subroutine quad_coefs_1D

  function quad_fun_1D(w0,f0,a,w)
    use nrtype
    implicit none
    real(dp) :: quad_fun_1D
    real(dp), intent(in) :: w0
    real(dp), intent(in) :: f0
    real(dp), intent(in) :: a
    real(dp), intent(in) :: w
    quad_fun_1D = f0 + a*(w-w0)**2
    return
  end function quad_fun_1D

  subroutine quad_coefs_2D(w0,f0,w1,f1,w2,f2,w3,f3,a,b,c)
    use nrtype
    implicit none
    real(dp), dimension(2), intent(in) :: w0
    real(dp), intent(in) :: f0
    real(dp), dimension(2), intent(in) :: w1
    real(dp), intent(in) :: f1
    real(dp), dimension(2), intent(in) :: w2
    real(dp), intent(in) :: f2
    real(dp), dimension(2), intent(in) :: w3
    real(dp), intent(in) :: f3
    real(dp), intent(out) :: a
    real(dp), intent(out) :: b
    real(dp), intent(out) :: c

    integer(i4b) :: info
    integer(i4b), dimension(3) :: ipiv
    real(dp), dimension(3,3) :: mat
    real(dp), dimension(3,1) :: vec

    ! set the right hand side
    vec(1,1) = f1-f0
    vec(2,1) = f2-f0
    vec(3,1) = f3-f0

    ! set the system matrix
    mat(1,1) = (w1(1)-w0(1))**2
    mat(1,2) = 2.0_dp*(w1(1)-w0(1))*(w1(2)-w0(2))
    mat(1,3) = (w1(2)-w0(2))**2
    mat(2,1) = (w2(1)-w0(1))**2
    mat(2,2) = 2.0_dp*(w2(1)-w0(1))*(w2(2)-w0(2))
    mat(2,3) = (w2(2)-w0(2))**2
    mat(3,1) = (w3(1)-w0(1))**2
    mat(3,2) = 2.0_dp*(w3(1)-w0(1))*(w3(2)-w0(2))
    mat(3,3) = (w3(2)-w0(2))**2


    ! solve the linear system
    call dgesv(3,1,mat,3,ipiv,vec,3,info)

    ! set the coefficients
    a = vec(1,1)
    b = vec(2,1)
    c = vec(3,1)
    
    return
  end subroutine quad_coefs_2D

  function quad_fun_2D(w0,f0,a,b,c,w)
    use nrtype
    implicit none
    real(dp) :: quad_fun_2D
    real(dp), dimension(2), intent(in) :: w0
    real(dp), intent(in) :: f0
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp), intent(in) :: c
    real(dp), dimension(2), intent(in) :: w
    quad_fun_2D = f0 + a*(w(1)-w0(1))**2                &
                     + 2.0_dp*b*(w(1)-w0(1))*(w(2)-w0(2)) &
                     + c*(w(2)-w0(2))**2                
    return
  end function quad_fun_2D

  function quad_fun_ND(n,w0,f0,aa,w)
    use nrtype
    implicit none
    real(dp) :: quad_fun_ND
    integer(i4b), intent(in) :: n
    real(dp), dimension(n), intent(in) :: w0
    real(dp), intent(in) :: f0
    real(dp), dimension(n,n), intent(in) :: aa
    real(dp), dimension(n), intent(in) :: w

    integer(i4b) :: i,j,k

    quad_fun_ND = f0 + dot_product(matmul(aa,w-w0),w-w0)

    return
  end function quad_fun_ND

  subroutine quad_coefs_ND(n,w0,f0,wa,fa,aa)
    use nrtype
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), dimension(n), intent(in) :: w0
    real(dp), intent(in) :: f0
    real(dp), dimension(n,n*(n+1)/2), intent(in) :: wa
    real(dp), dimension(n*(n+1)/2), intent(in) :: fa
    real(dp), dimension(n,n), intent(out) :: aa

    integer(i4b) :: i,j,k,l,info
    integer(i4b), dimension(n*(n+1)/2) :: ipiv
    real(dp), dimension(n*(n+1)/2,n*(n+1)/2) :: amat
    real(dp), dimension(n*(n+1)/2) :: aavec


    ! set the right hand side
    aavec(:) = fa(:) - f0

    ! loop over the different evaluation points
    do l = 1,n*(n+1)/2

       ! set the matrix elements
       k = 0
       do i = 1,n
          k = k+1
          amat(l,k) = (wa(i,l)-w0(i))**2
          do j = i+1,n
             k = k+1
             amat(l,k) = 2.0_dp*(wa(i,l)-w0(i))*(wa(j,l)-w0(j))
          end do
       end do
       
    end do

    ! solve the linear system
    call dgesv(n*(n+1)/2,1,amat,n*(n+1)/2,ipiv,aavec,n*(n+1)/2,info)
    if(info /= 0) stop 'quad_coefs_ND: problem solving linear system'

    ! package the coefficients in the symmetric matrix (upper triangle)
    k = 0
    do i = 1,n
       do j = i,n
          k = k+1
          aa(i,j) = aavec(k)
       end do
    end do

    ! fill in the lower triangle
    do i = 1,n
       do j = 1,i-1
          aa(i,j) = aa(j,i)
       end do
    end do

    
    return
  end subroutine quad_coefs_ND


  subroutine random_sample_SN(n,w)
    use nrtype
    use ran_state
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), dimension(n+1), intent(out) :: w

    integer(i4b) :: i
    real(dp) :: norm
    
    norm = 0.0_dp
    do i = 1,n+1
       call gasdev(w(i))
    end do
    norm = sqrt(dot_product(w,w))
    w = w/norm
    
    return
  end subroutine random_sample_SN


  subroutine point_correlation_S2(lmin,lmax,qq,th,cor)
    use nrtype
    use module_function
    implicit none
    integer(i4b), intent(in) :: lmin
    integer(i4b), intent(in) :: lmax
    real(dp), dimension(lmax+1), intent(in) :: qq
    real(dp), intent(in) :: th
    real(dp), intent(out) :: cor

    integer(i4b) :: l
    real(dp) :: fac
    real(dp), dimension(2) :: x,xp,xc

    
    cor = 0.0_dp
    do l = lmin,lmax

       fac = (1.0_dp + mu_S2_point**2*l*(l+1))**s_S2_point
       fac = 1.0_dp/fac

       call legendre(th,l,0,x,xp,xc)
       
       cor = cor + sqrt((2*l+1.0_dp)/fourpi_d)*fac*qq(l+1)*x(1)
       
    end do

    return
  end subroutine point_correlation_S2


  subroutine gaussian_confidence_set(n,alpha,ss)
    use nrtype
    use ran_state
    implicit none
    integer(i4b), intent(in) :: n
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: ss

    integer(i4b), parameter :: maxit = 100
    integer(i4b), parameter :: nsamp = 100000
    integer(i4b) :: isamp,i,it
    real(dp) :: nll,err,ss1,ss2,p,ran
    real(dp), dimension(nsamp) :: samp

    ! sample from the distribution
    do isamp = 1,nsamp
       nll = 0.0_dp
       do i = 1,n
          call gasdev(ran)
          nll = nll + 0.5_dp*ran**2
       end do
       samp(isamp) = nll     
    end do

    ! the value for s
    err = alpha*0.0001_dp
    ss1 = 0.0_dp
    ss2 = maxval(samp)
    do it = 1,maxit       
       ss = 0.5_dp*(ss1+ss2)
       p = ecdf(ss,samp)
       if(abs(1.0_dp-alpha-p) <= err) exit       
       if( p < 1.0_dp-alpha) ss1 = ss
       if( p > 1.0_dp-alpha) ss2 = ss
       if(it == maxit) stop 'gaussian_confidence_set: no convergence'
    end do
    

    return
  end subroutine gaussian_confidence_set

  
  function ecdf(x,samp)
    real(dp) :: ecdf
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: samp
    integer(i4b) :: i,n
    n = size(samp)
    ecdf = 0.0_dp
    do i = 1,n
       if(samp(i) <= x) ecdf = ecdf + 1.0_dp/n
    end do
    
    return
  end function ecdf

  
end module module_S2_point

