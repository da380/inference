module module_mesh

  use nrtype
  implicit none

  ! store the mesh parameters
  integer(i4b), parameter :: ngll = 5
  integer(i4b), save :: nspec,nglob
  integer(i4b), dimension(:,:), allocatable, save :: ibool

  real(dp), dimension(ngll) :: xigll,wgll
  real(dp), dimension(ngll,ngll) :: hprime
  real(dp), dimension(:,:), allocatable, save :: r_node
  real(dp), dimension(:,:), allocatable, save :: vel_node
  real(dp), dimension(:), allocatable, save :: jac_element

  
contains


  subroutine read_radial_velocity_model(vfile,vfile_out)
    use nrtype
    use module_interp
    character(len=*), intent(in) :: vfile
    character(len=*), intent(in), optional :: vfile_out

    logical(lgt) :: ltmp
    integer(i4b) :: io,ios,nknot,iknot,ispec,inode
    real(dp) :: rtmp,vel
    real(dp), dimension(:), allocatable :: r_knot
    real(dp), dimension(:,:), allocatable :: vel_knot
    
    ! open the model file
    inquire(file=trim(vfile),exist = ltmp)
    if(.not.ltmp) stop 'read_radial_velocity_model: velocity file does not exist'

    open(newunit = io,file=trim(vfile),action='read')

    ! work out number of knots
    nknot = 0
    do
       read(io,*,iostat =  ios) rtmp
       if(ios /= 0) exit
       nknot = nknot+1
    end do
    rewind(io)

    ! read in the model
    allocate(r_knot(nknot),vel_knot(nknot,2))
    do iknot = 1,nknot
       read(io,*) r_knot(iknot),vel_knot(iknot,1)
    end do    
    close(io)

    ! make the cubic splines!
    call spline(r_knot,vel_knot(:,1),1.0e40_dp,1.0e40_dp,vel_knot(:,2))




    ! evaluate the model on the mesh nodes
    allocate(vel_node(ngll,nspec))
    do ispec = 1,nspec
       do inode = 1,ngll
          rtmp = r_node(inode,ispec)
          vel = splint(r_knot,vel_knot(:,1),vel_knot(:,2),rtmp)
          vel_node(inode,ispec) = vel
       end do
    end do

    if(present(vfile_out)) then
       open(newunit = io,file=trim(vfile_out))
           do ispec = 1,nspec
       do inode = 1,ngll
          write(io,*) r_node(inode,ispec),vel_node(inode,ispec)
       end do
    end do
       close(io)       
    end if
    
    return
  end subroutine read_radial_velocity_model


  
  subroutine mesh_disk(drmax)
    use nrtype
    use module_GLL
    implicit none
    real(dp), intent(in) :: drmax

    integer(i4b) :: inode,jnode,ispec,k
    real(dp) :: r1,r2,dr


    ! set up the GLL points and weights  
    call zwgljd(xigll,wgll,ngll,0.0_dp,0.0_dp)
    if(mod(ngll,2) /= 0) xigll((ngll-1)/2+1) = 0.0_dp
    
    ! get the derivatives of the Lagrange polynomials at 
    ! the GLL points; recall that  hprime(i,j)=h'_{j}(xigll_{i}) 
    do jnode=1,ngll
       do inode=1,ngll
          hprime(inode,jnode) = lagrange_deriv_GLL(jnode-1,inode-1,xigll,ngll)
       end do
    end do
    
    ! work out the number of spectral elements
    nspec = 1.0_dp/drmax + 1
    dr = 1.0_dp/nspec

    ! allocate the mesh arrays
    allocate(r_node(ngll,nspec))
    allocate(jac_element(nspec))

      ! loop over the spectral elements and build up the mesh 
    r1 = 0.0_dp
    do ispec = 1,nspec
       
       ! set the end point
       r2 = r1 + dr
       
       ! loop over the nodes of the spectral element
       do inode = 1,ngll
          
        ! set the node position
          r_node(inode,ispec) = r1 + 0.5_dp*(xigll(inode)+1.0_dp)*(r2-r1)
          
       end do

       ! calculate the Jacobian for the element
       jac_element(ispec) = 0.5_dp*(r2-r1)
       
       ! update the start point
       r1 = r2
       
    end do

    allocate(ibool(ngll,nspec))
    k = 1
    do ispec = 1,nspec
       do inode = 1,ngll
          if(inode /= 1) k = k+1
          ibool(inode,ispec) = k
       end do
    end do
    
  ! number of global nodes
  nglob = maxval(ibool)
    
    return
  end subroutine mesh_disk

  
  


  subroutine radial_SEM_velocity(x,c,dcdx,ddcdxx)
    use nrtype
    use module_GLL
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out), optional :: c
    real(dp), dimension(:), intent(out), optional :: dcdx
    real(dp), dimension(:,:), intent(out), optional :: ddcdxx

    logical(lgt) :: bnd
    integer(i4b) :: inode,ispec,jnode
    real(dp) :: r,dr,xi,r1,r2
    real(dp), dimension(ngll) :: h,hp

    ! work out the value of xi
    bnd = .false.
    r = x(1)
    dr = r_node(ngll,nspec)/nspec
    ispec = r/dr + 1
    if(ispec < 1) then
       ispec = 1
       r = r_node(ngll,1)
    end if
    if(ispec > nspec) then
       ispec = nspec
       r = r_node(ngll,nspec)
       bnd = .true.
    end if
    r1 = r_node(1,ispec)
    r2 = r_node(ngll,ispec)
    xi = 2.0_dp*(r-r1)/(r2-r1)-1.0_dp

    ! get the lagrange polynomials at the radius
    call lagrange_any(xi,ngll,xigll,h,hp)

    if(present(c)) then

       c = 0.0_dp
       do inode = 1,ngll
          c = c + vel_node(inode,ispec)*h(inode)
       end do

    end if

    if(present(dcdx)) then

       
       dcdx(1) = 0.0_dp
       if(.not.bnd) then
          do inode = 1,ngll
             dcdx(1) = dcdx(1) + vel_node(inode,ispec) &
                               * hp(inode)/jac_element(ispec)
          end do
       end if
       dcdx(2) = 0.0_dp

       
    end if

    if(present(ddcdxx)) then
       stop 'radial_velocity: second derivatives not implemented yet!'       
    end if
    
    
    return
  end subroutine radial_SEM_velocity
  
end module module_mesh
