program sobolev_disk

  use nrtype
  use module_util
  use module_mesh
  use module_eig
  use module_ODE
  implicit none
  

  integer(i4b) :: info,i,io,j,inode,ispec,ia,k,m,nmax,nr,ir,nth,ith,n
  real(dp) :: drmax,c,dr,r,phi,p0,dth,th,fun,x,y
  real(dp), dimension(2) :: dcdx
  real(dp), dimension(4) :: z
  real(dp), dimension(:), allocatable :: eval
  real(dp), dimension(:,:), allocatable :: evec

  
  call get_integer(' m = ',m)
  call get_integer(' nmax = ',nmax)
  call get_integer(' n = ',n)
    
  ! build the mesh
  drmax = 0.05_dp*twopi_d/nmax
  call mesh_disk(drmax)

  ! solve the eigenvalue problem
  allocate(eval(nmax))
  allocate(evec(nglob,nmax))
  call scalar_disk_eig('D',m,nmax,eval,evec)


  open(newunit = io,file='eval.out')
  do i = 1,nmax
     write(io,*) i,eval(i)
  end do
  close(io)

  open(newunit = io,file='evec.out')
  do ispec = 1,nspec
     do inode = 1,ngll
        ia = ibool(inode,ispec)
        write(io,*) r_node(inode,ispec),evec(ia,:)
     end do
  end do
  close(io)

  nth = max(30*m,360)
  dth = twopi_d/(nth-1)
  open(newunit = io,file = 'evec2D.out')
  write(io,*) ngll*nspec,nth,0
  do ispec = 1,nspec
     do inode = 1,ngll
        r = r_node(inode,ispec)
        ia = ibool(inode,ispec)
        do ith = 1,nth
           th = (ith-1)*dth
           fun = evec(ia,n)*cos(m*th)
           x = r*cos(th)
           y = r*sin(th)
           write(io,*) x,y,fun
        end do
     end do
  end do
  close(io)


  





  
end program sobolev_disk
