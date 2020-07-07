module module_eig
  
  use nrtype
  implicit none

  
contains

  subroutine scalar_disk_eig(bc,m,nmax,eval,evec)
    use nrtype
    use module_util
    use module_mesh
    implicit none
    character(len=1), intent(in) :: bc
    integer(i4b), intent(in) :: m
    integer(i4b), intent(in) :: nmax
    real(dp), dimension(:), intent(out) :: eval
    real(dp), dimension(:,:), intent(out) :: evec
    

    integer(i4b) :: n,ispec,inode,jnode,knode, &
                    ia,ja,ka,info,io,i,j,kda,kdb,ldab,ldbb
    real(dp), dimension(:,:), allocatable :: a
    real(dp), dimension(:,:), allocatable :: b
    real(dp), dimension(:), allocatable   :: w
    real(dp), dimension(:,:), allocatable :: z  
    real(dp), dimension(:), allocatable   :: work

    if(bc /= 'D' .and. bc /= 'N') stop 'problem with boundary conditions'

    !=====================================================!
    !              set up the various arrays              !
    !=====================================================!
    
    ! set and store the band widths
    kda = ngll-1
    ldab = kda+1
    kdb = 0
    ldbb = kdb+1

    ! set the size of the linear system
    if(m == 0) then
       if(bc == 'D') then
          n = nglob-1
       else 
          n = nglob
       end if
    else
       if(bc == 'D') then
          n = nglob-2
       else 
          n = nglob-1
       end if
    end if

    ! allocate the arrays
    allocate(a(ldab,n))
    allocate(b(ldbb,n))
    allocate(w(n))
    allocate(z(n,n))
    allocate(work(3*n-2))

    !=========================================================!
    !                 form the system matrix                  !
    !=========================================================!
    
    ! initialise the system sysrices
    a = 0.0_dp
    b = 0.0_dp
    
    ! begin loop over the spectral elements
    do ispec = 1,nspec

       do inode = 1,ngll

          ! set the first index
          ia = ibool(inode,ispec)

          ! skip the final element do to BC
          if(bc == 'D') then
             if(ia == nglob) cycle
          end if
          
          ! adjust indices if needed
          if(m > 0) then
             if(ia == 1) cycle
             ia = ia-1
          end if

          ! add the non-derivative term to the b matrix
          ka = kdb+1
          b(ka,ia) = b(ka,ia) + r_node(inode,ispec) &
                              * wgll(inode)         &
                              * jac_element(ispec)

          if(m > 0) then
             ka = kda+1
             a(ka,ia) = a(ka,ia) + m*m/r_node(inode,ispec) &
                                 * wgll(inode)         &
                                 * jac_element(ispec)
             
          end if

          
          do jnode = 1,ngll

             ! set the second index
             ja = ibool(jnode,ispec)

             ! skip the final element do to BC
             if(bc == 'D') then
                if(ja == nglob) cycle
             end if
             
             ! adjust indices if needed
             if(m > 0) then
                if(ja == 1) cycle
                ja = ja-1
             end if
             
             ! skip the lower triangle
             if(ia > ja) cycle
             
             ! set the reduces column index
             ka = kda+1+ia-ja

             do knode = 1,ngll

                a(ka,ja) = a(ka,ja) + r_node(knode,ispec) &
                                    * hprime(knode,inode) &
                                    * hprime(knode,jnode) & 
                                    * wgll(knode)         &
                                    / jac_element(ispec)

                
             end do
             
          end do
       end do
       
    end do
    ! end loop over the spectral elements

    ! if m = 0 add a small number to B to make it definite
    ! this is a bit crude, but seems to work okay
    if(m == 0) then
       b(1,1) = 0.0001_dp*(minval(abs(b(1,2:n))))
    end if
    

    !=========================================================!
    !               solve the eigenvalue problem              !
    !=========================================================!

    ! call the LAPACK routine
    call dsbgv('V','U',n,kda,kdb,a,ldab,b,ldbb,w,z,n,work,info)
    if(info /= 0) stop 'gevp_disk_dirichlet: problem with diagonalisation'


    ! store the eigenvalues
    eval = w(1:nmax)

    ! store the eigenvectors
    do ispec = 1,nspec
       do inode = 1,ngll
          ia = ibool(inode,ispec)
          if(m == 0) then
             if(bc == 'D') then
                if(ia == nglob) then
                   evec(ia,:) = 0.0_dp
                else
                   evec(ia,:) = z(ia,1:nmax)
                end if
             else                
                evec(ia,:) = z(ia,1:nmax)
             end if
          else
             if(bc == 'D') then             
                if(ia == nglob .or. ia == 1) then
                   evec(ia,:) = 0.0_dp
                else
                   evec(ia,:) = z(ia-1,1:nmax)
                end if
             else
                if(ia == 1) then
                   evec(ia,:) = 0.0_dp
                else
                   evec(ia,:) = z(ia-1,1:nmax)
                end if                
             end if
          end if             
       end do
    end do

    
    
    return
  end subroutine scalar_disk_eig



  
end module module_eig
