
module pairwise
  use kinds
  use sim_tools
  use cell
  use system
  use potentials
!***************************************************************************************************
  implicit none
  public
!***************************************************************************************************



!***************************************************************************************************


  integer,parameter,public                             :: pair_ncontacts = 32
  integer,public                                       :: pair_count
  integer,dimension(:),allocatable,public              :: pair_lookupbond
  real(kind=dbl),dimension(:),allocatable,public       :: pair_first,pair_second,pair_third,pair_fourth
  real(kind=dbl),dimension(:,:),allocatable,public     :: pair_rij
  real(kind=dbl),public                                :: typical_contact_force

  integer,dimension(:),allocatable,public     :: pair_list_ni
  integer,dimension(:,:),allocatable,public   :: pair_list_ij
contains



subroutine init_pairwise_arrays()

  if( allocated(pair_lookupbond) .eqv. .false. ) then
    allocate(pair_lookupbond(0:npart*pair_ncontacts*2-1))
    allocate(pair_first(0:npart*pair_ncontacts-1),pair_second(0:npart*pair_ncontacts-1))
    allocate(pair_third(0:npart*pair_ncontacts-1),pair_fourth(0:npart*pair_ncontacts-1))
    allocate(pair_rij(0:ndim-1,0:npart*pair_ncontacts-1))
    allocate(pair_list_ni(0:npart-1),pair_list_ij(0:pair_ncontacts,0:npart-1))
  end if

  pair_lookupbond(:) = 0
  pair_first(:)   = 0._dbl
  pair_second(:)  = 0._dbl
  pair_third(:)   = 0._dbl
  pair_fourth(:)  = 0._dbl
  pair_rij(:,:)   = 0._dbl


end subroutine init_pairwise_arrays





subroutine create_chain_list()
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  integer                                                :: i,icel
  integer,dimension(ndim)                                :: id

  frc_cell%hoc_a(:) = 0

  do i = 1 , size(r,2)

    id(:) = floor( r(:,i) / frc_cell%diam(:) )
    icel  = frc_cell%n_xyz(1)*id(2) + id(1) + 1 
    if(ndim.eq.3)icel  = icel + frc_cell%n_xyz(1)*frc_cell%n_xyz(2)*id(3)

    frc_cell%list_a(i)   = icel

  end do

  do i = 1 , size(r,2)
    frc_cell%ll_a(i) = frc_cell%hoc_a(frc_cell%list_a(i))
    frc_cell%hoc_a(frc_cell%list_a(i)) = i
  end do


end subroutine create_chain_list



subroutine update_chain_list(id,old_cell,new_cell)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  integer,intent(in)             :: id,old_cell,new_cell

  integer        :: i,j,k,l

  j = frc_cell%hoc_a(old_cell)
  i = 1

  if(j.eq.id) then
    k = frc_cell%ll_a(j)
    frc_cell%hoc_a(old_cell) = k
  else 
    do while ( j .ne. 0 )
      i = i + 1
      k = j
      j = frc_cell%ll_a(j)
      if(j.eq.id)then
        j = 0
        l = frc_cell%ll_a(id)
      end if
    end do

    frc_cell%ll_a(k) = l
  end if

  frc_cell%ll_a(id) = frc_cell%hoc_a(new_cell)
  frc_cell%hoc_a(new_cell) = id
  frc_cell%list_a(i) = new_cell

end subroutine update_chain_list



subroutine compute_poly_everything()
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  integer                                                :: i,j,l,m,icel,jcel,ncel
  real(kind=dbl),dimension(ndim)                         :: drij,ri,fij
  real(kind=dbl)                                         :: rij2,dij,rcut2
  type(const_vec)                                        :: out_full


  allocate(out_full%c(4))
  out_full%c(:) = 0._dbl


  pair_count = 0
  pair_lookupbond(:) = 0 
  pair_first(:)  = 0._dbl
  pair_second(:) = 0._dbl
  pair_third(:)  = 0._dbl
  pair_fourth(:) = 0._dbl
  pair_rij(:,:)  = 0._dbl

  typical_grad = 0._dbl
  thermo_pre   = 0._dbl
  thermo_sigma = 0._dbl
  thermo_pot   = 0._dbl

  pair_list_ni(:) = 0
  pair_list_ij(:,:) = 0

  f(:,:) = 0._dbl

  call create_chain_list()

  do i = 1 , npart

    ri(:) = r(:,i)
    icel  = frc_cell%list_a(i)

    do ncel = 1 , frc_cell%nc

      jcel = frc_cell%maps( ( icel - 1 )*frc_cell%nc + ncel)

      if ( jcel .gt. 0 ) then

        j = frc_cell%hoc_a(jcel)

        do while ( j .ne. 0 )

          if ( j .gt. i ) then

            drij(:) = r(:,j) - ri(:)
            call apply_lepbc(drij,boxl2,boxl,box_delrx)
            rij2 = dot_product(drij,drij)

            dij   = get_dij_poly(d(i),d(j),dij_arg_poly)
            rcut2 = get_rcut_poly(dij,rc_arg_poly)

            if ( rij2 .lt. rcut2 ) then

              pair_lookupbond(2*pair_count)   = i - 1
              pair_lookupbond(2*pair_count+1) = j - 1
              pair_rij(:,pair_count) = drij(:)

              pair_list_ij(pair_list_ni(i-1),i-1) = j-1
              pair_list_ij(pair_list_ni(j-1),j-1) = i-1
              pair_list_ni(i-1) = pair_list_ni(i-1) + 1
              pair_list_ni(j-1) = pair_list_ni(j-1) + 1

              call get_full_poly(rij2,dij,phir_arg_poly,out_full)
              pair_first(pair_count)  = out_full%c(1)
              pair_second(pair_count) = out_full%c(2)
              pair_third(pair_count)  = out_full%c(3)
              pair_fourth(pair_count) = out_full%c(4)


              f(:,i) = f(:,i) + drij(:)*pair_first(pair_count)
              f(:,j) = f(:,j) - drij(:)*pair_first(pair_count)

              thermo_pre   = thermo_pre - pair_first(pair_count)*dot_product(drij,drij)
              thermo_sigma = thermo_sigma + pair_first(pair_count)*drij(1)*drij(2)
              thermo_pot   = thermo_pot + get_pot_poly(rij2,dij,pot_arg_poly)
    

              pair_count = pair_count + 1


            end if

          end if
          j = frc_cell%ll_a(j)
        end do

      end if

    end do

  end do

  deallocate(out_full%c)

  do i = 1 , npart
    typical_grad = typical_grad + dot_product(f(:,i),f(:,i))
  end do
  typical_grad = sqrt(typical_grad/real(npart,kind=dbl))

  thermo_pre = thermo_pre / (real(ndim,kind=dbl)*box_volume)
  thermo_sigma = thermo_sigma / box_volume


return
end subroutine compute_poly_everything





subroutine compute_hessian(hessian)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  real(kind=dbl),dimension(0:2,0:ndim*ndim*npart+ndim*ndim*ndim*pair_count-1),intent(out) :: hessian
  
  integer                                        :: i,j,m,l,alpha,beta,a,b,cpt_line
  real(kind=dbl)                                 :: thisHes
  real(kind=dbl),dimension(0:ndim*ndim*npart-1)  :: diagonal

  hessian(:,:) = 0._dbl
  diagonal(:)  = 0._dbl

  cpt_line = 0
  do l = 0 , pair_count - 1

    i = pair_lookupbond(2*l)
    j = pair_lookupbond(2*l+1)
    m = 0

    do alpha = 0 , ndim - 1
      do beta =  0 , ndim - 1
        a = ndim*i + alpha
        b = ndim*j + beta
        thisHes = -pair_second(l)*pair_rij(alpha,l)*pair_rij(beta,l) - pair_first(l)*kron_delta(alpha,beta)

        hessian(0,cpt_line) = a
        hessian(1,cpt_line) = b
        hessian(2,cpt_line) = thisHes
        cpt_line = cpt_line + 1
        hessian(0,cpt_line) = b
        hessian(1,cpt_line) = a
        hessian(2,cpt_line) = thisHes
        cpt_line = cpt_line + 1

        diagonal(ndim*ndim*i + m) = diagonal(ndim*ndim*i + m) - thisHes
        diagonal(ndim*ndim*j + m) = diagonal(ndim*ndim*j + m) - thisHes
        m = m + 1
      end do
    end do

  end do

  do i = 0 , npart -1
    m = 0
    do alpha = 0 , ndim - 1
      do beta =  0 , ndim - 1
        a = ndim*i + alpha
        b = ndim*i + beta

        hessian(0,cpt_line) = a
        hessian(1,cpt_line) = b
        hessian(2,cpt_line) = diagonal(ndim*ndim*i+m)
        cpt_line = cpt_line + 1

        m = m + 1
      end do
    end do
  end do


end subroutine compute_hessian



function kron_delta(a,b) result(d)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  integer, intent(in) :: a,b
  integer             :: d

  d = 0
  if(a.eq.b)d=1

end function kron_delta

end module pairwise



