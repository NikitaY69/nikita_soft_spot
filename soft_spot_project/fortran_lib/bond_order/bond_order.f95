
module bond_order
  use kinds
  use sim_tools
!***************************************************************************************************
  implicit none
  public
!***************************************************************************************************

contains



subroutine compute_2d_bop(r,list_ni,list_ij,boxl,delrx,bop)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  real(kind=dbl),dimension(:,:),intent(in)        :: r
  integer,dimension(:),intent(in)                 :: list_ni
  integer,dimension(:,:),intent(in)               :: list_ij
  real(kind=dbl),dimension(:),intent(in)          :: boxl
  real(kind=dbl),intent(in)                       :: delrx

  real(kind=dbl),dimension(size(r,2)),intent(out) :: bop


  integer                                         :: i,j,id_j
  real(kind=dbl),dimension(size(r,1))             :: drij,bop_ij,boxl2
  real(kind=dbl)                                  :: rij,cc,ss
  real(kind=dbl),dimension(size(r,1),size(r,2))   :: tmp_bop


  tmp_bop(:,:) = 0._dbl
  bop(:)       = 0.0_dbl
  boxl2(:)     = boxl(:)/2._dbl

  do i = 1 , size(r,2)
    do id_j = 1 , list_ni(i)

      j = list_ij(id_j,i)

      drij = r(:,j) - r(:,i)
      call apply_lepbc(drij,boxl2,boxl,delrx)
      rij = sqrt(dot_product(drij,drij))

      cc = drij(1) / rij
      ss = drij(2) / rij

      bop_ij(1) = cc**6 + 15.0_dbl*( (cc**2)*(ss**4) - (cc**4)*(ss**2) ) - ss**6
      bop_ij(2) = 6.0_dbl * ( ss*(cc**5) +cc*(ss**5) ) - 20.0_dbl*((cc*ss)**3)

      tmp_bop(:,i) = tmp_bop(:,i) + bop_ij(:)

    end do
    if(list_ni(i).gt.0)bop(i)=norm2(tmp_bop(:,i)/real(list_ni(i),kind=dbl))
  end do



end subroutine compute_2d_bop





subroutine compute_2d_theta(r,d,epsilon,list_ni,list_ij,boxl,delrx,theta)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  real(kind=dbl),dimension(:,:),intent(in)        :: r
  real(kind=dbl),dimension(:),intent(in)          :: d
  real(kind=dbl),intent(in)                       :: epsilon
  integer,dimension(:),intent(in)                 :: list_ni
  integer,dimension(:,:),intent(in)               :: list_ij
  real(kind=dbl),dimension(:),intent(in)          :: boxl
  real(kind=dbl),intent(in)                       :: delrx

  real(kind=dbl),dimension(size(r,2)),intent(out) :: theta


  integer                                         :: i,j,k,idj,idk,cpt,n_theta
  real(kind=dbl),dimension(size(r,1))             :: drij,drik,drjk,boxl2
  real(kind=dbl)                                  :: rij,rik,rjk,dij,dik,djk,theta1,theta2,diff_theta
  real(kind=dbl),dimension(:),allocatable         :: theta1_list,diff_theta_list
  integer,dimension(:),allocatable                :: id_theta1_list

  theta(:)       = 0.0_dbl
  boxl2(:)     = boxl(:)/2._dbl
  allocate(theta1_list(2),diff_theta_list(2),id_theta1_list(2))

  do i = 1 , size(r,2)
    cpt = 0
    n_theta = list_ni(i)*(list_ni(i)-1)/2
    deallocate(theta1_list,diff_theta_list,id_theta1_list)
    allocate(theta1_list(n_theta),diff_theta_list(n_theta),id_theta1_list(n_theta))


    do idj = 1 , list_ni(i)
      j = list_ij(idj,i)
      dij = 0.5_dbl*(d(i)+d(j))*(1._dbl-epsilon*abs(d(i)-d(j)))
      drij = r(:,j) - r(:,i)
      call apply_lepbc(drij,boxl2,boxl,delrx)
      rij = sqrt(dot_product(drij,drij))

      do idk = idj + 1 , list_ni(i)
        k = list_ij(idk,i)

        dik = 0.5_dbl*(d(i)+d(k))*(1._dbl-epsilon*abs(d(i)-d(k)))
        drik = r(:,k) - r(:,i)
        call apply_lepbc(drik,boxl2,boxl,delrx)
        rik = sqrt(dot_product(drik,drik))

        djk = 0.5_dbl*(d(j)+d(k))*(1._dbl-epsilon*abs(d(j)-d(k)))
        drjk = r(:,k) - r(:,j)
        call apply_lepbc(drjk,boxl2,boxl,delrx)
        rjk = sqrt(dot_product(drjk,drjk))

        theta1 = acos(dot_product(drij,drik)/(rij*rik))
        theta2 = acos((-djk*djk + dik*dik + dij*dij)/(2._dbl*dik*dij))

        diff_theta = theta1-theta2
        cpt = cpt + 1
        theta1_list(cpt) = theta1
        diff_theta_list(cpt) = diff_theta
        id_theta1_list(cpt) = cpt

      end do
    end do

    call bubble_sorting(theta1_list,id_theta1_list)

    do idj = 1 , list_ni(i)
      theta(i)=theta(i)+abs(diff_theta_list(id_theta1_list(idj)))
    end do
    if(list_ni(i).gt.0)theta(i)=theta(i)/list_ni(i)

  end do

end subroutine compute_2d_theta




subroutine bubble_sorting(a,id_a)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  real(kind=dbl),dimension(:),intent(inout) :: a
  integer,dimension(size(a)),intent(inout)  :: id_a

  real(kind=dbl) :: temp
  integer        :: i,j,temp_id
  logical        :: swapped


  do j = size(a)-1, 1, -1
    swapped = .false.
    do i = 1, j
      if (a(i) .gt. a(i+1)) then
        temp = a(i)
        temp_id = id_a(i)
        a(i) = a(i+1)
        a(i+1) = temp
        id_a(i) = id_a(i+1)
        id_a(i+1) = temp_id
        swapped = .true.
      end if
    end do
    if (.not. swapped) exit
  end do

end subroutine bubble_sorting




end module bond_order