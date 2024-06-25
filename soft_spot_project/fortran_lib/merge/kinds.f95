
module kinds
!***************************************************************************************************
! precision used in the code and some derived types
!***************************************************************************************************
  implicit none
  public
!***************************************************************************************************

  integer,parameter  :: dbl = 8   ! regular 8 byte real
  integer,parameter  :: sgl = 4   ! regular 4 byte real

  integer, parameter :: lint = selected_int_kind(15) ! regular 8 byte int

  type const_vec
    real(dbl),dimension(:),allocatable      :: c(:) ! vector of parameters
  end type const_vec

!***************************************************************************************************

contains


subroutine init_const_vec(dertype, n)
  type(const_vec),intent(inout) :: dertype
  integer,intent(in)            :: n
  allocate(dertype%c(n))
  dertype%c(:) = 0.0_dbl
end subroutine init_const_vec

subroutine destroy_const_vec(dertype)
  type(const_vec),intent(inout) :: dertype
  deallocate(dertype%c)
end subroutine destroy_const_vec




end module kinds