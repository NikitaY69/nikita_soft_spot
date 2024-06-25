
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

!***************************************************************************************************

contains



end module kinds