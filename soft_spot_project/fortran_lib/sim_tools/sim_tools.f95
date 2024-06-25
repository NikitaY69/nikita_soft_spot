
module sim_tools
!***************************************************************************************************
!
!
!
!***************************************************************************************************
  use kinds
!***************************************************************************************************
  implicit none
  public
!***************************************************************************************************

  contains




subroutine apply_lepbc(pos,box2,box,drx_shear)
  implicit none
!***************************************************************************************************
!
!
!
!***************************************************************************************************

  real(kind=dbl),intent(in)                    :: drx_shear
  real(kind=dbl),dimension(:),intent(in)       :: box2,box
  real(kind=dbl),dimension(:),intent(out)      :: pos

  integer                                      :: i
  integer                                      :: cory


  cory   = nint(pos(2)/box(2))
  pos(1) = pos(1) - real(cory,kind=dbl)*drx_shear


  do i = 1 , size(pos)

    if ( pos(i) .gt. box2(i) ) then
      pos(i) = pos(i) - box(i)
    elseif ( pos(i) .lt. -box2(i) ) then
      pos(i) = pos(i) + box(i)
    end if

  end do

return
end subroutine apply_lepbc



subroutine put_back_in_box(pos,box)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  real(kind=dbl),dimension(:),intent(in)       :: box
  real(kind=dbl),dimension(:),intent(out)      :: pos

  integer                                      :: i


  do i = 1 , size(pos)

    if ( pos(i) .gt. box(i) ) then
      pos(i) = pos(i) - box(i)
    elseif ( pos(i) .lt. 0._dbl ) then
      pos(i) = pos(i) + box(i)
    end if

  end do


return
end subroutine put_back_in_box




end module sim_tools
