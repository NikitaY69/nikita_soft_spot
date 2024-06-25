
module potentials
  use kinds
  use system
!***************************************************************************************************
  implicit none
  public
!***************************************************************************************************

  interface

    function func_type_a(q1,arg) result(return_value)
      use kinds
      real(kind=dbl),intent(in)    :: q1
      type(const_vec),intent(in)   :: arg
      real(kind=dbl)               :: return_value
    end function func_type_a

    function func_type_b(q1,q2,arg) result(return_value)
      use kinds
      real(kind=dbl),intent(in)    :: q1,q2
      type(const_vec),intent(in)   :: arg
      real(kind=dbl)               :: return_value
    end function func_type_b

    subroutine func_type_c(q1,q2,arg,argout)
      use kinds
      real(kind=dbl),intent(in)      :: q1,q2
      type(const_vec),intent(in)     :: arg
      type(const_vec),intent(inout)  :: argout
    end subroutine func_type_c


  end interface




  !poly force field
  procedure(func_type_b), pointer :: get_pot_poly
  procedure(func_type_b), pointer :: get_phir_poly
  procedure(func_type_a), pointer :: get_rcut_poly
  procedure(func_type_c), pointer :: get_full_poly
  procedure(func_type_b), pointer :: get_dij_poly
  type(const_vec)       :: pot_arg_poly
  type(const_vec)       :: phir_arg_poly
  type(const_vec)       :: rc_arg_poly
  type(const_vec)       :: dij_arg_poly



contains






subroutine init_poly_force_field(type,arg)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  character(len=50),intent(in)            :: type
  real(kind=dbl),dimension(:),intent(in)  :: arg



  write(*,*) '=> init polydisperse force field:'
  write(*,*) 'potential:',type
  write(*,*) 'parameter:',arg



  select case (type)

    case ( 'swap_ipl' ) ! repulsive soft spheres n=12
      get_pot_poly  => pot_swap_ipl
      get_phir_poly => phir_swap_ipl
      get_rcut_poly => rcut_swap_ipl
      get_dij_poly  => dij_swap_ipl
      get_full_poly => full_swap_ipl

      if ( allocated(pot_arg_poly%c) .eqv. .false. ) then
        allocate(pot_arg_poly%c(4),rc_arg_poly%c(1),dij_arg_poly%c(1),phir_arg_poly%c(3))
      end if
      pot_arg_poly%c(1) = arg(1) !v0
      pot_arg_poly%c(2) =  -28._dbl*((1._dbl/arg(2))**12.) !c_0
      pot_arg_poly%c(3) =  48._dbl*((1._dbl/arg(2))**14.)  !c_2
      pot_arg_poly%c(4) =  -21._dbl*((1._dbl/arg(2))**16.) !c_4

      phir_arg_poly%c(1) = arg(1) !v0
      phir_arg_poly%c(2) =  48._dbl*((1._dbl/arg(2))**14.)  !c_2
      phir_arg_poly%c(3) =  -21._dbl*((1._dbl/arg(2))**16.) !c_4

      rc_arg_poly%c(1)  = arg(2) ! rcutone
      dij_arg_poly%c(1) = arg(3) ! epsilon


    case default
      stop 'ERROR: your potential is undefined'

  end select

end subroutine init_poly_force_field





!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************



function dij_swap_ipl(di,dj,arg) result(dij2)
  implicit none
!***************************************************************************************************
!arg=epsilon
!***************************************************************************************************

  real(kind=dbl),intent(in)               :: di,dj
  type(const_vec),intent(in)              :: arg

  real(kind=dbl)                          :: dij2

  dij2 = 0.5_dbl*(di+dj)*(1._dbl-arg%c(1)*abs(di-dj))
  dij2 = dij2 * dij2

return
end function dij_swap_ipl



function rcut_swap_ipl(dij2,arg) result(rc2)
  implicit none
!***************************************************************************************************
!arg=rcutone
!***************************************************************************************************

  real(kind=dbl),intent(in)               :: dij2
  type(const_vec),intent(in)              :: arg

  real(kind=dbl)                          :: rc2

  rc2 = arg%c(1)*arg%c(1)*dij2

return
end function rcut_swap_ipl



function pot_swap_ipl(r2,dij2,arg) result(pot)
  implicit none
!***************************************************************************************************
! arg=v0,c0,c2,c4
!***************************************************************************************************

  real(kind=dbl),intent(in)               :: r2,dij2
  type(const_vec),intent(in)              :: arg

  real(kind=dbl)                          :: pot
  real(kind=dbl)                          :: sigr2,sigr4

  sigr2  = dij2 / r2
  sigr4  = sigr2*sigr2

  pot = arg%c(1)*sigr4*sigr4*sigr4 + arg%c(2) + arg%c(3)/sigr2 + arg%c(4)/sigr4


return
end function pot_swap_ipl



function phir_swap_ipl(r2,dij2,arg) result(phir)
  implicit none
!***************************************************************************************************
! arg=v0,c2,c4
!***************************************************************************************************

  real(kind=dbl),intent(in)               :: r2,dij2
  type(const_vec),intent(in)              :: arg

  real(kind=dbl)                          :: phir
  real(kind=dbl)                          :: sigr2,sigr4


  sigr2  = dij2 / r2
  sigr4  = sigr2*sigr2

  phir = -12._dbl*arg%c(1)*sigr4*sigr4*sigr4 + 2._dbl*arg%c(2)/sigr2 + 4._dbl*arg%c(3)/sigr4
  phir = phir / r2

return
end function phir_swap_ipl



subroutine full_swap_ipl(r2,dij2,arg,arg_out)
  implicit none
!***************************************************************************************************
! arg=v0,c2,c4
!***************************************************************************************************

  real(kind=dbl),intent(in)               :: r2,dij2
  type(const_vec),intent(in)              :: arg

  type(const_vec),intent(inout)           :: arg_out

  real(kind=dbl)                          :: sigr2,sigr6,sigr12,dij4

  sigr2   = dij2/r2
  sigr6   = sigr2*sigr2*sigr2
  sigr12  = sigr6*sigr6
  dij4    = dij2*dij2


  arg_out%c(1)  = ( 2._dbl*arg%c(2) + 4._dbl*arg%c(3)/sigr2  - arg%c(1)*12._dbl*sigr12*sigr2 )/dij2
  arg_out%c(2)  = (8._dbl*arg%c(3)+arg%c(1)*168._dbl*sigr12*sigr2*sigr2)/dij4
  arg_out%c(3)  = -arg%c(1)*2688._dbl*sigr12*sigr6/(dij4*dij2)
  arg_out%c(4)  = arg%c(1)*48384._dbl*sigr12*sigr6*sigr2/(dij4*dij4)


return
end subroutine full_swap_ipl


!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************



end module potentials






