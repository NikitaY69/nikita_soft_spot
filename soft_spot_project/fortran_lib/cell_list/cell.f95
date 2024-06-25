
module cell
  use kinds
  use system
!***************************************************************************************************
!
!***************************************************************************************************
  implicit none
  public
!***************************************************************************************************

  type cell_info
    integer                                    :: nt
    integer                                    :: nc
    integer,dimension(:),allocatable           :: n_xyz
    real(kind=dbl),dimension(:),allocatable    :: diam
    integer,dimension(:),allocatable           :: ll_a,hoc_a
    integer,dimension(:),allocatable           :: list_a
    integer,dimension(:),allocatable           :: maps
  end type cell_info

!***************************************************************************************************

  type(cell_info),public     :: frc_cell
  real(kind=dbl)             :: frc_cell_rc_min

!***************************************************************************************************


  contains



subroutine init_cell()
  implicit none
!***************************************************************************************************
! 
!***************************************************************************************************

  !write(*,*) '=> init cell list:',frc_cell_rc_min

  if ( allocated(frc_cell%n_xyz) .eqv. .true. ) then
    deallocate(frc_cell%n_xyz,frc_cell%diam)
    allocate(frc_cell%n_xyz(size(boxl)),frc_cell%diam(size(boxl)))
  else
    allocate(frc_cell%n_xyz(size(boxl)),frc_cell%diam(size(boxl)))
  end if
  frc_cell%n_xyz(:) = 0
  frc_cell%diam(:)  = 0.0_dbl

  frc_cell%n_xyz(:) = floor( boxl(:) / frc_cell_rc_min )
  if(frc_cell%n_xyz(1).lt.3)frc_cell%n_xyz(1) = 3
  if(frc_cell%n_xyz(2).lt.3)frc_cell%n_xyz(2) = 3
  if((ndim.eq.3) .and. (frc_cell%n_xyz(3).lt.3))frc_cell%n_xyz(3) = 3

  frc_cell%nt = product(frc_cell%n_xyz(:))

  !print*,frc_cell%n_xyz
  frc_cell%diam(:)  = boxl(:) / real(frc_cell%n_xyz(:),kind=dbl)


  frc_cell%nc = 10
  if(ndim.eq.3)frc_cell%nc = 30


  if ( allocated(frc_cell%hoc_a) .eqv. .true. ) then
    deallocate(frc_cell%hoc_a,frc_cell%maps)
    allocate(frc_cell%hoc_a(frc_cell%nt),frc_cell%maps(frc_cell%nc*frc_cell%nt))
  else
    allocate(frc_cell%hoc_a(frc_cell%nt),frc_cell%maps(frc_cell%nc*frc_cell%nt))
  end if

  frc_cell%hoc_a(:)  = 0
  frc_cell%maps(:)   = 0


  if ( allocated(frc_cell%ll_a) .eqv. .false. ) then
    allocate(frc_cell%ll_a(npart),frc_cell%list_a(npart))
  end if

  frc_cell%ll_a(:)   = 0
  frc_cell%list_a(:) = 0


  if ( ndim .eq. 2 ) then
    call construct_2dcell_map()
    call update_shear_2dcell_map()
  end if

  if ( ndim .eq. 3 ) then
    call construct_3dcell_map()
    call update_shear_3dcell_map()
  end if

  !write(*,*)'---------------------------------------------------------------------------'

return
end subroutine init_cell



function get_2d_icell(idx,idy,n_xyz) result(icell)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  integer,intent(in)              :: idx,idy
  integer,dimension(2),intent(in) :: n_xyz

  integer                         :: icell

  icell =  1 + mod(idx-1+n_xyz(1),n_xyz(1)) 
  icell = icell + mod(idy-1+n_xyz(2),n_xyz(2))*n_xyz(1) 

return
end function get_2d_icell



subroutine construct_2dcell_map()
  implicit none
!***************************************************************************************************
! 
!***************************************************************************************************

  integer                              :: i,j,imap

  frc_cell%maps(:)    = 0

  do i = 1 , frc_cell%n_xyz(1)

    do j = 1 , frc_cell%n_xyz(2) !- 1 ! without top and bottom cells

      imap = ( get_2d_icell(i,j,frc_cell%n_xyz) - 1 ) * frc_cell%nc ! shift in the array

      frc_cell%maps( imap + 1 )  = get_2d_icell(i,j,frc_cell%n_xyz)
      frc_cell%maps( imap + 2 )  = get_2d_icell(i+1,j,frc_cell%n_xyz)
      frc_cell%maps( imap + 3 )  = get_2d_icell(i-1,j,frc_cell%n_xyz)

      frc_cell%maps( imap + 4 )  = get_2d_icell(i,j+1,frc_cell%n_xyz)
      frc_cell%maps( imap + 5 )  = get_2d_icell(i+1,j+1,frc_cell%n_xyz)
      frc_cell%maps( imap + 6 )  = get_2d_icell(i-1,j+1,frc_cell%n_xyz)

      frc_cell%maps( imap + 7 )  = get_2d_icell(i,j-1,frc_cell%n_xyz)
      frc_cell%maps( imap + 8 )  = get_2d_icell(i+1,j-1,frc_cell%n_xyz)
      frc_cell%maps( imap + 9 )  = get_2d_icell(i-1,j-1,frc_cell%n_xyz)

    end do

  end do


return
end subroutine construct_2dcell_map



subroutine update_shear_2dcell_map()
  implicit none
!***************************************************************************************************
! 
!***************************************************************************************************

  integer                              :: i,j,imap,dx,jcel

  dx = floor( box_delrx / frc_cell%diam(1) ) ! flow in x direction

  j = 1 ! bottom boundary

  do i = 1 , frc_cell%n_xyz(1)

    imap = ( get_2d_icell(i,j,frc_cell%n_xyz) - 1 ) * frc_cell%nc ! shift in the array

    jcel = get_2d_icell(i+dx,j-1,frc_cell%n_xyz)
    frc_cell%maps( imap + 7 )  = jcel

    jcel = get_2d_icell(i+1+dx,j-1,frc_cell%n_xyz)
    frc_cell%maps( imap + 8 )  = jcel

    jcel = get_2d_icell(i-1+dx,j-1,frc_cell%n_xyz)
    frc_cell%maps( imap + 9 )  = jcel

    jcel = get_2d_icell(i+2+dx,j-1,frc_cell%n_xyz)
    frc_cell%maps( imap + 10 )  = jcel

  end do


  j = frc_cell%n_xyz(2) ! upper boundary

  do i = 1 , frc_cell%n_xyz(1)

    imap = ( get_2d_icell(i,j,frc_cell%n_xyz) - 1 ) * frc_cell%nc ! shift in the array

    jcel = get_2d_icell(i-dx,j+1,frc_cell%n_xyz)
    frc_cell%maps( imap + 4 )  = jcel

    jcel = get_2d_icell(i+1-dx,j+1,frc_cell%n_xyz)
    frc_cell%maps( imap + 5 )  = jcel

    jcel = get_2d_icell(i-1-dx,j+1,frc_cell%n_xyz)
    frc_cell%maps( imap + 6 )  = jcel

    jcel = get_2d_icell(i-2-dx,j+1,frc_cell%n_xyz)
    frc_cell%maps( imap + 10 )  = jcel


  end do



return
end subroutine update_shear_2dcell_map



function get_3d_icell(idx,idy,idz,n_xyz) result(icell)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  integer,intent(in)              :: idx,idy,idz
  integer,dimension(3),intent(in) :: n_xyz

  integer                         :: icell


  icell =  1 + mod(idx-1+n_xyz(1),n_xyz(1)) 
  icell = icell + mod(idy-1+n_xyz(2),n_xyz(2))*n_xyz(1) 
  icell = icell + mod(idz-1+n_xyz(3),n_xyz(3))*n_xyz(1)*n_xyz(2)


return
end function get_3d_icell



subroutine construct_3dcell_map()
  implicit none
!***************************************************************************************************
! 
!***************************************************************************************************

  integer                              :: i,j,k,imap

  frc_cell%maps(:)    = 0

  do i = 1 , frc_cell%n_xyz(1)

    do j = 2 , frc_cell%n_xyz(2) - 1 ! without top and bottom cells

      do k = 1 , frc_cell%n_xyz(3)

        imap = ( get_3d_icell(i,j,k,frc_cell%n_xyz) - 1 ) * frc_cell%nc ! shift in the array

        frc_cell%maps( imap + 1 )  = get_3d_icell(i,j,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 2 )  = get_3d_icell(i+1,j,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 3 )  = get_3d_icell(i-1,j,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 4 )  = get_3d_icell(i,j+1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 5 )  = get_3d_icell(i+1,j+1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 6 )  = get_3d_icell(i-1,j+1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 7 )  = get_3d_icell(i,j-1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 8 )  = get_3d_icell(i+1,j-1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 9 )  = get_3d_icell(i-1,j-1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 10 ) = get_3d_icell(i,j,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 11 ) = get_3d_icell(i+1,j,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 12 ) = get_3d_icell(i-1,j,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 13 ) = get_3d_icell(i,j+1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 14 ) = get_3d_icell(i+1,j+1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 15 ) = get_3d_icell(i-1,j+1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 16 ) = get_3d_icell(i,j-1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 17 ) = get_3d_icell(i+1,j-1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 18 ) = get_3d_icell(i-1,j-1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 19 ) = get_3d_icell(i,j,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 20 ) = get_3d_icell(i+1,j,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 21 ) = get_3d_icell(i-1,j,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 22 ) = get_3d_icell(i,j+1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 23 ) = get_3d_icell(i+1,j+1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 24 ) = get_3d_icell(i-1,j+1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 25 ) = get_3d_icell(i,j-1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 26 ) = get_3d_icell(i+1,j-1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 27 ) = get_3d_icell(i-1,j-1,k+1,frc_cell%n_xyz)


      end do

    end do

  end do



  ! only bottom
  j = 1
  do i = 1 , frc_cell%n_xyz(1)

      do k = 1 , frc_cell%n_xyz(3)

        imap = ( get_3d_icell(i,j,k,frc_cell%n_xyz) - 1 ) * frc_cell%nc ! shift in the array

        frc_cell%maps( imap + 1 )  = get_3d_icell(i,j,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 2 )  = get_3d_icell(i+1,j,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 3 )  = get_3d_icell(i-1,j,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 4 )  = get_3d_icell(i,j+1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 5 )  = get_3d_icell(i+1,j+1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 6 )  = get_3d_icell(i-1,j+1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 7 )  = get_3d_icell(i,j,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 8 )  = get_3d_icell(i+1,j,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 9 )  = get_3d_icell(i-1,j,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 10 ) = get_3d_icell(i,j+1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 11 ) = get_3d_icell(i+1,j+1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 12 ) = get_3d_icell(i-1,j+1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 13 ) = get_3d_icell(i,j,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 14 ) = get_3d_icell(i+1,j,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 15 ) = get_3d_icell(i-1,j,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 16 ) = get_3d_icell(i,j+1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 17 ) = get_3d_icell(i+1,j+1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 18 ) = get_3d_icell(i-1,j+1,k+1,frc_cell%n_xyz)

      end do

  end do


  ! only top
  j = frc_cell%n_xyz(2)
  do i = 1 , frc_cell%n_xyz(1)

      do k = 1 , frc_cell%n_xyz(3)

        imap = ( get_3d_icell(i,j,k,frc_cell%n_xyz) - 1 ) * frc_cell%nc ! shift in the array

        frc_cell%maps( imap + 1 )  = get_3d_icell(i,j,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 2 )  = get_3d_icell(i+1,j,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 3 )  = get_3d_icell(i-1,j,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 4 )  = get_3d_icell(i,j-1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 5 )  = get_3d_icell(i+1,j-1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 6 )  = get_3d_icell(i-1,j-1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 7 )  = get_3d_icell(i,j,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 8 )  = get_3d_icell(i+1,j,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 9 )  = get_3d_icell(i-1,j,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 10 ) = get_3d_icell(i,j-1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 11 ) = get_3d_icell(i+1,j-1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 12 ) = get_3d_icell(i-1,j-1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 13 ) = get_3d_icell(i,j,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 14 ) = get_3d_icell(i+1,j,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 15 ) = get_3d_icell(i-1,j,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 16 ) = get_3d_icell(i,j-1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 17 ) = get_3d_icell(i+1,j-1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 18 ) = get_3d_icell(i-1,j-1,k+1,frc_cell%n_xyz)

      end do

  end do



return
end subroutine construct_3dcell_map



subroutine update_shear_3dcell_map()
  implicit none
!***************************************************************************************************
! 
!***************************************************************************************************

  integer                              :: i,j,k,imap,dx,jcel




  dx = floor( box_delrx / frc_cell%diam(1) ) ! flow in x direction

   j = 1 ! bottom boundary

  do i = 1 , frc_cell%n_xyz(1)

      do k = 1 , frc_cell%n_xyz(3)

        imap = ( get_3d_icell(i,j,k,frc_cell%n_xyz) - 1 ) * frc_cell%nc ! shift in the array

        jcel = get_3d_icell(i+dx,j-1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 19 )  = jcel

        jcel = get_3d_icell(i+1+dx,j-1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 20 )  = jcel

        jcel = get_3d_icell(i-1+dx,j-1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 21 )  = jcel

        jcel = get_3d_icell(i+2+dx,j-1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 22 )  = jcel

        jcel = get_3d_icell(i+dx,j-1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 23 )  = jcel

        jcel = get_3d_icell(i+1+dx,j-1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 24 )  = jcel

        jcel = get_3d_icell(i-1+dx,j-1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 25 )  = jcel

        jcel = get_3d_icell(i+2+dx,j-1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 26 )  = jcel

        jcel = get_3d_icell(i+dx,j-1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 27 )  = jcel

        jcel = get_3d_icell(i+1+dx,j-1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 28 )  = jcel

        jcel = get_3d_icell(i-1+dx,j-1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 29 )  = jcel

        jcel = get_3d_icell(i+2+dx,j-1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 30 )  = jcel

      end do

  end do





  j = frc_cell%n_xyz(2) ! upper boundary


  do i = 1 , frc_cell%n_xyz(1)

      do k = 1 , frc_cell%n_xyz(3)


        imap = ( get_3d_icell(i,j,k,frc_cell%n_xyz) - 1 ) * frc_cell%nc ! shift in the array

        jcel = get_3d_icell(i-dx,j+1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 19 )  = jcel

        jcel = get_3d_icell(i+1-dx,j+1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 20 )  = jcel

        jcel = get_3d_icell(i-1-dx,j+1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 21 )  = jcel

        jcel = get_3d_icell(i-2-dx,j+1,k,frc_cell%n_xyz)
        frc_cell%maps( imap + 22 )  = jcel

        jcel = get_3d_icell(i-dx,j+1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 23 )  = jcel

        jcel = get_3d_icell(i+1-dx,j+1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 24 )  = jcel

        jcel = get_3d_icell(i-1-dx,j+1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 25 )  = jcel

        jcel = get_3d_icell(i-2-dx,j+1,k-1,frc_cell%n_xyz)
        frc_cell%maps( imap + 26 )  = jcel

        jcel = get_3d_icell(i-dx,j+1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 27 )  = jcel

        jcel = get_3d_icell(i+1-dx,j+1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 28 )  = jcel

        jcel = get_3d_icell(i-1-dx,j+1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 29 )  = jcel

        jcel = get_3d_icell(i-2-dx,j+1,k+1,frc_cell%n_xyz)
        frc_cell%maps( imap + 30 )  = jcel


      end do

  end do



return
end subroutine update_shear_3dcell_map




function get_particle_cell(pos) result(icel)
  implicit none
!***************************************************************************************************
!
!***************************************************************************************************

  real(kind=dbl),dimension(:),intent(in)     :: pos
  integer                                    :: icel
  integer,dimension(size(pos))               :: id

  id(:) = floor( pos(:) / frc_cell%diam(:) )
  icel  = frc_cell%n_xyz(1)*id(2) + id(1) + 1 
  if ( ndim .eq. 3 ) icel  = icel + frc_cell%n_xyz(1)*frc_cell%n_xyz(2)*id(3)

end function get_particle_cell



end module cell