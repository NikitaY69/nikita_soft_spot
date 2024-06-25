
module system
  ! use omp_lib
  use kinds
!***************************************************************************************************
!
!***************************************************************************************************
  implicit none
  public
!***************************************************************************************************


  ! openMP
  integer                                         :: nthreads

  ! dimensions
  integer                                         :: ndim
  integer                                         :: npart
  integer                                         :: ntypes
  real(kind=dbl)                                  :: inv_ndim
  real(kind=dbl)                                  :: real_npart

  ! box parameters
  real(kind=dbl)                                  :: box_volume
  real(kind=dbl),dimension(:),allocatable         :: boxl,boxl2
  real(kind=dbl)                                  :: box_strain,box_delrx
  real(kind=dbl)                                  :: box_rho


  ! particle arrays
  integer,dimension(:),allocatable                :: ptypes
  integer,dimension(:),allocatable                :: mtypes
  real(kind=dbl),dimension(:),allocatable         :: mass
  real(kind=dbl),dimension(:),allocatable         :: d
  real(kind=dbl),dimension(:,:),allocatable       :: r
  real(kind=dbl),dimension(:,:),allocatable       :: vel
  real(kind=dbl),dimension(:,:),allocatable       :: f
  real(kind=dbl),dimension(:,:,:),allocatable     :: s
  real(kind=dbl)                                  :: typical_grad,grad_scale


  ! current thermo observables
  real(kind=dbl)                                  :: thermo_temp
  real(kind=dbl)                                  :: thermo_pot,thermo_pre,thermo_sigma



!***************************************************************************************************

contains



subroutine init_system_arrays()
  write(*,*) '=> init system arrays'

! !$omp parallel shared(nthreads)
! !$omp master

!   nthreads = omp_get_num_threads()
!   write(*,*) 'Running OpenMP version using ',nthreads,' thread(s).'

! !$omp end master
! !$omp end parallel


  if( allocated(boxl) .eqv. .false. ) then
    allocate(boxl(ndim),boxl2(ndim))
    allocate(r(ndim,npart),f(ndim,npart),d(npart),vel(ndim,npart),mass(npart))
    allocate(ptypes(npart))
  end if

  inv_ndim    = 1._dbl/real(ndim,kind=dbl)
  boxl(:)     = 0._dbl
  boxl2(:)    = 0._dbl
  r(:,:)      = 0._dbl
  vel(:,:)    = 0._dbl
  f(:,:)      = 0._dbl
  d(:)        = 1._dbl
  ptypes(:)   = 1
  mass(:)     = 1._dbl

  real_npart = real(npart,kind=dbl)

end subroutine init_system_arrays






subroutine print_system_info()
  write(*,*)'---------------------------------------------------------------------------'
  write(*,*) '=> print system info:'
  write(*,*) 'd=',ndim,', N=',npart
  write(*,*) 'density      ',real(box_rho,kind=sgl)
  write(*,*) 'strain       ',real(box_strain,kind=sgl)
  write(*,*) 'box geometry ',real(boxl,kind=sgl)
  write(*,*)'---------------------------------------------------------------------------'
end subroutine print_system_info



subroutine switch_diameter(i,j)
  implicit none
  integer,intent(in) :: i,j
  real(kind=dbl)     :: temp
  temp = d(j)
  d(j) = d(i)
  d(i) = temp
end subroutine switch_diameter




end module system
