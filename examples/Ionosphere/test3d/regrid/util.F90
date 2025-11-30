module util_module

  use prec,only:rp

  implicit none

  interface ismember
    module procedure ismember_char,ismember_real_1d,ismember_real_2d
  endinterface ismember

  interface find_index
    module procedure find_index_char,find_index_real_1d,find_index_real_2d
  endinterface find_index

  interface unique
    module procedure unique_1d,unique_2d
  endinterface unique

  contains
!-----------------------------------------------------------------------
  elemental function isclose(a,b,rtol,atol) result(flag)
! check whether a and b are numerically close

    real(kind=rp),intent(in) :: a,b
    real(kind=rp),intent(in),optional :: rtol,atol
    logical :: flag

    real(kind=rp) :: reltol,abstol

    if (present(rtol)) then
      reltol = rtol
    else
      reltol = 1e-9_rp
    endif
    if (present(atol)) then
      abstol = atol
    else
      abstol = 0
    endif

    flag = abs(a-b) <= max(reltol * max(abs(a), abs(b)), abstol)

  endfunction isclose
!-----------------------------------------------------------------------
  pure function ismember_char(group,instance) result(flag)
! test if a string is contained in another string array

    character(len=*),dimension(:),intent(in) :: group
    character(len=*),intent(in) :: instance
    logical :: flag

    integer :: i

    flag = .false.
    do i = 1,size(group)
      if (len_trim(instance)>0 .and. len_trim(group(i))>0 .and. &
          trim(instance)==trim(group(i))) then
        flag = .true.
        exit
      endif
    enddo

  endfunction ismember_char
!-----------------------------------------------------------------------
  pure function ismember_real_1d(group,instance,rtol,atol) result(flag)

    real(kind=rp),dimension(:),intent(in) :: group
    real(kind=rp),intent(in) :: instance
    real(kind=rp),intent(in),optional :: rtol,atol
    logical :: flag

    integer :: i
    real(kind=rp) :: reltol,abstol

    if (present(rtol)) then
      reltol = rtol
    else
      reltol = 1e-9_rp
    endif
    if (present(atol)) then
      abstol = atol
    else
      abstol = 0
    endif

    flag = .false.
    do i = 1,size(group)
      if (isclose(instance,group(i),reltol,abstol)) then
        flag = .true.
        exit
      endif
    enddo

  endfunction ismember_real_1d
!-----------------------------------------------------------------------
  pure function ismember_real_2d(group,instance,rtol,atol) result(flag)

    real(kind=rp),dimension(:,:),intent(in) :: group
    real(kind=rp),dimension(size(group,1)),intent(in) :: instance
    real(kind=rp),intent(in),optional :: rtol,atol
    logical :: flag

    integer :: i
    real(kind=rp) :: reltol,abstol

    if (present(rtol)) then
      reltol = rtol
    else
      reltol = 1e-9_rp
    endif
    if (present(atol)) then
      abstol = atol
    else
      abstol = 0
    endif

    flag = .false.
    do i = 1,size(group,2)
      if (all(isclose(instance,group(:,i),reltol,abstol))) then
        flag = .true.
        exit
      endif
    enddo

  endfunction ismember_real_2d
!-----------------------------------------------------------------------
  pure function find_index_char(group,instance) result(idx)
! find the index of a string in another string array

    character(len=*),dimension(:),intent(in) :: group
    character(len=*),intent(in) :: instance
    integer :: idx

    do idx = 1,size(group)
      if (len_trim(instance)>0 .and. len_trim(group(idx))>0 .and. &
          trim(instance)==trim(group(idx))) return
    enddo
    idx = 0

  endfunction find_index_char
!-----------------------------------------------------------------------
  pure function find_index_real_1d(group,instance,rtol,atol) result(idx)

    real(kind=rp),dimension(:),intent(in) :: group
    real(kind=rp),intent(in) :: instance
    real(kind=rp),intent(in),optional :: rtol,atol
    integer :: idx

    real(kind=rp) :: reltol,abstol

    if (present(rtol)) then
      reltol = rtol
    else
      reltol = 1e-9_rp
    endif
    if (present(atol)) then
      abstol = atol
    else
      abstol = 0
    endif

    do idx = 1,size(group)
      if (isclose(instance,group(idx),reltol,abstol)) return
    enddo
    idx = 0

  endfunction find_index_real_1d
!-----------------------------------------------------------------------
  pure function find_index_real_2d(group,instance,rtol,atol) result(idx)

    real(kind=rp),dimension(:,:),intent(in) :: group
    real(kind=rp),dimension(size(group,1)),intent(in) :: instance
    real(kind=rp),intent(in),optional :: rtol,atol
    integer :: idx

    real(kind=rp) :: reltol,abstol

    if (present(rtol)) then
      reltol = rtol
    else
      reltol = 1e-9_rp
    endif
    if (present(atol)) then
      abstol = atol
    else
      abstol = 0
    endif

    do idx = 1,size(group,2)
      if (all(isclose(instance,group(:,idx),reltol,abstol))) return
    enddo
    idx = 0

  endfunction find_index_real_2d
!-----------------------------------------------------------------------
  pure subroutine unique_1d(x,n,uniq,rtol,atol)

    real(kind=rp),dimension(:),intent(in) :: x
    integer,intent(out) :: n
    real(kind=rp),dimension(size(x)),intent(out) :: uniq
    real(kind=rp),intent(in),optional :: rtol,atol

    integer :: i
    real(kind=rp) :: reltol,abstol

    if (present(rtol)) then
      reltol = rtol
    else
      reltol = 1e-9_rp
    endif
    if (present(atol)) then
      abstol = atol
    else
      abstol = 0
    endif

    n = 0
    do i = 1,size(x)
      if (.not. ismember_real_1d(uniq(1:n),x(i),reltol,abstol)) then
        n = n+1
        uniq(n) = x(i)
      endif
    enddo

  endsubroutine unique_1d
!-----------------------------------------------------------------------
  pure subroutine unique_2d(x,n,uniq,rtol,atol)

    real(kind=rp),dimension(:,:),intent(in) :: x
    integer,intent(out) :: n
    real(kind=rp),dimension(size(x,1),size(x,2)),intent(out) :: uniq
    real(kind=rp),intent(in),optional :: rtol,atol

    integer :: i
    real(kind=rp) :: reltol,abstol

    if (present(rtol)) then
      reltol = rtol
    else
      reltol = 1e-9_rp
    endif
    if (present(atol)) then
      abstol = atol
    else
      abstol = 0
    endif

    n = 0
    do i = 1,size(x,2)
      if (.not. ismember_real_2d(uniq(:,1:n),x(:,i),reltol,abstol)) then
        n = n+1
        uniq(:,n) = x(:,i)
      endif
    enddo

  endsubroutine unique_2d
!-----------------------------------------------------------------------
  pure function rotate_s2c(theta,phi,sph) result(cart)
! https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates

    real(kind=rp),intent(in) :: theta,phi
    real(kind=rp),dimension(3),intent(in) :: sph
    real(kind=rp),dimension(3) :: cart

    real(kind=rp),dimension(3,3) :: A

    A(1,:) = (/sin(theta)*cos(phi), cos(theta)*cos(phi),-sin(phi)/)
    A(2,:) = (/sin(theta)*sin(phi), cos(theta)*sin(phi), cos(phi)/)
    A(3,:) = (/cos(theta)         ,-sin(theta)         , 0.0_rp  /)
    cart = matmul(A,sph)

  endfunction rotate_s2c
!-----------------------------------------------------------------------
endmodule util_module
