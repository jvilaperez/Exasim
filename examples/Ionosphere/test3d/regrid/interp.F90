module interp_module

  use prec,only:rp

  implicit none

! 2D index increases faster in latitudes
! (stack along latitudes first then longitudes)

! 0-based indexing is used to simplify index calculation (no wrap needed)
! (nlat=3)
!  8  9 10 11
!  4  5  6  7
!  0  1  2  3 (nlon=4)

! when converting from 1D to 2D index
! a +1 is applied to convert from 0-based to 1-based indexing

! when converting from 2D to 1D index
! a -1 is applied to convert from 1-based to 0-based indexing

! for 2D interpolation from an orthogonal mesh to a non-orthogonal mesh,
! first find the 4-point bounding box of the non-orthogonal mesh in the orthogonal mesh,
! then apply inverse distance weighting based on the 4-point bounding box.
! a virtual pole of the orthogonal mesh is constructed based on the average of the first ring.

  contains
!-----------------------------------------------------------------------
  pure function interp2d(fp,box,wt) result(f)
! 2D interpolation for structured grid (general purpose)
! fp (input): values at input locations
! box (input): bounding box of the output location
! wt (input): weights of each points in the bounding box
! f (output): value at output location

    real(kind=rp),dimension(:,:),intent(in) :: fp
    integer,dimension(4),intent(in) :: box
    real(kind=rp),dimension(4),intent(in) :: wt
    real(kind=rp) :: f

    integer :: nlatp,nlonp,n,latidx,lonidx
    real(kind=rp) :: sp,np,fp0

    nlatp = size(fp,1)
    nlonp = size(fp,2)

! construct virtual poles
    sp = sum(fp(1,:))/nlonp
    np = sum(fp(nlatp,:))/nlonp

    f = 0
    do n = 1,4

! poles are double counted, but one of the weights is zero
      if (box(n) == -1) then
        fp0 = sp
      elseif (box(n) == -2) then
        fp0 = np
      else
        latidx = box(n)/nlonp+1
        lonidx = modulo(box(n),nlonp)+1
        fp0 = fp(latidx,lonidx)
      endif

      f = f+wt(n)*fp0
    enddo

  endfunction interp2d
!-----------------------------------------------------------------------
  pure function weight(lat,lon,latp,lonp,tol,box) result(wt)
! calculate weights for 2D interpolation on a sphere (inverse distance weighting)
! lat,lon (input): output latitude and longitude in radian
! latp,lonp (input): input latitudes and longitudes in radian
!   required to be monotonic (increasing or decreasing)
! tol (input, optional): tolerance of distance
!   below which two points are treated as overlapped
! box (input, optional): 4-point bounding box of the output location
! wt (output): weights of each points in the bounding box

! don't include periodic points in lonp

    real(kind=rp),intent(in) :: lat,lon
    real(kind=rp),dimension(:),intent(in) :: latp,lonp
    real(kind=rp),intent(in),optional :: tol
    integer,dimension(4),intent(in),optional :: box
    real(kind=rp),dimension(4) :: wt

    real(kind=rp),parameter :: pi = 4*atan(1.0_rp)
    integer :: nlatp,nlonp,n,latidx,lonidx
    integer,dimension(4) :: box_in
    real(kind=rp) :: tol_in,d
    real(kind=rp),dimension(4) :: invdist

    nlatp = size(latp)
    nlonp = size(lonp)

    if (present(tol)) then
      tol_in = tol
    else
      tol_in = 1e-9_rp
    endif

    if (present(box)) then
      box_in = box
    else
      box_in = bound(lat,lon,latp,lonp)
    endif

! initialize all weights to zero
    invdist = 0

! near poles, add pole to the bounding box
    if (box_in(3)==-1 .or. box_in(4)==-1 .or. &
        box_in(3)==-2 .or. box_in(4)==-2) then
      do n = 1,2
        latidx = box_in(n)/nlonp+1
        lonidx = modulo(box_in(n),nlonp)+1
        d = distance(lat,lon,latp(latidx),lonp(lonidx))
        if (d > tol_in) invdist(n) = 1/d
      enddo

      if (box_in(3)==-1 .or. box_in(4)==-1) then ! south pole
        d = distance(lat,lon,-pi/2,0.0_rp)
        if (d > tol_in) invdist(3) = 1/d
      else ! north pole
        d = distance(lat,lon,pi/2,0.0_rp)
        if (d > tol_in) invdist(3) = 1/d
      endif
    else ! internal latitudes
      do n = 1,4
        latidx = box_in(n)/nlonp+1
        lonidx = modulo(box_in(n),nlonp)+1
        d = distance(lat,lon,latp(latidx),lonp(lonidx))
        if (d > tol_in) invdist(n) = 1/d
      enddo
    endif

    wt = invdist/sum(invdist)

  endfunction weight
!-----------------------------------------------------------------------
  pure function bound(lat,lon,latp,lonp) result(box)
! given coordinates, find the bounding box of a mesh in mesh p

    real(kind=rp),intent(in) :: lat,lon
    real(kind=rp),dimension(:),intent(in) :: latp,lonp
    integer,dimension(4) :: box

    real(kind=rp),parameter :: pi = 4*atan(1.0_rp), period = pi*2
    integer :: nlatp,nlonp,ilat,ilon,i0,i1
    real(kind=rp) :: lonwrap

    nlatp = size(latp)
    nlonp = size(lonp)

    lonwrap = wrap(lon,lonp(1),period) ! wrap longitude
    ilon = find(lonp,lonwrap)

    if (ilon == 0) then
      i0 = nlonp
    else
      i0 = ilon
    endif
    if (ilon == nlonp) then
      i1 = 1
    else
      i1 = ilon+1
    endif

    ilat = find(latp,lat)

    if (ilat == 0) then ! bounded by the south pole
      box(1) = i0-1
      box(2) = i1-1
      box(3:4) = -1
    elseif (ilat == nlatp) then ! bounded by the north pole
      box(1) = (nlatp-1)*nlonp+i0-1
      box(2) = (nlatp-1)*nlonp+i1-1
      box(3:4) = -2
    else
      box(1) = (ilat-1)*nlonp+i0-1
      box(2) = (ilat-1)*nlonp+i1-1
      box(3) = ilat*nlonp+i1-1
      box(4) = ilat*nlonp+i0-1
    endif

  endfunction bound
!-----------------------------------------------------------------------
  pure function interp1d(x,xp,fp,period,left,right) result(f)
! 1D linear interpolation, linear extrapolation if out of bound
! x (input): output locations
! xp (input): input locations
!   required to be monotonic (increasing or decreasing)
! fp (input): values at input locations
! f (output): values at output locations
! period (input,optional): period of xp if it is a circular coordinate
! left,right (input,optional): extrapolation flags if out of bound

! if period is present, then xp range should be smaller than period

    real(kind=rp),dimension(:),intent(in) :: x,xp
    real(kind=rp),dimension(size(xp)),intent(in) :: fp
    real(kind=rp),intent(in),optional :: period
    logical,intent(in),optional :: left,right
    real(kind=rp),dimension(size(x)) :: f

    integer :: nx,nxp,i,ix,i0,i1
    real(kind=rp) :: xwrap,dxp

    nx = size(x)
    nxp = size(xp)

    if (present(period)) then
      do i = 1,nx
        xwrap = wrap(x(i),xp(1),period) ! wrap around
        ix = find(xp,xwrap)

! wrap around if the query location is out of bound
        if (ix==0 .or. ix==nxp) then
          i0 = nxp
          i1 = 1
          dxp = xp(1)+period-xp(nxp)
        else
          i0 = ix
          i1 = ix+1
          dxp = xp(ix+1)-xp(ix)
        endif

        f(i) = ((xp(i1)-x(i))*fp(i0) + (x(i)-xp(i0))*fp(i1)) / dxp
      enddo
    else
      do i = 1,nx
        ix = find(xp,x(i))

! extrapolate based on flags if the query location is out of bound
        if (ix == 0) then
          i0 = 1
          i1 = 1
          if (present(left)) then
            if (left) i1 = 2
          endif
        elseif (ix == nxp) then
          i0 = nxp
          i1 = nxp
          if (present(right)) then
            if (right) i0 = nxp-1
          endif
        else
          i0 = ix
          i1 = ix+1
        endif
        dxp = xp(i1)-xp(i0)

        if (i0 == i1) then
          f(i) = fp(i0)
        else
          f(i) = ((xp(i1)-x(i))*fp(i0) + (x(i)-xp(i0))*fp(i1)) / dxp
        endif
      enddo
    endif

  endfunction interp1d
!-----------------------------------------------------------------------
  pure function find(x,x0) result(i)
! 1D binary search
! x (input): value array, required to be monotonic
! x0 (input): value to be found
! i (output): the index of x0 in x
!   satisfying x(i)<=x0<x(i+1) (if x is increasing)
!   or         x(i)>=x0>x(i+1) (if x is decreasing)
!   0 and size(x) indicate x0 is out of the range of x

! The monotonicity of x is not enforced in this function.
! Caller should ensure x is in order (increasing or decreasing),
! otherwise the result is meaningless.

    real(kind=rp),dimension(:),intent(in) :: x
    real(kind=rp),intent(in) :: x0
    integer :: i

    integer :: nx,i0,i1

    nx = size(x)

! x is in increasing order
    if (x(1) < x(nx)) then
      if (x0 < x(1)) then
        i = 0
      elseif (x0 >= x(nx)) then
        i = nx
      else
        i0 = 1
        i1 = nx
        do while (i0+1 < i1)
          i = (i0+i1)/2
          if (x(i) <= x0) then
            i0 = i
          else
            i1 = i
          endif
        enddo
        i = (i0+i1)/2
      endif

! x is in decreasing order
    else
      if (x0 > x(1)) then
        i = 0
      elseif (x0 <= x(nx)) then
        i = nx
      else
        i0 = 1
        i1 = nx
        do while (i0+1 < i1)
          i = (i0+i1)/2
          if (x(i) >= x0) then
            i0 = i
          else
            i1 = i
          endif
        enddo
        i = (i0+i1)/2
      endif
    endif

  endfunction find
!-----------------------------------------------------------------------
  elemental function wrap(x,x0,period) result(xwrap)
! given x, return wrapped x within [x0, x0+period)

    real(kind=rp),intent(in) :: x,x0,period
    real(kind=rp) :: xwrap

    integer :: n

    if (x < x0) then
      n = floor((x0-x)/period)+1
      xwrap = x+period*n
    elseif (x >= x0+period) then
      n = floor((x-x0-period)/period)+1
      xwrap = x-period*n
    else
      xwrap = x
    endif

  endfunction wrap
!-----------------------------------------------------------------------
  elemental function distance(lat1,lon1,lat2,lon2) result(d)
! calculate the distance of two points on a sphere (great circle distance)
! lat1,lat2 (input): longitudes in radian
! lon1,lon2 (input): latitudes (not polar angle) in radian
! d (output): great circle distance in radian

    real(kind=rp),intent(in) :: lat1,lon1,lat2,lon2
    real(kind=rp) :: d

    real(kind=rp) :: c1,s1,c2,s2,dlon,clon,slon

    c1 = cos(lat1)
    s1 = sin(lat1)
    c2 = cos(lat2)
    s2 = sin(lat2)
    dlon = lon1-lon2
    clon = cos(dlon)
    slon = sin(dlon)
    d = atan2(sqrt((c2*slon)**2+(c1*s2-s1*c2*clon)**2),s1*s2+c1*c2*clon)

  endfunction distance
!-----------------------------------------------------------------------
endmodule interp_module
