program driver

  use prec,only:rp
  use mesh_module,only:nnode,nalt,nelement,ncell,nnode_per_cell, &
    alt,lat_int=>lat,lon_int=>lon,latr,lonr,elements, &
    read_coord,setup_nodes,setup_elements,reorder_in,reorder_out
  use input_module,only:read_data, &
    nlon,nlat,nlev,ntime,nvar,time, &
    lon_ext=>lon,lat_ext=>lat,z,var_ext=>variable
  use output_module,only:write_data
  use regrid_module,only:regrid, &
    init_esmf,init_ext,init_int,init_regrid
  use igrf_module,only:igrf14syn
  use interp_module,only:interp1d
  use util_module,only:rotate_s2c

  implicit none

  integer,parameter :: isv = 0, itype = 1, nscalar = 5, nvector = nvar-nscalar
  real(kind=8),parameter :: date = 2002 + 80/365.0_8
  real(kind=rp),parameter :: pi = 4*atan(1.0_rp), dtr = pi/180
  character(len=*),dimension(*),parameter :: &
    varname = (/'N ','T ','NU','P ','L ','H ','C ', &
      'VX','VY','VZ','EX','EY','EZ','UX','UY','UZ'/)
  integer :: k0,k1,inode,i,j,k,t,ivar,ivx,ivy,ivz,iex,iey,iez
  real(kind=8) :: colat,elong,h,bx,by,bz,f
  real(kind=rp) :: theta,phi
  real(kind=rp),dimension(3) :: v,vperp
  real(kind=rp),dimension(:,:),allocatable :: b2,vector
  real(kind=rp),dimension(:,:,:),allocatable :: b
  real(kind=rp),dimension(:,:,:,:),allocatable :: var_int
  real(kind=rp),dimension(:,:,:,:,:),allocatable :: var_ext_alt
  integer,dimension(:,:),allocatable :: cellidx_in,nodeidx_in,nodeidx_out,altidx_out

  call init_esmf

  call read_coord('coord.nc')
  call setup_nodes
  call setup_elements

  call read_data('/glade/work/haonan/ionmodel.nc')

  call init_ext(nlon,nlat,nalt,ntime,nvar,lon_ext,lat_ext)
  call init_int(nnode,nalt,ntime,nvar,nelement,lon_int,lat_int,elements)
  call init_regrid

  allocate(b2(nnode,nalt))
  allocate(b(nnode,nalt,3))
  do inode = 1,nnode
    colat = 90-lat_int(inode)
    if (lon_int(inode) < 0) then
      elong = lon_int(inode)+360
    else
      elong = lon_int(inode)
    endif
    theta = pi/2-latr(inode)
    phi = lonr(inode)
    do k = 1,nalt
      h = alt(k)
      call igrf14syn(isv,date,itype,h,colat,elong,bx,by,bz,f)
      b2(inode,k) = f**2
      b(inode,k,:) = rotate_s2c(theta,phi,real((/-bz,-bx,by/),kind=rp))
    enddo
  enddo

  allocate(var_ext_alt(nlon,nlat,nalt,ntime,nvar))
  do concurrent (i = 1:nlon, j = 1:nlat, t = 1:ntime) ! no extrapolation of ion density for now
    do k0 = 1,nalt
      if (alt(k0) >= z(i,j,1,t)) exit
    enddo
    do k1 = nalt,1,-1
      if (alt(k1) <= z(i,j,nlev,t)) exit
    enddo
    var_ext_alt(i,j,k0:k1,t,1) = interp1d(alt(k0:k1),z(i,j,:,t),var_ext(i,j,:,t,1))
    var_ext_alt(i,j,1:k0-1,t,1) = var_ext(i,j,1,t,1)
    var_ext_alt(i,j,k1+1:nalt,t,1) = var_ext(i,j,nlev,t,1)
  enddo
  do concurrent (i = 1:nlon, j = 1:nlat, t = 1:ntime, ivar = 2:nscalar)
    var_ext_alt(i,j,:,t,ivar) = exp(interp1d(alt,z(i,j,:,t),log(var_ext(i,j,:,t,ivar))))
  enddo

  allocate(vector(nalt,nvector))
  do concurrent (i = 1:nlon, j = 1:nlat, t = 1:ntime)
    do ivar = 1,nvector
      vector(:,ivar) = interp1d(alt,z(i,j,:,t),var_ext(i,j,:,t,nscalar+ivar))
    enddo
    theta = pi/2-lat_ext(j)*dtr
    phi = lon_ext(i)*dtr
    do concurrent (k = 1:nalt, ivar = 1:nvector/3)
      v = rotate_s2c(theta,phi,(/vector(k,ivar*3),-vector(k,ivar*3-1),vector(k,ivar*3-2)/))
      var_ext_alt(i,j,k,t,nscalar+ivar*3-2) = v(1)
      var_ext_alt(i,j,k,t,nscalar+ivar*3-1) = v(2)
      var_ext_alt(i,j,k,t,nscalar+ivar*3) = v(3)
    enddo
  enddo

  allocate(var_int(nnode,nalt,ntime,nvar))
  var_int = regrid(nnode,nlon,nlat,nalt,ntime,nvar,var_ext_alt)

! zero out the parallel component of drift velocity and electric field
! this might not be necessary
  do ivar = 1,nvar
    if (trim(varname(ivar)) == 'VX') ivx = ivar
    if (trim(varname(ivar)) == 'VY') ivy = ivar
    if (trim(varname(ivar)) == 'VZ') ivz = ivar
    if (trim(varname(ivar)) == 'EX') iex = ivar
    if (trim(varname(ivar)) == 'EY') iey = ivar
    if (trim(varname(ivar)) == 'EZ') iez = ivar
  enddo

  do concurrent (inode = 1:nnode, k = 1:nalt, t = 1:ntime)
    v = (/var_int(inode,k,t,ivx),var_int(inode,k,t,ivy),var_int(inode,k,t,ivz)/)
    vperp = v-dot_product(v,b(inode,k,:))*b(inode,k,:)/b2(inode,k)
    var_int(inode,k,t,ivx) = vperp(1)
    var_int(inode,k,t,ivy) = vperp(2)
    var_int(inode,k,t,ivz) = vperp(3)

    v = (/var_int(inode,k,t,iex),var_int(inode,k,t,iey),var_int(inode,k,t,iez)/)
    vperp = v-dot_product(v,b(inode,k,:))*b(inode,k,:)/b2(inode,k)
    var_int(inode,k,t,iex) = vperp(1)
    var_int(inode,k,t,iey) = vperp(2)
    var_int(inode,k,t,iez) = vperp(3)
  enddo

  allocate(cellidx_in(nnode,nalt))
  allocate(nodeidx_in(nnode,nalt))
  call reorder_in(cellidx_in,nodeidx_in)

  allocate(nodeidx_out(ncell,nnode_per_cell))
  allocate(altidx_out(ncell,nnode_per_cell))
  call reorder_out(nodeidx_out,altidx_out)

  call write_data(nnode,nalt,ntime,nvar,ncell,nnode_per_cell, &
    lon_int,lat_int,alt,time,b*1e-9_rp,var_int, &
    cellidx_in,nodeidx_in,nodeidx_out,altidx_out, &
    'data.nc',varname)

endprogram driver
