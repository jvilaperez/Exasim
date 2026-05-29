program driver

  use prec,only:rp
  use mesh_module,only:ncell,nnode_per_cell, &
    nnode_int=>nnode,nelement_int=>nelement,nalt_int=>nalt, &
    elements_int=>elements,alt_int=>alt, &
    lon_int=>lon,lat_int=>lat,lonr,latr, &
    read_coord,setup_nodes,setup_elements, &
    reorder_in,reorder_out
  use input_module,only:read_data,nvar, &
    nnode_ext=>nnode,nelement_ext=>nelement,nalt_ext=>nalt, &
    elements_ext=>elements,var_ext=>variable,alt_ext=>alt, &
    nodeGlon,nodeGlat,elementGlon,elementGlat
  use output_module,only:write_data
  use regrid_module,only:regrid, &
    init_esmf,init_ext,init_int,init_regrid
  use igrf_module,only:igrf14syn
  use interp_module,only:interp1d
  use util_module,only:rotate_s2c

  implicit none

  integer,parameter :: isv = 0, itype = 1, timeidx = 7, &
    nscalar = 30, nvector = nvar-nscalar
  real(kind=8),parameter :: date = 2002 + 80/365.0_8
  real(kind=rp),parameter :: pi = 4*atan(1.0_rp), dtr = pi/180
  character(len=*),dimension(*),parameter :: varname = &
    (/'OP      ','O2P     ','NP      ','N2P     ','NOP     ', &
      'NU_OP   ','NU_O2P  ','NU_NP   ','NU_N2P  ','NU_NOP  ', &
      'PROD_OP ','PROD_O2P','PROD_NP ','PROD_N2P','PROD_NOP', &
      'LOSS_OP ','LOSS_O2P','LOSS_NP ','LOSS_N2P','LOSS_NOP', &
      'TI      ','TE      ','TN      ','KAPPA_I ','KAPPA_E ', &
      'HEAT_I  ','HEAT_E  ','QEI     ','QIN     ','QEN     ', &
      'UI      ','VI      ','WI      ', &
      'EX      ','EY      ','EZ      ', &
      'UN      ','VN      ','WN      '/)
  integer :: inode,ielem,k,ivar,iex,iey,iez
  real(kind=8) :: colat,elong,h,bx,by,bz,f
  real(kind=rp) :: theta,phi
  real(kind=rp),dimension(3) :: e,eperp
  real(kind=rp),dimension(:,:),allocatable :: b2,vector
  real(kind=rp),dimension(:,:,:),allocatable :: b,var_ext2int,var_int
  integer,dimension(:,:),allocatable :: cellidx_in,nodeidx_in,nodeidx_out,altidx_out
  logical,dimension(:),allocatable :: positive

  call init_esmf

  call read_coord('coord.nc')
  call setup_nodes
  call setup_elements

  call read_data('/glade/work/haonan/ionmodel.nc',timeidx)

  call init_ext(nnode_ext,nelement_ext,nalt_int,nvar, &
    elements_ext,nodeGlon,nodeGlat,elementGlon,elementGlat)
  call init_int(nnode_int,nelement_int,nalt_int,nvar, &
    elements_int,lon_int,lat_int)
  call init_regrid

  allocate(b2(nnode_int,nalt_int))
  allocate(b(nnode_int,nalt_int,3))
  do inode = 1,nnode_int
    colat = 90-lat_int(inode)
    if (lon_int(inode) < 0) then
      elong = lon_int(inode)+360
    else
      elong = lon_int(inode)
    endif
    theta = pi/2-latr(inode)
    phi = lonr(inode)
    do k = 1,nalt_int
      h = alt_int(k)*1e-3_rp
      call igrf14syn(isv,date,itype,h,colat,elong,bx,by,bz,f)
      b2(inode,k) = f**2
      b(inode,k,:) = rotate_s2c(theta,phi,real((/-bz,-bx,by/),kind=rp))
    enddo
  enddo

  allocate(positive(nalt_ext))
  allocate(var_ext2int(nelement_ext,nalt_int,nvar))
  do concurrent (ielem = 1:nelement_ext, ivar = 1:nscalar)
    positive = var_ext(ielem,:,ivar) > 0
    var_ext2int(ielem,:,ivar) = exp(interp1d(alt_int,pack(alt_ext,positive),log(pack(var_ext(ielem,:,ivar),positive))))
  enddo

  allocate(vector(nalt_int,nvector))
  do ielem = 1,nelement_ext
    do ivar = 1,nvector
      vector(:,ivar) = interp1d(alt_int,alt_ext,var_ext(ielem,:,nscalar+ivar))
    enddo
    theta = pi/2-elementGlat(ielem)*dtr
    phi = elementGlon(ielem)*dtr
    do concurrent (k = 1:nalt_int, ivar = 1:nvector/3)
      e = rotate_s2c(theta,phi,(/vector(k,ivar*3),-vector(k,ivar*3-1),vector(k,ivar*3-2)/))
      var_ext2int(ielem,k,nscalar+ivar*3-2) = e(1)
      var_ext2int(ielem,k,nscalar+ivar*3-1) = e(2)
      var_ext2int(ielem,k,nscalar+ivar*3) = e(3)
    enddo
  enddo

  allocate(var_int(nnode_int,nalt_int,nvar))
  var_int = regrid(nelement_ext,nnode_int,nalt_int,nvar,var_ext2int)

! zero out the parallel component of drift velocity and electric field
! this might not be necessary
  do ivar = 1,nvar
    if (trim(varname(ivar)) == 'EX') iex = ivar
    if (trim(varname(ivar)) == 'EY') iey = ivar
    if (trim(varname(ivar)) == 'EZ') iez = ivar
  enddo

  do concurrent (inode = 1:nnode_int, k = 1:nalt_int)
    e = (/var_int(inode,k,iex),var_int(inode,k,iey),var_int(inode,k,iez)/)
    eperp = e-dot_product(e,b(inode,k,:))*b(inode,k,:)/b2(inode,k)
    var_int(inode,k,iex) = eperp(1)
    var_int(inode,k,iey) = eperp(2)
    var_int(inode,k,iez) = eperp(3)
  enddo

  allocate(cellidx_in(nnode_int,nalt_int))
  allocate(nodeidx_in(nnode_int,nalt_int))
  call reorder_in(cellidx_in,nodeidx_in)

  allocate(nodeidx_out(ncell,nnode_per_cell))
  allocate(altidx_out(ncell,nnode_per_cell))
  call reorder_out(nodeidx_out,altidx_out)

  call write_data(nnode_int,nalt_int,nvar,ncell,nnode_per_cell, &
    lon_int,lat_int,alt_int,b*1e-9_rp,var_int, &
    cellidx_in,nodeidx_in,nodeidx_out,altidx_out, &
    'data.nc',varname)

endprogram driver
