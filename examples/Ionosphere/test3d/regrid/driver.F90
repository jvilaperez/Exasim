program driver

  use prec,only:rp
  use mesh_module,only:ncell,nnode_per_cell, &
    nnode_int=>nnode,nelement_int=>nelement,nalt_int=>nalt, &
    elements_int=>elements,alt_int=>alt, &
    lon_int=>lon,lat_int=>lat,lonr,latr, &
    read_coord,setup_nodes,setup_elements, &
    reorder_in,reorder_out
  use input_module,only:read_data,nvar_ext=>nvar, &
    iop_ext=>iop,io2p_ext=>io2p,inp_ext=>inp,in2p_ext=>in2p,inop_ext=>inop, &
    iv1_ext=>iv1,iv2_ext=>iv2,iv3_ext=>iv3,ivp,iti_ext=>iti,ite_ext=>ite, &
    iprod_op_ext=>iprod_op,iprod_o2p_ext=>iprod_o2p, &
    iprod_np_ext=>iprod_np,iprod_n2p_ext=>iprod_n2p,iprod_nop_ext=>iprod_nop, &
    iloss_op_ext=>iloss_op,iloss_o2p_ext=>iloss_o2p, &
    iloss_np_ext=>iloss_np,iloss_n2p_ext=>iloss_n2p,iloss_nop_ext=>iloss_nop, &
    iped,iqsum,ie1_ext=>ie1,ie2_ext=>ie2,ie3_ext=>ie3, &
    itn_ext=>itn,iu1_ext=>iu1,iu2_ext=>iu2,iu3_ext=>iu3,io2,io1,ihe,in2, &
    nnode_ext=>nnode,nelement_ext=>nelement,nlev, &
    elements_ext=>elements,var_ext=>variable,alt_ext=>alt, &
    nodeGlon,nodeGlat,elementGlon,elementGlat
  use output_module,only:write_data,nvar_int=>nvar, &
    iop_int=>iop,io2p_int=>io2p,inp_int=>inp,in2p_int=>in2p,inop_int=>inop, &
    iv1_int=>iv1,iv2_int=>iv2,iv3_int=>iv3,iti_int=>iti,ite_int=>ite, &
    iprod_op_int=>iprod_op,iprod_o2p_int=>iprod_o2p, &
    iprod_np_int=>iprod_np,iprod_n2p_int=>iprod_n2p,iprod_nop_int=>iprod_nop, &
    iloss_op_int=>iloss_op,iloss_o2p_int=>iloss_o2p, &
    iloss_np_int=>iloss_np,iloss_n2p_int=>iloss_n2p,iloss_nop_int=>iloss_nop, &
    inu_op,inu_o2p,inu_np,inu_n2p,inu_nop,ie1_int=>ie1,ie2_int=>ie2,ie3_int=>ie3, &
    iheat_i,iheat_e,iqei,iqin,iqen,ikappa_i,ikappa_e, &
    itn_int=>itn,iu1_int=>iu1,iu2_int=>iu2,iu3_int=>iu3
  use regrid_module,only:regrid, &
    init_esmf,init_ext,init_int,init_regrid
  use igrf_module,only:igrf14syn
  use calculate_module,only:calculate
  use interp_module,only:interp1d
  use util_module,only:rotate_s2c

  implicit none

  integer,parameter :: isv = 0, itype = 1, timeidx = 5
  real(kind=8),parameter :: date = 2002 + 80/365.0_8
  real(kind=rp),parameter :: pi = 4*atan(1.0_rp), dtr = pi/180
  integer :: inode,ielem,k,ivar
  real(kind=8) :: colat,elong,h,bx,by,bz,f
  real(kind=rp) :: theta,phi
  real(kind=rp),dimension(3) :: v,vperp
  real(kind=rp),dimension(:,:,:),allocatable :: b,b_int,var_ext_alt,var_ext_mesh,var_int
  integer,dimension(:,:),allocatable :: cellidx_in,nodeidx_in,nodeidx_out,altidx_out
  real(kind=rp),dimension(:),allocatable :: z_ext,v_ext
  logical,dimension(:),allocatable :: positive

  call init_esmf

  call read_coord('/glade/work/haonan/coord.nc')
  call setup_nodes
  call setup_elements

  call read_data('/glade/work/haonan/ionmodel.nc',timeidx)

  call init_ext(nnode_ext,nelement_ext,nalt_int,nvar_ext, &
    elements_ext,nodeGlon,nodeGlat,elementGlon,elementGlat)
  call init_int(nnode_int,nelement_int,nalt_int,nvar_ext, &
    elements_int,lon_int,lat_int)
  call init_regrid

  allocate(z_ext(nlev))
  allocate(v_ext(nlev))
  allocate(positive(nlev))
  allocate(b(nnode_int,nalt_int,3))
  allocate(b_int(nnode_int,nalt_int,4))
  allocate(var_ext_alt(nelement_ext,nalt_int,nvar_ext))
  allocate(var_ext_mesh(nnode_int,nalt_int,nvar_ext))
  allocate(var_int(nnode_int,nalt_int,nvar_int))
  allocate(cellidx_in(nnode_int,nalt_int))
  allocate(nodeidx_in(nnode_int,nalt_int))
  allocate(nodeidx_out(ncell,nnode_per_cell))
  allocate(altidx_out(ncell,nnode_per_cell))

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
      b(inode,k,:) = rotate_s2c(theta,phi,real((/-bz,-bx,by/),kind=rp))
      b_int(inode,k,1:3) = b(inode,k,:)/f
      b_int(inode,k,4) = f*1e-9_rp
    enddo
  enddo

  do concurrent (ielem = 1:nelement_ext, ivar = 1:nvar_ext)
    z_ext = alt_ext(ielem,:)
    v_ext = var_ext(ielem,:,ivar)
    if (any((/iv1_ext,iv2_ext,iv3_ext,ivp, &
              ie1_ext,ie2_ext,ie3_ext, &
              iu1_ext,iu2_ext,iu3_ext/) == ivar)) then
      var_ext_alt(ielem,:,ivar) = interp1d(alt_int,z_ext,v_ext)
    else
      positive = v_ext > 0
      var_ext_alt(ielem,:,ivar) = exp(interp1d(alt_int, &
        pack(z_ext,positive),log(pack(v_ext,positive))))
    endif
  enddo

  do ielem = 1,nelement_ext
    theta = pi/2-elementGlat(ielem)*dtr
    phi = elementGlon(ielem)*dtr

    do k = 1,nalt_int
      v = rotate_s2c(theta,phi, &
        (/var_ext_alt(ielem,k,iv3_ext), &
          -var_ext_alt(ielem,k,iv2_ext), &
          var_ext_alt(ielem,k,iv1_ext)/))
      var_ext_alt(ielem,k,iv1_ext) = v(1)
      var_ext_alt(ielem,k,iv2_ext) = v(2)
      var_ext_alt(ielem,k,iv3_ext) = v(3)

      v = rotate_s2c(theta,phi, &
        (/var_ext_alt(ielem,k,ie3_ext), &
          -var_ext_alt(ielem,k,ie2_ext), &
          var_ext_alt(ielem,k,ie1_ext)/))
      var_ext_alt(ielem,k,ie1_ext) = v(1)
      var_ext_alt(ielem,k,ie2_ext) = v(2)
      var_ext_alt(ielem,k,ie3_ext) = v(3)

      v = rotate_s2c(theta,phi, &
        (/var_ext_alt(ielem,k,iu3_ext), &
          -var_ext_alt(ielem,k,iu2_ext), &
          var_ext_alt(ielem,k,iu1_ext)/))
      var_ext_alt(ielem,k,iu1_ext) = v(1)
      var_ext_alt(ielem,k,iu2_ext) = v(2)
      var_ext_alt(ielem,k,iu3_ext) = v(3)
    enddo
  enddo

  var_ext_mesh = regrid(nelement_ext,nnode_int,nalt_int,nvar_ext,var_ext_alt)

  do concurrent (inode = 1:nnode_int, k = 1:nalt_int)
    var_int(inode,k,iop_int) = var_ext_mesh(inode,k,iop_ext)*1e6_rp
    var_int(inode,k,io2p_int) = var_ext_mesh(inode,k,io2p_ext)*1e6_rp
    var_int(inode,k,inp_int) = var_ext_mesh(inode,k,inp_ext)*1e6_rp
    var_int(inode,k,in2p_int) = var_ext_mesh(inode,k,in2p_ext)*1e6_rp
    var_int(inode,k,inop_int) = var_ext_mesh(inode,k,inop_ext)*1e6_rp
    var_int(inode,k,iv1_int) = (var_ext_mesh(inode,k,iv1_ext)+ &
      var_ext_mesh(inode,k,ivp)*b_int(inode,k,1))*1e-2_rp
    var_int(inode,k,iv2_int) = (var_ext_mesh(inode,k,iv2_ext)+ &
      var_ext_mesh(inode,k,ivp)*b_int(inode,k,2))*1e-2_rp
    var_int(inode,k,iv3_int) = (var_ext_mesh(inode,k,iv3_ext)+ &
      var_ext_mesh(inode,k,ivp)*b_int(inode,k,3))*1e-2_rp
    var_int(inode,k,iti_int) = var_ext_mesh(inode,k,iti_ext)
    var_int(inode,k,ite_int) = var_ext_mesh(inode,k,ite_ext)
    var_int(inode,k,iprod_op_int) = var_ext_mesh(inode,k,iprod_op_ext)*1e6_rp
    var_int(inode,k,iprod_o2p_int) = var_ext_mesh(inode,k,iprod_o2p_ext)*1e6_rp
    var_int(inode,k,iprod_np_int) = var_ext_mesh(inode,k,iprod_np_ext)*1e6_rp
    var_int(inode,k,iprod_n2p_int) = var_ext_mesh(inode,k,iprod_n2p_ext)*1e6_rp
    var_int(inode,k,iprod_nop_int) = var_ext_mesh(inode,k,iprod_nop_ext)*1e6_rp
    var_int(inode,k,iloss_op_int) = var_ext_mesh(inode,k,iloss_op_ext)
    var_int(inode,k,iloss_o2p_int) = var_ext_mesh(inode,k,iloss_o2p_ext)
    var_int(inode,k,iloss_np_int) = var_ext_mesh(inode,k,iloss_np_ext)
    var_int(inode,k,iloss_n2p_int) = var_ext_mesh(inode,k,iloss_n2p_ext)
    var_int(inode,k,iloss_nop_int) = var_ext_mesh(inode,k,iloss_nop_ext)
    var_int(inode,k,itn_int) = var_ext_mesh(inode,k,itn_ext)
    var_int(inode,k,iu1_int) = var_ext_mesh(inode,k,iu1_ext)*1e-2_rp
    var_int(inode,k,iu2_int) = var_ext_mesh(inode,k,iu2_ext)*1e-2_rp
    var_int(inode,k,iu3_int) = var_ext_mesh(inode,k,iu3_ext)*1e-2_rp

! zero out parallel electric field
    v = (/var_ext_mesh(inode,k,ie1_ext), &
          var_ext_mesh(inode,k,ie2_ext), &
          var_ext_mesh(inode,k,ie3_ext)/)
    vperp = v-dot_product(v,b_int(inode,k,1:3))*b_int(inode,k,1:3)
    var_int(inode,k,ie1_int) = vperp(1)*1e2_rp
    var_int(inode,k,ie2_int) = vperp(2)*1e2_rp
    var_int(inode,k,ie3_int) = vperp(3)*1e2_rp

    call calculate(var_ext_mesh(inode,k,itn_ext), &
      var_ext_mesh(inode,k,iti_ext),var_ext_mesh(inode,k,ite_ext), &
      var_ext_mesh(inode,k,io2),var_ext_mesh(inode,k,io1), &
      var_ext_mesh(inode,k,ihe),var_ext_mesh(inode,k,in2), &
      var_ext_mesh(inode,k,iop_ext),var_ext_mesh(inode,k,io2p_ext), &
      var_ext_mesh(inode,k,inp_ext),var_ext_mesh(inode,k,in2p_ext), &
      var_ext_mesh(inode,k,inop_ext), &
      var_ext_mesh(inode,k,iped),var_ext_mesh(inode,k,iqsum), &
      var_int(inode,k,ie1_int),var_int(inode,k,ie2_int),var_int(inode,k,ie3_int), &
      var_ext_mesh(inode,k,iu1_ext)*1e-2_rp, &
      var_ext_mesh(inode,k,iu2_ext)*1e-2_rp, &
      var_ext_mesh(inode,k,iu3_ext)*1e-2_rp, &
      b(inode,k,1)*1e-9_rp,b(inode,k,2)*1e-9_rp,b(inode,k,3)*1e-9_rp, &
      var_int(inode,k,inu_op),var_int(inode,k,inu_o2p), &
      var_int(inode,k,inu_np),var_int(inode,k,inu_n2p), &
      var_int(inode,k,inu_nop), &
      var_int(inode,k,iheat_i),var_int(inode,k,iheat_e), &
      var_int(inode,k,iqei),var_int(inode,k,iqin),var_int(inode,k,iqen), &
      var_int(inode,k,ikappa_i),var_int(inode,k,ikappa_e))
  enddo

  call reorder_in(cellidx_in,nodeidx_in)
  call reorder_out(nodeidx_out,altidx_out)

  call write_data(nnode_int,nalt_int,ncell,nnode_per_cell, &
    lon_int,lat_int,alt_int,b_int,var_int, &
    cellidx_in,nodeidx_in,nodeidx_out,altidx_out, &
    '/glade/work/haonan/data.nc')

endprogram driver
