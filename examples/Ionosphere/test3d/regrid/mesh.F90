module mesh_module

  use prec,only:rp

  implicit none

  real(kind=rp),parameter :: rtol = 0, atol = 1e-6_rp, rtd = 45/atan(1.0_rp)

  integer :: ncell,nnode_per_cell,nsubcell,nnode_per_subcell,nalt,nnode,nelement
  real(kind=rp) :: re
  real(kind=rp),dimension(:),allocatable :: r,alt,lat,lon,latr,lonr
  real(kind=rp),dimension(:,:),allocatable :: x,y,z,nodes
  integer,dimension(:,:),allocatable :: conn,elements

  contains
!-----------------------------------------------------------------------
  subroutine read_coord(filename)

    use netcdf,only:nf90_open,nf90_inq_dimid,nf90_inquire_dimension, &
      nf90_inq_varid,nf90_get_var,nf90_close,nf90_noerr,nf90_nowrite

    character(len=*),intent(in) :: filename

    integer :: stat,ncid,dimid,varid
    real(kind=8) :: r8
    real(kind=8),dimension(:,:),allocatable :: values

    stat = nf90_open(trim(filename),nf90_nowrite,ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_open',stat)

    stat = nf90_inq_dimid(ncid,'cell',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=ncell)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    stat = nf90_inq_dimid(ncid,'node_per_cell',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nnode_per_cell)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    stat = nf90_inq_dimid(ncid,'subcell',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nsubcell)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    stat = nf90_inq_dimid(ncid,'node_per_subcell',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nnode_per_subcell)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    allocate(values(ncell,nnode_per_cell))
    allocate(x(ncell,nnode_per_cell))
    allocate(y(ncell,nnode_per_cell))
    allocate(z(ncell,nnode_per_cell))
    allocate(conn(nsubcell,nnode_per_subcell))

    stat = nf90_inq_varid(ncid,'r',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,r8)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    re = r8

    stat = nf90_inq_varid(ncid,'x',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    x = values

    stat = nf90_inq_varid(ncid,'y',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    y = values

    stat = nf90_inq_varid(ncid,'z',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    z = values

    stat = nf90_inq_varid(ncid,'conn',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,conn)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_close',stat)

    deallocate(values)

  endsubroutine read_coord
!-----------------------------------------------------------------------
  subroutine setup_nodes

    use util_module,only:isclose,unique

    integer :: nrep
    logical,dimension(ncell,nnode_per_cell) :: layer
    real(kind=rp),dimension(ncell,nnode_per_cell) :: rad
    real(kind=rp),dimension(ncell*nnode_per_cell) :: r_tmp
    real(kind=rp),dimension(:),allocatable :: coord_z
    real(kind=rp),dimension(:,:),allocatable :: nodes_repeat,nodes_tmp

    nrep = ncell*nnode_per_cell
    rad = sqrt(x**2+y**2+z**2)
    call unique(reshape(rad,(/nrep/)),nalt,r_tmp,rtol,atol)
    allocate(r(nalt))
    allocate(alt(nalt))
    r = r_tmp(1:nalt)
    alt = r-re

    layer = isclose(rad,r(1),rtol,atol)
    nrep = count(layer)
    allocate(nodes_repeat(3,nrep))
    allocate(nodes_tmp(3,nrep))
    nodes_repeat(1,:) = pack(x,layer)
    nodes_repeat(2,:) = pack(y,layer)
    nodes_repeat(3,:) = pack(z,layer)
    call unique(nodes_repeat,nnode,nodes_tmp,rtol,atol)
    allocate(nodes(3,nnode))
    nodes = nodes_tmp(:,1:nnode)

    allocate(coord_z(nnode))
    allocate(latr(nnode))
    allocate(lonr(nnode))
    allocate(lat(nnode))
    allocate(lon(nnode))
    coord_z = nodes(3,:)/r(1)
    where (abs(coord_z) > 1) coord_z = sign(1.0_rp,coord_z)
    latr = asin(coord_z)
    lonr = atan2(nodes(2,:),nodes(1,:))
    lat = latr*rtd
    lon = lonr*rtd

    deallocate(coord_z)
    deallocate(nodes_repeat)
    deallocate(nodes_tmp)

  endsubroutine setup_nodes
!-----------------------------------------------------------------------
  subroutine setup_elements

    use util_module,only:isclose,find_index

    integer,parameter :: maxelem = 10000
    integer :: icell,ilayer,i,isubcell,inode
    integer,dimension(9) :: node_per_cell
    integer,dimension(nnode_per_subcell,maxelem) :: elements_tmp

    nelement = 0
    do icell = 1,ncell
      do ilayer = 1,3
        if (isclose(r(1),norm2( &
            (/x(icell,(ilayer-1)*9+1), &
              y(icell,(ilayer-1)*9+1), &
              z(icell,(ilayer-1)*9+1)/)),rtol,atol)) then
          do i = 1,9
            node_per_cell(i) = find_index(nodes, &
              (/x(icell,(ilayer-1)*9+i), &
                y(icell,(ilayer-1)*9+i), &
                z(icell,(ilayer-1)*9+i)/),rtol,atol)
          enddo
          do isubcell = 1,nsubcell
            do inode = 1,nnode_per_subcell
              elements_tmp(inode,nelement+isubcell) = node_per_cell(conn(inode,isubcell))
            enddo
          enddo
          nelement = nelement+nsubcell
        endif
      enddo
    enddo

    allocate(elements(nnode_per_subcell,nelement))
    elements = elements_tmp(:,1:nelement)

  endsubroutine setup_elements
!-----------------------------------------------------------------------
  pure subroutine reorder_in(cellidx,nodeidx)
! reorder the CG mesh to the corresponding DG mesh (input)

    use util_module,only:isclose

    integer,dimension(nnode,nalt),intent(out) :: cellidx,nodeidx

    integer :: inode,k,icell,i
    real(kind=rp),dimension(ncell,nnode_per_cell) :: rad
    real(kind=rp),dimension(3,ncell,nnode_per_cell) :: coord_norm
    real(kind=rp),dimension(3,nnode) :: nodes_norm

    rad = sqrt(x**2+y**2+z**2)
    coord_norm(1,:,:) = x/rad
    coord_norm(2,:,:) = y/rad
    coord_norm(3,:,:) = z/rad
    nodes_norm = nodes/r(1)

    do concurrent (inode = 1:nnode, k = 1:nalt)
      loop: do icell = 1,ncell
        do i = 1,nnode_per_cell
          if (isclose(rad(icell,i),r(k),rtol,atol) .and. &
              all(isclose(coord_norm(:,icell,i),nodes_norm(:,inode),rtol,atol))) then
            cellidx(inode,k) = icell
            nodeidx(inode,k) = i
            exit loop
          endif
        enddo
      enddo loop
    enddo

  endsubroutine reorder_in
!-----------------------------------------------------------------------
  pure subroutine reorder_out(nodeidx,altidx)
! reorder the DG mesh to the corresponding CG mesh (output)

    use util_module,only:find_index

    integer,dimension(ncell,nnode_per_cell),intent(out) :: nodeidx,altidx

    integer :: icell,i
    real(kind=rp),dimension(ncell,nnode_per_cell) :: rad
    real(kind=rp),dimension(3,ncell,nnode_per_cell) :: coord_norm
    real(kind=rp),dimension(3,nnode) :: nodes_norm

    rad = sqrt(x**2+y**2+z**2)
    coord_norm(1,:,:) = x/rad
    coord_norm(2,:,:) = y/rad
    coord_norm(3,:,:) = z/rad
    nodes_norm = nodes/r(1)

    do concurrent (icell = 1:ncell, i = 1:nnode_per_cell)
      nodeidx(icell,i) = find_index(nodes_norm,coord_norm(:,icell,i),rtol,atol)
      altidx(icell,i) = find_index(r,rad(icell,i),rtol,atol)
    enddo

  endsubroutine reorder_out
!-----------------------------------------------------------------------
  subroutine handle_error(funcname,ncerr)

    use netcdf,only:nf90_strerror

    character(len=*),intent(in) :: funcname
    integer,intent(in) :: ncerr

    write(6,"('NetCDF error encountered: ',a,', when calling ',a)") &
      trim(nf90_strerror(ncerr)),trim(funcname)

  endsubroutine handle_error
!-----------------------------------------------------------------------
endmodule mesh_module
