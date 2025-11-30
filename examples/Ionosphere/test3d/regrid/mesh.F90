module mesh_module

  use prec,only:rp

  implicit none

  real(kind=rp),parameter :: rtol = 0, atol = 1e-6_rp, rtd = 45/atan(1.0_rp)

  integer :: ncell,nnode_per_cell,nalt,nnode,nelement
  real(kind=rp),dimension(:),allocatable :: alt,lat,lon,latr,lonr
  real(kind=rp),dimension(:,:),allocatable :: x,y,z,nodes
  integer,dimension(:,:),allocatable :: elements

  contains
!-----------------------------------------------------------------------
  subroutine read_coord(filename)

    use netcdf,only:nf90_open,nf90_inq_dimid,nf90_inquire_dimension, &
      nf90_inq_varid,nf90_get_var,nf90_close,nf90_noerr,nf90_nowrite

    character(len=*),intent(in) :: filename

    integer :: stat,ncid,dimid,varid
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

    allocate(values(ncell,nnode_per_cell))
    allocate(x(ncell,nnode_per_cell))
    allocate(y(ncell,nnode_per_cell))
    allocate(z(ncell,nnode_per_cell))

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

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_close',stat)

    deallocate(values)

  endsubroutine read_coord
!-----------------------------------------------------------------------
  subroutine setup_nodes

    use util_module,only:isclose,unique

    integer :: nrep
    logical,dimension(ncell,nnode_per_cell) :: layer
    real(kind=rp),dimension(ncell,nnode_per_cell) :: r
    real(kind=rp),dimension(ncell*nnode_per_cell) :: alt_tmp
    real(kind=rp),dimension(:),allocatable :: coord_z
    real(kind=rp),dimension(:,:),allocatable :: nodes_repeat,nodes_tmp

    nrep = ncell*nnode_per_cell
    r = sqrt(x**2+y**2+z**2)
    call unique(reshape(r,(/nrep/)),nalt,alt_tmp,rtol,atol)
    allocate(alt(nalt))
    alt = alt_tmp(1:nalt)

    layer = isclose(r,alt(1),rtol,atol)
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
    coord_z = nodes(3,:)/alt(1)
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
    integer :: icell,ilayer,i
    integer,dimension(9) :: node_per_cell
    integer,dimension(4,maxelem) :: elements_tmp

    nelement = 0
    do icell = 1,ncell
      do ilayer = 1,3
        if (isclose(alt(1),norm2( &
            (/x(icell,(ilayer-1)*9+1), &
              y(icell,(ilayer-1)*9+1), &
              z(icell,(ilayer-1)*9+1)/)),rtol,atol)) then
          do i = 1,9
            node_per_cell(i) = find_index(nodes, &
              (/x(icell,(ilayer-1)*9+i), &
                y(icell,(ilayer-1)*9+i), &
                z(icell,(ilayer-1)*9+i)/),rtol,atol)
          enddo
          elements_tmp(:,nelement+1) = (/node_per_cell(1),node_per_cell(2),node_per_cell(5),node_per_cell(4)/)
          elements_tmp(:,nelement+2) = (/node_per_cell(2),node_per_cell(3),node_per_cell(6),node_per_cell(5)/)
          elements_tmp(:,nelement+3) = (/node_per_cell(4),node_per_cell(5),node_per_cell(8),node_per_cell(7)/)
          elements_tmp(:,nelement+4) = (/node_per_cell(5),node_per_cell(6),node_per_cell(9),node_per_cell(8)/)
          nelement = nelement+4
        endif
      enddo
    enddo

    allocate(elements(4,nelement))
    elements = elements_tmp(:,1:nelement)

  endsubroutine setup_elements
!-----------------------------------------------------------------------
  pure subroutine reorder_in(cellidx,nodeidx)
! reorder the CG mesh to the corresponding DG mesh (input)

    use util_module,only:isclose

    integer,dimension(nnode,nalt),intent(out) :: cellidx,nodeidx

    integer :: inode,k,icell,i
    real(kind=rp),dimension(ncell,nnode_per_cell) :: r
    real(kind=rp),dimension(3,ncell,nnode_per_cell) :: coord_norm
    real(kind=rp),dimension(3,nnode) :: nodes_norm

    r = sqrt(x**2+y**2+z**2)
    coord_norm(1,:,:) = x/r
    coord_norm(2,:,:) = y/r
    coord_norm(3,:,:) = z/r
    nodes_norm = nodes/alt(1)

    do concurrent (inode = 1:nnode, k = 1:nalt)
      loop: do icell = 1,ncell
        do i = 1,nnode_per_cell
          if (isclose(r(icell,i),alt(k),rtol,atol) .and. &
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
    real(kind=rp),dimension(ncell,nnode_per_cell) :: r
    real(kind=rp),dimension(3,ncell,nnode_per_cell) :: coord_norm
    real(kind=rp),dimension(3,nnode) :: nodes_norm

    r = sqrt(x**2+y**2+z**2)
    coord_norm(1,:,:) = x/r
    coord_norm(2,:,:) = y/r
    coord_norm(3,:,:) = z/r
    nodes_norm = nodes/alt(1)

    do concurrent (icell = 1:ncell, i = 1:nnode_per_cell)
      nodeidx(icell,i) = find_index(nodes_norm,coord_norm(:,icell,i),rtol,atol)
      altidx(icell,i) = find_index(alt,r(icell,i),rtol,atol)
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
