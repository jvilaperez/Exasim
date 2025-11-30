module output_module

  implicit none

  contains
!-----------------------------------------------------------------------
  subroutine write_data(nnode,nalt,ntime,nvar,ncell,nnode_per_cell, &
    lon,lat,alt,time,b,variable,cellidx_in,nodeidx_in,nodeidx_out,altidx_out, &
    filename,varname)

    use prec,only:rp
    use netcdf,only:nf90_create,nf90_def_dim,nf90_def_var,nf90_enddef, &
      nf90_put_var,nf90_close,nf90_netcdf4,nf90_double,nf90_int,nf90_noerr

    integer,intent(in) :: nnode,nalt,ntime,nvar,ncell,nnode_per_cell
    real(kind=rp),dimension(nnode),intent(in) :: lon,lat
    real(kind=rp),dimension(nalt),intent(in) :: alt
    real(kind=rp),dimension(ntime),intent(in) :: time
    real(kind=rp),dimension(nnode,nalt,3),intent(in) :: b
    real(kind=rp),dimension(nnode,nalt,ntime,nvar),intent(in) :: variable
    integer,dimension(nnode,nalt),intent(in) :: cellidx_in,nodeidx_in
    integer,dimension(ncell,nnode_per_cell),intent(in) :: nodeidx_out,altidx_out
    character(len=*),intent(in) :: filename
    character(len=*),dimension(nvar),intent(in) :: varname

    integer :: stat,ncid,dimid_node,dimid_alt,dimid_time, &
      dimid_cell,dimid_node_per_cell, &
      varid_lon,varid_lat,varid_alt,varid_time,ivar, &
      varid_cellidx_in,varid_nodeidx_in, &
      varid_nodeidx_out,varid_altidx_out
    integer,dimension(3) :: varid_b
    integer,dimension(nvar) :: varid

    stat = nf90_create(trim(filename),nf90_netcdf4,ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_create',stat)

    stat = nf90_def_dim(ncid,'node',nnode,dimid_node)
    if (stat /= nf90_noerr) call handle_error('nf90_def_dim',stat)

    stat = nf90_def_dim(ncid,'alt',nalt,dimid_alt)
    if (stat /= nf90_noerr) call handle_error('nf90_def_dim',stat)

    stat = nf90_def_dim(ncid,'time',ntime,dimid_time)
    if (stat /= nf90_noerr) call handle_error('nf90_def_dim',stat)

    stat = nf90_def_dim(ncid,'cell',ncell,dimid_cell)
    if (stat /= nf90_noerr) call handle_error('nf90_def_dim',stat)

    stat = nf90_def_dim(ncid,'node_per_cell',nnode_per_cell,dimid_node_per_cell)
    if (stat /= nf90_noerr) call handle_error('nf90_def_dim',stat)

    stat = nf90_def_var(ncid,'lon',nf90_double,dimid_node,varid_lon)
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    stat = nf90_def_var(ncid,'lat',nf90_double,dimid_node,varid_lat)
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    stat = nf90_def_var(ncid,'alt',nf90_double,dimid_alt,varid_alt)
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    stat = nf90_def_var(ncid,'time',nf90_double,dimid_time,varid_time)
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    stat = nf90_def_var(ncid,'BX',nf90_double, &
      (/dimid_node,dimid_alt/),varid_b(1))
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    stat = nf90_def_var(ncid,'BY',nf90_double, &
      (/dimid_node,dimid_alt/),varid_b(2))
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    stat = nf90_def_var(ncid,'BZ',nf90_double, &
      (/dimid_node,dimid_alt/),varid_b(3))
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    do ivar = 1,nvar
      stat = nf90_def_var(ncid,trim(varname(ivar)),nf90_double, &
        (/dimid_node,dimid_alt,dimid_time/),varid(ivar))
      if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)
    enddo

    stat = nf90_def_var(ncid,'cellidx_in',nf90_int, &
      (/dimid_node,dimid_alt/),varid_cellidx_in)
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    stat = nf90_def_var(ncid,'nodeidx_in',nf90_int, &
      (/dimid_node,dimid_alt/),varid_nodeidx_in)
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    stat = nf90_def_var(ncid,'nodeidx_out',nf90_int, &
      (/dimid_cell,dimid_node_per_cell/),varid_nodeidx_out)
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    stat = nf90_def_var(ncid,'altidx_out',nf90_int, &
      (/dimid_cell,dimid_node_per_cell/),varid_altidx_out)
    if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)

    stat = nf90_enddef(ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_enddef',stat)

    stat = nf90_put_var(ncid,varid_lon,real(lon,kind=8))
    if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)

    stat = nf90_put_var(ncid,varid_lat,real(lat,kind=8))
    if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)

    stat = nf90_put_var(ncid,varid_alt,real(alt,kind=8))
    if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)

    stat = nf90_put_var(ncid,varid_time,real(time,kind=8))
    if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)

    do ivar = 1,3
      stat = nf90_put_var(ncid,varid_b(ivar),real(b(:,:,ivar),kind=8))
      if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)
    enddo

    do ivar = 1,nvar
      stat = nf90_put_var(ncid,varid(ivar),real(variable(:,:,:,ivar),kind=8))
      if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)
    enddo

    stat = nf90_put_var(ncid,varid_cellidx_in,cellidx_in)
    if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)

    stat = nf90_put_var(ncid,varid_nodeidx_in,nodeidx_in)
    if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)

    stat = nf90_put_var(ncid,varid_nodeidx_out,nodeidx_out)
    if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)

    stat = nf90_put_var(ncid,varid_altidx_out,altidx_out)
    if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_close',stat)

  endsubroutine write_data
!-----------------------------------------------------------------------
  subroutine handle_error(funcname,ncerr)

    use netcdf,only:nf90_strerror

    character(len=*),intent(in) :: funcname
    integer,intent(in) :: ncerr

    write(6,"('NetCDF error encountered: ',a,', when calling ',a)") &
      trim(nf90_strerror(ncerr)),trim(funcname)

  endsubroutine handle_error
!-----------------------------------------------------------------------
endmodule output_module
