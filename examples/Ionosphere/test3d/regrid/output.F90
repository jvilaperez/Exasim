module output_module

  use prec,only:rp

  implicit none

  character(len=*),dimension(*),parameter :: varname = &
    (/'OP      ','O2P     ','NP      ','N2P     ','NOP     ', &
      'V1      ','V2      ','V3      ','TI      ','TE      ', &
      'PROD_OP ','PROD_O2P','PROD_NP ','PROD_N2P','PROD_NOP', &
      'LOSS_OP ','LOSS_O2P','LOSS_NP ','LOSS_N2P','LOSS_NOP', &
      'NU_OP   ','NU_O2P  ','NU_NP   ','NU_N2P  ','NU_NOP  ', &
      'E1      ','E2      ','E3      ','HEAT_I  ','HEAT_E  ', &
      'QEI     ','QIN     ','QEN     ','KAPPA_I ','KAPPA_E ', &
      'TN      ','U1      ','U2      ','U3      '/)
  integer,parameter :: nvar = size(varname), &
    iop = 1, io2p = 2, inp = 3, in2p = 4, inop = 5, &
    iv1 = 6, iv2 = 7, iv3 = 8, iti = 9, ite = 10, &
    iprod_op = 11, iprod_o2p = 12, iprod_np = 13, iprod_n2p = 14, iprod_nop = 15, &
    iloss_op = 16, iloss_o2p = 17, iloss_np = 18, iloss_n2p = 19, iloss_nop = 20, &
    inu_op = 21, inu_o2p = 22, inu_np = 23, inu_n2p = 24, inu_nop = 25, &
    ie1 = 26, ie2 = 27, ie3 = 28, iheat_i = 29, iheat_e = 30, &
    iqei = 31, iqin = 32, iqen = 33, ikappa_i = 34, ikappa_e = 35, &
    itn = 36, iu1 = 37, iu2 = 38, iu3 = 39

  contains
!-----------------------------------------------------------------------
  subroutine write_data(nnode,nalt,ncell,nnode_per_cell, &
    lon,lat,alt,b,variable,cellidx_in,nodeidx_in,nodeidx_out,altidx_out, &
    filename)

    use netcdf,only:nf90_create,nf90_def_dim,nf90_def_var,nf90_enddef, &
      nf90_put_var,nf90_close,nf90_netcdf4,nf90_double,nf90_int,nf90_noerr

    integer,intent(in) :: nnode,nalt,ncell,nnode_per_cell
    real(kind=rp),dimension(nnode),intent(in) :: lon,lat
    real(kind=rp),dimension(nalt),intent(in) :: alt
    real(kind=rp),dimension(nnode,nalt,4),intent(in) :: b
    real(kind=rp),dimension(nnode,nalt,nvar),intent(in) :: variable
    integer,dimension(nnode,nalt),intent(in) :: cellidx_in,nodeidx_in
    integer,dimension(ncell,nnode_per_cell),intent(in) :: nodeidx_out,altidx_out
    character(len=*),intent(in) :: filename

    integer :: stat,ncid,dimid_node,dimid_alt, &
      dimid_cell,dimid_node_per_cell, &
      varid_lon,varid_lat,varid_alt,ivar, &
      varid_cellidx_in,varid_nodeidx_in, &
      varid_nodeidx_out,varid_altidx_out
    integer,dimension(4) :: varid_b
    integer,dimension(nvar) :: varid
    character(len=*),dimension(4),parameter :: bname = (/'b1','b2','b3','B '/)

    stat = nf90_create(trim(filename),nf90_netcdf4,ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_create',stat)

    stat = nf90_def_dim(ncid,'node',nnode,dimid_node)
    if (stat /= nf90_noerr) call handle_error('nf90_def_dim',stat)

    stat = nf90_def_dim(ncid,'alt',nalt,dimid_alt)
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

    do ivar = 1,4
      stat = nf90_def_var(ncid,trim(bname(ivar)),nf90_double, &
        (/dimid_node,dimid_alt/),varid_b(ivar))
      if (stat /= nf90_noerr) call handle_error('nf90_def_var',stat)
    enddo

    do ivar = 1,nvar
      stat = nf90_def_var(ncid,trim(varname(ivar)), &
        nf90_double,(/dimid_node,dimid_alt/),varid(ivar))
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

    do ivar = 1,4
      stat = nf90_put_var(ncid,varid_b(ivar),real(b(:,:,ivar),kind=8))
      if (stat /= nf90_noerr) call handle_error('nf90_put_var',stat)
    enddo

    do ivar = 1,nvar
      stat = nf90_put_var(ncid,varid(ivar),real(variable(:,:,ivar),kind=8))
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
