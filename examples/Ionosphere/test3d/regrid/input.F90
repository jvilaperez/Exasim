module input_module

  use prec,only:rp

  implicit none

  character(len=*),dimension(*),parameter :: varname = &
    (/'OP  ','TI  ','NE  ','NU  ','PROD','LOSS','HEAT','QEI ','QIN ', &
      'UI  ','VI  ','WI  ','EX  ','EY  ','EZ  ','UN_A','VN_A','WN_A'/)
  integer,parameter :: nvar = size(varname)

  integer :: nelement,nnode,nalt
  integer,dimension(:,:),allocatable :: elements
  real(kind=rp),dimension(:),allocatable :: alt, &
    nodeGlon,nodeGlat,elementGlon,elementGlat
  real(kind=rp),dimension(:,:,:),allocatable :: variable

  contains
!-----------------------------------------------------------------------
  subroutine read_data(filename,timeidx)

    use netcdf,only:nf90_open,nf90_inq_dimid,nf90_inquire_dimension, &
      nf90_inq_varid,nf90_get_var,nf90_close,nf90_noerr,nf90_nowrite

    character(len=*),intent(in) :: filename
    integer,intent(in) :: timeidx

    real(kind=rp),dimension(nvar),parameter :: varscale = &
      (/1e6_rp,1.0_rp,1e6_rp,1.0_rp,1e6_rp,1.0_rp,1e7_rp,1e7_rp,1e7_rp, &
        1e-2_rp,1e-2_rp,1e-2_rp,1.0_rp,1.0_rp,1.0_rp,1e-2_rp,1e-2_rp,1e-2_rp/)
    integer :: stat,ncid,dimid,dimlen,varid,ivar
    real(kind=8),dimension(:),allocatable :: values1d
    real(kind=4),dimension(:,:,:),allocatable :: values3d

    stat = nf90_open(trim(filename),nf90_nowrite,ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_open',stat)

    stat = nf90_inq_dimid(ncid,'node_per',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=dimlen)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    if (dimlen /= 3) then
      write(6,"('Input file does not have a consistent element type with the model')")
      stop 'Invalid input file'
    endif

    stat = nf90_inq_dimid(ncid,'element',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nelement)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    stat = nf90_inq_dimid(ncid,'node',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nnode)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    stat = nf90_inq_dimid(ncid,'alt',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nalt)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    allocate(values1d(max(nnode,nelement,nalt)))
    allocate(values3d(nelement,nalt,1))
    allocate(elements(3,nelement))
    allocate(nodeGlon(nnode))
    allocate(nodeGlat(nnode))
    allocate(elementGlon(nelement))
    allocate(elementGlat(nelement))
    allocate(alt(nalt))
    allocate(variable(nelement,nalt,nvar))

    stat = nf90_inq_varid(ncid,'elementConn',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,elements)
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    stat = nf90_inq_varid(ncid,'nodeGlon',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values1d(1:nnode))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    nodeGlon = values1d(1:nnode)

    stat = nf90_inq_varid(ncid,'nodeGlat',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values1d(1:nnode))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    nodeGlat = values1d(1:nnode)

    stat = nf90_inq_varid(ncid,'lon',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values1d(1:nelement))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    elementGlon = values1d(1:nelement)

    stat = nf90_inq_varid(ncid,'lat',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values1d(1:nelement))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    elementGlat = values1d(1:nelement)

    stat = nf90_inq_varid(ncid,'alt',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values1d(1:nalt))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    alt = values1d(1:nalt)/100

    do ivar = 1,nvar
      stat = nf90_inq_varid(ncid,trim(varname(ivar)),varid)
      if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

      stat = nf90_get_var(ncid,varid,values3d, &
        start=(/1,1,timeidx/),count=(/nelement,nalt,1/))
      if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

      variable(:,:,ivar) = values3d(:,:,1)*varscale(ivar)
    enddo

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_close',stat)

    deallocate(values1d)
    deallocate(values3d)

  endsubroutine read_data
!-----------------------------------------------------------------------
  subroutine handle_error(funcname,ncerr)

    use netcdf,only:nf90_strerror

    character(len=*),intent(in) :: funcname
    integer,intent(in) :: ncerr

    write(6,"('NetCDF error encountered: ',a,', when calling ',a)") &
      trim(nf90_strerror(ncerr)),trim(funcname)

    stop

  endsubroutine handle_error
!-----------------------------------------------------------------------
endmodule input_module
