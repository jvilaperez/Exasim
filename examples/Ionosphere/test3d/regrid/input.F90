module input_module

  use prec,only:rp

  implicit none

  character(len=*),dimension(*),parameter :: varname = &
    (/'OP      ','O2P     ','NP      ','N2P     ','NOP     ', &
      'UI_ExB  ','VI_ExB  ','WI_ExB  ','PARVEL  ','TI      ','TE      ', &
      'PROD_OP ','PROD_O2P','PROD_NP ','PROD_N2P','PROD_NOP', &
      'LOSS_OP ','LOSS_O2P','LOSS_NP ','LOSS_N2P','LOSS_NOP', &
      'PED     ','QSUM    ','EX      ','EY      ','EZ      ', &
      'TN      ','UN      ','VN      ','WN      ', &
      'O2_CM3  ','O1_CM3  ','HE_CM3  ','N2_CM3  '/)
  integer,parameter :: nvar = size(varname), &
    iop = 1, io2p = 2, inp = 3, in2p = 4, inop = 5, &
    iv1 = 6, iv2 = 7, iv3 = 8, ivp = 9, iti = 10, ite = 11, &
    iprod_op = 12, iprod_o2p = 13, &
    iprod_np = 14, iprod_n2p = 15, iprod_nop = 16, &
    iloss_op = 17, iloss_o2p = 18, &
    iloss_np = 19, iloss_n2p = 20, iloss_nop = 21, &
    iped = 22, iqsum = 23, ie1 = 24, ie2 = 25, ie3 = 26, &
    itn = 27, iu1 = 28, iu2 = 29, iu3 = 30, &
    io2 = 31, io1 = 32, ihe = 33, in2 = 34

  integer :: nelement,nnode,nlev
  integer,dimension(:,:),allocatable :: elements
  real(kind=rp),dimension(:),allocatable :: &
    nodeGlon,nodeGlat,elementGlon,elementGlat
  real(kind=rp),dimension(:,:),allocatable :: alt
  real(kind=rp),dimension(:,:,:),allocatable :: variable

  contains
!-----------------------------------------------------------------------
  subroutine read_data(filename,timeidx)

    use netcdf,only:nf90_open,nf90_inq_dimid,nf90_inquire_dimension, &
      nf90_inq_varid,nf90_get_var,nf90_close,nf90_noerr,nf90_nowrite

    character(len=*),intent(in) :: filename
    integer,intent(in) :: timeidx

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

    stat = nf90_inq_dimid(ncid,'lev',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nlev)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    nlev = nlev-1 ! exclude the filling level for now

    allocate(values1d(max(nnode,nelement,nlev)))
    allocate(values3d(nelement,nlev,1))
    allocate(elements(3,nelement))
    allocate(nodeGlon(nnode))
    allocate(nodeGlat(nnode))
    allocate(elementGlon(nelement))
    allocate(elementGlat(nelement))
    allocate(alt(nelement,nlev))
    allocate(variable(nelement,nlev,nvar))

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

    stat = nf90_inq_varid(ncid,'ZG',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values3d, &
      start=(/1,1,timeidx/),count=(/nelement,nlev,1/))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    alt = values3d(:,:,1)/100

    do ivar = 1,nvar
      stat = nf90_inq_varid(ncid,trim(varname(ivar)),varid)
      if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

      stat = nf90_get_var(ncid,varid,values3d, &
        start=(/1,1,timeidx/),count=(/nelement,nlev,1/))
      if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

      variable(:,:,ivar) = values3d(:,:,1)
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
