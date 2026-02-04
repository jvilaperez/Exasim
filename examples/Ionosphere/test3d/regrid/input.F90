module input_module

  use prec,only:rp

  implicit none

  character(len=*),dimension(*),parameter :: varname = &
    (/'OP          ','TI          ','NU_OP       ', &
      'OP_PROD     ','OP_LOSS_COEF', &
      'TI_HEAT     ','TI_COOL_COEF', &
      'UI_ExB      ','VI_ExB      ','WI_ExB      ', &
      'EEX         ','EEY         ','EEZ         ', &
      'UN          ','VN          ','WN          '/)
  integer,parameter :: nvar = size(varname)

  integer :: nlon,nlat,nlev,ntime
  real(kind=rp),dimension(:),allocatable :: lon,lat,lev,time
  real(kind=rp),dimension(:,:,:,:),allocatable :: z
  real(kind=rp),dimension(:,:,:,:,:),allocatable :: variable

  contains
!-----------------------------------------------------------------------
  subroutine read_data(filename)

    use util_module,only:rotate_s2c
    use netcdf,only:nf90_open,nf90_inq_dimid,nf90_inquire_dimension, &
      nf90_inq_varid,nf90_get_var,nf90_close,nf90_noerr,nf90_nowrite

    character(len=*),intent(in) :: filename

    real(kind=rp),dimension(nvar),parameter :: varscale = &
      (/1e6_rp,1.0_rp,1.0_rp, &
        1e6_rp,1.0_rp, &
        0.1_rp,1.0_rp, &
        1e-2_rp,1e-2_rp,1e-2_rp, &
        1e2_rp,1e2_rp,1e2_rp, &
        1e-2_rp,1e-2_rp,1e-2_rp/)
    integer :: stat,ncid,dimid,varid,ivar
    real(kind=8),dimension(:),allocatable :: values1d
    real(kind=4),dimension(:,:,:,:),allocatable :: values4d

    stat = nf90_open(trim(filename),nf90_nowrite,ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_open',stat)

    stat = nf90_inq_dimid(ncid,'lon',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nlon)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    stat = nf90_inq_dimid(ncid,'lat',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nlat)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    stat = nf90_inq_dimid(ncid,'lev',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nlev)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

! exlude the top level (filled)
    nlev = nlev-1

    stat = nf90_inq_dimid(ncid,'time',dimid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=ntime)
    if (stat /= nf90_noerr) call handle_error('nf90_inquire_dimension',stat)

    allocate(values1d(max(nlon,nlat,nlev,ntime)))
    allocate(values4d(nlon,nlat,nlev,ntime))
    allocate(lon(nlon))
    allocate(lat(nlat))
    allocate(lev(nlev))
    allocate(time(ntime))
    allocate(z(nlon,nlat,nlev,ntime))
    allocate(variable(nlon,nlat,nlev,ntime,nvar))

    stat = nf90_inq_varid(ncid,'lon',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values1d(1:nlon))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    lon = values1d(1:nlon)

    stat = nf90_inq_varid(ncid,'lat',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values1d(1:nlat))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    lat = values1d(1:nlat)

    stat = nf90_inq_varid(ncid,'lev',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values1d(1:nlev),count=(/nlev/))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    lev = values1d(1:nlev)

    stat = nf90_inq_varid(ncid,'time',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values1d(1:ntime))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

    time = values1d(1:ntime)

    stat = nf90_inq_varid(ncid,'Z',varid)
    if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,values4d,count=(/nlon,nlat,nlev,ntime/))
    if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

! fields can be defined at midpoints and interfaces
! for now, just discard those minor differences
    z = values4d*1e-5_rp

    do ivar = 1,nvar
      stat = nf90_inq_varid(ncid,trim(varname(ivar)),varid)
      if (stat /= nf90_noerr) call handle_error('nf90_inq_varid',stat)

      stat = nf90_get_var(ncid,varid,values4d,count=(/nlon,nlat,nlev,ntime/))
      if (stat /= nf90_noerr) call handle_error('nf90_get_var',stat)

      variable(:,:,:,:,ivar) = values4d*varscale(ivar)
    enddo

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) call handle_error('nf90_close',stat)

    deallocate(values1d)
    deallocate(values4d)

  endsubroutine read_data
!-----------------------------------------------------------------------
  subroutine handle_error(funcname,ncerr)

    use netcdf,only:nf90_strerror

    character(len=*),intent(in) :: funcname
    integer,intent(in) :: ncerr

    write(6,"('NetCDF error encountered: ',a,', when calling ',a)") &
      trim(nf90_strerror(ncerr)),trim(funcname)

  endsubroutine handle_error
!-----------------------------------------------------------------------
endmodule input_module
