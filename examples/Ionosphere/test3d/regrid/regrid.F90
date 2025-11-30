module regrid_module

  use prec,only:rp
  use ESMF,only:ESMF_Grid,ESMF_Mesh,ESMF_Field,ESMF_RouteHandle

  implicit none

  logical,parameter :: debug = .false.

  type(ESMF_Grid) :: extGrid
  type(ESMF_Mesh) :: intMesh
  type(ESMF_Field) :: extField,intField
  type(ESMF_RouteHandle) :: rh

  contains
!-----------------------------------------------------------------------
  subroutine init_esmf

    use ESMF,only:ESMF_Initialize,ESMF_SUCCESS

    integer :: rc

    call ESMF_Initialize(rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_Initialize')

  endsubroutine init_esmf
!-----------------------------------------------------------------------
  subroutine init_ext(nlon,nlat,nalt,ntime,nfield,lon,lat)

    use ESMF,only:ESMF_GridCreate1PeriDim,ESMF_GridAddCoord, &
      ESMF_GridGetCoord,ESMF_GridValidate, &
      ESMF_FieldCreate,ESMF_FieldValidate, &
      ESMF_KIND_R8,ESMF_INDEX_GLOBAL, &
      ESMF_TYPEKIND_R8,ESMF_SUCCESS

    integer,intent(in) :: nlon,nlat,nalt,ntime,nfield
    real(kind=rp),dimension(nlon),intent(in) :: lon
    real(kind=rp),dimension(nlat),intent(in) :: lat

    integer :: rc
    real(kind=ESMF_KIND_R8),dimension(:),pointer :: farrayPtr

    extGrid = ESMF_GridCreate1PeriDim( &
      countsPerDEDim1=(/nlon/),coordDep1=(/1/), &
      countsPerDEDim2=(/nlat/),coordDep2=(/2/), &
      indexflag=ESMF_INDEX_GLOBAL,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_GridCreate1PeriDim')

    call ESMF_GridAddCoord(grid=extGrid,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_GridAddCoord')

    call ESMF_GridGetCoord(grid=extGrid, &
      coordDim=1,farrayPtr=farrayPtr,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_GridGetCoord')

    farrayPtr = lon

    call ESMF_GridGetCoord(grid=extGrid, &
      coordDim=2,farrayPtr=farrayPtr,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_GridGetCoord')

    farrayPtr = lat

    call ESMF_GridValidate(grid=extGrid,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_GridValidate')

    extField = ESMF_FieldCreate(grid=extGrid, &
      typekind=ESMF_TYPEKIND_R8,indexflag=ESMF_INDEX_GLOBAL, &
      ungriddedLBound=(/1,1,1/),ungriddedUBound=(/nalt,ntime,nfield/),rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldCreate')

    call ESMF_FieldValidate(field=extField,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldValidate')

  endsubroutine init_ext
!-----------------------------------------------------------------------
  subroutine init_int(nnode,nalt,ntime,nfield,nelement,lon,lat,elements)

    use ESMF,only:ESMF_MeshCreate,ESMF_FieldCreate,ESMF_FieldValidate, &
      ESMF_KIND_R8,ESMF_MESHELEMTYPE_QUAD,ESMF_TYPEKIND_R8,ESMF_SUCCESS

    integer,intent(in) :: nnode,nalt,ntime,nfield,nelement
    real(kind=rp),dimension(nnode),intent(in) :: lon,lat
    integer,dimension(4,nelement),intent(in) :: elements

    integer :: ielem,inode,rc
    integer,dimension(nelement*4) :: elementConn
    real(kind=ESMF_KIND_R8),dimension(nnode*2) :: nodeCoords

    do ielem = 1,nelement
      do inode = 1,4
        elementConn((ielem-1)*4+inode) = elements(inode,ielem)
      enddo
    enddo

    do inode = 1,nnode
      nodeCoords(inode*2-1) = lon(inode)
      nodeCoords(inode*2) = lat(inode)
    enddo

    intMesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
      nodeIds=(/(inode,inode=1,nnode)/),nodeCoords=nodeCoords, &
      elementIds=(/(ielem,ielem=1,nelement)/), &
      elementTypes=(/(ESMF_MESHELEMTYPE_QUAD,ielem=1,nelement)/), &
      elementConn=elementConn,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_MeshCreate')

    intField = ESMF_FieldCreate(mesh=intMesh,typekind=ESMF_TYPEKIND_R8, &
      ungriddedLBound=(/1,1,1/),ungriddedUBound=(/nalt,ntime,nfield/),rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldCreate')

    call ESMF_FieldValidate(field=intField,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldValidate')

  endsubroutine init_int
!-----------------------------------------------------------------------
  subroutine init_regrid

    use ESMF,only:ESMF_FieldRegridStore,ESMF_FieldPrint, &
      ESMF_RouteHandlePrint,ESMF_EXTRAPMETHOD_NEAREST_IDAVG,ESMF_SUCCESS

    integer :: srcTermProcessing,pipelineDepth,rc

    srcTermProcessing = 0
    pipelineDepth = 16

    call ESMF_FieldRegridStore( &
      srcField=extField,dstField=intField, &
      extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG, &
      srcTermProcessing=srcTermProcessing, &
      pipelineDepth=pipelineDepth, &
      routehandle=rh,checkFlag=debug,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldRegridStore')

    if (debug) then
      call ESMF_FieldPrint(field=extField,rc=rc)
      if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldPrint')

      call ESMF_FieldPrint(field=intField,rc=rc)
      if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldPrint')

      call ESMF_RouteHandlePrint(routehandle=rh,rc=rc)
      if (rc /= ESMF_SUCCESS) call handle_error('ESMF_RouteHandlePrint')
    endif

  endsubroutine init_regrid
!-----------------------------------------------------------------------
  function regrid(nnode,nlon,nlat,nalt,ntime,nfield,fin) result(fout)

    use ESMF,only:ESMF_FieldGet,ESMF_FieldRegrid, &
      ESMF_KIND_R8,ESMF_TERMORDER_SRCSEQ,ESMF_SUCCESS

    integer,intent(in) :: nnode,nlon,nlat,nalt,ntime,nfield
    real(kind=rp),dimension(nlon,nlat,nalt,ntime,nfield),intent(in) :: fin
    real(kind=rp),dimension(nnode,nalt,ntime,nfield) :: fout

    integer :: rc
    real(kind=ESMF_KIND_R8),dimension(:,:,:,:),pointer :: fptr4d
    real(kind=ESMF_KIND_R8),dimension(:,:,:,:,:),pointer :: fptr5d

    call ESMF_FieldGet(field=extField,farrayPtr=fptr5d,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldGet')

    fptr5d = fin

    call ESMF_FieldRegrid(srcField=extField,dstField=intField, &
      routehandle=rh,termorderflag=ESMF_TERMORDER_SRCSEQ, &
      checkflag=debug,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldRegrid')

    call ESMF_FieldGet(field=intField,farrayPtr=fptr4d,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldGet')

    fout = fptr4d

  endfunction regrid
!-----------------------------------------------------------------------
  subroutine handle_error(funcname)

    use ESMF,only:ESMF_Finalize,ESMF_END_ABORT

    character(len=*),intent(in) :: funcname

    write(6,"('ESMF error in calling ',a,'. Finalizing...')") trim(funcname)
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  endsubroutine handle_error
!-----------------------------------------------------------------------
endmodule regrid_module
