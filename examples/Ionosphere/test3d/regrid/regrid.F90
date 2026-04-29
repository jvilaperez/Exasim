module regrid_module

  use prec,only:rp
  use ESMF,only:ESMF_Mesh,ESMF_Field,ESMF_RouteHandle

  implicit none

  logical,parameter :: debug = .false.

  type(ESMF_Mesh) :: extMesh,intMesh
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
  subroutine init_ext(nnode,nelement,nalt,nfield, &
    elements,nodeGlon,nodeGlat,elementGlon,elementGlat)

    use ESMF,only:ESMF_MeshCreate,ESMF_FieldCreate,ESMF_FieldValidate, &
      ESMF_KIND_R8,ESMF_MESHELEMTYPE_TRI,ESMF_TYPEKIND_R8, &
      ESMF_MESHLOC_ELEMENT,ESMF_SUCCESS

    integer,intent(in) :: nnode,nelement,nalt,nfield
    integer,dimension(3,nelement),intent(in) :: elements
    real(kind=rp),dimension(nnode),intent(in) :: nodeGlon,nodeGlat
    real(kind=rp),dimension(nelement),intent(in) :: elementGlon,elementGlat

    integer :: ielem,inode,rc
    integer,dimension(nelement*3) :: elementConn
    real(kind=ESMF_KIND_R8),dimension(nelement*2) :: elementCoords
    real(kind=ESMF_KIND_R8),dimension(nnode*2) :: nodeCoords

    do ielem = 1,nelement
      elementCoords(ielem*2-1) = elementGlon(ielem)
      elementCoords(ielem*2) = elementGlat(ielem)
      do inode = 1,3
        elementConn((ielem-1)*3+inode) = elements(inode,ielem)
      enddo
    enddo

    do inode = 1,nnode
      nodeCoords(inode*2-1) = nodeGlon(inode)
      nodeCoords(inode*2) = nodeGlat(inode)
    enddo

    extMesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
      nodeIds=(/(inode,inode=1,nnode)/),nodeCoords=nodeCoords, &
      elementIds=(/(ielem,ielem=1,nelement)/), &
      elementTypes=(/(ESMF_MESHELEMTYPE_TRI,ielem=1,nelement)/), &
      elementConn=elementConn,elementCoords=elementCoords,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_MeshCreate')

    extField = ESMF_FieldCreate(mesh=extMesh, &
      typekind=ESMF_TYPEKIND_R8,meshloc=ESMF_MESHLOC_ELEMENT, &
      ungriddedLBound=(/1,1/),ungriddedUBound=(/nalt,nfield/),rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldCreate')

    call ESMF_FieldValidate(field=extField,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldValidate')

  endsubroutine init_ext
!-----------------------------------------------------------------------
  subroutine init_int(nnode,nelement,nalt,nfield,elements,lon,lat)

    use ESMF,only:ESMF_MeshCreate,ESMF_FieldCreate,ESMF_FieldValidate, &
      ESMF_KIND_R8,ESMF_MESHELEMTYPE_QUAD,ESMF_TYPEKIND_R8,ESMF_SUCCESS

    integer,intent(in) :: nnode,nelement,nalt,nfield
    integer,dimension(4,nelement),intent(in) :: elements
    real(kind=rp),dimension(nnode),intent(in) :: lon,lat

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
      ungriddedLBound=(/1,1/),ungriddedUBound=(/nalt,nfield/),rc=rc)
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
  function regrid(nelement,nnode,nalt,nfield,fin) result(fout)

    use ESMF,only:ESMF_FieldGet,ESMF_FieldRegrid, &
      ESMF_KIND_R8,ESMF_TERMORDER_SRCSEQ,ESMF_SUCCESS

    integer,intent(in) :: nelement,nnode,nalt,nfield
    real(kind=rp),dimension(nelement,nalt,nfield),intent(in) :: fin
    real(kind=rp),dimension(nnode,nalt,nfield) :: fout

    integer :: rc
    real(kind=ESMF_KIND_R8),dimension(:,:,:),pointer :: farrayPtr

    call ESMF_FieldGet(field=extField,farrayPtr=farrayPtr,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldGet')

    farrayPtr = fin

    call ESMF_FieldRegrid(srcField=extField,dstField=intField, &
      routehandle=rh,termorderflag=ESMF_TERMORDER_SRCSEQ, &
      checkflag=debug,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldRegrid')

    call ESMF_FieldGet(field=intField,farrayPtr=farrayPtr,rc=rc)
    if (rc /= ESMF_SUCCESS) call handle_error('ESMF_FieldGet')

    fout = farrayPtr

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
