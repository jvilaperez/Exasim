#!/usr/bin/env python
# coding: utf-8


import numpy
import esmpy


def interp2d(variable, lat, lon, latitude, longitude):
    nlat = len(lat)
    nlon = len(lon)

    grid = esmpy.Grid(max_index=numpy.array([nlon, nlat]),
                      num_peri_dims=1,
                      periodic_dim=0,
                      pole_dim=1,
                      coord_sys=esmpy.CoordSys.SPH_DEG,
                      coord_typekind=esmpy.TypeKind.R8,
                      staggerloc=esmpy.StaggerLoc.CENTER,
                      pole_kind=[esmpy.PoleKind.NONE, esmpy.PoleKind.MONOPOLE])

    loncoord = grid.get_coords(coord_dim=0, staggerloc=esmpy.StaggerLoc.CENTER)
    for ilat in range(nlat):
        loncoord[:, ilat] = lon

    latcoord = grid.get_coords(coord_dim=1, staggerloc=esmpy.StaggerLoc.CENTER)
    for ilon in range(nlon):
        latcoord[ilon, :] = lat

    srcfield = esmpy.Field(grid=grid,
                           typekind=esmpy.TypeKind.R8,
                           staggerloc=esmpy.StaggerLoc.CENTER)
    srcfield.data[:] = variable.T

    locstream = esmpy.LocStream(location_count=len(latitude), coord_sys=esmpy.CoordSys.SPH_DEG)

    locstream['ESMF:Lon'] = longitude
    locstream['ESMF:Lat'] = latitude

    dstfield = esmpy.Field(grid=locstream,
                           typekind=esmpy.TypeKind.R8,
                           staggerloc=esmpy.StaggerLoc.CENTER)

    regrid = esmpy.Regrid(srcfield=srcfield,
                          dstfield=dstfield,
                          regrid_method=esmpy.RegridMethod.BILINEAR,
                          pole_method=esmpy.PoleMethod.NONE,
                          line_type=esmpy.LineType.GREAT_CIRCLE,
                          extrap_method=esmpy.ExtrapMethod.NEAREST_IDAVG)

    return dstfield.data
