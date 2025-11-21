import numpy as np
import os

import Preprocessing, Postprocessing, Gencode, Mesh

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.tri import Triangulation
from cartopy.crs import PlateCarree, NorthPolarStereo, SouthPolarStereo
from matplotlib.collections import LineCollection


def getgridelements(sdg2d,master):

    p2d,t2d = Postprocessing.getnodes(sdg2d[:,0:2,:],master['elemtype'],master['porder'])
    xgrid = p2d[:,t2d].transpose(1,0,2)

    mask_xgrid = (xgrid[:,0,:].min(axis=0) < 0) & (xgrid[:,0,:].max(axis=0) > 0)

    # Add 2*pi to negative nodes of those elements
    xgrid[:, 0, mask_xgrid] += (xgrid[:, 0, mask_xgrid] < 0) * 360
    lon0 = xgrid[:, 0, mask_xgrid] - (xgrid[:, 0, mask_xgrid]>0)*360
    lat0 = xgrid[:, 1, mask_xgrid]

    # append lon0,lat0 as new elements
    xgrid0 = np.array([lon0, lat0]).transpose(1,0,2)
    xgrid = np.append(xgrid, xgrid0, axis=2)

    return xgrid


def getslice(xdg,udg,master,ilev):
    """
    xdg: (npe, 3, ne) high-order DG mesh on cube-sphere
    udg: (npe, ne) high-order solution DG mesh on xdg
    master struct with solution information
    ilev: layer at which we plot solution, e.g. 0-nlev
    """    

    x = xdg[:,0,:]; y = xdg[:,1,:]; z = xdg[:,2,:]
    r = np.sqrt(x**2 + y**2 + z**2); lon = np.rad2deg(np.arctan2(y,x)); lat = np.rad2deg(np.arcsin(z/r))

    # Identify layer of elements
    tol = 1e-3  # tolerance for uniqueness
    altdg = np.round(r / tol) * tol
    alt1d = np.unique(altdg)
    altObj = alt1d[ilev]

    layerdg = (altdg == altObj)
    good_elem = np.where(np.any(layerdg, axis=0))[0]
    perm = np.unique(layerdg[:,good_elem],axis=1)
    good_elem = np.where(np.all(layerdg == perm[:, 0][:, None], axis=0))[0]

    # Extract mesh/field in such layer
    xdg2d = xdg[perm[:,0], :,:]
    xdg2d = xdg2d[:,:,good_elem]

    udg2d = udg[perm[:, 0],:]
    udg2d = udg2d[:,good_elem]

    # Create high-order grid for display
    sdg = np.zeros_like(xdg)
    sdg[:,0,:] = lon; sdg[:,1,:] = lat; sdg[:,2,:] = r
    sdg2d = sdg[perm[:, 0],:,:]
    sdg2d = sdg2d[:,:,good_elem] 
    xgrid = getgridelements(sdg2d,master)

    # Increase visualization resolution
    visorder = 3*master['porder']
    xpe0,telem0 = Preprocessing.masternodes(master['porder'],master['nd']-1,master['elemtype'])[0:2]
    xpe,telem = Preprocessing.masternodes(visorder,master['nd']-1,master['elemtype'])[0:2]
    visshape = Preprocessing.mkshape(master['porder'],xpe0,xpe,master['elemtype'])
    visshape = visshape[:,:,0].T

    npe2d, nd2d, nelem = xdg2d.shape
    xdgvis = np.matmul(visshape,np.reshape(xdg2d,(visshape.shape[1], nd2d*nelem), 'F'));
    xdgvis = np.reshape(xdgvis,(visshape.shape[0], nd2d, nelem),'F');

    udgvis = np.matmul(visshape,np.reshape(udg2d,(visshape.shape[1], nelem), 'F'));
    udgvis = np.reshape(udgvis,(visshape.shape[0], nelem),'F');


    # Create CG mesh and visualization field
    cgnodes, cgelcon, cgcells, celltype = Postprocessing.createcggrid(xdgvis,telem)[0:4]
    xcg2d = cgnodes[cgcells.astype(int),:].transpose(1,2,0)

    r2d = np.sqrt(xcg2d[:,0,:]**2 + xcg2d[:,1,:]**2 + xcg2d[:,2,:]**2)
    lon2d = np.arctan2(xcg2d[:,1,:], xcg2d[:,0,:])
    lat2d = np.arctan2(xcg2d[:,2,:], np.sqrt(xcg2d[:,0,:]**2 + xcg2d[:,1,:]**2))

    # Correcting visualization at transition
    mask_elem = (lon2d.min(axis=0) < 0) & (lon2d.max(axis=0) > 0)

    # Add 2*pi to negative nodes of those elements
    lon2d[:, mask_elem] += (lon2d[:, mask_elem] < 0) * 2*np.pi
    # append new element with negative longitudes
    lon2d = np.append(lon2d, lon2d[:, mask_elem] - (lon2d[:, mask_elem]>0)*2*np.pi, axis=1)
    lat2d = np.append(lat2d, lat2d[:, mask_elem], axis=1)

    xcg2d = np.zeros((lon2d.shape[0],2,lon2d.shape[1]))
    xcg2d[:,0,:] = np.rad2deg(lon2d)
    xcg2d[:,1,:] = np.rad2deg(lat2d)

    ucg2d = udgvis[telem-1,:].transpose(1,0,2).reshape(telem.shape[1],-1,order= 'F')
    ucg2d = np.append(ucg2d, ucg2d[:, mask_elem], axis=1)

    return xcg2d, ucg2d, altObj, xgrid


def plot_mesh_edges(xgrid, ax, color='k', lw=0.05):
    """
    xgrid: (npe, 2, ne) array of element vertex coordinates (lon,lat)
    ax : matplotlib axes
    Draws only the edges of the mesh.
    """
    lines = []
    npe, _, ne = xgrid.shape

    for e in range(ne):
        poly = xgrid[:, :, e]    # shape (npe,2)

        # Loop around edges
        for i in range(npe):
            p1 = poly[i]
            p2 = poly[(i+1) % npe]
            lines.append([p1, p2])

    lc = LineCollection(lines, colors=color, linewidths=lw)
    ax.add_collection(lc)



def plotslice(xcg2d,ucg2d, xgrid=None, coastlines=True, coastcolor='k', coastlw=0.5, color='w', lw=0.05):
    ne = xcg2d.shape[2] # Extract number of elements

    # Collect all vertices and triangulate each polygon
    all_x = []
    all_y = []
    all_vals = []
    triangles = []
    offset = 0

    for e in range(ne):
        poly = xcg2d[:,:, e]       # shape (npe, 2)
        x = poly[:, 0]
        y = poly[:, 1]
        val = ucg2d[:, e]
        
        # Triangulate the polygon (simple fan triangulation from first vertex)
        n = poly.shape[0]
        for i in range(1, n-1):
            triangles.append([offset, offset+i, offset+i+1])
        
        all_x.extend(x)
        all_y.extend(y)
        all_vals.extend(val)
        offset += n

    all_x = np.array(all_x)
    all_y = np.array(all_y)
    all_vals = np.array(all_vals)
    triangles = np.array(triangles)

    # Create Triangulation
    triang = Triangulation(all_x, all_y, triangles)

    # Normalize colors
    cmap = cm.viridis
    norm = plt.Normalize(vmin=all_vals.min(), vmax=all_vals.max())

    # Plot using tripcolor with linear interpolation
    fig, ax = plt.subplots(figsize=(10, 6),subplot_kw=dict(projection=PlateCarree()))
    # ax.set_global()  # Sets extent to [-180, 180, -90, 90]
    tpc = ax.tripcolor(triang, all_vals, shading='gouraud', cmap=cmap, norm=norm,transform=PlateCarree())

    if coastlines:
        ax.coastlines(color=coastcolor, linewidth=coastlw)

    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=PlateCarree())
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=PlateCarree())

    ax.set_aspect('equal')
    ax.set_xlim([-180, 180])
    ax.set_ylim([-90, 90])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # ax.set_title('Field visualization on unstructured quad mesh')

    # Add colorbar
    fig.colorbar(tpc, ax=ax, shrink=0.7)
    if xgrid is not None:
        plot_mesh_edges(xgrid, ax, color, lw)

    return fig, ax