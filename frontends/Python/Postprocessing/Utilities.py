# import external modules
import numpy as np
import os

import Preprocessing, Postprocessing, Gencode, Mesh


def retrieveData(inpath,outpath):

    pde = Preprocessing.readapp(inpath + "/app.bin")
    master = Postprocessing.getmaster(inpath+"/master.bin")
    mpiprocs = pde['ndims'][0]
    mesh = Postprocessing.readmeshmpi(inpath + "/mesh",mpiprocs)
    _, sol = Postprocessing.readsolmpi(outpath + "/outudg",mpiprocs)
    xdg, udg, vdg = Postprocessing.readxdgsolmpi(inpath + "/sol", mesh["elempart"])
    

    return pde, master, mesh, xdg, udg, vdg, sol


