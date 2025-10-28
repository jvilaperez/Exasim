import numpy as np
from pathlib import Path

class MeshGroup:
    def __init__(self, combined):
        self.combined = combined

    def __getitem__(self, key):
        if isinstance(key, int):
            return {attr: vals[key] for attr, vals in self.combined.items()}
        elif isinstance(key, str):
            return self.combined[key]
        else:
            raise TypeError("Key must be int (rank) or str (attribute)")

    def __len__(self):
        first_attr = next(iter(self.combined.values()))
        return len(first_attr)
    
    def __repr__(self):
        keys = list(self.combined.keys())
        n_parts = len(self)
        return f"<mesh struct: attributes={keys}, n_partitions={n_parts}>"



def readmeshstruct(filemesh):
    mesh = {'nsize' : []};

    tm = np.fromfile(open(filemesh, "r"), dtype=np.float64);
    tm = np.int_(tm);

    sz = np.int_(tm[0]);
    k1 = 1;
    k2 = k1+(sz);
    mesh['nsize'] = np.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][0];
    mesh['ndims'] = np.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][1];
    mesh['facecon'] = np.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][2];
    mesh['eblks'] = np.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][3];
    mesh['fblks'] = np.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][4];
    mesh['nbsd'] = np.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][5];
    mesh['elemsend'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][6];
    mesh['elemrecv'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][7];
    mesh['elemsendpts'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][8];
    mesh['elemrecvpts'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][9];
    mesh['elempart'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][10];
    mesh['elempartpts'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][11];
    mesh['cgelcon'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][12];
    mesh['rowent2elem'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][13];
    mesh['cgent2dgent'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][14];
    mesh['colent2elem'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][15];
    mesh['rowe2f1'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][16];
    mesh['cole2f1'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][17];
    mesh['ent2ind1'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][18];
    mesh['rowe2f2'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][19];
    mesh['cole2f2'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][20];
    mesh['ent2ind2'] = tm[k1:k2];

    return mesh


def readmeshmpi(base, nprocs):

    if nprocs == 1:  
        filesol = f"{base}.bin"
        if not Path(filesol).exists():
            raise FileNotFoundError(f"Cannot open file: {filesol}")
        mesh_single = readmeshstruct(filesol)

        mesh = {key: [val] for key, val in mesh_single.items()}
        return MeshGroup(mesh)

    meshes = []
    mesh = {}
    for i in range(1, nprocs + 1):
        filesol = f"{base}{i}.bin"
        if not Path(filesol).exists():
            raise FileNotFoundError(f"Cannot open file: {filesol}")
        meshi = readmeshstruct(filesol)
        meshes.append(meshi)

    keys = meshes[0].keys()
    mesh = {key: [m[key] for m in meshes] for key in keys}

    return MeshGroup(mesh)