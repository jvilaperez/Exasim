import os
import numpy as np

def readsolpart(base, elempart, nsteps=1, stepoffsets=0):
    nproc = len(elempart)
    k = np.zeros((nproc, 4), dtype=float)

    # --- First pass: read headers and metadata
    for i in range(nproc):
        fname = f"{base}_np{i}.bin"
        if not os.path.exists(fname):
            raise FileNotFoundError(f"Cannot open file: {fname}")

        with open(fname, "rb") as f:
            hdr = np.fromfile(f, dtype=np.float64, count=3)
            if hdr.size != 3:
                raise ValueError(f"Header read failed for {fname}")

            n1, n2, n3 = hdr
            N = int(n1 * n2 * n3)

            L = os.path.getsize(fname) / 8  # number of doubles
            timesteps = (L - 3) / N
            k[i, :] = [n1, n2, n3, timesteps]

    # --- Consistency checks
    if not np.allclose(k[:, 0], k[0, 0]):
        raise ValueError("First dimension in binary files do not match.")
    if not np.allclose(k[:, 1], k[0, 1]):
        raise ValueError("Second dimension in binary files do not match.")
    if not np.allclose(k[:, 3], k[0, 3]):
        raise ValueError("Number of timesteps in binary files do not match.")

    # --- If user requested sol output
    if nsteps is not None and nsteps > 0:
        n1, n2 = int(k[0, 0]), int(k[0, 1])
        ne = int(np.sum(k[:, 2]))
        sol = np.zeros((n1, n2, ne, nsteps), dtype=np.float64)

        offset_bytes = int(stepoffsets * n1 * n2 * 8)

        for n in range(nproc):
            fname = f"{base}_np{n}.bin"
            with open(fname, "rb") as f:
                hdr = np.fromfile(f, dtype=np.float64, count=3)
                if hdr.size != 3:
                    raise ValueError(f"Header read failed for {fname}")

                n1, n2, n3 = map(int, hdr)
                N = n1 * n2 * n3

                if stepoffsets > 0:
                    f.seek(stepoffsets * N * 8, os.SEEK_CUR)

                elempartK = np.array(elempart[n][:n3], dtype=int)
                # if elempart.min() == 1:
                #     elempart = elempart - 1
                
                for i in range(nsteps):
                    tm = np.fromfile(f, dtype=np.float64, count=N)
                    if tm.size != N:
                        raise ValueError(f"Unexpected end of file in {fname}")
                    sol[:, :, elempartK, i] = tm.reshape((n1, n2, n3), order='F')

        return k, sol

    return k, None