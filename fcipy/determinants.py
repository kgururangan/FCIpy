import numpy as np

def read_determinants(file):

    with open(file, "r") as f:
        for i, line in enumerate(f.readlines()):
            if i == 0: # first line stores N_int, something, ndet
                N_int, _, ndet = [int(x) for x in line.split()]
                det = np.zeros(N_int * 2 * ndet, dtype=np.int64, order="F")
                continue
            # load the numbers into a flat array
            det[i - 1] = np.int64(line.split()[0])

    # Reshape the array in Fortran order (F order is important!)
    det = np.reshape(det, (N_int, 2, ndet), order="F")
    return det
