import struct
from profileclass import Profile


class DoseMatrix:
    def __init__(self, numHist=1, weight=1,
                 nx=1, ny=1, nz=1,
                 XPosition=0.0, YPosition=0.0, ZPosition=0.0,
                 xPixelSize=1.0, yPixelSize=1.0, zPixelSize=1.0):
        self.numHist = numHist  # double, 8 bytes
        self.weight = weight  # float, 4 bytes
        self.nx = nx  # int, 4 bytes
        self.ny = ny  # int, 4 bytes
        self.nz = nz  # int, 4 bytes
        self.XPosition = XPosition  # float, 4 bytes
        self.YPosition = YPosition  # float, 4 bytes
        self.ZPosition = ZPosition  # float, 4 bytes
        self.xPixelSize = xPixelSize  # float, 4 bytes
        self.yPixelSize = yPixelSize  # float, 4 bytes
        self.zPixelSize = zPixelSize  # float, 4 bytes
        self.dosArr = [[[0.0 for _ in range(nz)] for _ in range(ny)] for _ in range(nx)]  # dose matrix
        self.uncArr = [[[0.0 for _ in range(nz)] for _ in range(ny)] for _ in range(nx)]  # uncertainty matrix
        self.isUnc = False

    def load_from_bindose_file(self, file_path, ):
        with open(file_path, 'rb') as file:
            # Reading header information
            self.numHist = struct.unpack('d', file.read(8))[0]
            self.weight = struct.unpack('f', file.read(4))[0]
            self.nx = struct.unpack('i', file.read(4))[0]
            self.ny = struct.unpack('i', file.read(4))[0]
            self.nz = struct.unpack('i', file.read(4))[0]
            self.XPosition = struct.unpack('f', file.read(4))[0]
            self.YPosition = struct.unpack('f', file.read(4))[0]
            self.ZPosition = struct.unpack('f', file.read(4))[0]
            self.xPixelSize = struct.unpack('f', file.read(4))[0]
            self.yPixelSize = struct.unpack('f', file.read(4))[0]
            self.zPixelSize = struct.unpack('f', file.read(4))[0]
            self.isUnc = True

            total_voxels = self.nx * self.ny * self.nz

            # Initializing dose and uncertainty arrays
            self.dosArr = [0.0] * total_voxels
            self.uncArr = [0.0] * total_voxels

            # Reading dose values
            for i in range(total_voxels):
                self.dosArr[i] = struct.unpack('f', file.read(4))[0]

            # Reading uncertainty values
            for i in range(total_voxels):
                self.uncArr[i] = struct.unpack('f', file.read(4))[0]

    def save_to_file(self, file_path):
        # Implement file writing logic here
        pass

    def dose_matrix(self, ix, iy, iz):
        index = ix + iy * self.nx + iz * self.nx * self.ny
        return self.dosArr[index]

    def unc_matrix(self, ix, iy, iz):
        index = ix + iy * self.nx + iz * self.nx * self.ny
        return self.uncArr[index]

    def value_from_position(self, x, y, z, mode='D'):
        if mode == 'D':
            arr = self.dosArr
        elif mode == 'U':
            arr = self.uncArr
        else:
            print('Incorrect mode. Use D for dose and U for uncertainty.')
            return None

        # Convert position to voxel index
        ix = int((x - self.XPosition) / self.xPixelSize)
        iy = int((y - self.YPosition) / self.yPixelSize)
        iz = int((z - self.ZPosition) / self.zPixelSize)

        # Calculate distances from point to voxel boundaries
        dx = (x - self.XPosition) / self.xPixelSize - ix
        dy = (y - self.YPosition) / self.yPixelSize - iy
        dz = (z - self.ZPosition) / self.zPixelSize - iz

        # Ensure indices are within bounds
        ix = max(0, min(self.nx - 2, ix))
        iy = max(0, min(self.ny - 2, iy))
        iz = max(0, min(self.nz - 2, iz))

        # Tri-linear interpolation
        c00 = arr[ix + iy * self.nx + iz * self.nx * self.ny] * (1 - dx) + \
              arr[(ix + 1) + iy * self.nx + iz * self.nx * self.ny] * dx
        c01 = arr[ix + iy * self.nx + (iz + 1) * self.nx * self.ny] * (1 - dx) + \
              arr[(ix + 1) + iy * self.nx + (iz + 1) * self.nx * self.ny] * dx
        c10 = arr[ix + (iy + 1) * self.nx + iz * self.nx * self.ny] * (1 - dx) + \
              arr[(ix + 1) + (iy + 1) * self.nx + iz * self.nx * self.ny] * dx
        c11 = arr[ix + (iy + 1) * self.nx + (iz + 1) * self.nx * self.ny] * (1 - dx) + \
              arr[(ix + 1) + (iy + 1) * self.nx + (iz + 1) * self.nx * self.ny] * dx

        c0 = c00 * (1 - dy) + c10 * dy
        c1 = c01 * (1 - dy) + c11 * dy

        c = c0 * (1 - dz) + c1 * dz

        return c

    def get_dose_profile(self, idim, off1, off2):
        if idim == 0:
            nbin = self.nx
            binsize = self.xPixelSize
            origen = self.XPosition
        elif idim == 1:
            nbin = self.ny
            binsize = self.yPixelSize
            origen = self.YPosition
        elif idim == 2:
            nbin = self.nz
            binsize = self.zPixelSize
            origen = self.ZPosition
        else:
            raise ValueError("Invalid dimension")

        valores = [0.0] * nbin
        for i in range(nbin):
            x = origen + i * binsize
            if idim == 0:
                valores[i] = self.value_from_position(x, off1, off2)
            elif idim == 1:
                valores[i] = self.value_from_position(off1, x, off2)
            elif idim == 2:
                valores[i] = self.value_from_position(off1, off2, x)

        # Assuming DoseProfile is a class you have defined elsewhere
        dp = Profile(nbin=nbin, binsize=binsize, origin=origen, values=valores)

        return dp
