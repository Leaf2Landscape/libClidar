# libClidar

A collection of C programs and libraries for processing lidar data, including voxelisation of terrestrial and airborne laser scanning (TLS/ALS) point clouds, simulation of large-footprint full-waveform lidar (e.g. GEDI), and estimation of plant area index (PAI) and 3D vegetation structure.

Developed by Steven Hancock. Released under the GNU General Public License v3.

## Components

| Directory | Description | Key tools |
|-----------|-------------|-----------|
| **libclidar** | Core C library for lidar I/O and voxel operations | (library only) |
| **gedisimulator** | Simulates large-footprint full-waveform lidar from ALS data | `gediRat`, `gediMetric`, `mapLidar`, `lasPoints`, `collocateWaves`, `addNoiseHDF` |
| **voxelate** | Creates voxel maps from full-waveform ALS or TLS data | `voxelate`, `voxelTLS`, `dumpTLScloud`, `voxDecimate` |
| **voxel_lidar** | PAI mapping and gap-fraction calibration from ALS+TLS | `optimiseTLS`, `rxpVox` |
| **rxptobinary** | Converts Riegl .rxp scan files to binary or ASCII | `readRXP` |
| **salcamake** | Processes data from the SALCA dual-wavelength TLS | `hedgehog`, `salcaCal`, `tiffSALCA` |

## Dependencies

- **System libraries:** GSL, HDF5, GDAL, libgeotiff, libtiff, proj
- **Hancock C-tools:** https://bitbucket.org/StevenHancock/tools
- **cmpfit / Minpack:** https://www.physics.wisc.edu/~craigm/idl/cmpfit.html
- **Riegl RiVLib** (for `readRXP` only): proprietary, not included

## Building with Dev Container

A `.devcontainer` configuration is provided for building all components. It uses Fedora 43 and installs all required dependencies automatically.

### Prerequisites

- Docker or Podman
- VS Code with the Dev Containers extension (or any devcontainer-compatible tool)
- Riegl RiVLib and RDBLib archives placed in `riegl_libs/`:
  - `rivlib-2_6_0-x86_64-linux-gcc11.zip`
  - `rdblib-2.4.1-x86_64-linux.tar.gz`

### Usage

1. Open the project in VS Code and select **Reopen in Container**, or build manually:
   ```bash
   cd .devcontainer
   docker build -t libclidar -f Dockerfile ..
   docker run -it -v $(pwd)/..:/workspaces/libClidar libclidar bash
   ```

2. Inside the container, run the build script:
   ```bash
   bash .devcontainer/build-all.sh
   ```
   Compiled binaries are placed in `bin/`.

### Environment variables

The container sets these automatically. If building outside the container, you must define them:

| Variable | Description |
|----------|-------------|
| `HANCOCKTOOLS_ROOT` | Path to Hancock C-tools |
| `CMPFIT_ROOT` | Path to cmpfit/Minpack |
| `LIBCLIDAR_ROOT` | Path to `libclidar/` |
| `VOXELATE_ROOT` | Path to `voxelate/` |
| `GEDIRAT_ROOT` | Path to `gedisimulator/` |
| `GSL_ROOT` | Path to GSL library |
| `HDF5_LIB` | Path to HDF5 installation |
| `RiVLib_DIR` | Path to RiVLib (for readRXP) |

## Example: TLS voxelisation pipeline

```bash
# 1. Convert Riegl RXP to binary format (applying SOP transform)
readRXP -input scan.rxp -output scan.bin -trans scan.DAT

# 2. Create voxel map (1m resolution)
voxelTLS -input scan.bin -res 1 1 1 -output voxels.dat

# 3. Export voxel point cloud
dumpTLScloud -input voxels.dat -output cloud.asc
```

## References

- Hancock, S., et al. (2019). The GEDI Simulator: A Large-Footprint Waveform Lidar Simulator for Calibration and Validation of Spaceborne Missions. *Earth and Space Science*, 6, 294-310.

## License

GNU General Public License v3.0 -- see individual source file headers for details.
