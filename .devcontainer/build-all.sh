#!/usr/bin/env bash
#
# Build all libClidar sub-projects inside the devcontainer.
# Run from the workspace root: /workspaces/libClidar
#
set -uo pipefail

WORKSPACE="/workspaces/libClidar"
BINDIR="$WORKSPACE/bin"
mkdir -p "$BINDIR"

export ARCH=${ARCH:-$(uname -m)}
export HANCOCKTOOLS_ROOT=${HANCOCKTOOLS_ROOT:-/opt/tools}
export CMPFIT_ROOT=${CMPFIT_ROOT:-/opt/cmpfit}
export LIBCLIDAR_ROOT=${LIBCLIDAR_ROOT:-$WORKSPACE/libclidar}
export VOXELATE_ROOT=${VOXELATE_ROOT:-$WORKSPACE/voxelate}
export GEDIRAT_ROOT=${GEDIRAT_ROOT:-$WORKSPACE/gedisimulator}
export GSL_ROOT=${GSL_ROOT:-/usr/lib64}
export HDF5_LIB=${HDF5_LIB:-/usr}
export RiVLib_DIR=${RiVLib_DIR:-/opt/rivlib/rivlib-2_6_0-x86_64-linux-gcc11}

ok=0
fail=0
report=""

build_make_target() {
    local dir="$1"
    local target="$2"
    local name
    name=$(basename "$dir")

    echo "──────────────────────────────────────────────"
    echo "  Building $name :: $target"
    echo "──────────────────────────────────────────────"

    pushd "$dir" > /dev/null
    if make THIS="$target" 2>&1; then
        cp "$target" "$BINDIR/" 2>/dev/null || true
        echo "  ✓ $target"
        report+="  OK   $name/$target\n"
        ((ok++)) || true
    else
        echo "  ✗ $target FAILED"
        report+="  FAIL $name/$target\n"
        ((fail++)) || true
    fi
    popd > /dev/null
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. libclidar  (object files only — compiled as part of dependent projects)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════"
echo "  [1/6] libclidar — shared library objects"
echo "══════════════════════════════════════════════"
echo "  (compiled automatically by dependent projects)"

# ─────────────────────────────────────────────────────────────────────────────
# 2. gedisimulator
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════"
echo "  [2/6] gedisimulator"
echo "══════════════════════════════════════════════"
cd "$WORKSPACE/gedisimulator" && make clean 2>/dev/null || true

for prog in gediRat gediMetric mapLidar collocateWaves lasPoints \
            addNoiseHDF lgw2hdf fitTXpulse clipLidar gediAscii2hdf; do
    build_make_target "$WORKSPACE/gedisimulator" "$prog"
done

# ─────────────────────────────────────────────────────────────────────────────
# 3. salcamake
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════"
echo "  [3/6] salcamake"
echo "══════════════════════════════════════════════"
cd "$WORKSPACE/salcamake" && make clean 2>/dev/null || true

for prog in hedgehog; do
    build_make_target "$WORKSPACE/salcamake" "$prog"
done

# cal_salca sub-programs
cd "$WORKSPACE/salcamake/cal_salca" && make clean 2>/dev/null || true
for prog in salcaCal; do
    build_make_target "$WORKSPACE/salcamake/cal_salca" "$prog"
done

# tiffSALCA
cd "$WORKSPACE/salcamake/tiffSALCA" && make clean 2>/dev/null || true
for prog in tiffSALCA; do
    build_make_target "$WORKSPACE/salcamake/tiffSALCA" "$prog"
done

# ─────────────────────────────────────────────────────────────────────────────
# 4. voxelate
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════"
echo "  [4/6] voxelate"
echo "══════════════════════════════════════════════"
cd "$WORKSPACE/voxelate" && make clean 2>/dev/null || true

for prog in voxelate voxelTLS voxDecimate dumpTiffAscii dumpTLScloud; do
    build_make_target "$WORKSPACE/voxelate" "$prog"
done

# ─────────────────────────────────────────────────────────────────────────────
# 5. voxel_lidar
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════"
echo "  [5/6] voxel_lidar"
echo "══════════════════════════════════════════════"
cd "$WORKSPACE/voxel_lidar" && make clean 2>/dev/null || true

for prog in compareWaves optimiseTLS; do
    build_make_target "$WORKSPACE/voxel_lidar" "$prog"
done

# ─────────────────────────────────────────────────────────────────────────────
# 6. rxptobinary  (CMake build — needs RiVLib)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════"
echo "  [6/6] rxptobinary (CMake)"
echo "══════════════════════════════════════════════"

RXPDIR="$WORKSPACE/rxptobinary"
RXPBUILD="$RXPDIR/build"
mkdir -p "$RXPBUILD"

if cmake -S "$RXPDIR" -B "$RXPBUILD" \
       -DRiVLib_DIR="$RiVLib_DIR" \
       -DCMAKE_INSTALL_PREFIX="$BINDIR" 2>&1 \
   && cmake --build "$RXPBUILD" 2>&1; then
    cp "$RXPBUILD/readRXP" "$BINDIR/" 2>/dev/null || true
    echo "  ✓ readRXP"
    report+="  OK   rxptobinary/readRXP\n"
    ((ok++)) || true
else
    echo "  ✗ readRXP FAILED"
    report+="  FAIL rxptobinary/readRXP\n"
    ((fail++)) || true
fi

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════"
echo "  Build summary"
echo "══════════════════════════════════════════════"
echo -e "$report"
echo "  Succeeded: $ok   Failed: $fail"
echo ""
echo "  Binaries installed to: $BINDIR"
echo "══════════════════════════════════════════════"
