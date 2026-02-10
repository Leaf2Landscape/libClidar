#!/usr/bin/env python3

import warnings
import numpy as np
import argparse
from sys import exit


###############################
# Combines output of voxelTLS #
# from multiple scans         #
#                             #
# svenhancock@gmail.com       #
# November 2018               #
###############################

# Column indices (0-based) in voxelTLS output
# x y z gap0 inHit0 inMiss0 gapTo0 hits0 miss0 PAIbel0 volSamp0 volTrans0 PAIrad0 meanRefl0 meanZen0 meanRange0
COL_X, COL_Y, COL_Z = 0, 1, 2
COL_GAP, COL_GAPTO, COL_HITS = 3, 6, 7
COL_PAIB, COL_VOLSAMP, COL_VOLTRANS, COL_PAIRAD = 9, 10, 11, 12
DATA_COLS = [COL_GAP, COL_GAPTO, COL_HITS, COL_PAIB, COL_VOLSAMP, COL_VOLTRANS, COL_PAIRAD]
# Indices into the extracted data array (after selecting DATA_COLS)
I_GAP, I_GAPTO, I_HITS, I_PAIB, I_VOLSAMP, I_VOLTRANS, I_PAIRAD = range(7)

NCOLS_PER_SCAN = 13  # number of data columns per scan in voxelTLS output

HEADER = ("x y z 4_mGapTo 5_maxGapTo 6_mGapIn 7_mPAIb 8_mPAIrad "
          "9_totBeam 10_bPAIb 11_maxVol 12_maxVolTrans "
          "13_mPAIb_volTrans 14_mPAIb_nBeams")


def getCmdArgs():
    p = argparse.ArgumentParser(
        description="Combines voxel PAI estimates from multiple scans")
    p.add_argument("--input", dest="inNamen", default=" ",
                   type=str, help="Input of combined scans")
    p.add_argument("--inList", dest="listNamen", default=" ",
                   type=str, help="Input list of single scans")
    p.add_argument("--output", dest="outNamen", default="meanVox.vox",
                   type=str, help="Output filename (default: meanVox.vox)")
    p.add_argument("--minGap", dest="minGap", type=float, default=0.0,
                   help="Minimum trustable gap fraction (default: 0)")
    return p.parse_args()


def readFile(path):
    """Read a voxelTLS ASCII file, return raw numpy array (skip header)."""
    return np.loadtxt(path, skiprows=1, dtype=np.float32)


def readSingleInput(path):
    """Read a combined multi-scan file. Returns pos (N,3) and data (N, nVars, nScans)."""
    raw = readFile(path)
    pos = raw[:, :3].astype(np.float64)
    nDataCols = raw.shape[1] - 3
    nScans = nDataCols // NCOLS_PER_SCAN
    # Extract only the columns we need for each scan
    scans = []
    for s in range(nScans):
        offset = 3 + s * NCOLS_PER_SCAN
        cols = [offset + (c - 3) for c in DATA_COLS]
        scans.append(raw[:, cols])
    # Stack: (N, nVars, nScans)
    data = np.stack(scans, axis=2)
    return pos, data


def readMultipleInputs(listPath):
    """Read individual scan files listed in a text file. Returns pos and data."""
    with open(listPath) as f:
        files = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    nScans = len(files)

    pos = None
    scans = []
    for i, path in enumerate(files):
        print(f"\rReading scan: {i+1}/{nScans}       ", end="")
        raw = readFile(path)
        if pos is None:
            pos = raw[:, :3].astype(np.float64)
        scans.append(raw[:, DATA_COLS])

    # Stack: (N, nVars, nScans)
    data = np.stack(scans, axis=2)
    return pos, data


def wavg(values, weights):
    """Weighted average along scan axis (axis=1), returning -1 where weights sum to 0."""
    wsum = np.nansum(weights, axis=1)
    good = wsum > 0
    out = np.full(values.shape[0], -1.0, dtype=np.float32)
    out[good] = np.nansum(values[good] * weights[good], axis=1) / wsum[good]
    return out


def combineData(pos, data, minGap):
    """Combine per-scan voxel data into weighted averages.
    data shape: (N, nVars, nScans). Uses NaN masking instead of np.ma."""
    nVox = data.shape[0]

    # Replace nodata with NaN — mask where gapTo < minGap
    gapTo = data[:, I_GAPTO, :].copy()  # (N, nScans)
    bad = gapTo < minGap
    gapTo[bad] = np.nan

    gapIn = data[:, I_GAP, :].copy()
    gapIn[bad] = np.nan
    paib = data[:, I_PAIB, :].copy()
    paib[bad] = np.nan
    pairad = data[:, I_PAIRAD, :].copy()
    pairad[bad] = np.nan
    nbeam = data[:, I_HITS, :].copy()
    nbeam[bad] = np.nan
    volSamp = data[:, I_VOLSAMP, :].copy()
    volSamp[bad] = np.nan
    volTrans = data[:, I_VOLTRANS, :].copy()
    volTrans[bad] = np.nan

    all_nan = np.all(bad, axis=1)

    # Allocate output: (N, 11)
    out = np.full((nVox, 11), -1.0, dtype=np.float32)

    # Weighted averages using gapTo as weight
    w = gapTo.copy()
    w[np.isnan(w)] = 0.0
    out[:, 0] = wavg(gapTo, w)                    # mGapTo
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out[:, 1] = np.where(all_nan, -1.0, np.nanmax(gapTo, axis=1))  # maxGapTo
    out[:, 2] = wavg(gapIn, w)                     # mGapIn
    out[:, 3] = wavg(paib, w)                      # mPAIb
    out[:, 4] = wavg(pairad, w)                    # mPAIrad

    # Total beams (sum, not weighted average)
    nbeam_clean = nbeam.copy()
    nbeam_clean[np.isnan(nbeam_clean)] = 0.0
    out[:, 5] = np.sum(nbeam_clean, axis=1)        # totBeam

    # Beland's best-view PAI (from scan with max gapTo)
    gapTo_for_argmax = gapTo.copy()
    gapTo_for_argmax[np.isnan(gapTo_for_argmax)] = -1.0
    best = np.argmax(gapTo_for_argmax, axis=1)     # (N,)
    rows = np.arange(nVox)
    bPAIb = paib[rows, best]
    bPAIb[np.isnan(bPAIb)] = -1.0
    out[:, 6] = bPAIb                              # bPAIb
    vs = volSamp[rows, best]
    vt = volTrans[rows, best]
    with np.errstate(divide='ignore', invalid='ignore'):
        out[:, 7] = np.where((vt > 0) & ~np.isnan(vt), vs / vt, -1.0)  # maxVol

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out[:, 8] = np.where(all_nan, -1.0, np.nanmax(volTrans, axis=1))  # maxVolTrans

    # PAIb weighted by volTrans
    badVT = data[:, I_VOLTRANS, :] < minGap
    paib_vt = data[:, I_PAIB, :].copy()
    paib_vt[badVT] = np.nan
    wVT = data[:, I_VOLTRANS, :].copy()
    wVT[badVT] = 0.0
    out[:, 9] = wavg(paib_vt, wVT)                # mPAIb_volTrans

    # PAIb weighted by nBeam
    badNB = data[:, I_HITS, :] < minGap
    paib_nb = data[:, I_PAIB, :].copy()
    paib_nb[badNB] = np.nan
    wNB = data[:, I_HITS, :].copy()
    wNB[badNB] = 0.0
    out[:, 10] = wavg(paib_nb, wNB)               # mPAIb_nBeams

    return out


if __name__ == '__main__':
    cmdargs = getCmdArgs()
    inNamen = cmdargs.inNamen
    listNamen = cmdargs.listNamen
    minGap = cmdargs.minGap
    outNamen = cmdargs.outNamen

    # Read data
    if inNamen != " ":
        pos, data = readSingleInput(inNamen)
    elif listNamen != " ":
        pos, data = readMultipleInputs(listNamen)
    else:
        print("No input specified")
        exit(1)

    print(f"\n{data.shape[0]} voxels, {data.shape[2]} scans")

    # Combine
    print("Processing...", end="", flush=True)
    out = combineData(pos, data, minGap)
    print(" done")

    # Write output
    print("Writing...", end="", flush=True)
    result = np.column_stack([pos, out])
    # totBeam column (index 8 = pos 3 cols + out col 5) as integer
    fmt = ['%.4f'] * 3 + ['%.6f'] * 5 + ['%d'] + ['%.6f'] * 5
    np.savetxt(outNamen, result, header=HEADER, comments='', fmt=fmt, delimiter=' ')
    print(f" {outNamen}")
