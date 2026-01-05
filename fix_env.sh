#!/bin/bash
# ------------------------------------------------------------------
# Fix GLIBCXX and library path issues for dd-env
# ------------------------------------------------------------------

# Activate the conda environment
source ~/anaconda3/etc/profile.d/conda.sh   # adjust path if needed
conda activate dd-env

# Ensure env libraries are preferred over system ones
export _OLD_LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Prevent interference from user site packages
export PYTHONNOUSERSITE=1

# Force-load the conda libstdc++ if system version is too old
if [ -f "${CONDA_PREFIX}/lib/libstdc++.so.6" ]; then
  export _OLD_LD_PRELOAD="${LD_PRELOAD:-}"
  export LD_PRELOAD="${CONDA_PREFIX}/lib/libstdc++.so.6${LD_PRELOAD:+:$LD_PRELOAD}"
fi

echo "dd-env is now active"
