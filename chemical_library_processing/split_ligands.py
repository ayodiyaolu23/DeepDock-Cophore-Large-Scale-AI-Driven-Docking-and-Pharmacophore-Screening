#!/usr/bin/env python3
import os
import shutil
from math import ceil

# === USER INPUT ===
source_folder = "/home/projects/scratch/diyaolua/directory/test_1_ligands_pdbqt"   # change this
num_folders = 4                                 # number of groups to split into
base_output = os.path.join(source_folder, "test_split_1")  # optional parent folder

# === SCRIPT START ===
os.makedirs(base_output, exist_ok=True)

# Get all files (ignore directories)
all_files = [f for f in os.listdir(source_folder)
             if os.path.isfile(os.path.join(source_folder, f))]

# Sort files to ensure consistent order
all_files.sort()
n_files = len(all_files)
if n_files == 0:
    raise ValueError("No files found in source folder.")

# Calculate chunk size
chunk_size = ceil(n_files / num_folders)

# Split into groups and copy/move
for i in range(num_folders):
    start = i * chunk_size
    end = min(start + chunk_size, n_files)
    group_files = all_files[start:end]

    # Create destination subfolder
    dest_folder = os.path.join(base_output, f"group_{i+1}")
    os.makedirs(dest_folder, exist_ok=True)

    # Move or copy each file
    for f in group_files:
        src_path = os.path.join(source_folder, f)
        dst_path = os.path.join(dest_folder, f)
        shutil.move(src_path, dst_path)  # use move; change to shutil.copy if you want copies

    print(f"Group {i+1}: moved {len(group_files)} files")

print("Done. Files distributed equally into 6 folders.")

