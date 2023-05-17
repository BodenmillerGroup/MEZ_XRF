# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 16:09:11 2022

@author: WinAdm
"""

import os
import shutil
from pathlib import Path
import h5py
import pandas as pd


#%%
parent_dir = Path(r"C:\Users\MerrickS\OneDrive\Work\2_UZH\Papers\1_MEZ_XRF")
in_dir = parent_dir / "paper" / "Additional material" / "XRF_RawData_forpaper" / "XRF_RawData_forpaper" / "h5_files"
out_dir = in_dir.parent / "all_SABER_hdf_files"
out_dir.mkdir(exist_ok=True, parents=True)

#%%

h5_paths = list(in_dir.glob("*/*/*.h5"))

for hdf in h5_paths:
    out_path = out_dir / hdf.name
    print(out_path)
    shutil.copyfile(hdf, out_path)
    
h5_outnames = [i.stem for i in list(out_dir.glob('*.h5'))]
    
#%% check have scan name in metadata csv

metadata_path = Path(r"C:/Users/MerrickS/OneDrive/Work/2_UZH/Papers/1_MEZ_XRF/data/raw/xrf/xrf_scan_metadata.csv")

df_scan_metadata = pd.read_csv(metadata_path)

metadata_scan_names = df_scan_metadata['scan_name'][df_scan_metadata['beamtime'] == 'LS3011'].to_list()

for scan in metadata_scan_names:   
    check_path = out_dir / f'{scan}.h5'
    if check_path.exists():
        pass
    else:
        print(check_path)
        
#%% make hdf_paths links for scan metadata
metadata_path = Path(r"C:/Users/MerrickS/OneDrive/Work/2_UZH/Papers/1_MEZ_XRF/data/raw/xrf/xrf_scan_metadata.csv")

df_scan_metadata = pd.read_csv(metadata_path)

metadata_scan_names = df_scan_metadata['scan_name']


scan_dir = Path(r'C:\Users\MerrickS\OneDrive\Work\2_UZH\Papers\1_MEZ_XRF\data\raw\xrf')

"scans/axo_thin_film_c00/axo_thin_film_c00_0005/axo_thin_film_c00_0005.h5",
"."

hdf_paths = []
for scan in metadata_scan_names:   
    matched_hdfs = list(scan_dir.glob(f'*/*/*/*/*{scan}.h5'))
    if len(matched_hdfs) == 0:
        matched_hdfs = list(scan_dir.glob(f'*/*/*/*{scan}.h5'))
    if len(matched_hdfs) == 0:
        matched_hdfs = list(scan_dir.glob(f'*/*/*{scan}.h5'))
    if len(matched_hdfs) == 0:
        matched_hdfs = list(scan_dir.glob(f'*/*{scan}.h5'))
        
    if len(matched_hdfs) == 0:
        print('no match for', scan)
        
    if len(matched_hdfs) > 1:
        print('multiple matches for', scan)
        
    else:
        print(scan, matched_hdfs[0].name)
        
    # hdf_relative_path = str(matched_hdfs[0]).split(str(scan_dir))[-1].replace('\\', '/')
    hdf_relative_path = matched_hdfs[0]
        
    hdf_paths.append(hdf_relative_path)
    
df_scan_metadata['hdf_filepath'] = hdf_paths

df_scan_metadata.to_csv(metadata_path)




    # scans/GelatinStandards/12_nothing/blank_env_10Hz_0001.h5
    
    # print(scan)
    
    
    

#     if str(scan) in h5_outnames:
#         pass
#     else:
#         missing.append(scan)
        
# df_missing = pd.Series(missing)