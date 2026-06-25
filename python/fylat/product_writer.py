"""HDF5 product writer — Python replacement for Fortran io_module.f90 write functions.
Writes CLM, CLA products with correct dataset names, attributes, and data types.
"""

import os
from typing import Optional, Tuple

import h5py
import numpy as np


def write_cloud_mask(output_path: str, cm_bitarray: np.ndarray, qa_bitarray: np.ndarray):
    """Write L2 Cloud Mask product (FY3D_MERSI_ORBT_L2_CLM).

    Args:
        output_path: Full path to output HDF5 file.
        cm_bitarray: Cloud mask bits, shape (nElem, nLine, 6) uint8.
        qa_bitarray: QA bits, shape (nElem, nLine, 10) uint8.
    """
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    nElem, nLine, _ = cm_bitarray.shape

    with h5py.File(output_path, 'w') as f:
        # Cloud_Mask dataset
        ds_cm = f.create_dataset('Cloud_Mask', data=cm_bitarray.transpose(1,2,0).astype(np.uint8),
                                 chunks=(100, 100, 6), compression='gzip',
                                 compression_opts=5)
        ds_cm.attrs['units'] = ''
        ds_cm.attrs['valid_range'] = np.array([0, 255], dtype=np.int32)
        ds_cm.attrs['_FillValue'] = np.int32(0)
        ds_cm.attrs['Intercept'] = np.float32(0.0)
        ds_cm.attrs['Slope'] = np.float32(1.0)
        ds_cm.attrs['long_name'] = 'fylat MERSI_II Cloud Mask'
        ds_cm.attrs['band_name'] = ''

        # Quality_Assurance dataset
        ds_qa = f.create_dataset('Quality_Assurance', data=qa_bitarray.transpose(1,2,0).astype(np.uint8),
                                 chunks=(100, 100, 10), compression='gzip',
                                 compression_opts=5)
        ds_qa.attrs['units'] = ''
        ds_qa.attrs['valid_range'] = np.array([0, 255], dtype=np.int32)
        ds_qa.attrs['_FillValue'] = np.int32(0)
        ds_qa.attrs['Intercept'] = np.float32(0.0)
        ds_qa.attrs['Slope'] = np.float32(1.0)
        ds_qa.attrs['long_name'] = 'fylat MERSI_II Cloud Mask Quality Assurance'
        ds_qa.attrs['band_name'] = ''


def write_cloud_amount(output_path: str, cloud_amount: np.ndarray,
                       cloud_amount_qa: np.ndarray,
                       lon_5km: np.ndarray, lat_5km: np.ndarray):
    """Write L2 Cloud Amount product (5km resolution).

    Args:
        output_path: Full path to output HDF5 file.
        cloud_amount: Cloud fraction 0-100, shape (ix_5km, iy_5km) int16.
        cloud_amount_qa: QA flags, shape (ix_5km, iy_5km) uint8.
        lon_5km: Longitude grid, shape (ix_5km, iy_5km) float32 degrees.
        lat_5km: Latitude grid, shape (ix_5km, iy_5km) float32 degrees.
    """
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)

    with h5py.File(output_path, 'w') as f:
        # Cloud Fraction
        ds_cf = f.create_dataset('5-min granule Cloud Fraction',
                                 data=cloud_amount.astype(np.int16),
                                 chunks=True, compression='gzip', compression_opts=5)
        ds_cf.attrs['units'] = 'none'
        ds_cf.attrs['valid_range'] = np.array([0, 100], dtype=np.int32)
        ds_cf.attrs['_FillValue'] = np.int32(-999)
        ds_cf.attrs['Intercept'] = np.float32(0.0)
        ds_cf.attrs['Slope'] = np.float32(1.0)

        # QA flags
        ds_qa = f.create_dataset('5-min granule Cloud Fraction QA_flags',
                                 data=cloud_amount_qa.astype(np.uint8),
                                 chunks=True, compression='gzip', compression_opts=5)
        ds_qa.attrs['units'] = 'none'
        ds_qa.attrs['valid_range'] = np.array([0, 2], dtype=np.int32)
        ds_qa.attrs['_FillValue'] = np.int32(0)

        # Latitude (stored as scaled int16: degrees * 100)
        lat_scaled = (lat_5km * 100.0).astype(np.int16)
        ds_lat = f.create_dataset('Latitude', data=lat_scaled,
                                  chunks=True, compression='gzip', compression_opts=5)
        ds_lat.attrs['units'] = 'degree'
        ds_lat.attrs['valid_range'] = np.array([-9000, 9000], dtype=np.int32)
        ds_lat.attrs['_FillValue'] = np.int32(65535)
        ds_lat.attrs['Intercept'] = np.float32(0.0)
        ds_lat.attrs['Slope'] = np.float32(0.01)

        # Longitude (stored as scaled int16: degrees * 100)
        lon_scaled = (lon_5km * 100.0).astype(np.int16)
        ds_lon = f.create_dataset('Longitude', data=lon_scaled,
                                  chunks=True, compression='gzip', compression_opts=5)
        ds_lon.attrs['units'] = 'degree'
        ds_lon.attrs['valid_range'] = np.array([-18000, 18000], dtype=np.int32)
        ds_lon.attrs['_FillValue'] = np.int32(65535)
        ds_lon.attrs['Intercept'] = np.float32(0.0)
        ds_lon.attrs['Slope'] = np.float32(0.01)


def write_intermediate(output_path: str,
                       ref_vis: np.ndarray, tbb_ir: np.ndarray,
                       snow_mask: np.ndarray, eco: np.ndarray,
                       precip_water: np.ndarray, sfctmp: np.ndarray,
                       **kwargs):
    """Write intermediate diagnostic HDF5 file (matching Fortran format)."""
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)

    with h5py.File(output_path, 'w') as f:
        f.create_dataset('ref_vis', data=ref_vis.astype(np.float32),
                         chunks=(100, 100, min(19, ref_vis.shape[2])),
                         compression='gzip', compression_opts=5)
        f.create_dataset('tbb_ir', data=tbb_ir.astype(np.float32),
                         chunks=(100, 100, 6), compression='gzip', compression_opts=5)
        f.create_dataset('snow_mask', data=snow_mask.astype(np.int8),
                         chunks=(100, 100), compression='gzip', compression_opts=5)
        f.create_dataset('eco', data=eco.astype(np.int8),
                         chunks=(100, 100), compression='gzip', compression_opts=5)
        if precip_water is not None:
            f.create_dataset('precip_water', data=precip_water.astype(np.float32),
                             chunks=(100, 100), compression='gzip', compression_opts=5)
        if sfctmp is not None:
            f.create_dataset('sfctmp', data=sfctmp.astype(np.float32),
                             chunks=(100, 100), compression='gzip', compression_opts=5)
