import os
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from MassLynxRawDataReader import DataReader

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def to_mspix(fname, width, height):
    mspix_folder = Path(fname).with_suffix(".msPix")
    os.mkdir(mspix_folder) # will fail if folder already exists
    reader = DataReader(fname)

    axis = reader.fast_get_data_scans(reader.scanrange, None)[:,0]
    # np.savetxt(os.path.join(mspix_folder, "axis.csv"), axis, fmt='%.6f', delimiter=',')
    pixel_channels = np.zeros(height * width, dtype=np.uint32)
    total_ion = np.zeros_like(pixel_channels)
    spectrum = np.zeros_like(axis)
    # u32 numbers selected and hard-coded because a single data file is being converted.
    # More general conversion software would allow for other numeric formats depending on data
    # density and needs 
    with open(os.path.join(mspix_folder, "indices.u32"), "ab") as ind_file:
        with open(os.path.join(mspix_folder, "intensities.u32"), "ab") as int_file:
            for i, s in enumerate(reader.scans):
                mzs, ints = reader.readerMS.ReadScan(reader.function, s)
                ints = np.array(ints, dtype=np.uint32)
                mzs = np.array(mzs)
                mzs = mzs[ints > 0]
                ints = ints[ints > 0]
                indices = np.array(np.searchsorted(axis, mzs), dtype=np.uint32)
                assert(indices.size == ints.size)
                indices.tofile(ind_file)
                ints.tofile(int_file)
                pixel_channels[i] = indices.size
                total_ion[i] = np.sum(ints)
                np.add.at(spectrum, indices, ints)
    total_ion.tofile(os.path.join(mspix_folder, "pixel_intensities.u32"))
    pixel_channels.tofile(os.path.join(mspix_folder, "pixel_channels.u32"))
    # total_ion = np.reshape(total_ion, (height, width))
    # pixel_channels = np.reshape(pixel_channels, (height, width))
    # ti_img = Image.fromarray(total_ion)
    # pc_img = Image.fromarray(pixel_channels)
    # ti_img.save(os.path.join(mspix_folder, "total_ion.tiff"), compression="tiff_adobe_deflate")
    # pc_img.save(os.path.join(mspix_folder, "pixel_channels.tiff"), compression="tiff_adobe_deflate")
    
    # np.savetxt(os.path.join(mspix_folder, "spectrum.csv"), spectrum, fmt='%.6f', delimiter=',')

    with open(os.path.join(mspix_folder, "metadata.json"), 'w') as f:
        reader.get_stats()
        # This metadata is for an example - it is likely that not all metadata in the Waters
        # raw file is extracted or preserved here.
        metadata = {
            "msPix_version": "1.0.0",
            "file_name": fname,
            "times": np.around(reader.times, decimals=3).tolist(),
            "image_width_pixels": width,
            "image_height_pixels": height,
            "spectral_channels": list(axis),
            "spectral_intensities": list(spectrum),
        }

            # 'msPix_version': "1.0.0",
            # 'spectral_channels': self.spectral_channels.tolist(),
            # 'spectral_channels_dtype': str(self.spectral_channels.dtype),
            # 'spectral_intensities': self.spectral_intensities.tolist(),
            # 'spectral_intensities_dtype': str(self.spectral_intensities.dtype),
            # 'image_width_pixels': int(self.image_width_pixels),
            # 'image_height_pixels':int(self.image_height_pixels),
            # 'imzML_metadata': p.metadata.pretty(),

        header_items = reader.reader.GetHeaderItems([i for i in range(36)])
        metadata |= {f'Header Item {k}': v for k, v in enumerate(header_items)}
        metadata |= {k: v for k, v in zip(reader.stat_names, reader.stat_vals)}
        f.write(json.dumps(metadata)) 

raw_file = '20240912_Kidney workshop Analyte 1.raw'
width, height = 446, 302 # these are known apriori for the raw_file
to_mspix(raw_file, width, height)



