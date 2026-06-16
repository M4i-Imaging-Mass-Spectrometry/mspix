# -*- coding: utf-8 -*-

# Copyright 2025 Ian Anthony and Kasper Krijnen
#
# Licensed under the MIT license found at http://opensource.org/licenses/MIT .
# May not be used, copied, or modified except in compliance with the MIT license.

import os, json
import numpy as np
from os.path import join

class Mspix:
    '''
    Generate an Mspix object for reading and writing .mspix files.

    It is not recommended to modify fields of the Mspix object directly.
    The following exceptions apply:
    
    - Adding entries to the ``mspix.meta`` dictionary.

    Notes
    -----
    Any modifications are not saved to disk unless the ``.save`` method is called.

    '''
    def __init__(self, filename: str, mz_decimals=None, threshold=0.0):
        '''
        Initialize the object using different methods depending on the input path.

        The initialization behavior depends on the type of ``filename`` provided:

        1. **Directory** containing the seven files that constitute an .mspix dataset  
            - The ``intensities`` and ``indices`` files are memory-mapped and accessed with
                numpy.memmap().  
            - All other files are loaded fully into RAM.
        2. **.mspix HDF5 file**  
            - Requires the additional ``h5py`` library.
            - The ``intensities`` and ``indices`` files are accessed as HDF5 hyperslabs and not
                loaded into RAM.
            - All other files are loaded fully into RAM.
        3. **.imzML file** with a corresponding .ibd file in the same directory  
            - Requires the additional ``pyimzml`` library.
            - Everything is loaded into RAM
  
        Parameters
        ----------
        filename : str
            Path to the .imzML file, .mspix directory, or .mspix HDF5 file.
        mz_decimals : int, optional
            Decimal precision used to round m/z values when loading data from .imzML files.
            - If None (recommended for continuous-profile or continuous-centroid data), m/z values are left unchanged.
            - If an integer, m/z values are rounded to that number of decimals.
            This can be used as a workaround to align unaligned m/z bins of different spectra in processed-profile or processed-centroid data.
        threshold : float, optional
            Intensity threshold cutoff used when loading data from .imzML files. Values less than or equal to this threshold are excluded. 
            Defaults to 0.0, excluding zero-intensity values.

        Returns
        -------
        self : Mspix
            The initialized object instance.

        The following internal variables are created during initialization:
        ``self.file`` : str
            The initialization parameter ``filename``.
        ``self.image_width_pixels`` : int  and ``self.image_height_pixels`` : int  
            Number of columns and rows, respectively.
        ``self.size`` : int
            Total number of pixels (``width * height``).
        ``self.contents`` : list or None  
            Either a list of files in the dataset directory or ``None`` for .imzML and .mspix inputs.
        ``self.spectral_channels`` : ndarray 64-bit floats
            Common mass (m/z) or time axis for all pixels.
        ``self.spectral_intensities`` : ndarray  64-bit floats
            Total summed intensity spectrum across all pixels.
        ``self.indices`` : ndarray
            Indices into the spectral axis.
        ``self.intensities`` : ndarray
            Intensities corresponding to ``self.indices``.
        ``self.pixel_intensities`` : ndarray (2D)  
            Total ion current (TIC) per pixel.
        ``self.pixel_channels`` : ndarray (2D)  
            Number of spectral entries per pixel.
        ``self.offsets`` : ndarray  
            Table of pixel offsets into the ``indices`` and ``intensities`` files.

        Notes
        Any modifications are not saved to disk unless the ``.save`` method is called.
        - ``pyimzML`` must be installed to use this function on a .imzML file.
        '''         

        self.file = filename
        if os.path.isdir(self.file):
            self.__init__mspix_loose()
        elif self.file.lower().endswith(".mspix"):
            self.__init__mspix_packed()
        elif self.file.lower().endswith(".imzml"):
            self.__init__imzml(mz_decimals, threshold)
        else:
            raise ValueError(f"Error: file {filename} not recognized as valid input file.")
        self._update_internal_state()

    def __init__imzml(self, mz_decimals, threshold):
        from pyimzml.ImzMLParser import ImzMLParser         
        p = ImzMLParser(self.file)
        self.contents = None
        coordinates = np.array(p.coordinates) - 1
        
        mzs, ints = p.getspectrum(0)
        mz_dtype, int_dtype = mzs.dtype, ints.dtype

        self.spectral_channels = np.unique(
            np.concatenate([
                mzs[ints > threshold]
                for i in range(len(p.coordinates))
                if (mzs := p.getspectrum(i)[0], ints := p.getspectrum(i)[1])[0].size > 0
            ])
        )

        if mz_decimals != None:
            self.spectral_channels = np.round(self.spectral_channels, decimals=mz_decimals)
        
        self.spectral_channels = np.sort(np.unique(self.spectral_channels))
        indices_dtype = self._size_to_dtype(self.spectral_channels.size)
        self.image_width_pixels, self.image_height_pixels = np.max(coordinates[:,0]) + 1, np.max(coordinates[:,1]) + 1
        self.pixel_intensities = np.zeros((self.image_height_pixels, self.image_width_pixels), dtype=int_dtype)
        self.pixel_channels = np.zeros((self.image_height_pixels, self.image_width_pixels), dtype=indices_dtype)
        self.intensities, self.indices = [], []
        for i, coords in enumerate(p.coordinates):
            mzs, ints = p.getspectrum(i)
            mzs = mzs[ints > threshold]
            if mz_decimals != None:
                mzs = np.round(mzs, decimals=mz_decimals)   
            ints = ints[ints > threshold]
            
            if len(ints) == 0:
                self.pixel_channels[coords[1] - 1, coords[0] - 1] = 1
                self.pixel_intensities[coords[1] - 1, coords[0] - 1] = 0
                self.intensities.append(0)
                self.indices.append(0)
                continue
            
            indices = np.searchsorted(self.spectral_channels, mzs)   
            unique_indices = np.unique(indices)
            final_indices = np.arange(unique_indices.size)[np.searchsorted(unique_indices, indices)] 
            new_ints = self._index_add(unique_indices, final_indices, ints)
            self.intensities.append(new_ints)
            self.indices.append(unique_indices)
            self.pixel_intensities[coords[1] - 1, coords[0] - 1] = np.sum(new_ints)
            self.pixel_channels[coords[1] - 1, coords[0] - 1] = new_ints.size
        
        self.intensities = np.array(np.concatenate(self.intensities), dtype=int_dtype)
        self.indices = np.array(np.concatenate(self.indices), dtype=indices_dtype)
        self.spectral_intensities = self._index_add(self.spectral_channels, self.indices, self.intensities)
        self.meta = {
            'mspix_version': "1.0.5",
            'spectral_channels': self.spectral_channels.tolist(),
            'spectral_intensities': self.spectral_intensities.tolist(),
            'image_width_pixels': int(self.image_width_pixels),
            'image_height_pixels':int(self.image_height_pixels),
            'imzML_metadata': p.metadata.pretty(),
        }

    def __init__mspix_loose(self):
        self.contents = os.listdir(self.file)
        missing_files = [
            f for f in ["metadata", "pixel_intensities", "pixel_channels", "indices", "intensities"]
            if f not in [os.path.splitext(x)[0] for x in self.contents]
        ]
        if missing_files:
            raise FileNotFoundError(f"Missing required file(s): {', '.join(missing_files)}")
        with open(join(self.file, "metadata.json"), 'r') as f:
            self.meta = json.load(f)
        self.spectral_channels = np.array(self.meta['spectral_channels'], dtype=np.float64)
        self.spectral_intensities = np.array(self.meta['spectral_intensities'], dtype=np.float64)
        height, width = self.meta['image_height_pixels'], self.meta['image_width_pixels']
        self.pixel_intensities = np.fromfile(join(self.file, self._get_fname('pixel_intensities')), dtype=self._get_fsize('pixel_intensities')).reshape(height, width)
        self.pixel_channels = np.fromfile(join(self.file, self._get_fname('pixel_channels')), dtype=self._get_fsize('pixel_channels')).reshape(height, width)
        self.indices = self._mmap_content('indices')
        self.intensities = self._mmap_content('intensities')

    def __init__mspix_packed(self):
        import h5py
        self.contents = None
        self.hdf5file = h5py.File(self.file, 'r')
        self.meta = json.loads(self.hdf5file['meta'][()].decode('utf-8'))
        self.spectral_channels = np.array(self.meta['spectral_channels'], dtype=np.float64)
        self.spectral_intensities = np.array(self.meta['spectral_intensities'], dtype=np.float64)
        self.pixel_intensities = np.array(self.hdf5file['pixel_intensities'][()])
        self.pixel_channels = np.array(self.hdf5file['pixel_channels'][()])
        self.indices = self.hdf5file['indices']
        self.intensities = self.hdf5file['intensities']
    
    def _update_internal_state(self):
        self.image_height_pixels = self.pixel_channels.shape[0]
        self.image_width_pixels = self.pixel_channels.shape[1]
        self.size = self.pixel_channels.size
        self.offsets = np.roll(np.append(np.cumsum(self.pixel_channels.flatten(), dtype=np.uint64), 0), 1) 
        
    def _mmap_content(self, f):
        return np.memmap(join(self.file, self._get_fname(f)), dtype=self._get_fsize(f), mode='r')

    def _get_fname(self, file_name):
        return [f for f in self.contents if f.startswith(file_name)][0]

    def _get_fsize(self, file_name):
        size = f"{os.path.splitext(self._get_fname(file_name))[-1][1:]}"
        if (size.lower().startswith('u') or size.lower().startswith('i')) and 'int' not in size:
            size = size.replace('i', 'int').replace('u', 'uint')
        return  size
    
    @staticmethod
    def _size_to_dtype(maximum):
        return np.uint8 if maximum < 2**8 else np.uint16 if maximum < 2**16 else np.uint32

    @staticmethod
    def _index_add(array, indices, intensities):
        new = np.zeros_like(array, dtype = intensities.dtype) 
        np.add.at(new, indices, intensities)
        return new
        
    def get_pixels(self, pixel_indices):
        '''
        Takes an array of pixel indices or a single pixel index and returns the summed set of index
        and intensity values.

        Can be used to obtain the integrated mass spectrum of a selected area
        within the mass image.

        Parameters
        ----------
        pixel_indices : int or array-like of int
            Single pixel index or an array of pixel indices to integrate.

        Returns
        -------
        indices : ndarray
            Index values into ``self.spectral_channels`` for the selected pixels.
        intensities : ndarray
            Summed intensity values for the selected pixels.
            Each intensity corresponds to the spectral_channel indicated by the matching element in
            ``indices``.
        '''
        starts, ends = self.offsets[pixel_indices], self.offsets[pixel_indices + 1]
        if type(pixel_indices) == int or type(pixel_indices) == float or pixel_indices.size == 1:
            selected = np.arange(starts, ends, dtype=np.uint64)
        else:
            selected = np.concatenate([np.arange(s, e, dtype=np.uint64) for s, e in zip(starts, ends)])
        if selected.size == 0:
            return np.array([0]), np.array([0])
        indices, intensities = self.indices[selected], self.intensities[selected]
        integrated_indices = np.array(np.sort(np.unique(indices)))
        indices = np.arange(integrated_indices.size)[np.searchsorted(integrated_indices, indices)]
        return integrated_indices, self._index_add(integrated_indices, indices, intensities)
          
    def ion_image(self, low_index, high_index, batch_size=None):
        '''
        Construct an ion image by accumulating intensities within a specified index range.

        The function uses chunked processing to avoid loading all data into memory. Only intensities
        with indices within ``[low_index, high_index]`` are accumulated. Chunked processing is not
        supported for mspix objects initiated from .imzML files. 

        Parameters
        ----------
        low_index : int
            Lower bound of the index range (inclusive) into ``self.spectral_channels``.
        high_index : int
            Upper bound of the index range (inclusive) into ``self.spectral_channels``.
        batch_size : int, optional
            Number of rows to process per chunck. If ``None``, defaults to 32.
            Set to 0 to disable chunked processing and load all data into memory.

        Returns
        -------
        image : ndarray
            2D array where each element represents the summed intensity of all
            spectral channels within ``[low_index, high_index]`` for the corresponding pixel.
           
        '''
        if batch_size is None:
            batch_size = 32
        if batch_size == 0:
            self.int_to_pix = np.repeat(np.arange(self.size), self.pixel_channels.flatten())
            mask_ints = np.where((self.indices >= low_index) & (self.indices <= high_index), self.intensities, 0) 
            return np.reshape(self._index_add(self.offsets[:-1], self.int_to_pix, mask_ints),self.pixel_intensities.shape)
        
        if batch_size > 0:
            height, width = self.image_height_pixels, self.image_width_pixels
            image = np.zeros((height, width), dtype=self.pixel_intensities.dtype)
            for i in range(0, height, batch_size):
                start_row = i
                rows_in_batch = min(batch_size, height - i)
                row_offset = i * width
                start = int(self.offsets[row_offset])
                end = int(self.offsets[row_offset + rows_in_batch * width])
                length = end - start
                if length == 0:
                    continue
                if isinstance(self.intensities, np.memmap):
                    indices_batch = np.memmap(
                        self.indices.filename, dtype=self.indices.dtype, mode='r',
                        offset=self.indices.dtype.itemsize * start, shape=(length,)
                    )
                    intensities_batch = np.memmap(
                        self.intensities.filename, dtype=self.intensities.dtype, mode='r',
                        offset=self.intensities.dtype.itemsize * start, shape=(length,)
                    )
                elif hasattr(self.intensities, 'file') and 'h5py' in str(type(self.intensities)):
                    indices_batch = self.indices[start:end]
                    intensities_batch = self.intensities[start:end]
                masked_ints = np.where((indices_batch >= low_index) & (indices_batch <= high_index), intensities_batch, 0)
                pixel_channels_batch = self.pixel_channels[i : i + rows_in_batch, :].reshape(-1)
                int_to_pix = np.repeat(np.arange(rows_in_batch * width), pixel_channels_batch)
                offsets_batch = self.offsets[row_offset : row_offset + rows_in_batch * width]
                flat_result = self._index_add(offsets_batch, int_to_pix, masked_ints)
                image[i : i + rows_in_batch, :] = flat_result.reshape(rows_in_batch, width)
                del indices_batch
                del intensities_batch
            return image

    def ion_image_spectral_channel(self, spectral_channel, tolerance=0.5, batch_size=None):
        '''
        Generate an ion image for a specified spectral channel within a given tolerance.

        This function wraps ``ion_image`` to produce an ion image centered on a target
        spectral channel. The tolerance defines the half-width of the window around
        ``spectral_channel``. For example, a tolerance of 0.5 creates a window one
        channel wide centered on the target value.

        Parameters
        ----------
        spectral_channel : float
            Target spectral channel value (e.g., m/z or time-of-flight).
        tolerance : float, optional
            Half-width of the spectral window centered on ``spectral_channel``.
            Default is 0.5, corresponding to a total width of 1 channel.
        batch_size : int, optional
            Number of pixel chunks to process at once. If ``None``, defaults to 32.
            Passed directly to ``ion_image``.

        Returns
        -------
        image : ndarray
            2D array where each element represents the summed intensity of all
            spectral channels within the range
            ``[spectral_channel - tolerance, spectral_channel + tolerance]``
            for the corresponding pixel.

        Notes
        -----
        ``spectral_channel`` values typically represent m/z, but may represent
        times-of-flight or other quantities depending on the file type.
        '''
        return self.ion_image(
            self.spectral_channel_to_index(spectral_channel - tolerance),
            self.spectral_channel_to_index(spectral_channel + tolerance),
            batch_size=batch_size
        )

    def spectral_channel_to_index(self, spectral_channel):
        '''Converts a m/z value to the nearest index.'''
        return np.searchsorted(self.spectral_channels, spectral_channel)

    def index_to_spectral_channel(self, index):
        '''Converts a mass-axis index to the corresponding m/z value'''
        return self.spectral_channels[index]
    
    def pixel_index_to_row_col(self, pixel_index):
        '''Converts a pixel index to the corresponding row and column'''
        return pixel_index // self.image_width_pixels, pixel_index % self.image_width_pixels

    def row_col_to_pixel_index(self, row, col):
        '''Converts a row and column to the corresponding pixel index'''
        return row * self.image_width_pixels + col

    def save_imzml(self, fname, mode='processed', spec_type='profile'):
        '''
        Save the Mspix dataset as an .imzML file and corresponding .ibd file.
        The method exports the current Mspix data to the specified location in either *processed* or
         *continuous* mode and spectrum type *profile* or *centroid*, using the ``pyimzML`` library.

        Parameters
        ----------
        fname : str
            Output filename (without extension) for the .imzML and .ibd files.
        mode : {'processed', 'continuous'}, optional
            ``'processed'``: Each pixel stores its own mass (m/z) axis and
            corresponding intensity values.
            ``'continuous'``: All pixels share a common mass (m/z) axis, and each pixel stores
             intensity values for every position along that axis.
            Default is ``'processed'``.
        spec_type : {'profile', 'centroid'}, optional
            Type of spectrum to record in the output file. Default is ``'profile'``.

        Notes
        -----
        - ``pyimzML`` must be installed to use this function.
        '''
        from pyimzml.ImzMLWriter import ImzMLWriter
        
        mz_dtype = np.dtype(self.spectral_channels.dtype).type
        intensity_dtype = np.dtype(self.spectral_intensities.dtype).type

        if mode == 'processed':
            with ImzMLWriter(fname, mz_dtype=mz_dtype, intensity_dtype=intensity_dtype, mode=mode, spec_type=spec_type) as f:
                for i in range(self.pixel_intensities.size):
                    indices, intensities = self.get_pixels(np.array(i))
                    mzs = np.array(self.spectral_channels[indices], dtype=mz_dtype)
                    row, col = self.pixel_index_to_row_col(i)
                    if mzs.size == 0:
                        mzs = np.zeros(1, dtype=mz_dtype)
                        intensities = np.zeros_like(mzs)            
                    f.addSpectrum(mzs, np.array(intensities, dtype=intensity_dtype), (col + 1, row + 1))

        if mode == 'continuous':
            with ImzMLWriter(fname, mz_dtype, intensity_dtype, mode=mode, spec_type=spec_type) as f:
                mzs = self.spectral_channels
                for i in range(self.pixel_intensities.size):
                    intensities = np.zeros(len(mzs), dtype=self.pixel_intensities.dtype)
                    indices, pix_int = self.get_pixels(np.array(i))
                    intensities[indices] = pix_int
                    row, col = self.pixel_index_to_row_col(i)
                    f.addSpectrum(mzs, np.array(intensities, dtype=intensity_dtype), (col + 1, row + 1))


    def save(self, fname, format='loose', compression=None, compression_opts=None):
        '''
        Save the Mspix dataset to disk in either loose or packed format.

        The method stores the current Mspix data to the specified location, either
        as a directory containing individual component files (*loose* format) or as
        a single ``.mspix`` (HDF5) file (*packed* format).

        Parameters
        ----------
        fname : str
            Output path for saving the dataset. Can be a directory name (for loose
            format) or a filename ending in ``.mspix`` (for packed format).

        format : {'loose', 'packed'}, optional
            Storage format:
            
            ``'loose'``: Save all component files (metadata, indices, intensities,
            and derived data) to a directory.

            ``'packed'``: Save all components into a single ``.mspix`` file.

            Default is ``'loose'``.

        compression : {None, 'gzip', 'lzf'}, optional (only for packed format)
            
            ``'gzip'``: Uses the lossless gzip compression filter available in HDF5 for the 
            packed format. Good compression moderate speed.
            
            ``'lzf'``: Uses the lossless lzf compression filter available in HDF5 for the packed
            format. Low to moderate compression, very fast.
        
        compression_opts : int
            Sets the compression level for gzip with an integer between 0-9, default is 4.  
        
        Notes
        -----
        The *packed* format requires the ``h5py`` library.
        '''
        self.meta["spectral_intensities"] = [round(x, 2) for x in self.meta["spectral_intensities"]]
        if format == 'loose':
            self._save_loose(fname)
        elif format == 'packed':
            self._save_packed(fname, compression, compression_opts)
        else:
            raise ValueError(f"Unknown format '{format}'. Expected  'loose' or 'packed'.")

    def _save_loose(self, fname): 
        assert not os.path.exists(fname), f"Error: {fname} already exists. Choose another name."
        os.mkdir(fname)
        self.indices.tofile(join(fname, f"indices.{self.indices.dtype}"))
        self.intensities.tofile(join(fname, f"intensities.{self.intensities.dtype}"))
        self.pixel_channels.tofile(join(fname, f'pixel_channels.{self.pixel_channels.dtype}'))
        self.pixel_intensities.tofile(join(fname, f'pixel_intensities.{self.pixel_intensities.dtype}'))
        with open(join(fname, "metadata.json"), 'w') as f:
            json.dump(self.meta, f, indent=0)

    def _save_packed(self, fname, compression , compression_opts):
        import h5py
        assert not os.path.exists(fname), f"Error: {fname} already exists. Choose another name."
        with h5py.File(fname + ".mspix", 'w') as file :
            file.create_dataset('meta', data=json.dumps(self.meta, indent=0))
            file.create_dataset("indices", shape=self.indices.shape, dtype=self.indices.dtype, data=self.indices, 
                compression=compression, compression_opts=compression_opts)
            file.create_dataset("intensities", shape=self.intensities.shape, dtype=self.intensities.dtype, data=self.intensities,
                compression=compression, compression_opts=compression_opts)
            file.create_dataset("pixel_intensities", shape=self.pixel_intensities.shape, dtype=self.pixel_intensities.dtype,
                data=self.pixel_intensities)
            file.create_dataset("pixel_channels", shape=self.pixel_channels.shape, dtype=self.pixel_channels.dtype,
                data=self.pixel_channels)
