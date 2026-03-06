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
    def __init__(self, filename: str):
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
            self.__init__imzml()
        else:
            raise ValueError(f"Error: file {filename} not recognized as valid input file.")
        self._update_internal_state()

    def __init__imzml(self):
        from pyimzml.ImzMLParser import ImzMLParser         
        p = ImzMLParser(self.file)
        self.contents = None
        coordinates = np.array(p.coordinates) - 1
        mzs, ints = p.getspectrum(0)
        mz_dtype, int_dtype = mzs.dtype, ints.dtype
        self.spectral_channels = np.unique(
            np.concatenate([
                mzs[ints > 0]
                for i in range(len(p.coordinates))
                if (mzs := p.getspectrum(i)[0], ints := p.getspectrum(i)[1])[0].size > 0
            ])
        )
        indices_dtype = self._size_to_dtype(self.spectral_channels.size)
        self.image_width_pixels, self.image_height_pixels = np.max(coordinates[:,0]) + 1, np.max(coordinates[:,1]) + 1
        self.pixel_intensities = np.zeros((self.image_height_pixels, self.image_width_pixels), dtype=int_dtype)
        self.pixel_channels = np.zeros((self.image_height_pixels, self.image_width_pixels), dtype=indices_dtype)
        self.intensities, self.indices = [], []
        for i, coords in enumerate(p.coordinates):
            mzs, ints = p.getspectrum(i)
            mzs = mzs[ints > 0] 
            ints = ints[ints > 0]
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
            'mspix_version': "1.0.0",
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


    def save(self, fname, format='loose'):
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
            - ``'loose'``: Save all component files (metadata, indices, intensities,
            and derived data) to a directory.
            - ``'packed'``: Save all components into a single ``.mspix`` file.
            Default is ``'loose'``.

        Notes
        -----
        The *packed* format requires the ``h5py`` library.
        '''
        self.meta["spectral_intensities"] = [round(x, 2) for x in self.meta["spectral_intensities"]]
        if format == 'loose':
            self._save_loose(fname)
        elif format == 'packed':
            self._save_packed(fname)
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

    def _save_packed(self, fname, compression=None , compression_opts=None):
        import h5py
        assert not os.path.exists(fname), f"Error: {fname} already exists. Choose another name."
        with h5py.File(fname + ".mspix", 'w') as file :
            file.create_dataset('meta', data=json.dumps(self.meta, indent=0))
            file.create_dataset("indices", shape=self.indices.shape, dtype=self.indices.dtype, data=self.indices)
            file.create_dataset("intensities", shape=self.intensities.shape, dtype=self.intensities.dtype, data=self.intensities)
            file.create_dataset("pixel_intensities", shape=self.pixel_intensities.shape, dtype=self.pixel_intensities.dtype, data=self.pixel_intensities, compression=compression, compression_opts=compression_opts)
            file.create_dataset("pixel_channels", shape=self.pixel_channels.shape, dtype=self.pixel_channels.dtype, data=self.pixel_channels, compression=compression, compression_opts=compression_opts)

def specification():
    specification = '''
    1. Design philosophy
        mspix files are designed to be an efficient storage and analysis format for mass spectrometry
        imaging data produced using time-of-flight mass spectrometetry, with an emphasis on being
        highly space and access-time efficient when working with large numbers (>5 million) of
        pixels.

        The main goals of this file format may be summarized as:
            A. Store large pixel-count mass images in a space and access-efficient manner
                mspix files are designed to accomodate ~1 billion pixel mass images and are
                recommended for mass images exceeding 5 million pixels, mass images with a single
                mass axis, and/or targeted mass images.
            B. Provide easy inspection and simple, obvious, formats for the loose files.
                The intention behind this decision is to enable ease of tooling and libraries to
                be written by others.
            C. Allow for user-added data in the form of additional files, axes, and metadata
                This is done using a folder/HDF5-based structure.

    2. File or folder structure
        A mspix file may be either a 'loose' collection of files in a folder or a 'packed',
        HDF5-based file with a similar structure. The advantages of the loose representation is that
        it:
            A. Allows for multiple tools present to inspect the files
            B. Enables easy memory mapping larger files
            C. Enables ease of modificaiton to individual composite files
            D. May be easier to implement in some programming languages and environments that do not
                have extensive HDF5 support
        The advantages of the packed representation are that it facilitates file storage and sharing
        by ensuring related data remain together

    3. Loose representation file structure
        The loose file is a folder that contains the following files:
        metadata.json
            A JSON-compliant, UTF-8 encoded text file of all metadata for the image. At this time, 
            only the following entries are strictly required in this file:
                1. a semantically-versioned "mspix_version", which is "1.0.0". This version number
                 allows for backwards-compatible changes to the file-format.
                2. the width of the image in pixels with key "image_width_pixels".
                3. the height of the image in pixels with key "image_height_pixels".
                4. the list of spectral (mass) channels with key "spectral_channels".
                    This list of spectral channels is what the "indices" file indexes into. The
                    largest index in the indices file should be 1 less than the length of this
                    list of spectral channels. These spectral channels will most often be m/z 
                    values, however they may also be times-of-flight or other more abstract
                    values to enable easier or more consistent data calibration and processing.
                        In cases where these values are not m/z values, it is advised to provide
                    a list of m/z values as an additional metadata list. 
                5. the corresponding list of intensities for each spectral channel with key 
                    "spectral_intensities". The reason for the spectral_intensities entry is
                    to enable quick viewing of the overall spectrum of the mspix data to facilitate
                    some common operations such as spectral normalization and peak picking. This
                    data is redundantly stored and could be generated from the "indices" and 
                    "intensities" files.
                
                Entries 2, 3, 4, and 5 consist of numbers or lists of nubmers. All these numbers
                must all be interpretable as 64-bit IEEE floating point numbers. 

            In addition to the "mspix_version", it is advisable to add metadata such as:
                A. the physical dimensions of pixels
                B. the overall image's physical dimensions
                C. instrument settings that were used for the acquisiton (e.g., instrument polarity)
                D. image creation and construction parameters
                E. Checksums of the constituent files
                F. Relevant sample data (e.g., matrix or sample preparation used)
                G. Small annotations
                H. Alternative axis information (e.g., recalibrated mass axes, m/z labels, etc.)

            An example of such a file might be:
                {
                    "mspix_version": "1.0.0",
                    "image_width_pixels": 2412,
                    "image_height_pixels": 2411,
                    "spectral_channels": [1.000, 1.005, 1.010, ... 999.990, 999.995, 1000.000 ],
                    "spectral_intensities": [ 217, 151, 307, .... 58, 30, 176],
                    "original_file": "C:\\Users\\user\\Desktop\\example.tpx3",
                    "width_mm": 2.0025000934838317,
                    "height_mm": 2.001499891281128,
                    "acquisition_time_ns": 173700783473,
                    "tof_pulse_count": 3358385,
                    ...
                }
            with truncations indicated by an ellipsis in this example.

           
        pixel_intensities.xXX
            A little-endian encoded binary file consisting of either 8, 16, 32, or 64-bit numbers
            with the total number of numbers equal to the "image_width_pixels" multiplied by the
            "image_height_pixels". The ".xXX" file extension is a placeholder for indicating the
            numeric type and size of the file. Valid extensions are ".u8", ".u16", ".u32", ".u64",
            ".f32", and ".f64" where "u" indicates an unsigned integer, and "f" indicates an IEEE
            floating point number. The bits of the numbers are indicated with the 8, 16, etc. 

            This file is similar to the "spectral_intensities" entry in the "metadata.json" file in
            that its purpose is to facilitate some common operations such as image normalization and
            easy region-of-interest selection. This data is redundantly stored and could be
            generated from other files. 

        pixel_channels.uXX
            A little-endian encoded binary file consisting of either 8, 16, 32, or 64-bit numbers
            with the total number of numbers equal to the "image_width_pixels" multiplied by the
            "image_height_pixels". The ".uXX" file extension is a placeholder for indicating the
            numeric type and size of the file. Valid extensions are ".u8", ".u16", ".u32", and 
            ".u64" where "u" indicates an unsigned integer. The bits of the numbers are indicated
            with the 8, 16, etc. 
            
            Each "pixel_channels.uXX" pixel's value is the number of nonzero spectral channels in
            that pixel. E.g., if the mspix pixel has 15 nonzero m/z channels then the corresponding
            pixel_channels.uXX value would be 15.
            
            This file facilitates creation of an offset table (via cumulative summation of each
            pixel's channels) for access of the indices and intensities in the
            corresponding "intensities.xXX" and "indices.uXX" files specified below.

        intensities.xXX
            A little-endian encoded binary file of unsigned integers or floating point values. The
            size and type of the binary numbers is specified by the .xXX extension, where 'x' is
            either 'u' for 'unsigned integer' or 'f' for 'floating point value'. 'XX' represents the
            number of bits for the unsigned integer or floating point value and can be either 8, 16,
            32, or 64 if 'x' is 'u' or 32 or 64 if 'x' is 'f'. (e.g., "intensities.u8" and
            "intensities.f64" are valid names)

            Each value in the binary file is the intensity of a channel (i.e., m/z value) within a
            pixel. Which pixel is specified by the offset table created from the
            "pixel_channels.uXX" file. Which channel is specified by the corresponding unsigned
            index value in the "indices.uXX" file. Mass channels with intensities of 0 may be
            included but are recommended to be omitted to save space and, potentially, processing
            time.
            
            The number of intensities in the "intensities.uXX" and indices in the "indices.uXX" file
            must be identical.
        
        indices.uXX
            A little-endian encoded binary file of unsigned integers. The size of the integers
            is specified by the .uXX extension, where XX represents the number of bits for the
            unsigned integer and can be either 8, 16, 32, or 64, depending on the size of the axis
            (i.e., the number of entries in the axis.csv file).
            
            Each unsigned integer is the index of a nonzero intensity mass/time channel within a
            mspix pixel. The corresponding intensity of the indexed mass/time channel is found in
            the "intensities.uXX" file.
        
        The above files may be supplemented with additional (optional) files. These files may be
            convienient to group with the .mspix data for storage and specific program
            implementations. 
            Examples of such optional files are:
            
                A. Another UTF-8 based file that has the same number of entries as the mass axis
                    and time axis, but that contains molecular, fragment, or elemental
                    identifications of particular ions.
                B. A stored set of indices indicating different methods of centroiding that 
                    could be possible to reduce the data.
                C. A corresponding optical light microscopy image of the area imaged

    4. Packed representation file structure
        The packed file representation is an HDF5 file that contains equivalent information to the
        loose files as HDF5 datasets. These files translated to datasets are:
        
        metadata
            A 1x1-sized dataset with a single UTF-8-encoded, JSON-compliant string that is
            equivalent to the data of the loose-representation "metadata.json" file specified in 3.
        
        pixel_intensities
            A 2-dimensional dataset of one of the compatible datatypes corresponding to the 
            "pixel_intensities.xXX" loose-representation file specified in 3. The numeric type
            does not need to be otherwise specified as it is known by the HDF5 dataset. The
            2-dimensional dataset is natural as the pixel_intensities represent the overall image
            data. The width and height of the 2-dimensional dataset must correspond to the specified
            "image_width_pixels" and "image_height_pixels", respectively.
        
        pixel_channels
            A 2-dimensional dataset of one of the compatible datatypes corresponding to the 
            "pixel_channels.uXX" loose-representation file specified in 3. The numeric type
            does not need to be otherwise specified as it is known by the HDF5 dataset. The
            2-dimensional dataset is natural as the pixel_channels represents the overall image
            data. The width and height of the 2-dimensional dataset must correspond to the specified
            "image_width_pixels" and "image_height_pixels", respectively.

        indices
            A 1-dimensional dataset of one of the compatible datatypes corresponding to the
            "indices.uXX" loose-representation file specified in section 3. The numeric type does
            not need to be otherwise specified as it is known by the HDF5 dataset.

        intensities
            A 1-dimensional dataset of one of the compatible datatypes corresponding to the
            "intensities.xXX" loose-representation file specified in section 3. The numeric type
            does not need to be otherwise specified as it is known by the HDF5 dataset.
        
    5. Image orientation and shape
        All mspix images start with pixel at index 0, row: 0, column: 0 where row 0 is at the top
        of the image and column 0 is at the left of the image. The image then proceeds in the
        "x" direction (along the row 0, incrementing column; row-major) until it is at the right
        border, the next pixel is then at the far-left in row 1.

        The layout of mspix images is constant and independant of the stage or beam motion during
        acquisition (e.g. serpentine, typewriter, up-down, or left-right). If information
        regarding the stage motion during acquisition is desired, it may be provided in the
        metadata.json file or as a separate file within the folder/zip file.

    6. Limitations and restrictions
        The mspix file format takes advantage of some restrictions, limitations, and constraints 
        that are common in current mass spectrometry imaging experiments. These constraints do
        not apply to many mass spectrometry images, which means that mspix is not suitable for
        losslessly representing such images. Some considerations are that:
        
            A. The data are "sparse" and do not contain digitizer noise. This means that the
                "processed" and "continuous" options of, for exmaple, the .imzML file format do not
                apply as all mspix files are "processed". Similarly, "event-based" data that may be
                encountered with some TOF-SIMS data must be converted to "processed" data for mspix
                conversion.
            B. mspix images are rectangular. It is possible to store non-rectangular imaged areas by
                "cropping" or deleting pixels outside of a defined area. Represenging these empty
                "corner" pixels can be very space efficient (taking only a few bytes each).
                Images with pixels that are not sampled on a rectangular grid cannot be represented
                without re-binning the pixels (and thus some data loss) using mspix.
            C. The data must be represented with a single 1-dimensional axis. This excludes some
                multi-tandem-MS and ion mobility data as some data acquired with these modes are not
                able to be represented efficiently with a single axis. MS2 data with a common second
                axis may be able to be represented by extending the primary mass axis or providing
                additional metadata that tag specific pixels or axis channels with fragmentation
                information; however such representations are left unspecified.
            D. Only a single, 2-dimensional image will be contained in the mspix file. This means
                that 3-d images or depth-profiling images are not supported in the 1.0.0 version
                of mspix.
    '''
    print(specification)
    return
