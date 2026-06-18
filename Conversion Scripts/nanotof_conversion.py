


# -*- coding: utf-8 -*-

# Copyright 2023 Ian Anthony
#
# Licensed under the MIT license found at http://opensource.org/licenses/MIT .
# May not be used, copied, or modified except in compliance with the MIT license.

import os, json, PIL, zipfile, shutil, struct
import matplotlib.pyplot as plt
import numpy as np
from os.path import join

class msPix:
    def __init__(self, filename: str):
        self.file = filename
        self.__init__nanotof()
        self._update_internal_state()

    def __init__nanotof(self):
        self.meta, fsize = {}, os.path.getsize(self.file)
        with open((self.file), 'rb') as f:
            self._add_nanotof_header_to_meta(f)
            f.seek(int(self.meta['HeaderSize']), 0)
            pix = int(self.meta['ImagePixels'])
            self.rows = int( self.meta['Number Of Tiles Y'] * pix )
            self.cols = int( self.meta['Number Of Tiles X'] * pix )
            mass_data = [{} for i in range(self.rows * self.cols)]
            start, stop = self.meta['StartFlightTime'] * 1000, self.meta['StopFlightTime'] * 1000
            increment = self.meta['SpecBinIncr'] * 1000
            self.time_axis = np.linspace(start, stop, int((stop-start)/increment) + 1)
            mode = 1 if 'Inactive' in self.meta['MSMS Active'] else 8
            xtile, ytile, count = 0, 0, 0
            while True:
                if count % 1000 == 0:
                    print(f'{f.tell() / fsize * 100:.3f}% of file read')
                count += 1
                try:
                    block_id = struct.unpack('H', f.read(2))[0]
                    f.seek(10, 1)
                    record_size = struct.unpack('L', f.read(4))[0]
                    if block_id == mode:
                        data = np.frombuffer(f.read(record_size), dtype=np.uint64)
                        xs, ys, inds = self._parse_nanotof_data(data, xtile, ytile)
                        xs = xs[inds <= stop]
                        ys = ys[inds <= stop]
                        inds = inds[inds <= stop]
                        for x, y, i in zip(xs, ys, inds):
                            mass_data[(y * self.cols) + x][i] = mass_data[(y * self.cols) + x].get(i, 0) + 1
                    elif block_id == 6:
                        xtile = struct.unpack('L', f.read(4))[0] * pix
                        ytile = struct.unpack('L', f.read(4))[0] * pix
                    else:
                        f.seek(int(record_size), 1)
                except (EOFError, struct.error) as e:
                    print('Finished reading file')
                    self.sizes = np.zeros((self.rows, self.cols), dtype=np.uint16)
                    self.total = np.zeros((self.rows, self.cols), dtype=np.uint16)
                    self.mass_axis = self.time_axis 
                    indices, ints = [], []
                    for i, pixel in enumerate(mass_data):
                        if i % 1000 == 0:
                            print(f"processed {i} of {len(mass_data)} pixels")
                        row, col = self.pixel_index_to_row_col(i)
                        keys, values = zip(*pixel.items())
                        mass_data[i] = None
                        indices.extend(keys)
                        ints.extend(values)
                        self.total[row, col] = sum(values)
                        self.sizes[row, col] = len(keys)
                    self.indices = np.array(indices, dtype=self._size_to_dtype(np.max(indices)))
                    self.intensities = np.array(ints, dtype=self._size_to_dtype(np.max(ints)))
                    self.spectrum = np.zeros_like(self.time_axis)
                    np.add.at(self.spectrum, self.indices, self.intensities)
                    
                    self.meta['msPix_version'] = "1.0.0"
                    self.meta['spectral_channels'] = self.mass_axis.tolist()
                    self.meta['spectral_channels_dtype'] = str(self.mass_axis.dtype)
                    self.meta['spectral_intensities'] = self.spectrum.tolist()
                    self.meta['spectral_intensities_dtype'] = str(self.spectrum.dtype)
                    self.meta['image_width_pixels'] = int(self.rows)
                    self.meta['image_height_pixels'] = int(self.cols)
                    print('Finished parsing nanotof file')
                    return

    @staticmethod
    def _parse_nanotof_data(data, xtile, ytile):
        if data.size <= 0:
            return xtile, ytile, 0
        channels, pixind = data & 0x7FFFFFFF, (data >> 58) == 0
        channels = np.array(channels[pixind] / 1000, dtype=int)
        xs = np.array((data & 0x00_07FF_0000_0000) >> 32, dtype=int)[pixind] + xtile # 11 bits
        ys = np.array((data & 0x3F_F800_0000_0000) >> 43, dtype=int)[pixind] + ytile # 11 bits
        return xs, ys, channels       

    def _add_nanotof_header_to_meta(self, f):
        while (line := f.readline().strip()) != b'EOFH':
            if b'=' not in line and b':' not in line:
                continue
            a, b = line.split(b'=', 1) if b'=' in line else line.split(b':', 1)
            a, b = a.decode(), b.decode()
            try:
                b = [float(c.strip(")").strip("(")) for c in b.split(',')] if ',' in b else float(b)
            except Exception:
                b = b.strip().strip(')'.strip('(').strip('[').strip(']'))
            self.meta[a] = b

    def _update_internal_state(self):
        self.rows, self.cols, self.size = self.sizes.shape[0], self.sizes.shape[1], self.sizes.size
        self.offsets = np.roll(np.append(np.cumsum(self.sizes.flatten()), 0), 1)
        self.int_to_pix = np.repeat(np.arange(self.size), self.sizes.flatten())

    @staticmethod
    def _size_to_dtype(maximum):
        return np.uint8 if maximum < 2**8 else np.uint16 if maximum < 2**16 else np.uint32

    @staticmethod
    def _index_add(array, indices, intensities):
        new = np.zeros_like(array)
        np.add.at(new, indices, intensities)
        return new
 
    def pixel_index_to_row_col(self, pixel_index):
        '''Converts a pixel index to the corresponding row and column'''
        return pixel_index // self.cols, pixel_index % self.cols

    def _save_loose(self, fname): 
        '''Saves mspix data to fname as a zipped file (zipped=True) or a folder (zipped=False)'''
        assert not os.path.exists(fname), f"Error: {fname} already exists. Choose another name."
        os.mkdir(fname)
        self.indices.tofile(join(fname, f"indices.{self.indices.dtype}"))
        self.intensities.tofile(join(fname, f"intensities.{self.intensities.dtype}"))
        self.sizes.tofile(join(fname, f'pixel_channels.{self.sizes.dtype}'))
        self.total.tofile(join(fname, f'pixel_intensities.{self.total.dtype}'))
        with open(join(fname, "metadata.json"), 'w') as f:
            json.dump(self.meta, f, indent=0)


os.chdir(os.path.dirname(os.path.abspath(__file__)))
mspix = msPix("example_nanotof.raw")
mspix._save_loose("example_nanotof.mspix")
