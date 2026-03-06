# mspix File Format 
#
**mspix** files are designed to be an efficient storage and analysis format for mass spectrometry imaging data produced using time-of-flight mass spectrometetry, with an emphasis on being highly space and access-time efficient when working with large numbers (>5 million) of pixels.


### The main goals of this file format may be summarized as:
- Store large pixel-count mass images in a space and access-efficient manner mspix files are designed to accomodate ~1 billion pixel mass images and are recommended for mass images exceeding 5 million pixels, mass images with a single mass axis, and/or targeted mass images.


- Provide easy inspection and simple, obvious, formats for the loose files.
                The intention behind this decision is to enable ease of tooling and libraries to
                be written by others.
                
- Allow for user-added data in the form of additional files, axes, and metadata
                This is done using a folder/HDF5-based structure.

---

## File or Folder Structure

A `.mspix` file may be either a 'loose' collection of files in a folder or a 'packed',
        HDF5-based file with a similar structure. The advantages of the loose representation is that
        it:
- Allows for multiple tools present to inspect the files
- Enables easy memory mapping larger files
- Enables ease of modificaiton to individual composite files
- May be easier to implement in some programming languages and environments that do not have extensive HDF5 support

The advantages of the packed representation are that it facilitates file storage and sharing by ensuring related data remain together

---

## Loose Representation File Structure

The loose file is a folder that contains the following files:

### `metadata.json`

A JSON-compliant, UTF-8 encoded text file of all metadata for the image. At this time, only the following entries are strictly required in this file:

1. a semantically-versioned `"mspix_version"`, which is `"1.0.0"`. This version number allows for backwards-compatible changes to the file-format.  
2. the width of the image in pixels with key `"image_width_pixels"`.  
3. the height of the image in pixels with key `"image_height_pixels"`.  
4. the list of spectral (mass) channels with key `"spectral_channels"`.  
   This list of spectral channels is what the `"indices"` file indexes into. The largest index in the indices file should be 1 less than the length of this list of spectral channels. These spectral channels will most often be m/z values, however they may also be times-of-flight or other more abstract values to enable easier or more consistent data calibration and processing.  
   In cases where these values are not m/z values, it is advised to provide a list of m/z values as an additional metadata list.  
5. the corresponding list of intensities for each spectral channel with key `"spectral_intensities"`.  
   The reason for the `spectral_intensities` entry is to enable quick viewing of the overall spectrum of the mspix data to facilitate some common operations such as spectral normalization and peak picking. This data is redundantly stored and could be generated from the `"indices"` and `"intensities"` files.

Entries 2, 3, 4, and 5 consist of numbers or lists of numbers. All these numbers must all be interpretable as 64-bit IEEE floating point numbers.

In addition to the `"mspix_version"`, it is advisable to add metadata such as:

A. the physical dimensions of pixels  
B. the overall image's physical dimensions  
C. instrument settings that were used for the acquisition (e.g., instrument polarity)  
D. image creation and construction parameters  
E. Checksums of the constituent files  
F. Relevant sample data (e.g., matrix or sample preparation used)  
G. Small annotations  
H. Alternative axis information (e.g., recalibrated mass axes, m/z labels, etc.)

An example of such a file might be:

```text
{
    "mspix_version": "1.0.0",
    "image_width_pixels": 2412,
    "image_height_pixels": 2411,
    "spectral_channels": [1.000, 1.005, 1.010, ... 999.990, 999.995, 1000.000 ],
    "spectral_intensities": [ 217, 151, 307, ... 58, 30, 176],
    "original_file": "C:\\Users\\user\\Desktop\\example.tpx3",
    "width_mm": 2.0025000934838317,
    "height_mm": 2.001499891281128,
    "acquisition_time_ns": 173700783473,
    "tof_pulse_count": 3358385,
    ...,
}
```

with truncations indicated by an ellipsis in this example.

### `pixel_intensities.xXX`

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

### `pixel_channels.uXX`

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

### `intensities.xXX`
A little-endian encoded binary file of unsigned integers or floating point values. The size and type of the binary numbers is specified by the .xXX extension, where 'x' is
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

### `indices.uXX`
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

- Another UTF-8 based file that has the same number of entries as the mass axis and time axis, but that contains molecular, fragment, or elemental identifications of particular ions.
- A stored set of indices indicating different methods of centroiding that could be possible to reduce the data.
- A corresponding optical light microscopy image of the area imaged

## Packed representation file structure
The packed file representation is an HDF5 file that contains equivalent information to the
loose files as HDF5 datasets. These files translated to datasets are:

### `metadata`
A 1x1-sized dataset with a single UTF-8-encoded, JSON-compliant string that is equivalent to the data of the loose-representation "metadata.json" file specified in 3.

### `pixel_intensities`
A 2-dimensional dataset of one of the compatible datatypes corresponding to the 
            "pixel_intensities.xXX" loose-representation file specified in 3. The numeric type
            does not need to be otherwise specified as it is known by the HDF5 dataset. The
            2-dimensional dataset is natural as the pixel_intensities represent the overall image
            data. The width and height of the 2-dimensional dataset must correspond to the specified
            "image_width_pixels" and "image_height_pixels", respectively.
            
### `pixel_channels`
A 2-dimensional dataset of one of the compatible datatypes corresponding to the 
            "pixel_channels.uXX" loose-representation file specified in 3. The numeric type
            does not need to be otherwise specified as it is known by the HDF5 dataset. The
            2-dimensional dataset is natural as the pixel_channels represents the overall image
            data. The width and height of the 2-dimensional dataset must correspond to the specified
            "image_width_pixels" and "image_height_pixels", respectively.

### `indices`
A 1-dimensional dataset of one of the compatible datatypes corresponding to the
            "indices.uXX" loose-representation file specified in section 3. The numeric type does
            not need to be otherwise specified as it is known by the HDF5 dataset.

### `intensities`
A 1-dimensional dataset of one of the compatible datatypes corresponding to the
            "intensities.xXX" loose-representation file specified in section 3. The numeric type
            does not need to be otherwise specified as it is known by the HDF5 dataset.

##Image orientation and shape
All mspix images start with pixel at index 0, row: 0, column: 0 where row 0 is at the top
        of the image and column 0 is at the left of the image. The image then proceeds in the
        "x" direction (along the row 0, incrementing column; row-major) until it is at the right
        border, the next pixel is then at the far-left in row 1.

The layout of mspix images is constant and independant of the stage or beam motion during
        acquisition (e.g. serpentine, typewriter, up-down, or left-right). If information
        regarding the stage motion during acquisition is desired, it may be provided in the
        metadata.json file or as a separate file within the folder/zip file.

## Limitations and restrictions

The mspix file format takes advantage of some restrictions, limitations, and constraints 
that are common in current mass spectrometry imaging experiments. These constraints do
not apply to many mass spectrometry images, which means that mspix is not suitable for
losslessly representing such images. Some considerations are that:

- The data are "sparse" and do not contain digitizer noise. This means that the
   "processed" and "continuous" options of, for exmaple, the .imzML file format do not
   apply as all mspix files are "processed". Similarly, "event-based" data that may be
   encountered with some TOF-SIMS data must be converted to "processed" data for mspix
   conversion.

- mspix images are rectangular. It is possible to store non-rectangular imaged areas by
   "cropping" or deleting pixels outside of a defined area. Represenging these empty
   "corner" pixels can be very space efficient (taking only a few bytes each).
   Images with pixels that are not sampled on a rectangular grid cannot be represented
   without re-binning the pixels (and thus some data loss) using mspix.

- The data must be represented with a single 1-dimensional axis. This excludes some
   multi-tandem-MS and ion mobility data as some data acquired with these modes are not
   able to be represented efficiently with a single axis. MS2 data with a common second
   axis may be able to be represented by extending the primary mass axis or providing
   additional metadata that tag specific pixels or axis channels with fragmentation
   information; however such representations are left unspecified.

- Only a single, 2-dimensional image will be contained in the mspix file. This means
   that 3-d images or depth-profiling images are not supported in the 1.0.0 version
   of mspix.
