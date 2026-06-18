# mspix File Format 

**mspix** files are designed to be an efficient storage and analysis format for mass spectrometry imaging data produced using time-of-flight mass spectrometetry, with an emphasis on being highly space and access-time efficient when working with large numbers (>5 million) of pixels.

## Main goals of this file format:
- Store large pixel-count mass images in a space and access-efficient manner mspix files are designed to accomodate ~1 billion pixel mass images and are recommended for mass images exceeding 5-10 million pixels, mass images with a single mass axis, and/or targeted mass images.

- Provide easy inspection and simple, obvious, formats for the loose files.
                The intention behind this decision is to enable ease of tooling and libraries to
                be written by others.
                
- Allow for user-added data in the form of additional files, axes, and metadata
                This is done using a folder/HDF5-based structure.

- Only a single, 2-dimensional image will be contained in the mspix file. This means
   that 3-d images or depth-profiling images are not supported in the 1.0.5 version
   of mspix.
