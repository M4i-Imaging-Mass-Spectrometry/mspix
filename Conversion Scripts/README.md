# Scripts to convert formats to mspix

## nanotof_conversion.py

This script converts a simple NanoTOF (TOF-SIMS) acquired dataset to the mspix format. It implements a portion of mspix itself, which could be avoided by using the online mspix.py version. The easily-observed (text-based) NanoTOF metadata is simply copy/pasted into the mspix metadata.json file.

## WatersRawToMsPix.py

This script coverts a single Waters .raw file into the mspix format. It was custom-written for the single file and hard-codes the image dimensions. It simply pastes each recorded spectrum into another column. This script requires Waters MassLynxRawDataReader API (not included).

## Composed Dataset

This folder contains three sub-items:

1. A Rust tpx3 Conversion folder
   This conversion folder contains a file "mspix.rs" which is one module of a larger .tpx3 file analysis platform. This platform is not published at this time, however, an earlier version of this platform was published and is available on GitHub here: https://github.com/M4i-Imaging-Mass-Spectrometry/fast-mass-microscopy-example . The "mspix.rs" file provides a rough overview of the conversion process and should be easy to follow even without the associated data structures that are being read from. 
2. A Rust Merge Script folder
   This contains a full Rust "main.rs" and associated "cargo.toml" file that was used to merge seven tpx3 "passes" into a single (partial) mspix image.
3. A post_merge_cleanup_script.py . This Python file was used to finalize the mspix image and perform some tasks that were easier to accomplish in Python than Rust.

Together, these three represent the path that the multi-step fast mass microscopy image for the "Composed Surface" dataset had for conversion. A more streamlined, single-step conversion program is intended to be published at a later date. 