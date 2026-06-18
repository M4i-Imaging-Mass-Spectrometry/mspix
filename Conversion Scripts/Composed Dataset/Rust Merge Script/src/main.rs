use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write, BufWriter};
use std::path::PathBuf;

use anyhow::Result;
use byteorder::{LittleEndian, ReadBytesExt};
use nohash_hasher::IntMap;
use serde::Deserialize;

const MAX_WIDTH: usize = 42967;
const MAX_HEIGHT: usize = 26450;
const MASS_AXIS_LENGTH: usize = 3862;

#[derive(Debug, Deserialize)]
struct Metadata {
    image_width_pixels: usize,
    image_height_pixels: usize,
}

#[derive(Debug)]
struct Mspix {
    meta: Metadata,
    curr_offset: usize,
    pixel_channels: File,
    indices: File,
    intensities: File,
    index_offset: usize,
}

fn main() -> Result<()> {
    let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
    let dirs: Vec<_> = std::fs::read_dir(&cwd).unwrap()
        .filter_map(Result::ok)
        .filter(|entry| entry.path().is_dir() && entry.path().extension().and_then(|s| s.to_str()) == Some("mspix"))
        .collect();
    dbg!(&dirs);

    let mut mspix_files = vec![];
    for dir in dirs.iter() {
        dbg!(dir);
        let dir = dir.path();
        let meta: Metadata = serde_json::from_reader(File::open(dir.join("metadata.json")).unwrap()).unwrap();
        let pixel_channels = File::open(dir.join("pixel_channels.u16")).unwrap();
        let indices = File::open(dir.join("indices.u16")).unwrap();
        let intensities = File::open(dir.join("intensities.u16")).unwrap();
        let mspix = Mspix {
            meta,
            pixel_channels,
            indices,
            intensities,
            curr_offset: 0,
            index_offset: 0,
        };
        mspix_files.push(mspix);
    }

    std::fs::create_dir_all(cwd.join("merged")).unwrap();

    // std::process::exit(0);
    let mut merged_indices = BufWriter::new(File::create("merged/indices.u16").unwrap());
    let mut merged_intensities = BufWriter::new(File::create("merged/intensities.u16").unwrap());
    let mut merged_pixel_channels = BufWriter::new(File::create("merged/pixel_channels.u16").unwrap());

    let mut idx_buf = vec![0u8; 2 * MASS_AXIS_LENGTH];
    let mut int_buf = vec![0u8; 2 * MASS_AXIS_LENGTH];
    let mut rdr_idx;
    let mut rdr_int;

    let mut merged_map: IntMap<u16, u16> = IntMap::default();

    for y in 0..MAX_HEIGHT {
        if y % 100 == 0 {
            dbg!(y);
        }
        for x in 0..MAX_WIDTH {
            merged_map.clear();

            for (i, mspix_file) in mspix_files.iter_mut().enumerate() {
                if y >= mspix_file.meta.image_height_pixels || x >= mspix_file.meta.image_width_pixels {
                    continue;
                }
                
                mspix_file.pixel_channels.seek(SeekFrom::Start(mspix_file.curr_offset as u64 * 2)).unwrap();
                let count = mspix_file.pixel_channels.read_u16::<LittleEndian>().unwrap() as usize;
                if x == 0 && y % 100 == 0 {
                    println!("    {} {} {} {} {}", i, count, mspix_file.curr_offset, mspix_file.index_offset, mspix_file.indices.metadata().unwrap().len());
                }
                mspix_file.curr_offset += 1;
                if count > MASS_AXIS_LENGTH {
                    eprintln!("Skipping pixel ({}, {}) in file {}: unexpected count {}", x, y, i, count);
                    continue;
                } else if count == 0 {
                    continue;
                }
                let file_len = mspix_file.indices.metadata()?.len() as usize;
                let bytes_to_read = count * 2;
                if mspix_file.index_offset + bytes_to_read > file_len {
                    eprintln!("Skipping pixel ({}, {}) in file {}: trying to read {} bytes at offset {}, but file length is {}", 
                        x, y, i, bytes_to_read, mspix_file.index_offset, file_len);
                    continue;
                }

                mspix_file.indices.seek(SeekFrom::Start(mspix_file.index_offset as u64)).unwrap();
                mspix_file.intensities.seek(SeekFrom::Start(mspix_file.index_offset as u64)).unwrap();

                mspix_file.indices.read_exact(&mut idx_buf[..count * 2]).unwrap();
                mspix_file.intensities.read_exact(&mut int_buf[..count * 2]).unwrap();
                mspix_file.index_offset += count * 2;
                rdr_idx = std::io::Cursor::new(&idx_buf[..count * 2]);
                rdr_int = std::io::Cursor::new(&int_buf[..count * 2]);

                for _ in 0..count {
                    let index = rdr_idx.read_u16::<LittleEndian>().unwrap();
                    let intensity = rdr_int.read_u16::<LittleEndian>().unwrap();
                    *merged_map.entry(index).or_insert(0) += intensity as u16;
                }
            }

            merged_pixel_channels.write_all(&(merged_map.len() as u16).to_le_bytes()).unwrap();
            for (&index, &intensity) in merged_map.iter() {
                merged_indices.write_all(&(index as u16).to_le_bytes()).unwrap();
                merged_intensities.write_all(&intensity.to_le_bytes()).unwrap();
            }
        }
    }

    Ok(())
}
