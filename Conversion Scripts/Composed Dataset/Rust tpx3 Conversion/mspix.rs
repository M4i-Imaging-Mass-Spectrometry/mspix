use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use crate::ErrStr;
use crate::MSG;
use crate::interface;
use crate::interface::Interface;
use crate::mass;
use crate::mass::time_to_mass;
use crate::reader::TPX3Reader;
use crate::tpx::Coord;
use crate::writer::{convert_file_bitsize, write_to_bytes};

const fn num_to_bitsize(n: usize) -> usize {
    match n {
        0..=255 => 8,
        256..=65535 => 16,
        _ => 32,
    }
}

#[derive(Debug, serde::Serialize)]
pub struct MspixMaker {
    // mandatory fields
    pub mspix_version:         String,
    pub image_width_pixels:    u64,
    pub image_height_pixels:   u64,
    // additional fields
    #[serde(skip)]
    pub mspix_file:            PathBuf,
    #[serde(skip)]
    pub pix_ind_bitsize:       usize,
    #[serde(skip)]
    pub pix_int_bitsize:       usize,
    #[serde(skip)]
    pub offsets:               Vec<u32>,
    pub form:                  interface::Interface,
    pub original_file:         String,
    pub image_width_mm:        f64,
    pub image_height_mm:       f64,
    pub acquisition_time_ps:   i64,
    pub tof_pulse_count:       u64,
    pub highest_mz_value:      i64,
    pub spectral_channels_len: usize,
    pub spectral_channels:     Vec<f64>,
    pub spectral_intensities:  Vec<f64>,
    pub mass_channels:         Vec<f64>,
}

impl MspixMaker {
    pub fn new(form: &Interface, mspix_file: PathBuf) -> Result<MspixMaker, String> {
        let spectral_channels_len = ((form.tof_len / form.bin_time) + 1) as usize;
        let spectral_channels: Vec<f64> =
            (0..spectral_channels_len).map(|i| (form.bin_time as usize * i) as f64 / 1E3).collect();
        let mass_channels: Vec<f64> =
            spectral_channels.iter().map(|&t| time_to_mass(t / 1E3, form.cal)).collect();

        Ok(MspixMaker {
            // mandatory fields
            mspix_version: "1.0.0".to_string(),
            spectral_channels,
            spectral_intensities: vec![0.0; spectral_channels_len],
            image_width_pixels: form.img.cols() as u64,
            image_height_pixels: form.img.rows() as u64,
            // additional fields
            form: form.clone(),
            spectral_channels_len,
            mspix_file,
            pix_ind_bitsize: 0, // size in bits of the index blob
            pix_int_bitsize: 0, // size in bits of the intensity blob
            offsets: vec![],
            original_file: form.img.file.display().to_string(),
            image_width_mm: form.img.width,
            image_height_mm: form.img.height,
            highest_mz_value: 0,
            acquisition_time_ps: 0,
            tof_pulse_count: 0,
            mass_channels,
        })
    }

    pub fn make_mspix(&mut self, coords: &[Coord], deads: &[u16]) -> Result<(), String> {
        let img = self.form.img.clone();
        let (cols, rows, fov, ppmm) = (img.cols(), img.rows(), img.fov, img.pixels_per_mm);
        let (sin, cos) = (img.angle.sin(), img.angle.cos());

        // we make the (relatively safe) assumption that no individual pixel will exceed 4 bytes
        let tic_img = img.buf(coords, deads);
        let tic_sums: Vec<u64> =
            tic_img.chunks(cols).map(|r| r.iter().map(|x| *x as u64).sum()).collect();
        self.pix_int_bitsize = num_to_bitsize(*tic_img.iter().max().unwrap() as usize);
        self.pix_ind_bitsize = num_to_bitsize(self.spectral_channels_len);
        let ints = format!("intensities.u{}", self.pix_int_bitsize);
        let inds = format!("indices.u{}", self.pix_ind_bitsize);
        let mut int_file = File::create_new(self.mspix_file.join(&ints)).s()?;
        let mut ind_file = File::create(self.mspix_file.join(inds)).s()?;
        let mut image = vec![Some(vec![Pixel::default(); cols]); rows];
        let mut pix_sums = vec![0; tic_sums.len()];
        let mut row = 0;
        MSG.set("Extraction progress: 0%");
        for (pulse, coord) in TPX3Reader::new(&img.file).zip(coords) {
            self.tof_pulse_count += 1;
            self.acquisition_time_ps = pulse.time;
            for hit in pulse.hits.iter().filter(|h| h.size > 1 || !h.is_dead(deads)) {
                let (icol, irow) = coord.index(hit, sin, cos, fov, ppmm);
                if irow < rows && icol < cols {
                    let tof_ps = (hit.toa - pulse.time).rem_euclid(self.form.tof_len); // positive modulo
                    let tof: i64 = mass::twc(tof_ps as f64, hit.tot as f64, self.form.twc) as i64;
                    let intensity = if img.unit_intensity { 1 } else { hit.tot };
                    assert!(intensity != 0, "{:?}", hit);
                    let index = (tof / self.form.bin_time).max(0); // all negatives times to index 0
                    if let Some(Some(current_row)) = image.get_mut(irow) {
                        let pixel = &mut current_row[icol];
                        pixel.indices.push(index);
                        pixel.intensities.push(intensity.into());
                        pix_sums[irow] += intensity as u64;
                        self.spectral_intensities[index as usize] += intensity as f64;
                    }
                }
            }
            while row < tic_sums.len() && pix_sums[row] >= tic_sums[row] {
                let status = if pix_sums[row] != tic_sums[row] { "Warning!" } else { "Ok!" };
                println!("Writing row {}; {}: {} {}", row, status, tic_sums[row], pix_sums[row]);
                let percent = (row as f64 / rows as f64) * 100.0;
                MSG.set(&format!("Extraction progress: {:.2}%", percent));
                for pixel in image[row].take().unwrap() {
                    self.write_pixel(pixel, &mut ind_file, &mut int_file);
                }
                row += 1;
            }
        }
        self.cleanup_and_finish(&ints, &tic_img);
        Ok(())
    }

    fn write_pixel(&mut self, mut pixel: Pixel, ind: &mut File, int: &mut File) {
        pixel.sort();
        let (indices, ints): (Vec<u32>, Vec<u32>) = pixel.make_indices_intensities();
        if !ints.is_empty() {
            debug_assert!(*ints.iter().min().unwrap() != 0);
        }
        write_to_bytes(&indices, self.pix_ind_bitsize, ind);
        write_to_bytes(&ints, self.pix_int_bitsize, int);
        self.offsets.push(indices.len() as u32);
        self.highest_mz_value =
            self.highest_mz_value.max((*ints.iter().max().unwrap_or(&0)).into());
    }

    fn cleanup_and_finish(&mut self, input: &str, tic_img: &[u32]) -> Result<(), String> {
        let pix_int_bitsize = num_to_bitsize(self.highest_mz_value as usize);
        assert!(self.pix_int_bitsize >= pix_int_bitsize);
        if self.pix_int_bitsize > pix_int_bitsize {
            // we can shrink int file
            convert_file_bitsize(
                self.pix_int_bitsize,
                pix_int_bitsize,
                &self.mspix_file.join(input),
                &self.mspix_file.join(&format!("intensities.u{}", pix_int_bitsize)),
            );
            self.pix_int_bitsize = pix_int_bitsize;
        }
        // bitsize to store maximum pixel intensity; should always be an unsigned integer
        let pixel_intensities_file = format!("pixel_intensities.u{}", self.pix_int_bitsize);
        let mut pixel_ints_file = File::create(self.mspix_file.join(pixel_intensities_file)).s()?;
        write_to_bytes(tic_img, self.pix_int_bitsize, &mut pixel_ints_file);
        // bitsize to store pixel indices; should always be an unsigned integer
        let pixel_channels_file = format!("pixel_channels.u{}", self.pix_ind_bitsize);
        let mut pixel_inds_file = File::create(self.mspix_file.join(pixel_channels_file)).s()?;
        write_to_bytes(&self.offsets, self.pix_ind_bitsize, &mut pixel_inds_file);

        let mut metadata = File::create(self.mspix_file.join("metadata.json")).s()?;
        metadata.write_all(serde_json::to_string(&self).unwrap().as_bytes()).unwrap();
        Ok(())
    }
}


#[derive(Clone, Default)]
pub struct Pixel {
    indices:     Vec<i64>, // index/time
    intensities: Vec<u32>, // intensity
}

impl Pixel {
    pub fn make_indices_intensities(&self) -> (Vec<u32>, Vec<u32>) {
        let (mut indices, mut intensities, mut prev) = (vec![], vec![], 0);
        for (&index, &intensity) in self.indices.iter().zip(&self.intensities) {
            if index > prev || index == 0 {
                indices.push(index as u32);
                intensities.push(intensity);
            } else if !indices.is_empty() {
                *intensities.last_mut().unwrap() += intensity;
            }
            if index >= 0 {
                prev = index;
            }
        }
        (indices, intensities)
    }

    pub fn sort(&mut self) {
        // this is so dumb. Why is there no function to sort two vectors based on one???
        let mut sort_indices = (0..self.indices.len()).collect::<Vec<_>>();
        sort_indices.sort_by_key(|&i| &self.indices[i]);
        self.indices = sort_indices.iter().map(|&i| self.indices[i]).collect::<Vec<_>>();
        self.intensities = sort_indices.iter().map(|&i| self.intensities[i]).collect::<Vec<_>>();
    }
}
