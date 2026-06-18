import os
import json
import numpy as np
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict

MAX_WIDTH = 42967
MAX_HEIGHT = 26450
NUM_PIXELS = MAX_WIDTH * MAX_HEIGHT
MASS_AXIS_LENGTH = 3862  # Hard-coded as specified

def load_metadata(folder):
    with open(Path(folder) / "metadata.json") as f:
        return json.load(f)

def compute_pixel_offsets(pixel_channels):
    return np.concatenate(([0], np.cumsum(pixel_channels[:-1])))

def process_datasets(input_dirs, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # Dictionary to collect pixel data as {pixel_idx: {index: intensity}}
    pixel_data = defaultdict(lambda: defaultdict(float))

    print("Processing datasets...")
    for input_dir in tqdm(input_dirs):
        print(f"processing {input_dir}")
        input_dir = Path(input_dir)
        meta = load_metadata(input_dir)
        width, height = meta["image_width_pixels"], meta["image_height_pixels"]
        assert width <= MAX_WIDTH and height <= MAX_HEIGHT

        num_pixels = width * height

        # Load memory-mapped binary files
        pixel_channels = np.memmap(input_dir / "pixel_channels.u16", dtype=np.uint16, mode='r', shape=(num_pixels,))
        indices = np.memmap(input_dir / "indices.u16", dtype=np.uint16, mode='r')
        intensities = np.memmap(input_dir / "intensities.u16", dtype=np.uint16, mode='r')

        offsets = compute_pixel_offsets(pixel_channels)

        for y in range(height):
            if y % 2000 == 0:
                print(f"row {y}")
            for x in range(width):
                pixel_idx = y * width + x
                global_idx = y * MAX_WIDTH + x
                count = pixel_channels[pixel_idx]
                if count == 0:
                    continue

                start = offsets[pixel_idx]
                end = start + count

                inds = indices[start:end]
                ints = intensities[start:end]

                accum = pixel_data[global_idx]
                for i, intensity in zip(inds, ints):
                    accum[i] += intensity

    print("Flattening merged data...")
    new_pixel_channels = np.zeros(NUM_PIXELS, dtype=np.uint32)
    total_entries = sum(len(v) for v in pixel_data.values())
    merged_indices = np.zeros(total_entries, dtype=np.uint32)
    merged_intensities = np.zeros(total_entries, dtype=np.uint32)

    offset = 0
    for pixel_idx in tqdm(range(NUM_PIXELS)):
        if pixel_idx in pixel_data:
            items = sorted(pixel_data[pixel_idx].items())
            indices_chunk, intensities_chunk = zip(*items)
            size = len(indices_chunk)
            new_pixel_channels[pixel_idx] = size
            merged_indices[offset:offset+size] = indices_chunk
            merged_intensities[offset:offset+size] = intensities_chunk
            offset += size

    print("Writing merged files...")
    new_pixel_channels.tofile(output_dir / "pixel_channels.u32")
    merged_indices[:offset].tofile(output_dir / "indices.u32")
    merged_intensities[:offset].tofile(output_dir / "intensities.u32")

    merged_meta = {
        "width": MAX_WIDTH,
        "height": MAX_HEIGHT,
        "mass_axis_length": MASS_AXIS_LENGTH,
    }
    with open(output_dir / "metadata.json", "w") as f:
        json.dump(merged_meta, f, indent=4)

    print(f"Merge complete. Output written to {output_dir}")

os.chdir(os.path.dirname(os.path.abspath(__file__)))
input_folders = [f for f in os.listdir() if f.endswith(".mspix")]
process_datasets(input_folders, "lemon_42_5x26_merged")
