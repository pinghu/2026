#!/usr/bin/env python3
"""Crop a rectangular region from all TIFF files in the current directory,
then assemble the cropped frames into an animated GIF ordered by the frame
number embedded in each filename (the digits after 'Overlay_').

Usage:
    python3 crop_images.py X1 Y1 X2 Y2 [FPS]

    FPS  Frames per second for the output GIF (default: 5)

Example:
    python3 crop_images.py 1851 3356 3207 3971
    python3 crop_images.py 1851 3356 3207 3971 10
"""

import sys
import re
import glob
from pathlib import Path
from PIL import Image


def crop_images(x1: int, y1: int, x2: int, y2: int, fps: float = 5) -> None:
    left  = min(x1, x2)
    upper = min(y1, y2)
    right = max(x1, x2)
    lower = max(y1, y2)

    tiff_files = sorted(glob.glob("*.tiff") + glob.glob("*.tif"))
    if not tiff_files:
        print("No TIFF files found in current directory.")
        return

    cropped_paths: list[tuple[int, str]] = []

    for tiff_path in tiff_files:
        stem = Path(tiff_path).stem
        out_name = f"{stem}_{x1}_{y1}_{x2}_{y2}.png"
        with Image.open(tiff_path) as img:
            cropped = img.crop((left, upper, right, lower))
            cropped.save(out_name)
        print(f"  Saved: {out_name}")

        # Extract frame number from stem (digits right after 'Overlay_')
        match = re.search(r"Overlay_(\d+)", stem, re.IGNORECASE)
        frame_num = int(match.group(1)) if match else 0
        cropped_paths.append((frame_num, out_name))

    print(f"\nDone. Cropped {len(tiff_files)} file(s).")

    # ---------- assemble animated GIF ----------
    cropped_paths.sort(key=lambda t: t[0])
    frames = [Image.open(p).convert("RGB") for _, p in cropped_paths]

    gif_name = f"Overlay_{x1}_{y1}_{x2}_{y2}.gif"
    duration_ms = int(1000 / fps)          # ms per frame

    frames[0].save(
        gif_name,
        save_all=True,
        append_images=frames[1:],
        loop=0,                            # loop forever
        duration=duration_ms,
        optimize=False,
    )
    print(f"Movie saved: {gif_name}  ({len(frames)} frames @ {fps} fps)")


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        print("Usage: python3 crop_images.py X1 Y1 X2 Y2 [FPS]")
        print("Example: python3 crop_images.py 1851 3356 3207 3971")
        print("         python3 crop_images.py 1851 3356 3207 3971 10")
        sys.exit(1)

    try:
        X1, Y1, X2, Y2 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        FPS = float(sys.argv[5]) if len(sys.argv) == 6 else 5.0
    except ValueError:
        print("Error: Coordinates must be integers; FPS must be a number.")
        sys.exit(1)

    crop_images(X1, Y1, X2, Y2, FPS)
