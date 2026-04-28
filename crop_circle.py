#!/usr/bin/env python3
"""Crop a circular region from all TIFF files in the current directory,
then assemble the cropped frames into an animated GIF ordered by the frame
number embedded in each filename (the digits after 'Overlay_').

The circle is defined by center O(X1, Y1) and a point B(X2, Y2) on its edge.
The radius is the Euclidean distance OB. Output images are square PNGs with
the circle visible and the area outside the circle transparent.

Usage:
    python3 crop_circle.py OX OY BX BY [FPS]

    OX, OY  Centre of the circle
    BX, BY  Any point on the circle edge — radius = distance(O, B)
    FPS     Frames per second for the output GIF (default: 5)

Example:
    python3 crop_circle.py 2612 1693 2227 1278
    python3 crop_circle.py 2612 1693 2227 1278 10
"""

import sys
import re
import math
import glob
from pathlib import Path
from PIL import Image, ImageDraw


def crop_circle(ox: int, oy: int, bx: int, by: int, fps: float = 5) -> None:
    # Radius = Euclidean distance from centre O to edge point B
    radius = math.sqrt((bx - ox) ** 2 + (by - oy) ** 2)
    r = int(math.ceil(radius))

    # Bounding box of the circle in the source image coordinate space
    left  = ox - r
    upper = oy - r
    right = ox + r
    lower = oy + r

    tiff_files = sorted(glob.glob("*.tiff") + glob.glob("*.tif"))
    if not tiff_files:
        print("No TIFF files found in current directory.")
        return

    cropped_paths: list[tuple[int, str]] = []

    for tiff_path in tiff_files:
        stem = Path(tiff_path).stem
        out_name = f"{stem}_O_{ox}_{oy}_B_{bx}_{by}.png"

        with Image.open(tiff_path) as img:
            img_w, img_h = img.size

            # Clamp the bounding box to the actual image dimensions
            crop_left  = max(left,  0)
            crop_upper = max(upper, 0)
            crop_right = min(right, img_w)
            crop_lower = min(lower, img_h)
            region = img.crop((crop_left, crop_upper, crop_right, crop_lower))

            # Create a full (2r × 2r) RGBA canvas and paste the region at
            # the correct offset so the circle stays centred
            diameter = r * 2
            canvas = Image.new("RGBA", (diameter, diameter), (0, 0, 0, 0))
            paste_x = crop_left - left
            paste_y = crop_upper - upper
            canvas.paste(region.convert("RGBA"), (paste_x, paste_y))

            # Build a circular alpha mask — pixels outside the circle become
            # fully transparent (alpha = 0)
            mask = Image.new("L", (diameter, diameter), 0)
            ImageDraw.Draw(mask).ellipse((0, 0, diameter - 1, diameter - 1), fill=255)
            canvas.putalpha(mask)

            canvas.save(out_name)

        print(f"  Saved: {out_name}")

        # Extract the frame number from the stem (digits right after 'Overlay_')
        match = re.search(r"Overlay_(\d+)", stem, re.IGNORECASE)
        frame_num = int(match.group(1)) if match else 0
        cropped_paths.append((frame_num, out_name))

    print(f"\nDone. Cropped {len(tiff_files)} file(s) with radius={radius:.1f}px.")

    # ---------- assemble animated GIF ----------
    # Sort frames by the numeric index so the animation plays in order
    cropped_paths.sort(key=lambda t: t[0])

    # GIF does not support full RGBA; convert to RGB (transparent areas → black)
    frames = [Image.open(p).convert("RGB") for _, p in cropped_paths]

    gif_name = f"Overlay_O_{ox}_{oy}_B_{bx}_{by}.gif"
    duration_ms = int(1000 / fps)   # Pillow expects duration in milliseconds

    frames[0].save(
        gif_name,
        save_all=True,
        append_images=frames[1:],
        loop=0,             # 0 = loop forever
        duration=duration_ms,
        optimize=False,
    )
    print(f"Movie saved: {gif_name}  ({len(frames)} frames @ {fps} fps)")


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        print("Usage: python3 crop_circle.py OX OY BX BY [FPS]")
        print("Example: python3 crop_circle.py 2612 1693 2227 1278")
        print("         python3 crop_circle.py 2612 1693 2227 1278 10")
        sys.exit(1)

    try:
        OX, OY, BX, BY = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        FPS = float(sys.argv[5]) if len(sys.argv) == 6 else 5.0
    except ValueError:
        print("Error: Coordinates must be integers; FPS must be a number.")
        sys.exit(1)

    crop_circle(OX, OY, BX, BY, FPS)
