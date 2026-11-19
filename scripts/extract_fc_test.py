import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path.home() / "cyberneuro"
REPO = ROOT / "ABIDEII-Social-Motivation-FC"
PROJECT = ROOT / "ABIDEII_project"
MASK_DIR = ROOT / "abide2_audit/dcan_small"

sys.path.insert(0, str(REPO / "scripts"))
from dcan_io import load_mask

config = json.loads(
    (REPO / "config/roi_specification_v1.json").read_text()
)

subject = "28752"

files = list((PROJECT / "func").glob(f"sub-{subject}_*.ptseries.nii"))
assert len(files) == 1, f"Expected one CIFTI file, found {len(files)}"

img = nib.load(str(files[0]))
data = np.asarray(img.dataobj)
names = list(img.header.get_axis(1).name)
index = {name: i for i, name in enumerate(names)}

mask_files = list(MASK_DIR.glob(f"sub-{subject}_*motion_mask.mat"))
assert len(mask_files) == 1, f"Expected one mask, found {len(mask_files)}"

mm = load_mask(str(mask_files[0]))
keep = None

for entry in mm[2]:
    f = entry[2]
    threshold = round(float(f["FD_threshold"][2][0]), 2)
    if threshold == 0.3:
        removed = np.asarray(f["frame_removal"][2], dtype=int)
        keep = removed == 0
        break

assert keep is not None, "0.3 mm mask not found"
assert len(keep) == data.shape[0], "Mask and CIFTI lengths differ"

usable = int(keep.sum())
print("Subject:", subject)
print("CIFTI shape:", data.shape)
print("Retained frames:", usable)

tr = float(img.header.get_axis(0).step)
usable_minutes = usable * tr / 60

print("Usable duration (minutes):", round(usable_minutes, 2))

assert usable_minutes >= 4, "Excluded: fewer than 4 usable minutes"

results = {}

for pathway, (roi1, roi2) in config["pathways"].items():
    z_values = []

    for hemi in ("L", "R"):
        name1 = config["rois"][roi1][hemi]
        name2 = config["rois"][roi2][hemi]

        x = data[keep, index[name1]]
        y = data[keep, index[name2]]

        valid = np.isfinite(x) & np.isfinite(y)
        assert valid.sum() >= 4, "Insufficient valid paired frames"

        r = float(np.corrcoef(x[valid], y[valid])[0, 1])
        assert np.isfinite(r), "Non-finite correlation"

        z = float(np.arctanh(np.clip(r, -0.999999, 0.999999)))
        z_values.append(z)

        print(f"{pathway} {hemi}: r={r:.4f}, z={z:.4f}")

    results[pathway] = float(np.mean(z_values))

print("\nBilateral Fisher-z connectivity:")
print(json.dumps(results, indent=2))
