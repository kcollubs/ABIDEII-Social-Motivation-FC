import csv
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

output_dir = REPO / "results"
output_dir.mkdir(exist_ok=True)

rows = []

for path in sorted((PROJECT / "func").glob("*.ptseries.nii")):
    subject = path.name.split("_")[0].replace("sub-", "")
    site = path.name.split("_")[1].replace("ses-", "").split("1")[0]

    row = {
        "participant_id": subject,
        "site": site,
        "status": "",
        "usable_minutes": "",
        "FC_NAc_rACC": "",
        "FC_AI_a24pr": ""
    }

    try:
        img = nib.load(str(path))
        data = np.asarray(img.dataobj)
        tr = float(img.header.get_axis(0).step)

        names = list(img.header.get_axis(1).name)
        index = {name: i for i, name in enumerate(names)}

        masks = list(MASK_DIR.glob(f"sub-{subject}_*motion_mask.mat"))
        if len(masks) != 1:
            raise ValueError(f"Expected one mask, found {len(masks)}")

        mm = load_mask(str(masks[0]))
        keep = None

        for entry in mm[2]:
            f = entry[2]
            threshold = round(float(f["FD_threshold"][2][0]), 2)

            if threshold == 0.3:
                removed = np.asarray(f["frame_removal"][2], dtype=int)
                keep = removed == 0
                break

        if keep is None or len(keep) != data.shape[0]:
            raise ValueError("Missing or mismatched censoring mask")


        if site == "OHSU":
            if data.shape[0] != 360:
                raise ValueError("Unexpected OHSU frame count")
            segments = [(0, 120), (120, 240), (240, 360)]
        else:
            segments = [(0, data.shape[0])]

        eligible = []
        for a, b in segments:
            indices = np.flatnonzero(keep[a:b]) + a
            if site != "OHSU" or len(indices) >= 50:
                eligible.append(indices)

        retained = sum(len(indices) for indices in eligible)
        usable_minutes = float(retained * tr / 60)
        row["usable_minutes"] = round(usable_minutes, 3)

        if usable_minutes < 4:
            row["status"] = "excluded_short_duration"
            rows.append(row)
            continue

        results = {}

        for pathway, (roi1, roi2) in config["pathways"].items():
            hemisphere_z = []

            for hemi in ("L", "R"):
                name1 = config["rois"][roi1][hemi]
                name2 = config["rois"][roi2][hemi]

                weighted_z = []
                weights = []

                for indices in eligible:
                    x = data[indices, index[name1]]
                    y = data[indices, index[name2]]
                    valid = np.isfinite(x) & np.isfinite(y)
                    n = int(valid.sum())

                    if n < 4:
                        raise ValueError("Insufficient valid paired frames")

                    r = float(np.corrcoef(x[valid], y[valid])[0, 1])
                    if not np.isfinite(r):
                        raise ValueError("Non-finite correlation")

                    z = float(np.arctanh(np.clip(r, -0.999999, 0.999999)))
                    weighted_z.append(z)
                    weights.append(n - 3)

                hemisphere_z.append(
                    float(np.average(weighted_z, weights=weights))
                )

            results[pathway] = float(np.mean(hemisphere_z))

        row["FC_NAc_rACC"] = results["NAc_rACC"]
        row["FC_AI_a24pr"] = results["AI_a24pr"]
        row["status"] = "included"

    except Exception as exc:
        row["status"] = "error"
        row["error"] = str(exc)

    rows.append(row)

out = output_dir / "fc_extraction_preliminary.csv"

fields = [
    "participant_id", "site", "status",
    "usable_minutes", "FC_NAc_rACC", "FC_AI_a24pr", "error"
]

with out.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

from collections import Counter
print("Saved:", out)
print("Participants:", len(rows))
print("Status counts:", dict(Counter(r["status"] for r in rows)))
