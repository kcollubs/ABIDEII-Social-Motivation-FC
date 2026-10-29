import os, re, glob, urllib.request, json, struct, sys
BASE="https://fcp-indi.s3.amazonaws.com/data/Projects/ABIDE2/Derivatives/DCAN"
PROJ="/nas/longleaf/home/sameliza/cyberneuro/ABIDEII_project"
out="/nas/longleaf/home/sameliza/cyberneuro/abide2_audit/dcan_small"
os.makedirs(out, exist_ok=True)
def get(url, rng=None):
    req=urllib.request.Request(url)
    if rng: req.add_header("Range", f"bytes=0-{rng}")
    with urllib.request.urlopen(req, timeout=120) as r: return r.read()
subs=[]
for f in sorted(glob.glob(f"{PROJ}/func/*.ptseries.nii")):
    m=re.match(r"(sub-\d+)_(ses-[A-Za-z0-9]+baseline)_", os.path.basename(f))
    subs.append((m.group(1), m.group(2), os.path.basename(f)))
log=[]
for sub, ses, fname in subs:
    prefix=f"{sub}/{ses}/func/{sub}_{ses}_task-rest"
    row={"sub":sub,"ses":ses}
    try:
        mm=get(f"{BASE}/{prefix}_desc-filtered_motion_mask.mat")
        open(f"{out}/{sub}_{ses}_motion_mask.mat","wb").write(mm); row["mask_bytes"]=len(mm)
    except Exception as e: row["mask_err"]=str(e)[:120]
    try:
        h=get(f"{BASE}/{prefix}_run-1_bold_timeseries.dtseries.nii", rng=1023)
        open(f"{out}/{sub}_{ses}_raw_hdr.bin","wb").write(h)
        dim=struct.unpack("<8q", h[16:80]); row["raw_nt"]=dim[6]; row["raw_ncol"]=dim[5]
        row["raw_step_s"]=struct.unpack("<8d", h[80:144])[1]
    except Exception as e: row["raw_err"]=str(e)[:120]
    log.append(row)
json.dump(log, open(f"{out}/fetch_log.json","w"), indent=1)
ok=sum(1 for r in log if "mask_bytes" in r); okr=sum(1 for r in log if "raw_nt" in r)
print("subjects", len(log), "masks", ok, "raw headers", okr)
errs=[r for r in log if "mask_err" in r or "raw_err" in r]
print("errors:", len(errs)); print(errs[:3])
