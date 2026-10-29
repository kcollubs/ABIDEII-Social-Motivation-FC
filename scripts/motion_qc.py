"""Motion QC summary for the ABIDE II ASD project (read-only on the project folder).

Inputs (all new or read-only):
  - ABIDEII_project phenotypic CSV and func/ file names (read only)
  - fmriprep_confounds/   downloaded desc-confounds_regressors.tsv (fMRIPrep)
  - mriqc_json/           downloaded MRIQC bold JSON (separate motion estimate)
  - dcan_small/           DCAN motion masks (*.mat) downloaded earlier
Outputs: motion_qc/ (new folder)
"""
import csv
import glob
import json
import math
import os
import re

import numpy as np
import pandas as pd

from dcan_io import load_mask

PROJ = '/nas/longleaf/home/sameliza/cyberneuro/ABIDEII_project'
ROOT = '/nas/longleaf/home/sameliza/cyberneuro/abide2_audit'
OUT = os.path.join(ROOT, 'motion_qc')
os.makedirs(OUT, exist_ok=True)
NA = float('nan')
THR_MAIN = 0.3
THR_SENS = (0.2, 0.5)


def isnan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


def mean(xs):
    xs = [v for v in xs if not isnan(v)]
    return sum(xs) / len(xs) if xs else NA


def r(x, nd):
    return '' if isnan(x) else round(x, nd)


# ---------- phenotype and imaging inventory (read only) ----------
pheno = {}
with open(os.path.join(PROJ, 'phenotypic', 'ASD Project Data(Phenotypic Data).csv'),
          encoding='utf-8-sig') as fh:
    for row in csv.DictReader(fh):
        if row['SUB_ID'].strip():
            pheno[row['SUB_ID'].strip()] = row

subs = {}
for f in sorted(glob.glob(os.path.join(PROJ, 'func', '*.ptseries.nii'))):
    m = re.search(r'sub-(\d+)_(ses-[A-Za-z0-9]+baseline)_', os.path.basename(f))
    subs[m.group(1)] = m.group(2)


# ---------- readers ----------
def dcan(sub):
    ses = subs[sub]
    mm = load_mask(os.path.join(ROOT, 'dcan_small', f'sub-{sub}_{ses}_motion_mask.mat'))
    grid = {}
    for e in mm[2]:
        f = e[2]
        t = round(float(f['FD_threshold'][2][0]), 2)
        grid[t] = dict(removed=[int(x) for x in f['frame_removal'][2]],
                       remaining=int(f['remaining_frame_count'][2][0]),
                       total=int(f['total_frame_count'][2][0]),
                       tr=float(f['epi_TR'][2][0]),
                       skip=int(f['skip'][2][0]))
    for t, g in grid.items():  # internal consistency of the mask
        assert g['total'] == len(g['removed']), (sub, t, 'length')
        assert g['total'] - sum(g['removed']) == g['remaining'], (sub, t, 'count')
    return grid


def fmriprep(sub):
    files = glob.glob(os.path.join(ROOT, 'fmriprep_confounds', f'sub-{sub}_*_desc-confounds_regressors.tsv'))
    runs = []
    for fp in files:
        m = re.search(r'_run-(\d+)_', os.path.basename(fp))
        run = int(m.group(1)) if m else 1
        vals = []
        with open(fp) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                v = row['framewise_displacement'].strip()
                vals.append(NA if v in ('n/a', '', 'nan', 'NaN') else float(v))
        runs.append(dict(run=run, file=os.path.basename(fp), vals=vals))
    runs.sort(key=lambda d: d['run'])
    return runs


def find_key(obj, key):
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for v in obj.values():
            hit = find_key(v, key)
            if hit is not None:
                return hit
    elif isinstance(obj, list):
        for v in obj:
            hit = find_key(v, key)
            if hit is not None:
                return hit
    return None


def mriqc(sub):
    files = sorted(glob.glob(os.path.join(ROOT, 'mriqc_json', f'sub-{sub}_*_bold.json')))
    vals = []
    for fp in files:
        js = json.load(open(fp))
        v = find_key(js, 'fd_mean')
        if v is not None:
            vals.append(float(v))
    return dict(n_files=len(files), fd=mean(vals) if vals else NA, values=vals)


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    return float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 else NA


# ---------- main loop ----------
summary, align_rows, ohsu_rows, mriqc_pairs, notes_all = [], [], [], [], []
frame_pairs_bins, frame_pairs_fd = [], []
keys_checked = False

for sub in sorted(subs):
    ses = subs[sub]
    site = pheno[sub]['SITE_ID']
    dx = pheno[sub]['DX_GROUP']
    grid = dcan(sub)
    total = grid[THR_MAIN]['total']
    tr = grid[THR_MAIN]['tr']

    runs = fmriprep(sub)
    fm_all = [v for run in runs for v in run['vals']]
    run_lens = [len(run['vals']) for run in runs]
    n_rows = len(fm_all)
    n_first_na = sum(1 for run in runs if run['vals'] and isnan(run['vals'][0]))
    n_other_na = sum(1 for run in runs for v in run['vals'][1:] if isnan(v))
    valid = [v for v in fm_all if not isnan(v)]
    has_fm = len(runs) > 0
    mean_fd = mean(fm_all) if has_fm else NA
    aligned = has_fm and n_rows == total

    note = []
    if not has_fm:
        note.append('no fMRIPrep confounds file on public S3')
    if len(runs) > 1:
        note.append(f'{len(runs)} runs ({",".join(str(r_["run"]) for r_ in runs)}) concatenated in mask')
    if has_fm and not aligned:
        note.append(f'fMRIPrep rows ({n_rows}) != mask frames ({total})')

    rem03 = grid[THR_MAIN]['removed']
    if aligned:
        hi = [(not isnan(fm_all[i]) and fm_all[i] > THR_MAIN) for i in range(total)]
        both = sum(1 for i in range(total) if rem03[i] == 1 and hi[i])
        dc_only = sum(1 for i in range(total) if rem03[i] == 1 and not hi[i])
        fm_only = sum(1 for i in range(total) if rem03[i] == 0 and hi[i])
        n_hi = sum(hi)
        bins = []
        for i in range(total):
            b = None
            for t in sorted(grid):
                if grid[t]['removed'][i] == 0:
                    b = t
                    break
            bins.append(b)
        pairs = [(bins[i], fm_all[i]) for i in range(total)
                 if bins[i] is not None and not isnan(fm_all[i])]
        rho = NA
        if len(pairs) > 2:
            rho = pearson(pd.Series([p[0] for p in pairs]).rank().values, pd.Series([p[1] for p in pairs]).rank().values)
        frame_pairs_bins += [p[0] for p in pairs]
        frame_pairs_fd += [p[1] for p in pairs]
        # run-start checks (OHSU: runs 2 and 3 start at frames 120, 240)
        starts = [sum(run_lens[:k]) for k in range(len(run_lens))]
        start_info = []
        for k, s in enumerate(starts):
            start_info.append(dict(run_index=k + 1, first_frame=s,
                                   fd_first_frame_na=isnan(fm_all[s]),
                                   dcan_removed_first3_at_0p3=sum(rem03[s:s + 3])))
        align_rows.append(dict(participant_id=sub, site=site, total_frames=total,
                               fmriprep_rows=n_rows, aligned=aligned,
                               frames_fd_gt_0p3=n_hi, dcan_removed_0p3=sum(rem03),
                               removed_and_fd_gt_0p3=both,
                               removed_but_fd_le_0p3_or_na=dc_only,
                               fd_gt_0p3_but_kept=fm_only,
                               spearman_dcan_bin_vs_fmriprep_fd=r(rho, 3),
                               run_starts=json.dumps(start_info)))
        if len(runs) > 1:
            for k, run in enumerate(runs):
                s = starts[k]
                seg = slice(s, s + len(run['vals']))
                kept = [1 - x for x in rem03[seg]]
                ohsu_rows.append(dict(participant_id=sub, site=site, run=run['run'],
                                      run_file=run['file'], run_frames=len(run['vals']),
                                      fmriprep_first_frame_na=isnan(run['vals'][0]),
                                      fmriprep_valid_fd=sum(1 for v in run['vals'] if not isnan(v)),
                                      run_mean_fd=r(mean(run['vals']), 3),
                                      dcan_removed_0p3=len(run['vals']) - sum(kept),
                                      retained_0p3=sum(kept),
                                      usable_min_0p3=round(sum(kept) * tr / 60, 2)))
    elif has_fm:
        note.append('frame-level comparison skipped')

    mq = mriqc(sub)
    if mq['n_files'] and not keys_checked:
        sample = json.load(open(sorted(glob.glob(os.path.join(ROOT, 'mriqc_json', '*_bold.json')))[0]))
        print('MRIQC JSON top-level keys (sample):', sorted(sample.keys())[:12], '... fd_mean found:',
              find_key(sample, 'fd_mean') is not None)
        keys_checked = True
    if has_fm and not isnan(mq['fd']):
        mriqc_pairs.append(dict(participant_id=sub, site=site, fmriprep_mean_fd=mean_fd,
                                mriqc_fd_mean=mq['fd'], mriqc_runs=mq['n_files']))

    ret03 = grid[THR_MAIN]['remaining']
    summary.append(dict(
        participant_id=sub, subject_label=f'sub-{sub}', session=ses, site=site, DX_GROUP=dx,
        TR_s=tr, total_frames=total, total_duration_min=round(total * tr / 60, 2),
        fmriprep_runs=len(runs), fmriprep_valid_fd_frames=len(valid),
        MEAN_FD_mm=r(mean_fd, 3),
        MEAN_FD_source=('fMRIPrep desc-confounds_regressors (pooled valid frames)' if has_fm
                        else 'unavailable (no per-frame FD found)'),
        retained_frames_0p3=ret03, censored_frames_0p3=total - ret03,
        pct_censored_0p3=round(100 * (total - ret03) / total, 1),
        usable_duration_min_0p3=round(ret03 * tr / 60, 2),
        retained_frames_0p2=grid[0.2]['remaining'],
        usable_duration_min_0p2=round(grid[0.2]['remaining'] * tr / 60, 2),
        retained_frames_0p5=grid[0.5]['remaining'],
        usable_duration_min_0p5=round(grid[0.5]['remaining'] * tr / 60, 2),
        MRIQC_fd_mean_separate_method=r(mq['fd'], 3), MRIQC_bold_files=mq['n_files'],
        mask_frames_match_fmriprep=aligned if has_fm else '',
        notes='; '.join(note)))

    notes_all.append(dict(participant_id=sub, fd_rows=n_rows, first_row_na=n_first_na,
                          other_na=n_other_na, valid=len(valid)))

# ---------- write tables ----------
sum_df = pd.DataFrame(summary)
sum_df.to_csv(os.path.join(OUT, 'fd_summary_table.tsv'), sep='\t', index=False)
pd.DataFrame(align_rows).to_csv(os.path.join(OUT, 'mask_alignment_checks.tsv'), sep='\t', index=False)
pd.DataFrame(ohsu_rows).to_csv(os.path.join(OUT, 'ohsu_runs.tsv'), sep='\t', index=False)
pd.DataFrame(notes_all).to_csv(os.path.join(OUT, 'fd_missing_value_log.tsv'), sep='\t', index=False)
mp = pd.DataFrame(mriqc_pairs)
mp.to_csv(os.path.join(OUT, 'mriqc_vs_fmriprep_overlap.tsv'), sep='\t', index=False)

# ---------- validation numbers ----------
val = {}
have = sum_df[sum_df.MEAN_FD_mm != '']
val['n_subjects'] = len(sum_df)
val['n_with_MEAN_FD'] = len(have)
val['n_without_MEAN_FD'] = len(sum_df) - len(have)
val['n_with_MRIQC_fd'] = int((sum_df.MRIQC_fd_mean_separate_method != '').sum())
val['n_MRIQC_only_no_fMRIPrep'] = int(((sum_df.MEAN_FD_mm == '') & (sum_df.MRIQC_fd_mean_separate_method != '')).sum())
val['n_no_FD_source'] = int(((sum_df.MEAN_FD_mm == '') & (sum_df.MRIQC_fd_mean_separate_method == '')).sum())
val['fd_first_row_na_total'] = int(sum(n['first_row_na'] for n in notes_all))
val['fd_other_na_total'] = int(sum(n['other_na'] for n in notes_all))
val['all_masks_internally_consistent'] = True
val['alignment_aligned_subjects'] = int(sum(1 for a in align_rows if a['aligned']))
val['alignment_fraction_removed_that_are_fd_gt_0p3'] = round(
    sum(a['removed_and_fd_gt_0p3'] for a in align_rows) / max(1, sum(a['dcan_removed_0p3'] for a in align_rows)), 3)
val['alignment_fraction_fd_gt_0p3_removed'] = round(
    sum(a['removed_and_fd_gt_0p3'] for a in align_rows) / max(1, sum(a['frames_fd_gt_0p3'] for a in align_rows)), 3)
val['alignment_pooled_spearman_bin_vs_fd'] = round(pearson(pd.Series(frame_pairs_bins).rank().values, pd.Series(frame_pairs_fd).rank().values), 3)
if len(mp):
    x, y = mp.mriqc_fd_mean.values, mp.fmriprep_mean_fd.values
    d = x - y
    val['mriqc_vs_fmriprep'] = dict(
        n=len(mp), pearson_r=round(pearson(x, y), 3),
        mean_bias_mriqc_minus_fmriprep=round(float(d.mean()), 3),
        sd_diff=round(float(d.std(ddof=1)), 3),
        loa_low=round(float(d.mean() - 1.96 * d.std(ddof=1)), 3),
        loa_high=round(float(d.mean() + 1.96 * d.std(ddof=1)), 3),
        median_abs_diff=round(float(np.median(np.abs(d))), 3))
json.dump(val, open(os.path.join(OUT, 'validation_summary.json'), 'w'), indent=1)

print(json.dumps(val, indent=1))
print('\nSITE SUMMARY')
g = sum_df.groupby('site')
rows = []
for s, d in g:
    m = pd.to_numeric(d.MEAN_FD_mm, errors='coerce').dropna()
    rows.append(dict(site=s, n=len(d), TR=d.TR_s.iloc[0], frames=f"{d.total_frames.min()}-{d.total_frames.max()}",
                     with_MEAN_FD=len(m),
                     MEAN_FD_median=round(m.median(), 3) if len(m) else '',
                     MEAN_FD_min=round(m.min(), 3) if len(m) else '',
                     MEAN_FD_max=round(m.max(), 3) if len(m) else '',
                     retained03_median=int(d.retained_frames_0p3.median()),
                     usable_min03_median=round(d.usable_duration_min_0p3.median(), 2),
                     usable_min03_min=round(d.usable_duration_min_0p3.min(), 2),
                     n_usable_lt_5min_03=int((d.usable_duration_min_0p3 < 5).sum()),
                     n_usable_lt_4min_03=int((d.usable_duration_min_0p3 < 4).sum())))
print(pd.DataFrame(rows).to_string(index=False))
