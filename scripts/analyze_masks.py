import sys, glob, os, re, json, struct, csv, collections
sys.path.insert(0,'/tmp')
exec(open('/dev/stdin').read()) if False else None
import zlib
TYPES={1:('b',1),2:('B',1),3:('h',2),4:('H',2),5:('i',4),6:('I',4),7:('f',4),9:('d',8),12:('q',8),13:('Q',8)}
def elems(buf):
    i=0
    while i < len(buf):
        t,nb=struct.unpack('<II',buf[i:i+8])
        if t>>16:
            nb=t>>16; t=t&0xffff; body=buf[i+4:i+4+nb]; i+=8
        else:
            body=buf[i+8:i+8+nb]; i+=8+((nb+7)//8)*8
        if t==15:
            for x in elems(zlib.decompress(body)): yield x
        else: yield (t,body)
def parse_matrix(body):
    subs=list(elems(body)); cls=subs[0][1][0]
    dims=struct.unpack('<%di'%(len(subs[1][1])//4), subs[1][1])
    if cls==1: return ('cell',dims,[parse_matrix(b) for t,b in subs[3:] if t==14])
    if cls==2:
        nlen=struct.unpack('<i',subs[3][1][:4])[0]; names=subs[4][1]
        fn=[names[i:i+nlen].rstrip(b'\x00').decode() for i in range(0,len(names),nlen)]
        return ('struct',dims,{nm:parse_matrix(subs[5+i][1]) for i,nm in enumerate(fn)})
    dt,raw=subs[3]
    if dt in TYPES:
        fmt,sz=TYPES[dt]; n=len(raw)//sz; vals=list(struct.unpack('<%d%s'%(n,fmt), raw[:n*sz]))
    elif dt==2: vals=raw.decode(errors='ignore')
    else: vals=raw
    return ('num',dims,vals)
def load_mask(path):
    data=open(path,'rb').read(); top={}
    for t,body in elems(data[128:]):
        if t==14:
            subs=list(elems(body)); top[subs[2][1].decode(errors='ignore')]=parse_matrix(body)
    return top['motion_data']
def flat(x): return x[2][0][2] if False else None
# ---- ptseries lengths
import nibabel as nib
PROJ='/nas/longleaf/home/sameliza/cyberneuro/ABIDEII_project'
site={r['SUB_ID']:r['SITE_ID'] for r in csv.DictReader(open(f'{PROJ}/phenotypic/ASD Project Data(Phenotypic Data).csv',encoding='utf-8-sig')) if r['SUB_ID'].strip()}
rows=[]; sweep_issues=[]
for p in sorted(glob.glob(f'{PROJ}/func/*.ptseries.nii')):
    m=re.search(r'(sub-\d+)_(ses-[A-Za-z0-9]+baseline)_',os.path.basename(p)); sub,ses=m.group(1),m.group(2)
    sid=sub.replace('sub-','')
    pt_nt=nib.load(p).shape[0]
    raw=open(f'dcan_small/{sub}_{ses}_raw_hdr.bin','rb').read()
    raw_nt=struct.unpack('<8q',raw[16:80])[6]
    mm=load_mask(f'dcan_small/{sub}_{ses}_motion_mask.mat')
    entries=mm[2]
    thr=[]; rem=[]; fr=[]; mfd=[]; fmt=[]; total=None; skip=None; tr=None
    for e in entries:
        f=e[2]
        thr.append(float(f['FD_threshold'][2][0]));
        fr.append(f['frame_removal'][2])
        rem.append(int(f['remaining_frame_count'][2][0])); mfd.append(f['remaining_frame_mean_FD'][2][0])
        total=int(f['total_frame_count'][2][0]); skip=int(f['skip'][2][0]); tr=float(f['epi_TR'][2][0]); fmt.append(f['format_string'][2])
    n=total
    # sanity: removed count == total - remaining
    bad=[i for i,(v,r) in enumerate(zip(fr,rem)) if (n - sum(v))!=r]
    if bad: sweep_issues.append((sub,'removal-count mismatch at',bad[:3]))
    # per-frame FD bin: first threshold where frame kept (removal==0)
    fdbin=[None]*n
    for k,(t,v) in enumerate(zip(thr,fr)):
        for i,x in enumerate(v):
            if x==0 and fdbin[i] is None: fdbin[i]=t
    max_t=max(thr)
    frames_gt_max=sum(1 for x in fdbin if x is None)
    matched=[t for t,r in zip(thr,rem) if r==pt_nt]
    rows.append(dict(sub=sub,ses=ses,site=site.get(sid,'?'),raw_nt=raw_nt,mask_total=n,ptseries_nt=pt_nt,
        skip=skip,tr=tr,n_thr=len(thr),thr_min=min(thr),thr_max=max_t,
        thresholds_matching_ptseries_len=matched[:5],
        fd_bins=fdbin, frames_above_max_thr=frames_gt_max, mfd_at_0_2=None))
json.dump(rows, open('mask_summary.json','w'))
print('subjects',len(rows))
print('sweep issues:',len(sweep_issues), sweep_issues[:3])
print('TR (mask):', collections.Counter((r['site'],r['tr']) for r in rows))
print('skip (mask):', collections.Counter((r['site'],r['skip']) for r in rows))
print('thr grid sample:', [ (r['thr_min'],r['thr_max'],r['n_thr']) for r in rows[:1]])
print('raw_nt vs mask_total vs ptseries_nt by site:')
c=collections.Counter((r['site'],r['raw_nt'],r['mask_total'],r['ptseries_nt']) for r in rows)
for k,v in sorted(c.items()): print('  site=%s raw=%d mask=%d pt=%d  n=%d'%(k[0],k[1],k[2],k[3],v))
print('ptseries length matches a threshold for N subjects:', sum(1 for r in rows if r['thresholds_matching_ptseries_len']))
