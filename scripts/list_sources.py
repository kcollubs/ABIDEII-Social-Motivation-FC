import urllib.request, urllib.parse, re, os, json, glob
S3="https://fcp-indi.s3.amazonaws.com/"
def keys(prefix, delim=False):
    out=[]; tok=None; pre=[]
    while True:
        q={"prefix":prefix,"list-type":"2"}
        if delim: q["delimiter"]="/"
        if tok: q["continuation-token"]=tok
        url=S3+"?"+urllib.parse.urlencode(q)
        x=urllib.request.urlopen(url,timeout=180).read().decode()
        out+=re.findall(r"<Key>([^<]+)</Key>",x)
        pre+=re.findall(r"<CommonPrefixes><Prefix>([^<]+)</Prefix>",x)
        m=re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>",x)
        if not m: break
        tok=m.group(1)
    return out, pre
PROJ="/nas/longleaf/home/sameliza/cyberneuro/ABIDEII_project"
subs={}
for f in glob.glob(PROJ+"/func/*.ptseries.nii"):
    m=re.search(r"(sub-\d+)_(ses-[A-Za-z0-9]+baseline)_",os.path.basename(f))
    subs[m.group(1)]=m.group(2)
res={}
for sub in sorted(subs):
    ses=subs[sub]
    dk,_=keys(f"data/Projects/ABIDE2/Derivatives/DCAN/{sub}/{ses}/func/")
    motion=[k.split('/')[-1] for k in dk if 'motion.tsv' in k or 'motion_mask' in k]
    fk,_=keys(f"data/Projects/ABIDE2/Outputs/fmriprep/fmriprep/{sub}/")
    conf=[k for k in fk if k.endswith('desc-confounds_regressors.tsv')]
    dn,_=keys(f"data/Projects/ABIDE2/Outputs/denoise/denoise/{sub}/")
    res[sub]=dict(ses=ses, dcan_motion_files=motion, fmriprep_confounds=conf,
                  denoise_keys=[k.split('/')[-1] for k in dn][:20])
json.dump(res,open("source_inventory.json","w"),indent=1)
n_fm=sum(1 for v in res.values() if v['fmriprep_confounds'])
n_dcanmot=sum(1 for v in res.values() if v['dcan_motion_files'] and any('motion.tsv' in x for x in v['dcan_motion_files']))
n_den=sum(1 for v in res.values() if v['denoise_keys'])
print("subjects",len(res),"| fMRIPrep confounds present:",n_fm,"| DCAN *motion.tsv present:",n_dcanmot,"| denoise outputs present:",n_den)
multi={k:len(v['fmriprep_confounds']) for k,v in res.items() if len(v['fmriprep_confounds'])!=1 and v['fmriprep_confounds']}
print("subjects with !=1 confounds file:",multi)
print("example:",[ (k,v['fmriprep_confounds']) for k,v in res.items() if len(v['fmriprep_confounds'])>1][:2])
