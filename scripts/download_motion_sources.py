import json, urllib.request, urllib.parse, re, os
S3="https://fcp-indi.s3.amazonaws.com/"
OUT="fmriprep_confounds"
def keys(prefix):
    out=[]; tok=None
    while True:
        q={"prefix":prefix,"list-type":"2"}
        if tok: q["continuation-token"]=tok
        x=urllib.request.urlopen(S3+"?"+urllib.parse.urlencode(q),timeout=180).read().decode()
        out+=re.findall(r"<Key>([^<]+)</Key>",x)
        m=re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>",x)
        if not m: break
        tok=m.group(1)
    return out
inv=json.load(open('source_inventory.json'))
manifest=[]
for sub in sorted(inv):
    for url in inv[sub]['fmriprep_confounds']:
        name=url.split('/')[-1]
        dst=os.path.join(OUT,name)
        if not os.path.exists(dst):
            data=urllib.request.urlopen(S3+url,timeout=300).read()
            open(dst,'wb').write(data)
        manifest.append(dict(sub=sub,file=name,bytes=os.path.getsize(dst),source='fMRIPrep desc-confounds_regressors.tsv'))
    mq=[k for k in keys(f"data/Projects/ABIDE2/Outputs/mriqc/mriqc/{sub}/") if k.endswith('_bold.json')]
    for url in mq:
        name=url.split('/')[-1]; dst=os.path.join('mriqc_json',name)
        if not os.path.exists(dst):
            open(dst,'wb').write(urllib.request.urlopen(S3+url,timeout=300).read())
        manifest.append(dict(sub=sub,file=name,bytes=os.path.getsize(dst),source='MRIQC bold JSON'))
json.dump(manifest,open('download_manifest.json','w'),indent=1)
fm=[m for m in manifest if m['source'].startswith('fMRIPrep')]
mq=[m for m in manifest if m['source'].startswith('MRIQC')]
print('fMRIPrep confounds files:',len(fm),' subjects:',len(set(m['sub'] for m in fm)),' MB:',round(sum(m['bytes'] for m in fm)/1e6,1))
print('MRIQC bold JSON files:',len(mq),' subjects:',len(set(m['sub'] for m in mq)))
