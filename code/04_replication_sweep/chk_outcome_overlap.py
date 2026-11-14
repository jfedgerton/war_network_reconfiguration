# chk_outcome_overlap.py  [2026-11-14]
# Outcome-overlap check used to assign pool_status "outcome_overlap" (see DECISIONS_20261114.md).
# For each job "id|<data file relative to data/raw/replication_sweep/>|<outcome variable>", compares the outcome's
# event years with COW dyadic MID 4.02: share of event years that are MID onset years vs share that are MID years.
# An outcome codes continuing dispute years when the onset share is well below the MID-year share (ratio < ~0.8).
# Usage (from the repository root): python3 code/04_replication_sweep/chk_outcome_overlap.py "101|work/101/data.dta|dispute1"
import pandas as pd, numpy as np, sys, os
R=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','data','raw','replication_sweep','')
cow=pd.read_csv(R+'work_p4/X090/dymid/dyadic_mid_4.02.csv',usecols=['statea','stateb','strtyr','year','disno'])
cow['a']=cow[['statea','stateb']].min(1); cow['b']=cow[['statea','stateb']].max(1)
ons=set(zip(cow.a,cow.b,cow.strtyr)); inv=set(zip(cow.a,cow.b,cow.year))
CAND=[('statea','stateb'),('ccode1','ccode2'),('state','target'),('ccodea','ccodeb'),('cowcode1','cowcode2'),('ccode_a','ccode_b'),('side_a','side_b'),('challenger','target'),('ccode_1','ccode_2'),('cowa','cowb'),('statenum1','statenum2'),('chal','tgt')]
def load(f,cols=None):
    p=R+f
    if f.endswith('.dta'):
        try: d=pd.read_stata(p,convert_categoricals=False,convert_dates=False)
        except ValueError:
            import pyreadstat; d,_=pyreadstat.read_dta(p)
        return d
    if f.endswith('.rds'):
        import pyreadr; return list(pyreadr.read_r(p).values())[0]
    if f.endswith('.tab'): return pd.read_csv(p,sep='\t')
    return pd.read_csv(p,sep=None,engine='python')
JOBS=[l.split('|') for l in sys.argv[1:]]
for jid,f,dv in JOBS:
    try:
        d=load(f)
        keys=next((c for c in CAND if c[0] in d.columns and c[1] in d.columns),None)
        if 'Year' in d.columns and 'year' not in d.columns: d=d.rename(columns={'Year':'year'})
        if keys is None and 'dyad' in d.columns and 'dyadid' not in d.columns: d['dyadid']=d['dyad']
        if keys is None and 'dyadid' in d.columns:
            d=d[d.dyadid.notna()].copy(); dd=d.dyadid.astype('int64'); d['_a']=dd//1000; d['_b']=dd%1000; keys=('_a','_b')
        if keys is None and 'ccode' in d.columns and 'year' in d.columns:
            d=d[d[dv].notna()].copy(); d=d.sort_values(['ccode','year']).drop_duplicates(['ccode','year'])
            d['prev']=d.groupby('ccode')[dv].shift(1); d['py']=d.groupby('ccode').year.shift(1)
            o=d[(d[dv]==1)&(d.py==d.year-1)]; one=d[d[dv]==1]
            so=set(zip(cow.statea,cow.strtyr))|set(zip(cow.stateb,cow.strtyr)); si=set(zip(cow.statea,cow.year))|set(zip(cow.stateb,cow.year))
            r=[jid,dv,'MONADIC',f'n={len(d)}',f'cont={(o.prev==1).mean():.2f}']
            for lab,S in [('ons',so),('inv',si)]:
                for k in (0,1): r.append(f'{lab}t+{k}={(pd.Series(list(zip(one.ccode.astype(int),one.year.astype(int)+k))).isin(S)).mean():.2f}')
            print(' '.join(r)); continue
        if dv not in d.columns: print(jid,'NO DV',dv,[c for c in d.columns][:30]); continue
        y='year' if 'year' in d.columns else None
        d=d[d[dv].notna()].copy()
        if keys: d=d.dropna(subset=list(keys)+([y] if y else [])); d[dv]=d[dv].astype(float)
        out=[jid,dv,f'n={len(d)}',f'ones={int((d[dv]==1).sum())}']
        if keys and y:
            u=d[list(keys)].min(1).astype(int).astype(str)+'_'+d[list(keys)].max(1).astype(int).astype(str) if True else None
            d['u']=u; d=d.sort_values(['u',y]).drop_duplicates(['u',y])
            d['prev']=d.groupby('u')[dv].shift(1); d['py']=d.groupby('u')[y].shift(1)
            o=d[(d[dv]==1)&(d.py==d[y]-1)]
            out.append(f'cont={(o.prev==1).mean():.2f}')
            one=d[d[dv]==1]; a=one[list(keys)].min(1).astype(int); b=one[list(keys)].max(1).astype(int); t=one[y].astype(int)
            for lab,S in [('ons',ons),('inv',inv)]:
                for k in (0,1):
                    out.append(f'{lab}t+{k}={np.mean([(x,z,w+k) in S for x,z,w in zip(a,b,t)]):.2f}')
        else: out.append(f'keys={keys} year={y}')
        print(' '.join(out))
    except Exception as e: print(jid,'ERR',repr(e)[:150])
