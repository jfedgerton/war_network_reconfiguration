# apply_audit2_20261110.py -- adjudicated fixes from two independent replication audits (2026-11-10)
# Inputs: registry.csv, study_flow.csv (backups in tmp/backup_20261110_audit2/). Run from code/04_replication_sweep/.
# Rules (90_summarise.R / 94_study_flow_table.R): one headline row per id (study_headline_term Y);
# pool_status non-empty iff headline; main iff in_headline Y; pooled rows must be significant at alpha_pub.
import csv, io

def read(fn):
    raw = open(fn, 'rb').read()
    bom = raw.startswith(b'\xef\xbb\xbf')
    eol = '\r\n' if b'\r\n' in raw else '\n'
    rows = list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))))
    return rows, list(rows[0].keys()), bom, eol

def write(fn, rows, fields, bom, eol):
    out = io.StringIO(); w = csv.DictWriter(out, fieldnames=fields, lineterminator=eol)
    w.writeheader(); w.writerows(rows)
    b = out.getvalue().encode('utf-8')
    open(fn, 'wb').write((b'\xef\xbb\xbf' + b) if bom else b)

reg, RF, rbom, reol = read('registry.csv')
R = {r['row_id']: r for r in reg}
TAG = '[2026-11-10 audit2] '

def note(r, txt):
    r['methods_note'] = (TAG + txt + (' | ' + r['methods_note'] if r['methods_note'] else ''))

# 1. promoted headline rows from the first 2026-11-10 fix: fill study-level author_history, interaction_hand
for new, old, ih in [('25', '19', 'Y'), ('57', '26', 'N'), ('116', '175', 'N'), ('168', '88', 'N')]:
    assert R[new]['study_headline_term'] == 'Y' and R[new]['author_history'] == ''
    R[new]['author_history'] = R[old]['author_history']; R[new]['interaction_hand'] = ih

# 2. headline term switches (both audits agree unless marked judgment)
def switch(old, new, why, family=None):
    o, n = R[old], R[new]
    assert o['id'] == n['id'] and o['study_headline_term'] == 'Y' and o['pool_status'] == 'main', (old, new)
    for k in ('author_history', 'alpha_pub'):
        n[k] = o[k]
    n['interaction_hand'] = n['interaction_hand'] or 'N'
    o['study_headline_term'] = 'N'; o['in_headline'] = 'N'; o['pool_status'] = ''
    n['study_headline_term'] = 'Y'; n['in_headline'] = 'Y'; n['pool_status'] = 'main'
    if family: n['family'] = family
    note(n, why)

switch('217', '194', 'headline is the dictatorial peace (single-party), Table 1; democ is a control')
switch('48', '68', 'bothdem not in published tables; Table 1 Model 4 lowest democracy -0.030 (0.007)')
switch('72', '38', 'pooled term was leader tenure; democracy targeting term defpol4 (model without tenure)')
switch('22', '151', 'abstract: democracy reduces MIDs consistently, capitalism does not; demlo; sun2cati (UN affinity) was misfiled', family='DEMOCRACY')
switch('222', '224', 'all MIDs is the primary outcome; p .078 meets the paper\'s .10 threshold')
switch('156', '73', 'PolityJoint product (Table 3) is the primary democracy measure')
switch('232', '206', 'lead claim is the developmental peace (judgment)', family='CAPITALISM')

# 3. family corrections, same term
R['240']['family'] = 'DEMOCRACY'; note(R['240'], 'joint-democracy term; alliance terms are null in the article -> DEMOCRACY')
R['56']['family'] = 'CAPITALISM'; note(R['56'], 'contract-intensive economy is a CAPITALISM term')

# 4. N08: article claim is that only joint democracy survives; pool that term -> not significant -> leaves pool
switch('257', '264', 'article claim: only joint democracy survives once dependencies are modelled', family='DEMOCRACY')
R['264']['in_headline'] = 'N'; R['264']['pool_status'] = 'null_common_sample'

# 5. out of scope / duplicate (headline row leaves the pool)
OOS = {'267': 'N01 term built from the MID network itself (two-star); methods paper',
       '263': 'N02 term built from the MID network itself (k-star); methods paper',
       '277': 'P9_012 intervention in civil wars',
       '292': 'P9_036 support for a civil-conflict party',
       '280': 'P9_016 joining ongoing wars',
       '283': 'P9_023 MID duration outcome',
       '320': 'P9_065 crisis termination hazard',
       '297': 'P9_041 event-data (IDEA) outcome',
       '140': '104 unit is the territorial claim (cross-section), not dyad/state-year (judgment)',
       '275': 'P9_008 refugee flows are not a listed interstate tie (judgment)',
       '87':  'X088 headline is a sequence of steps to war, not one coefficient'}
for rid, why in OOS.items():
    r = R[rid]; assert r['study_headline_term'] == 'Y' and r['pool_status'] == 'main', rid
    r['in_headline'] = 'N'; r['pool_status'] = 'out_of_scope'; note(r, 'out of scope: ' + why)
for rid, why in {'183': 'C7 qualitative paper; model duplicates C1 (Gartzke 2007 Model 5, same data)',
                 '137': '99 same article and model as C6 (Gartzke & Hewitt 2010); headline is the capitalist peace term in C6'}.items():
    r = R[rid]; assert r['study_headline_term'] == 'Y' and r['pool_status'] == 'main', rid
    r['in_headline'] = 'N'; r['pool_status'] = 'duplicate'; note(r, 'duplicate: ' + why)

# 5b. earlier 2026-11-10 exclusions left pool_status blank on headline rows (90_summarise.R would fail)
for rid, st in [('99', 'duplicate'), ('252', 'out_of_scope')]:
    assert R[rid]['study_headline_term'] == 'Y' and R[rid]['pool_status'] == '' and R[rid]['in_headline'] == 'N'
    R[rid]['pool_status'] = st

# 5c. P9_078: article uses GEE AR(1) (fn 12); script switched from pooled logit (rebuild still inside 0.5-2x rule)
assert R['328']['id'] == 'P9_078'
R['328']['estimator'] = 'GEE logit, AR(1) working corr. by calendar year, robust SE by dyad (article Table 4, fn 12)'
note(R['328'], 'estimator switched from pooled logit to GEE AR(1) to match the article; full-sample rebuild -0.0352 vs published -0.0201 (1.75x)')

# 6. bibliographic / published_b fills (minor)
for rid, b in [('133', '14.597'), ('35', '-0.004'), ('52', '-0.0130')]:
    assert R[rid]['published_b'] == ''; R[rid]['published_b'] = b
R['182']['short_cite'] = 'Bak 2018'; R['182']['year'] = '2018'
R['308']['short_cite'] = 'Smith & Vreeland 2023'

# checks
hr = [r for r in reg if r['study_headline_term'] == 'Y']
ids = [r['id'] for r in hr]; assert len(ids) == len(set(ids)), 'one headline row per id'
assert all((r['pool_status'] != '') == (r['study_headline_term'] == 'Y') for r in reg)
assert all((r['pool_status'] == 'main') == (r['in_headline'] == 'Y') for r in hr)
assert all(r['author_history'] in ('A_dyadic', 'A_monadic', 'B_duration', 'C_none') for r in hr if r['pool_status'] == 'main')
assert all(r['interaction_hand'] in ('Y', 'N') for r in hr if r['pool_status'] == 'main')
write('registry.csv', reg, RF, rbom, reol)
main = [r for r in hr if r['pool_status'] == 'main']
from collections import Counter
print('registry pooled rows', len(main), Counter(r['family'] for r in main))

# ---------------- study_flow.csv ----------------
sf, SF, sbom, seol = read('study_flow.csv')
def sfrow(sid):
    m = [s for s in sf if sid in s['sweep_id'].split(';') or s['sweep_id'] == sid]
    assert len(m) == 1, (sid, len(m)); return m[0]
def setsf(sid, status, why, family=None):
    s = sfrow(sid); s['final_status'] = status
    s['reason'] = TAG + why + (' | ' + s['reason'] if s['reason'] else '')
    if family: s['family'] = family
for sid, fam in [('X085', 'DEMOCRACY'), ('X043', 'CAPITALISM'), ('P8_031', 'CAPITALISM'), ('C2', 'DEMOCRACY')]:
    s = sfrow(sid); s['family'] = fam; s['reason'] = TAG + 'family corrected to ' + fam + ' | ' + s['reason']
setsf('N08', 'excluded_not_significant_baseline', 'headline joint democracy term not significant at baseline (p .46)', 'DEMOCRACY')
for sid, rid in [('N01', '267'), ('N02', '263'), ('P9_012', '277'), ('P9_036', '292'), ('P9_016', '280'), ('P9_023', '283'),
                 ('P9_065', '320'), ('P9_041', '297'), ('104', '140'), ('P9_008', '275'), ('X088', '87')]:
    setsf(sid, 'excluded_out_of_scope', OOS[rid])
setsf('C7', 'duplicate', 'qualitative paper; model duplicates C1')
s = sfrow('103'); s['year'] = '2018'
write('study_flow.csv', sf, SF, sbom, seol)
print('study_flow pooled', sum(s['final_status'] == 'pooled' for s in sf),
      Counter(s['family'] for s in sf if s['final_status'] == 'pooled'))
