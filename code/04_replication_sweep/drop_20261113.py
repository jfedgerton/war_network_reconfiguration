# [2026-11-13] PI decision: drop out-of-rule studies from the pool. Usage: python3 drop_20261113.py P9_046 N07 [P8_041]
import csv, sys, shutil, os
ids = sys.argv[1:]
RULES = {'P9_046': ('out_of_scope', 'excluded_out_of_scope', 'event-data outcome (high-conflict event counts), outside the inclusion rule'),
         'N07':    ('methods_reanalysis', 'excluded_out_of_scope', 'methods paper without a significant published headline estimate'),
         'P8_041': ('not_comparable', 'excluded_pi_decision', 'pooled term (joint male-suffrage democracy, do-file) is not the printed headline coefficient of Table 5 Model 1'),
         'X060':   ('out_of_scope', 'excluded_out_of_scope', 'composite threat-tie outcome (rivalry, fatal MID, crisis, sanctions, proxy support; establish or maintain) and no printed coefficient (Figure 2 only)'),
         'X079':   ('out_of_scope', 'excluded_out_of_scope', 'outcome is IMI territorial use-of-force initiation, not MID or war'),
         'N03':    ('methods_reanalysis', 'excluded_out_of_scope', 'methods paper; GWESF- is not the headline term and the outcome is any-MID tie presence, not onset'),
         'X077':   ('not_comparable', 'excluded_pi_decision', 'pooled term is the power x preference interaction and preference is built from dyadic conflict-cooperation history'),
         '50':     ('methods_reanalysis', 'excluded_out_of_scope', 'measurement paper (identity between trade-dependence and trade-share measures); pooled term is one half of a collinear pair'),
         'P9_078': ('out_of_scope', 'excluded_out_of_scope', 'pooled column is ICB crisis involvement (ICB-only outcome)'),
         'P8_037': ('out_of_scope', 'excluded_out_of_scope', 'outcome is ICB crisis initiation (ICB-only)'),
         'C6':     ('out_of_scope', 'excluded_out_of_scope', 'outcome is ICB crisis onset (appendix replication, ICB-only)'),
         'X025':   ('out_of_scope', 'excluded_out_of_scope', 'outcome is ICB crisis initiation; pooled term from the base model'),
         'P9_057': ('out_of_scope', 'excluded_out_of_scope', 'pooled term is in the target-resistance equation only (resistance outcome)'),
         'P9_045': ('out_of_scope', 'excluded_out_of_scope', 'pooled term is in the reciprocation equation (p < .10); initiation-equation term is null'),
         'X022':   ('out_of_scope', 'excluded_out_of_scope', 'outcome is militarized compellent threats, not MID or war'),
         '52':     ('not_comparable', 'excluded_pi_decision', 'no single headline coefficient: multinomial contrasts have mixed signs; Model 1 is the base model and the Model 2 term (0->use) runs against the hypothesis and does not reproduce in the rebuild')}
assert ids and all(i in RULES for i in ids)
os.makedirs('../../tmp/backup_add20261113/pre_drop4', exist_ok=True)
for f in ('registry.csv', 'study_flow.csv', '95_study_table.R'):
    d = '../../tmp/backup_add20261113/pre_drop4/' + f
    if not os.path.exists(d): shutil.copy2(f, d)
rows = list(csv.reader(open('registry.csv', newline=''))); H = rows[0]; ix = {k: i for i, k in enumerate(H)}; n = 0
for r in rows[1:]:
    if r[ix['id']] in ids and r[ix['pool_status']] == 'main':
        r[ix['pool_status']] = RULES[r[ix['id']]][0]; r[ix['in_headline']] = 'N'
        r[ix['methods_note']] = ('[drop 2026-11-13, PI] ' + RULES[r[ix['id']]][2] + '. ' + r[ix['methods_note']]).strip(); n += 1
assert n == len(ids), n
csv.writer(open('registry.csv', 'w', newline=''), lineterminator='\n').writerows(rows)
rows = list(csv.reader(open('study_flow.csv', newline=''))); H = rows[0]; ix = {k: i for i, k in enumerate(H)}; n = 0
for r in rows[1:]:
    sid = r[ix['sweep_id']].split(';')
    for i in ids:
        if i in sid and r[ix['final_status']] == 'pooled':
            r[ix['final_status']] = RULES[i][1]; r[ix['reason']] = '2026-11-13 PI drop: ' + RULES[i][2] + ' | ' + r[ix['reason']]; n += 1
assert n == len(ids), n
csv.writer(open('study_flow.csv', 'w', newline=''), lineterminator='\n').writerows(rows)
reg = list(csv.DictReader(open('registry.csv'))); npool = sum(x['pool_status'] == 'main' for x in reg)
s = open('95_study_table.R').read()
import re
s = re.sub(r'stopifnot\(nrow\(pool\) == \d+,', 'stopifnot(nrow(pool) == %d,' % npool, s)
s = re.sub(r'stopifnot\(n_rows == \d+\)', 'stopifnot(n_rows == %d)' % npool, s)
open('95_study_table.R', 'w').write(s)
print('dropped', ids, '-> pooled', npool)
