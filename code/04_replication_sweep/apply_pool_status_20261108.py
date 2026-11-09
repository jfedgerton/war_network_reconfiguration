# apply_pool_status_20261108.py -- one-time registry edit from the 2026-11-08 replication audit
# (claude/replication_audit_findings.md). Run from code/04_replication_sweep. Idempotent: refuses to run twice.
# Adds registry column pool_status (headline rows only) and published_check; applies the audit exclusions;
# recodes Benson 2005 alpha_pub to 0.10; relabels Bell & Long 2016 (X079) to the pooled model; updates study_flow.csv.
import csv, io, sys
D = "2026-11-08"

def read(fn):
    raw = open(fn, newline="").read()
    return list(csv.reader(io.StringIO(raw))), ("\r\n" if "\r\n" in raw else "\n")

def write(fn, rows, eol):
    with open(fn, "w", newline="") as f:
        csv.writer(f, lineterminator=eol).writerows(rows)

# ---------------- registry.csv ----------------
rows, eol = read("registry.csv")
H = rows[0]
if "pool_status" in H:
    sys.exit("registry.csv already has pool_status; nothing done")
H += ["pool_status", "published_check"]
R = [dict(zip(H, r + ["", ""])) for r in rows[1:]]
hd = {r["id"]: r for r in R if r["study_headline_term"] == "Y"}
assert len(hd) == sum(r["study_headline_term"] == "Y" for r in R)   # one headline row per id

def excl(i, status, note):
    r = hd[i]
    assert r["in_headline"] == "Y", i
    r["in_headline"] = "N"
    r["pool_status"] = status
    r["methods_note"] = (f"EXCLUDED from pooling {D}: {note} " + r["methods_note"]).strip()

excl("4", "not_reproduced", "not reproduced (full-sample strategic interdependence 1.80 vs published -0.009, not significant); specification rebuilt from the article text.")
excl("78", "rebuild_outside_rule", "protocol-B rebuild fails the round-2 match rule (baseline 0.41x published, N about twice the article's).")
excl("20", "rebuild_outside_rule", "rebuild fails the round-2 match rule (baseline ratio 0.47); GEE unstable once history is added.")
excl("X039", "rebuild_outside_rule", "rebuild fails the round-2 match rule (baseline ratio 2.09).")
excl("X017", "rebuild_outside_rule", "rebuild fails the round-2 match rule (baseline ratio 0.47; full sample 0.74x published).")
excl("P8_032", "author_null", "author argues a null (published 0.348, SE 0.270, reproduced exactly on the full sample); reported separately.")

# Benson 2005: the article treats p = .08 as significant (text calls p = 0.082 significant)
hd["78"]["alpha_pub"] = "0.10"

# Bell & Long 2016 (X079): pooled row is Table 1 Model 6 (territorial); the stored published value was our Model 2 reproduction
x = hd["X079"]
assert x["headline_model"].startswith("Table 1 Model 2")
x["headline_model"] = "Table 1 Model 6 (territorial issues)"
x["published_b"] = ""; x["published_se"] = ""
x["methods_note"] = (f"[{D}] headline_model relabeled from Model 2 to Model 6 (the pooled term); the published value previously "
                     "stored here was our Model 2 reproduction and was removed (published Model 6 value not yet entered). " + x["methods_note"]).strip()

# published values that are our own reproductions, not read from the article
for i in ["14", "32", "X079"]:
    hd[i]["published_check"] = "unverified"

# pool_status for every other headline row
for i, r in hd.items():
    if r["pool_status"]:
        continue
    n = r["methods_note"]
    if r["in_headline"] == "Y":
        r["pool_status"] = "main"
    elif "MISMATCH" in n or "fails the round-2 match rule" in n:
        r["pool_status"] = "rebuild_outside_rule"
    elif "not significant at baseline" in n:
        r["pool_status"] = "null_common_sample"
    elif "duplicate" in n:
        r["pool_status"] = "duplicate"
    elif "coup" in n or "out of scope" in n:
        r["pool_status"] = "out_of_scope"
    elif "Methods paper" in n:
        r["pool_status"] = "methods_reanalysis"
    elif i in ("X005", "P9_017"):
        r["pool_status"] = "not_comparable"
    else:
        sys.exit(f"unclassified headline row {i}: {n[:120]}")

from collections import Counter
cnt = Counter(r["pool_status"] for r in hd.values())
print(dict(cnt))
assert all((r["pool_status"] == "main") == (r["in_headline"] == "Y") for r in hd.values())
assert all(r["pool_status"] == "" for r in R if r["study_headline_term"] != "Y")
OUT_REG = ("registry.csv", [H] + [[r[k] for k in H] for r in R], eol)   # written at the end, after every check

# ---------------- alpha_pub_coding.csv ----------------
rows, eol = read("alpha_pub_coding.csv"); H2 = rows[0]
A = [dict(zip(H2, r)) for r in rows[1:]]
b = [a for a in A if a["id"] == "78"]; assert len(b) == 1
b[0]["alpha_pub"] = "0.10"; b[0]["alpha_used"] = "0.10"
b[0]["decision_note"] = (b[0]["decision_note"] + " " if b[0]["decision_note"] else "") + \
    f"{D}: recoded to 0.10; the text treats p = .078 and p = .082 as significant."
OUT_ALPHA = ("alpha_pub_coding.csv", [H2] + [[a[k] for k in H2] for a in A], eol)

# ---------------- study_flow.csv ----------------
rows, eol = read("study_flow.csv"); H3 = rows[0]
F = [dict(zip(H3, r)) for r in rows[1:]]
new = {"4": ("not_reproduced", "did not reproduce the article (audit " + D + ")"),
       "78": ("rebuilt_no_match", "rebuild fails the round-2 match rule (audit " + D + ")"),
       "20": ("rebuilt_no_match", "rebuild fails the round-2 match rule (audit " + D + ")"),
       "X039": ("rebuilt_no_match", "rebuild fails the round-2 match rule (audit " + D + ")"),
       "X017": ("rebuilt_no_match", "rebuild fails the round-2 match rule (audit " + D + ")"),
       "P8_032": ("author_null", "author argues a null; reported separately (audit " + D + ")")}
for i, (st, why) in new.items():
    m = [f for f in F if i in f["sweep_id"].split(";")]
    assert len(m) == 1 and m[0]["final_status"] == "pooled", i
    m[0]["final_status"] = st
    m[0]["reason"] = why + "; previously: " + m[0]["reason"]
OUT_FLOW = ("study_flow.csv", [H3] + [[f[k] for k in H3] for f in F], eol)
for o in (OUT_REG, OUT_ALPHA, OUT_FLOW):
    write(*o)
print("done")
