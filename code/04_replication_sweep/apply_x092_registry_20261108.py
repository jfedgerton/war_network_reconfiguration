# apply_x092_registry_20261108.py -- run from code/04_replication_sweep AFTER apply_pool_status_20261108.py.
# Updates the registry estimator text for X092 (calendar-spaced AR(1), reverted 2026-11-08). Refuses to run twice.
import csv, io, sys
raw = open("registry.csv", newline="").read(); eol = "\r\n" if "\r\n" in raw else "\n"
rows = list(csv.reader(io.StringIO(raw))); H = rows[0]; R = [dict(zip(H, r)) for r in rows[1:]]
OLD = "GEE logit, AR(1) working corr. with Stata force (order-based lags), semi-robust SE by dyad (xtgee ... corr(ar1) force robust)"
NEW = "GEE logit, AR(1) working corr. spaced by calendar year, semi-robust SE by dyad (author: xtgee ... corr(ar1) force robust)"
x = [r for r in R if r["id"] == "X092"]
if all(r["estimator"] == NEW for r in x): sys.exit("X092 already updated; nothing done")
assert len(x) == 3 and all(r["estimator"] == OLD for r in x)
for r in x:
    r["estimator"] = NEW
    r["methods_note"] = ("[2026-11-08] estimator reverted from Stata-force order-based AR(1) to calendar-spaced AR(1), which "
                         "reproduces the published estimates. " + r["methods_note"]).strip()
with open("registry.csv", "w", newline="") as f:
    csv.writer(f, lineterminator=eol).writerows([H] + [[r[k] for k in H] for r in R])
print("X092 registry rows updated")
