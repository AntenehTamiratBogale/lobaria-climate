"""Check the supplied companion datasets using only Python's standard library."""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def read(name):
    with (ROOT / 'data' / name).open(encoding='latin1', newline='') as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)

columns, env = read('data_env_prot.csv')
assert len(env) == 42 and len(columns) == 237
assert len({r['sites_names'] for r in env}) == 42
assert len({r['sites'] for r in env}) == 14
for row in env:
    assert all(math.isfinite(float(row[c])) and float(row[c]) >= 0 for c in columns[20:])
_, all_proteins = read('all_protein_limma_output.csv')
_, selected = read('DEPs_limma_output.csv')
assert len(all_proteins) == 217 and len(selected) == 73
significant = [r for r in all_proteins if float(r['adj.P.Val']) < .05]
assert len(significant) == 73
# Display names differ slightly between supplied tables; compare numerical records.
assert sorted((float(r['logFC']), float(r['adj.P.Val'])) for r in selected) == sorted((float(r['logFC']), float(r['adj.P.Val'])) for r in significant)
a_cols, abundance = read('tax_abundance_data.csv')
t_cols, taxonomy = read('tax_class_data.csv')
_, metadata = read('meta_data.csv')
assert len(abundance) == len(taxonomy) == 58
assert [r[a_cols[0]] for r in abundance] == [r[t_cols[0]] for r in taxonomy]
assert set(a_cols[1:]) == {r['Sites'] for r in metadata}
for name, count in [('algae_genotypes.csv',16), ('fungi_genotypes.csv',10)]:
    header, records = read(name)
    assert len(records) == 46 and len(header) == count
_, components = read('data_permanova_output.csv')
for category in {r['catagory'] for r in components}:
    assert abs(sum(float(r['rSquare']) for r in components if r['catagory']==category)-1) < .0001
print('PASS: dataset dimensions, identifiers, abundances, differential-protein selection, taxonomy alignment, genotype tables, and PERMANOVA component totals')
