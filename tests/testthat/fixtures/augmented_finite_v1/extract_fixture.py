#!/usr/bin/env python3
"""Deterministic compact extraction of the accepted private augmented reference.

Copies the accepted independent tables byte for byte and writes one derived
table per source directory, unit_if_analytic.csv, holding the exact original
strings of (fixture, coordinate, unit_id, analytic) for every coordinate except
the per-unit prediction.* rows and the raw M./Y. moment rows. No number is parsed,
recomputed or reformatted. The finite-difference and complex-step columns of
the source stay in the accepted host artifact (pinned in manifest.json).

usage: extract_fixture.py ACCEPTED_PRIVATE_REFERENCE_DIR OUTPUT_DIR
"""
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path

COPIED = ('values.csv', 'components.csv', 'covariance.csv', 'predictions.csv',
          'projection.csv', 'projection_if.csv', 'atoms.csv')
EXCLUDED_PREFIXES = ('prediction.', 'M.', 'Y.')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    src, out = Path(sys.argv[1]), Path(sys.argv[2])
    manifest = {'rule': ('byte copies of the accepted tables; unit_if_analytic.csv keeps the exact source '
                         'strings of fixture,coordinate,unit_id,analytic for coordinates not starting with '
                         + ', '.join(EXCLUDED_PREFIXES)),
                'sources': {}, 'files': {}}
    for sub, target in (('', 'main'), ('selection_only_v_v1', 'selection_only_v')):
        sdir = src / sub if sub else src
        tdir = out / target
        tdir.mkdir(parents=True, exist_ok=True)
        for name in COPIED:
            source = sdir / ('input/atoms.csv' if (name == 'atoms.csv' and not sub) else name)
            shutil.copyfile(source, tdir / name)
            manifest['sources'][f'{target}/{name}'] = {'source': str(source.relative_to(src)), 'sha256': sha256(source)}
        unit = sdir / 'unit_if.csv'
        manifest['sources'][f'{target}/unit_if.csv'] = {'source': str(unit.relative_to(src)), 'sha256': sha256(unit)}
        rows = 0
        with unit.open(newline='') as fin, (tdir / 'unit_if_analytic.csv').open('w', newline='') as fout:
            reader = csv.reader(fin)
            header = next(reader)
            assert header[:4] == ['fixture', 'coordinate', 'unit_id', 'analytic'], header
            fout.write('fixture,coordinate,unit_id,analytic\n')
            for rec in reader:
                if rec[1].startswith(EXCLUDED_PREFIXES):
                    continue
                fout.write(','.join(rec[:4]) + '\n')
                rows += 1
        manifest['sources'][f'{target}/unit_if_analytic.csv'] = {'derived_from': f'{target}/unit_if.csv', 'rows': rows}
    for path in sorted(list((out / 'main').rglob('*')) + list((out / 'selection_only_v').rglob('*'))):
        if path.is_file():
            rel = str(path.relative_to(out))
            data = path.read_bytes()
            manifest['files'][rel] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
                                      'md5': hashlib.md5(data).hexdigest()}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=1, sort_keys=True) + '\n')
    with (out / 'manifest.csv').open('w', newline='') as f:
        f.write('file,bytes,md5,sha256\n')
        for rel, meta in sorted(manifest['files'].items()):
            f.write(f"{rel},{meta['bytes']},{meta['md5']},{meta['sha256']}\n")


if __name__ == '__main__':
    main()
