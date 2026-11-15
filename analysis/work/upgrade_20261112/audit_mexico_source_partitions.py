"""Stream annual source keys and bind cached admissions to source partitions.

Identifiers remain in memory. Persist aggregate QA counts. Exclude blank
lines and exact SQL-export count trailers, rejecting other malformed rows.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import hashlib
import json
import re
import sys

sys.path.insert(0, 'F:/analysis/climate_biliary_20261102/work/pylibs')
import duckdb

BASE = Path(__file__).resolve().parents[2]
ROOT = Path('D:/数据库集合/_transfer/mx_stage')
OUT = BASE / 'outputs/upgrade_20261112'
CACHE = Path('F:/analysis/climate_biliary_20261102/work/mexico_n20_id_proced.parquet')
FOOTER = re.compile(r'^\((\d+) row\(s\) affected\)$')
ISO_DATE = re.compile(r'^\d{4}-\d{2}-\d{2}(?:[ T].*)?$')


def sha(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def state_key(value):
    value = str(value or '').strip()
    return str(int(value)) if value.isdigit() else value


def signature(identifier, raw_date, state, diagnosis):
    return (str(identifier or '').strip(), str(raw_date or '').strip()[:10],
            state_key(state), str(diagnosis or '').strip().upper())


def main():
    db = duckdb.connect()
    cached = db.execute('SELECT ID,INGRE,ENTIDAD,AFECPRIN FROM read_parquet(?)', [str(CACHE)]).fetchall()
    targets = Counter(signature(*row) for row in cached)
    strict = {key for key in targets if key[3].startswith('N20')
              and '2018-01-01' <= key[1] <= '2022-12-31'
              and key[2].isdigit() and 1 <= int(key[2]) <= 32}
    if sum(targets.values()) != 40907 or len(strict) != 32417:
        raise ValueError('Cached cohort definition/count needs review')
    bindings = defaultdict(list)
    summaries = []
    source_manifest = []
    for year in range(2018, 2023):
        path = ROOT / (f'mx{year}_egresos.' + ('csv' if year == 2018 else 'txt'))
        delimiter = ',' if year in (2018, 2019) else '|'
        seen = set()
        count = missing = blanks = footers = iso = 0
        footer_total = None
        with path.open('r', encoding='utf-8-sig', newline='') as source:
            reader = csv.reader(source, delimiter=delimiter, strict=True)
            header = next(reader)
            indexes = {name: header.index(name) for name in ['ID','INGRE','ENTIDAD','AFECPRIN']}
            for row in reader:
                if not row or (len(row) == 1 and not row[0].strip()):
                    blanks += 1
                    continue
                trailer = FOOTER.fullmatch(row[0].strip()) if len(row) == 1 else None
                if trailer:
                    if footers:
                        raise ValueError(f'Duplicate source trailer in {year}')
                    footer_total = int(trailer.group(1))
                    footers += 1
                    continue
                if len(row) != len(header) or footers:
                    raise ValueError(f'Unexpected source structure in {year} at line {reader.line_num}')
                identifier = row[indexes['ID']].strip()
                count += 1
                missing += int(not identifier)
                if identifier:
                    seen.add(identifier)
                raw_date = row[indexes['INGRE']].strip()
                iso += int(bool(ISO_DATE.fullmatch(raw_date)))
                key = signature(identifier, raw_date, row[indexes['ENTIDAD']], row[indexes['AFECPRIN']])
                if key in targets:
                    bindings[key].append(year)
        if footer_total is not None and footer_total != count:
            raise ValueError(f'SQL-export trailer total mismatch in {year}')
        entry = {'year': year, 'file': path.name, 'records': count, 'unique_ids': len(seen),
                 'missing_ids': missing, 'duplicate_rows_beyond_unique': count - missing - len(seen),
                 'blank_lines_excluded': blanks, 'sql_export_trailers_excluded': footers,
                 'trailer_record_total': footer_total, 'ISO_INGRE_records': iso}
        summaries.append(entry)
        source_manifest.append({'file': str(path), 'bytes': path.stat().st_size, 'sha256': sha(path)})
        print(json.dumps(entry), flush=True)
    unmatched = sum(targets[key] for key in targets if not bindings[key])
    ambiguous = sum(targets[key] for key in targets if len(bindings[key]) > 1)
    strict_unmatched = sum(key not in bindings for key in strict)
    strict_ambiguous = sum(len(bindings[key]) > 1 for key in strict)
    cross = Counter()
    for key in strict:
        years = bindings[key]
        if len(years) == 1 and years[0] != int(key[1][:4]):
            cross[(years[0], int(key[1][:4]))] += 1
    partition = {
        'cached_rows': len(cached), 'unmatched_cached_rows': unmatched, 'ambiguous_cached_rows': ambiguous,
        'strict_valid_state_admissions': len(strict), 'strict_unmatched': strict_unmatched,
        'strict_ambiguous': strict_ambiguous,
        'strict_admission_year_differs_from_source_partition': sum(cross.values()),
        'cross_year_counts': [{'source_file_year': source_year, 'admission_year': admission_year,
                               'admissions': count} for (source_year, admission_year), count in sorted(cross.items())],
        'matching_signature': 'ID + ISO admission date + state + principal diagnosis, in memory only',
        'existing_link_rule': 'procedure source file selected by admission date year'}
    keys_passed = all(row['missing_ids'] == 0 and row['duplicate_rows_beyond_unique'] == 0 for row in summaries)
    partitions_passed = strict_unmatched == 0 and strict_ambiguous == 0 and not cross
    report = {'passed': keys_passed and partitions_passed,
              'global_annual_ID_uniqueness_passed': keys_passed,
              'existing_procedure_partition_rule_passed': partitions_passed, 'years': summaries,
              'source_partition_audit': partition, 'sources': source_manifest, 'cache_sha256': sha(CACHE)}
    with (OUT / 'mexico_specific_global_key_audit.csv').open('w', encoding='utf-8-sig', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    (OUT / 'mexico_specific_global_key_audit.json').write_text(
        json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    print('SOURCE_KEY_AND_PARTITION_AUDIT', json.dumps(partition), flush=True)
    if not report['passed']:
        raise ValueError('Procedure linkage source-year rule requires correction before release')
    print('ANNUAL_ADMISSION_KEY_AND_PARTITION_VERIFIED', flush=True)


if __name__ == '__main__':
    main()
