"""
Pull the GenBank metadata for all 162 accessions in one go.

This answers the sampling-time question: for each sequence it retrieves the
collection date, the country and isolation source, any free-text note, and the
publication the submitters linked to the record. The source publication is what
actually states whether a sample was taken before treatment or after failure.

Run it anywhere with internet access. No API key is needed for 162 records,
but if you have one, set NCBI_API_KEY below and the rate limit rises.

    pip install pandas requests
    python fetch_genbank_metadata.py

Put this script in the same folder as master_162.csv. It writes
genbank_metadata.csv next to itself. Nothing else needs editing.
"""
import csv
import re
import time
import sys

import pandas as pd
import requests

SRC = 'master_162.csv'         # put this script in the same folder as the file
OUT = 'genbank_metadata.csv'
EMAIL = 'hoanghuy310592@gmail.com'   # NCBI asks for a contact address
NCBI_API_KEY = ''                  # optional
BATCH = 50

BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

# ------------------------------------------------------------------ accessions
if SRC.endswith('.xlsx'):
    df = pd.read_excel(SRC)
    col = [c for c in df.columns if 'ccession' in c][0]
else:
    df = pd.read_csv(SRC)
    col = [c for c in df.columns if 'ccession' in c.lower() or c.lower() == 'accession'][0]
accessions = [str(a).strip() for a in df[col].dropna().unique()]
print(f'{len(accessions)} accessions read from {SRC}')

# ------------------------------------------------------------------ fetch
records = []
for i in range(0, len(accessions), BATCH):
    chunk = accessions[i:i + BATCH]
    params = {'db': 'nuccore', 'id': ','.join(chunk),
              'rettype': 'gb', 'retmode': 'text', 'email': EMAIL}
    if NCBI_API_KEY:
        params['api_key'] = NCBI_API_KEY
    r = requests.get(BASE, params=params, timeout=120)
    r.raise_for_status()
    records.append(r.text)
    print(f'  fetched {i + len(chunk)}/{len(accessions)}')
    time.sleep(0.4)
text = '\n'.join(records)

# ------------------------------------------------------------------ parse
def field(block, pattern, flags=0):
    m = re.search(pattern, block, flags)
    return m.group(1).strip().replace('\n', ' ') if m else ''


rows = []
for block in text.split('\n//\n'):
    if 'LOCUS' not in block:
        continue
    acc = field(block, r'^VERSION\s+(\S+)', re.M) or field(block, r'^ACCESSION\s+(\S+)', re.M)

    # every REFERENCE in the record, with its PubMed id if present
    refs = []
    for rb in re.split(r'\nREFERENCE', block)[1:]:
        title = field(rb, r'TITLE\s+(.*?)\n\s{2,}JOURNAL', re.S)
        journal = field(rb, r'JOURNAL\s+(.*?)(?:\n\s{2,}PUBMED|\n\s{0,2}[A-Z]{4,})', re.S)
        pmid = field(rb, r'PUBMED\s+(\d+)')
        refs.append((re.sub(r'\s+', ' ', title), re.sub(r'\s+', ' ', journal), pmid))
    titles = ' | '.join(t for t, j, p in refs if t and 'Direct Submission' not in t)
    journals = ' | '.join(j for t, j, p in refs if j and not j.startswith('Submitted'))
    pmids = ','.join(p for t, j, p in refs if p)
    submitted = field(block, r'JOURNAL\s+Submitted\s+\(([^)]+)\)')

    rows.append({
        'accession': acc,
        'collection_date': field(block, r'/collection_date="([^"]+)"'),
        'country': field(block, r'/(?:country|geo_loc_name)="([^"]+)"'),
        'isolate': field(block, r'/isolate="([^"]+)"'),
        'isolation_source': field(block, r'/isolation_source="([^"]+)"'),
        'host': field(block, r'/host="([^"]+)"'),
        'note': field(block, r'/note="([^"]+)"', re.S),
        'submitted': submitted,
        'reference_title': titles,
        'reference_journal': journals,
        'pubmed': pmids,
    })

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print(f'\nwrote {OUT}: {len(out)} records')

# ------------------------------------------------------------------ summary
missing = set(a.split('.')[0] for a in accessions) - set(
    str(a).split('.')[0] for a in out['accession'])
if missing:
    print('NOT RETURNED:', sorted(missing))

print('\ndistinct source publications:')
for (t, j, p), n in out.groupby(
        ['reference_title', 'reference_journal', 'pubmed']).size().items():
    print(f'  n={n:3d}  PMID {p or "-"}  {t[:90]}')

print('\nrecords whose note mentions treatment timing:')
pat = re.compile(r'post[- ]?treat|after treat|relapse|failure|baseline|pre[- ]?treat|na[iï]ve|retreat',
                 re.I)
hits = out[out['note'].fillna('').str.contains(pat) |
           out['isolation_source'].fillna('').str.contains(pat) |
           out['isolate'].fillna('').str.contains(pat)]
print(hits[['accession', 'isolate', 'isolation_source', 'note']].to_string(index=False)
      if len(hits) else '  none - the timing is not in the record and must come from the papers')

print('\ncollection dates present for', int(out['collection_date'].astype(bool).sum()),
      'of', len(out), 'records')
