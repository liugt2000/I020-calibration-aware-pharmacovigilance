#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reference verification v3 — corrected author parsing, Crossref 429 retry."""
import re, json, time, urllib.request, urllib.parse

MANUSCRIPT = '../I020_manuscript_v1.0.md'

def get_json(url, timeout=30, retries=3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'ref-verifier/3.0 (mailto:research@example.com)'})
            return json.load(urllib.request.urlopen(req, timeout=timeout))
        except urllib.error.HTTPError as ex:
            if ex.code == 429 and attempt < retries - 1:
                time.sleep(2.0 * (attempt + 1)); continue
            raise
        except Exception:
            if attempt < retries - 1:
                time.sleep(1.0); continue
            raise

txt = open(MANUSCRIPT, encoding='utf-8').read()
refsec = txt.split('## 7. References')[1].split('## 8.')[0]
entries = {}
for line in refsec.split('\n'):
    m = re.match(r'^(\d+)\.\s+(.*)', line.strip())
    if not m: continue
    num = int(m.group(1)); raw = m.group(2)
    pmid = re.search(r'PMID\s+(\d+)', raw)
    doi_m = re.search(r'doi:(\S+)', raw)
    doi_val = doi_m.group(1).rstrip('.') if doi_m else None
    # author segment = text before first '. ' that follows author list; take up to first period-space
    author_seg = raw.split('. ')[0]
    entries[num] = {'raw': raw, 'pmid': pmid.group(1) if pmid else None,
                    'doi': doi_val, 'author_seg': author_seg}

pmids = [e['pmid'] for e in entries.values() if e['pmid']]
pm_data = {}
for i in range(0, len(pmids), 50):
    batch = pmids[i:i+50]
    u = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=' + ','.join(batch)
    d = get_json(u)['result']
    for pid in d.get('uids', []):
        r = d[pid]
        doi = next((a['value'] for a in r.get('articleids', []) if a.get('idtype') == 'doi'), '')
        pm_data[pid] = {
            'title': r.get('title', ''), 'journal': r.get('fulljournalname', ''),
            'journal_iso': r.get('source', ''), 'volume': r.get('volume', ''),
            'issue': r.get('issue', ''), 'pages': r.get('pages', ''),
            'pubdate': r.get('pubdate', ''), 'doi': doi,
            'authors': [a['name'] for a in r.get('authors', [])],
        }
    time.sleep(0.5)

def norm(s): return re.sub(r'[^a-z0-9]', '', (s or '').lower())

results = []
for num in sorted(entries):
    e = entries[num]; status = 'OK'; issues = []
    p = pm_data.get(e['pmid'], {}) if e['pmid'] else {}
    if e['pmid'] and not p:
        status = 'CRITICAL'; issues.append('PMID not in PubMed')
    # DOI vs PubMed authoritative
    if p.get('doi') and e['doi'] and norm(p['doi']) != norm(e['doi']):
        status = 'CRITICAL'
        issues.append(f"DOI MISMATCH ms={e['doi']} pubmed={p['doi']}")
    elif e['doi'] and not p.get('doi'):
        # no PMID doi in pubmed; try crossref resolvability
        pass
    # first author: PubMed surname should appear in manuscript author segment
    if p.get('authors'):
        psurname = p['authors'][0].split()[-1]
        if norm(psurname) not in norm(e['author_seg']):
            status = 'CRITICAL'
            issues.append(f"AUTHOR MISMATCH pubmed_surname={psurname} ms_author_seg={e['author_seg'][:40]}")
    # title core check: first 30 normalized chars of PubMed title in raw
    if p.get('title'):
        if norm(p['title'])[:30] not in norm(e['raw']):
            if status != 'CRITICAL': status = 'WARNING'
            issues.append('TITLE fragment not matched')
    # page/volume check
    if p.get('pages') and p.get('volume'):
        vol = re.search(r'\b'+re.escape(p['volume'])+r'\b', e['raw'])
        if not vol:
            if status != 'CRITICAL': status = 'WARNING'
            issues.append(f"VOLUME {p['volume']} not in manuscript entry")
    results.append({'num': num, 'status': status, 'issues': issues,
                    'manuscript': e, 'pubmed': p or None})

out = {'summary': {}, 'results': results}
from collections import Counter
out['summary'] = dict(Counter(r['status'] for r in results))
json.dump(out, open('ref_verification_v3.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=2)

print('Summary:', out['summary'])
for r in results:
    if r['status'] != 'OK':
        print(r['num'], r['status'], '|', '; '.join(r['issues']))
print('\n=== DOI mismatches detail ===')
for r in results:
    for iss in r['issues']:
        if 'DOI MISMATCH' in iss:
            pm = r['pubmed'] or {}
            print(f"[{r['num']}] MS DOI: {r['manuscript']['doi']}")
            print(f"      PubMed DOI: {pm.get('doi')}")
            print(f"      PubMed title: {pm.get('title')}")
