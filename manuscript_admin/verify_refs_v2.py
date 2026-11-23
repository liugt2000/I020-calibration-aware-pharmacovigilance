#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reference verification per nature-ref-verifier workflow.
Step 2.0: DOI resolvability via Crossref; Step 2/3: PubMed authoritative field comparison.
Reads the 33 refs from I020_manuscript_v1.0.md, emits JSON + markdown report.
"""
import re, json, time, urllib.request, urllib.parse

MANUSCRIPT = '../I020_manuscript_v1.0.md'

def get_json(url, timeout=30):
    req = urllib.request.Request(url, headers={'User-Agent': 'ref-verifier/2.0 (mailto:research@example.com)'})
    return json.load(urllib.request.urlopen(req, timeout=timeout))

# ---------- parse refs from manuscript ----------
txt = open(MANUSCRIPT, encoding='utf-8').read()
refsec = txt.split('## 7. References')[1].split('## 8.')[0]
entries = {}
for line in refsec.split('\n'):
    m = re.match(r'^(\d+)\.\s+(.*)', line.strip())
    if not m: continue
    num = int(m.group(1)); raw = m.group(2)
    pmid = re.search(r'PMID\s+(\d+)', raw)
    doi = re.search(r'doi:([^\s.]+(?:\.[^\s]+)*)', raw)
    # cleaner doi: capture up to whitespace or trailing period
    doi_m = re.search(r'doi:(\S+)', raw)
    doi_val = doi_m.group(1).rstrip('.') if doi_m else None
    first_author = raw.split('.')[0].strip()
    entries[num] = {
        'raw': raw,
        'pmid': pmid.group(1) if pmid else None,
        'doi': doi_val,
        'first_author': first_author,
    }

print('Parsed refs:', len(entries))

# ---------- PubMed authoritative metadata ----------
pmids = [e['pmid'] for e in entries.values() if e['pmid']]
pm_data = {}
for i in range(0, len(pmids), 50):
    batch = pmids[i:i+50]
    u = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=' + ','.join(batch)
    try:
        d = get_json(u)['result']
        for pid in d.get('uids', []):
            r = d[pid]
            doi = ''
            for aid in r.get('articleids', []):
                if aid.get('idtype') == 'doi':
                    doi = aid.get('value', '')
            auth = r.get('authors', [])
            pm_data[pid] = {
                'title': r.get('title', ''),
                'journal': r.get('fulljournalname', ''),
                'journal_iso': r.get('source', ''),
                'volume': r.get('volume', ''),
                'issue': r.get('issue', ''),
                'pages': r.get('pages', ''),
                'pubdate': r.get('pubdate', ''),
                'doi': doi,
                'authors': [a['name'] for a in auth],
            }
    except Exception as ex:
        print('PubMed batch error:', ex)
    time.sleep(0.4)

# ---------- Crossref DOI resolvability ----------
def crossref_check(doi):
    if not doi: return {'status': 'no-doi'}
    u = 'https://api.crossref.org/works/' + urllib.parse.quote(doi, safe='')
    try:
        m = get_json(u, timeout=20)['message']
        return {
            'status': 'resolved',
            'title': (m.get('title') or [''])[0],
            'container': (m.get('container-title') or [''])[0],
            'volume': m.get('volume', ''),
            'page': m.get('page', ''),
            'year': (m.get('published', {}).get('date-parts') or [['']])[0][0],
        }
    except Exception as ex:
        code = getattr(ex, 'code', None)
        return {'status': 'FAILED', 'error': str(ex)[:80]}

# ---------- field comparison ----------
def norm(s):
    if not s: return ''
    return re.sub(r'[^a-z0-9]', '', s.lower())

results = []
for num in sorted(entries):
    e = entries[num]
    status = 'OK'
    issues = []
    p = pm_data.get(e['pmid'], {}) if e['pmid'] else {}
    # PMID resolution
    if e['pmid'] and not p:
        status = 'CRITICAL'; issues.append('PMID not found in PubMed')
    # DOI crossref
    cr = crossref_check(e['doi']) if e['doi'] else {'status': 'no-doi'}
    if cr['status'] == 'FAILED':
        if status != 'CRITICAL': status = 'WARNING'
        issues.append('DOI not resolvable in Crossref: ' + cr.get('error', ''))
    # DOI vs PubMed authoritative DOI
    if p.get('doi') and e['doi'] and norm(p['doi']) != norm(e['doi']):
        status = 'CRITICAL'
        issues.append(f"DOI mismatch: manuscript={e['doi']} PubMed={p['doi']}")
    # title core-word check (manuscript title = substring after first author period, rough)
    if p.get('title'):
        mtitle_match = norm(p['title'])[:40] in norm(e['raw'])
        if not mtitle_match:
            # try Crossref title
            if cr.get('title') and norm(cr['title'])[:40] in norm(e['raw']):
                pass
            else:
                if status != 'CRITICAL': status = 'WARNING'
                issues.append('Title fragment not matched in PubMed record; verify manually')
    # first author check
    if p.get('authors'):
        pfirst = p['authors'][0].split()[-1]  # surname
        mfirst = e['first_author'].split()[-1]
        if norm(pfirst) != norm(mfirst):
            status = 'CRITICAL'
            issues.append(f"First author mismatch: manuscript={mfirst} PubMed={pfirst}")
    results.append({
        'num': num,
        'status': status,
        'issues': issues,
        'manuscript': e,
        'pubmed': p or None,
        'crossref': cr,
    })
    time.sleep(0.15)

out = {'summary': {}, 'results': results}
from collections import Counter
out['summary'] = dict(Counter(r['status'] for r in results))
json.dump(out, open('ref_verification_v2.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=2)

# ---------- markdown report ----------
lines = ['# Reference verification report (v2)', '', 'Summary: ' + str(out['summary']), '']
lines.append('| # | Status | Issues | PMID | PubMed title | DOI (manuscript / PubMed) |')
lines.append('|---|--------|--------|------|--------------|---------------------------|')
for r in results:
    pm = r['pubmed'] or {}
    iss = '; '.join(r['issues']) if r['issues'] else '-'
    doi_pm = pm.get('doi', '-')
    doi_ms = r['manuscript'].get('doi', '-')
    title = (pm.get('title', '')[:60] + '…') if pm.get('title', '') else '-'
    lines.append(f"| {r['num']} | {r['status']} | {iss} | {r['manuscript'].get('pmid') or '-'} | {title} | {doi_ms} / {doi_pm} |")
open('ref_verification_v2.md', 'w', encoding='utf-8').write('\n'.join(lines))
print('Summary:', out['summary'])
for r in results:
    if r['status'] != 'OK':
        print(r['num'], r['status'], '; '.join(r['issues']))
