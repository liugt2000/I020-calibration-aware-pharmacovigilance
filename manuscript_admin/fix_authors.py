#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Complete author lists to ICMJE/JTM convention (first 6 + et al.) using verified PubMed data.
Mapping new->old number is fixed from the renumbering step.
"""
import re, json

MANUSCRIPT = '../I020_manuscript_v1.0.md'
txt = open(MANUSCRIPT, encoding='utf-8').read()
head, rest = txt.split('## 7. References', 1)
ref_block, tail = rest.split('## 8.', 1)

# old->new as produced by renumber_refs.py (first appearance order)
old2new = {13:1,14:2,15:3,16:4,1:5,2:6,3:7,4:8,24:9,5:10,6:11,7:12,8:13,25:14,
           9:15,10:16,26:17,27:18,28:19,18:20,19:21,20:22,21:23,22:24,30:25,31:26,
           29:27,17:28,33:29,11:30,12:31,23:32,32:33}
new2old = {v:k for k,v in old2new.items()}

d = json.load(open('ref_verification_v3.json', encoding='utf-8'))
old_data = {r['num']: (r['pubmed'] or {}) for r in d['results']}

def norm(s): return re.sub(r'[^a-z0-9]', '', (s or '').lower())

def fmt_authors(authors):
    if not authors: return None
    if len(authors) <= 6:
        return ', '.join(authors) + '. '
    return ', '.join(authors[:6]) + ', et al. '

new_lines = []
changed = 0
for line in ref_block.strip().split('\n'):
    m = re.match(r'^(\d+)\.\s+(.*)', line.strip())
    if not m:
        new_lines.append(line); continue
    num = int(m.group(1)); content = m.group(2).strip()
    old = new2old.get(num)
    p = old_data.get(old, {})
    authors = p.get('authors', [])
    ptitle = p.get('title', '')
    if not authors or not ptitle:
        print(f'#{num}(old{old}): no data, keep'); new_lines.append(f'{num}. {content}'); continue
    # locate title in content (case-insensitive, allow bullet dashes)
    idx = norm(content).find(norm(ptitle)[:45])
    if idx == -1:
        print(f'#{num}(old{old}): TITLE NOT FOUND -> keep as-is'); new_lines.append(f'{num}. {content}'); continue
    # find corresponding char boundary in original content: scan to match normalized offset
    count = 0; charidx = 0
    while charidx < len(content) and count < idx:
        if re.match(r'[a-z0-9]', content[charidx], re.I): count += 1
        charidx += 1
    title_part = content[charidx:].lstrip('. ').strip()
    newcontent = fmt_authors(authors) + title_part
    new_lines.append(f'{num}. {newcontent}')
    changed += 1

out = head + '## 7. References\n\n' + '\n'.join(new_lines) + '\n\n## 8.' + tail
open(MANUSCRIPT, 'w', encoding='utf-8').write(out)
print(f'\nChanged {changed} entries.')
