#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Renumber in-text citations to strict Vancouver order (first-appearance order),
reorder the reference list accordingly, and convert to JTM (BMC) Vancouver format.
"""
import re

MANUSCRIPT = '../I020_manuscript_v1.0.md'
txt = open(MANUSCRIPT, encoding='utf-8').read()

head, rest = txt.split('## 7. References', 1)
ref_block, tail = rest.split('## 8.', 1)
refsec = ref_block  # between '## 7. References' and '## 8.'

# ---------- 1. parse reference entries (old numbers) ----------
entries = {}
for line in refsec.split('\n'):
    m = re.match(r'^(\d+)\.\s+(.*)', line.strip())
    if m:
        entries[int(m.group(1))] = m.group(2).strip()

# apply known corrections to reference TEXT before renumbering
# [11] DOI张冠李戴: 10.1016/j.jclinepi.2021.07.021 -> .014
entries[11] = entries[11].replace('10.1016/j.jclinepi.2021.07.021', '10.1016/j.jclinepi.2021.07.014')
# [10] full READUS-PV E&E title (already full in reference list, verify)
# no change needed if already full

# ---------- 2. collect first-appearance order of old numbers in head ----------
def expand_brackets(s):
    """expand one bracket content like '6,7' or '1-4' or '13' into list of ints"""
    out = []
    for part in s.replace('\u2013', '-').split(','):
        part = part.strip()
        if '-' in part:
            a, b = part.split('-')
            out.extend(range(int(a), int(b) + 1))
        elif part.isdigit():
            out.append(int(part))
    return out

first_order = []
for m in re.finditer(r'\[([0-9][0-9,\u2013\-\s]*)\]', head):
    for n in expand_brackets(m.group(1)):
        if n not in first_order:
            first_order.append(n)

# sanity: all 33 refs appear
assert set(first_order) == set(entries), (set(entries) - set(first_order), set(first_order) - set(entries))
assert len(first_order) == len(entries) == 33

# old -> new map
old2new = {old: i + 1 for i, old in enumerate(first_order)}
new2old = {v: k for k, v in old2new.items()}

# ---------- 3. renumber in-text citations in head ----------
def renumber_bracket(m):
    nums = expand_brackets(m.group(1))
    newnums = sorted(old2new[n] for n in nums)
    # compress consecutive runs into ranges
    parts = []
    i = 0
    while i < len(newnums):
        j = i
        while j + 1 < len(newnums) and newnums[j + 1] == newnums[j] + 1:
            j += 1
        if j - i >= 2:  # 3+ consecutive -> range
            parts.append(f'{newnums[i]}\u2013{newnums[j]}')
        elif j - i == 1:
            parts.append(str(newnums[i])); parts.append(str(newnums[j]))
        else:
            parts.append(str(newnums[i]))
        i = j + 1
    return '[' + ','.join(parts) + ']'

head_new = re.sub(r'\[([0-9][0-9,\u2013\-\s]*)\]', renumber_bracket, head)

# ---------- 4. build new JTM Vancouver reference list ----------
def to_jtm_vancouver(raw):
    """Convert current reference text to JTM/BMC Vancouver style.
    Keep it faithful to verified metadata; normalize italics and punctuation.
    """
    s = raw
    # remove markdown italics on journal names (Vancouver no italics)
    s = s.replace('*', '')
    # remove PMID tokens (not part of JTM Vancouver; kept in verification log)
    s = re.sub(r'\s*PMID\s+\d+\.?', '', s)
    # normalize doi: no trailing period handled by adding one at end
    s = re.sub(r'doi:\s*', 'doi:', s)
    s = s.strip()
    if not s.endswith('.'):
        s += '.'
    return s

new_lines = []
for new in range(1, 34):
    old = new2old[new]
    new_lines.append(f'{new}. {to_jtm_vancouver(entries[old])}')

new_refsec = '\n\n' + '\n'.join(new_lines) + '\n\n'

# ---------- 5. write back ----------
out = head_new + '## 7. References' + new_refsec + '## 8.' + tail
open(MANUSCRIPT, 'w', encoding='utf-8').write(out)

print('=== old -> new mapping (first appearance order) ===')
for old in first_order:
    print(f'  old {old:2d} -> new {old2new[old]:2d}')
print('\nRenumbered. New first-appearance order should be 1..33.')
# verify
body2 = out.split('## 7. References')[0]
seq2 = []
for m in re.finditer(r'\[([0-9][0-9,\u2013\-\s]*)\]', body2):
    for n in expand_brackets(m.group(1)):
        if n not in seq2: seq2.append(n)
print('New first-appearance order:', seq2)
print('Strict Vancouver 1..33:', seq2 == list(range(1, 34)))
