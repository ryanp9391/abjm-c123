# read-only offline sanity check of the Wolfram code cells: bracket balance and newline-at-depth-0 statement separation
import json, re, sys
path = sys.argv[1]
d = json.load(open(path, encoding='utf-8'))
for i, c in enumerate(d['cells']):
    if c.get('kind') != 2:
        continue
    src = c.get('source') or c.get('value') or ''
    if isinstance(src, list):
        src = ''.join(src)
    s = re.sub(r'"(?:\\.|[^"\\])*"', '""', src)
    s = re.sub(r'\(\*.*?\*\)', '', s, flags=re.S)
    st = []
    ok = True
    bad_lines = []
    depth_at_line_end = []
    t = s
    k = 0
    line = 1
    while k < len(t):
        two = t[k:k + 2]
        ch = t[k]
        if two == '[[':
            st.append('[['); k += 2; continue
        if two == ']]' and st and st[-1] == '[[':
            st.pop(); k += 2; continue
        if ch in '([{':
            st.append(ch)
        elif ch in ')]}':
            want = {')': '(', ']': '[', '}': '{'}[ch]
            if not st or st[-1] != want:
                ok = False
                break
            st.pop()
        elif ch == '\n':
            # a newline at depth 0 must follow a ';' (statement separator) unless the line is blank
            prev = t[:k].rstrip(' ')
            if not st and prev and not prev.endswith(';'):
                bad_lines.append(line)
            line += 1
        k += 1
    print('cell', i + 1, 'brackets', 'ok' if ok and not st else 'PROBLEM %s' % st[-3:], '| depth-0 newlines without ; at lines', bad_lines)
