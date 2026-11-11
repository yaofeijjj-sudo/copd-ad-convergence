import fitz,re,sys
STR=r'(?:\((?:\\.|[^\\)])*\)|<[0-9A-Fa-f\s]*>)'
ARR=re.compile(r'\[\s*((?:'+STR+r'\s*(?:-?\d*\.?\d+\s*)?)+)\]\s*TJ',re.S)
def merge(m):
    body=m.group(1)
    nums=[float(n) for n in re.findall(r'(?<=[)>])\s*(-?\d*\.?\d+)',body)]
    if any(abs(n)>1 for n in nums): return m.group(0)  # real spacing: keep
    parts=re.findall(STR,body,re.S)
    if all(p.startswith('(') for p in parts): return '('+''.join(p[1:-1] for p in parts)+') Tj'
    if all(p.startswith('<') for p in parts): return '<'+''.join(re.sub(r'\s','',p[1:-1]) for p in parts)+'> Tj'
    return m.group(0)
src,dst=sys.argv[1],sys.argv[2]
d=fitz.open(src)
for page in d:
    for x in page.get_contents():
        c=d.xref_stream(x).decode('latin1')
        n=len(ARR.findall(c)); c2=ARR.sub(merge,c)
        d.update_stream(x,c2.encode('latin1'))
        print('TJ before',n,'after',len(ARR.findall(c2)))
d.save(dst,garbage=0,deflate=True)
