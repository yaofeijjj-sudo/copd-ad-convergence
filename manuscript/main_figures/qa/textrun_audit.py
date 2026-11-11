import fitz, re, sys, json, collections
SHOW = re.compile(r"(\((?:\\.|[^\\)])*\)|<[0-9A-Fa-f]*>)\s*Tj|\]\s*TJ")
res = {}
for i in range(1, 9):
    d = fitz.open(f"figures/source_pdf/Fig{i}.pdf"); p = d[0]
    fonts = sorted({(f[3].split("+")[-1], f[2]) for f in p.get_fonts()})
    xo = [x for x in range(1, d.xref_length()) if d.xref_get_key(x, "Subtype")[1] == "/Form"]
    c = "".join(d.xref_stream(x).decode("latin-1") for x in p.get_contents())
    blocks = re.findall(r"\bBT\b(.*?)\bET\b", c, re.S)
    multi = []
    for b in blocks:
        n = len(SHOW.findall(b))
        if n > 1:
            txt = "".join(re.findall(r"\(((?:\\.|[^\\)])*)\)\s*Tj", b))
            multi.append((n, len(re.findall(r"Tf", b)), txt))
    sizes = collections.Counter(round(s["size"], 1) for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"] if s["text"].strip())
    res[f"Fig{i}"] = dict(fonts=fonts, form_xobjects=len(xo), text_objects=len(blocks), TJ_left=len(re.findall(r"\]\s*TJ", c)),
                          multi_show_objects=multi, sizes=dict(sizes))
    print(f"Fig{i}: fonts={fonts} Type3={'Type3' in str(fonts)} formXObj={len(xo)} textobjs={len(blocks)} TJ={res[f'Fig{i}']['TJ_left']} sizes={dict(sizes)}")
    for m in multi: print("   multi-run:", m)
json.dump(res, open("code/main_figures/qa/text_runs_audit.json", "w"), ensure_ascii=False, indent=1)
