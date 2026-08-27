#!/usr/bin/env python3
"""
scripts/latex_to_pandoc.py — rewrite the Frontiers manuscript into LaTeX that
pandoc can read faithfully.

Called by scripts/make_docx.sh; see that script for why this step exists at all.

The guiding rule here is that nothing may silently disappear. Every transform
either produces equivalent text or is reported. A similarity check run against a
document that quietly dropped a paragraph would give a meaningless answer, so at
the end this asserts that the body survived roughly intact.

Usage:
    python3 scripts/latex_to_pandoc.py IN.tex IN.aux OUT.tex
"""

import os
import re
import sys


# ------------------------------------------------------------------------------
# Cross-references, resolved from the .aux LaTeX already wrote
# ------------------------------------------------------------------------------
def read_labels(aux_path):
    """label -> printed number, e.g. 'fig:crossover' -> '7'."""
    labels = {}
    with open(aux_path, encoding="utf8", errors="replace") as fh:
        for line in fh:
            m = re.match(r"\\newlabel\{([^}]*)\}\{\{+([^}{]*)\}", line)
            if m:
                labels[m.group(1)] = m.group(2)
    return labels


# ------------------------------------------------------------------------------
# siunitx
#
# Expanded here rather than left to pandoc because the manuscript declares a
# custom \bp unit that pandoc rejects, and because a wrong unit in a numeric
# result is worse than a formatting nit.
# ------------------------------------------------------------------------------
UNITS = {
    "second": "s", "milli\\second": "ms", "micro\\second": "\u00b5s",
    "gibi\\byte": "GiB", "mebi\\byte": "MiB", "kibi\\byte": "KiB",
    "giga\\byte": "GB", "mega\\byte": "MB",
    "bp": "bp", "percent": "%", "hour": "h", "minute": "min",
    "mega\\hertz": "MHz", "giga\\hertz": "GHz", "watt": "W",
}


def unit_text(raw):
    key = raw.strip().replace("\\per\\second", "/s")
    key = re.sub(r"^\\", "", key)
    key = key.replace("\\", "\\", 1)
    lookup = raw.strip().lstrip("\\")
    for k, v in UNITS.items():
        if lookup == k.lstrip("\\") or raw.strip() == "\\" + k:
            return v
    # Fall back to the bare macro name, which is right for simple units.
    return re.sub(r"\\", "", raw).strip()


def group_digits(n):
    """
    Apply the manuscript's own \\sisetup{group-minimum-digits=4}.

    Without this, \\num{15248} renders as "15248" in the .docx while the PDF says
    "15 248". The separator is a thin space, matching siunitx, so the two
    documents read identically rather than merely equivalently.
    """
    n = n.strip()
    m = re.fullmatch(r"([+-]?)(\d+)(\.\d+)?", n)
    if not m:
        return n
    sign, whole, frac = m.group(1), m.group(2), m.group(3) or ""
    if len(whole) >= 4:
        parts = []
        while len(whole) > 3:
            parts.insert(0, whole[-3:])
            whole = whole[:-3]
        parts.insert(0, whole)
        whole = "\u2009".join(parts)
    return sign + whole + frac


def prose_words(t):
    """
    Count words of real prose, ignoring markup.

    Naively counting every alphabetic token makes this check useless: expanding
    \\SI{3.11}{\\gibi\\byte} to "3.11 GiB" turns two macro-name tokens into one
    word, so a faithful conversion looks like 8% data loss. Control sequences are
    stripped from both sides first, and only words of four letters or more are
    counted, which unit abbreviations (s, ms, GiB) cannot perturb.
    """
    # Reference keys are not prose. \label{sec:crossover} would otherwise count
    # the word "crossover", and vanish once labels are stripped, making a correct
    # conversion look lossy. Drop these macros together with their arguments.
    t = re.sub(r"\\(?:label|ref|eqref|autoref|cite[a-zA-Z]*|includegraphics"
               r"|bibliography|bibliographystyle|input|include|usepackage"
               r"|documentclass|setcounter|def|newcommand)"
               r"(?:\[[^\]]*\])?(?:\{[^{}]*\})+", " ", t)
    t = re.sub(r"\\[a-zA-Z@]+", " ", t)
    t = re.sub(r"[{}$~^_&%#]", " ", t)
    return len(re.findall(r"[A-Za-z][A-Za-z'-]{3,}", t))


def braced(s, i):
    """Return (content, index_after) for a {...} group starting at s[i] == '{'."""
    assert s[i] == "{"
    depth, j = 0, i
    while j < len(s):
        if s[j] == "{":
            depth += 1
        elif s[j] == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
        j += 1
    raise ValueError("unbalanced brace")


def expand_macro(s, name, arity, render):
    """Replace every \\name{..}x arity with render(*args), innermost-safe."""
    out, i = [], 0
    tok = "\\" + name
    while True:
        k = s.find(tok, i)
        # Reject a longer macro that merely starts with this name (\num vs \numlist).
        if k == -1:
            out.append(s[i:])
            break
        after = k + len(tok)
        if after < len(s) and (s[after].isalpha()):
            out.append(s[i:after])
            i = after
            continue
        out.append(s[i:k])
        j, args = after, []
        try:
            for _ in range(arity):
                while j < len(s) and s[j] in " \n":
                    j += 1
                a, j = braced(s, j)
                args.append(a)
        except (AssertionError, ValueError, IndexError):
            out.append(tok)
            i = after
            continue
        out.append(render(*args))
        i = j
    return "".join(out)


def expand_siunitx(s):
    s = expand_macro(s, "SIrange", 3,
                     lambda a, b, u: (f"{group_digits(a)}\u2013{group_digits(b)}"
                                      f"\u2009{unit_text(u)}").rstrip("\u2009"))
    s = expand_macro(s, "numrange", 2,
                     lambda a, b: f"{group_digits(a)}\u2013{group_digits(b)}")
    s = expand_macro(s, "numlist", 1,
                     lambda a: ", ".join(x.strip() for x in a.split(";")))
    s = expand_macro(s, "SI", 2,
                     lambda n, u: f"{group_digits(n)}\u2009{unit_text(u)}")
    s = expand_macro(s, "si", 1, lambda u: unit_text(u))
    s = expand_macro(s, "num", 1, group_digits)
    return s


# ------------------------------------------------------------------------------
def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__.strip())
    src, aux, dst = sys.argv[1:4]

    s = open(src, encoding="utf8").read()
    labels = read_labels(aux)

    # --- title and byline, which live in class macros pandoc cannot evaluate ---
    m = re.search(r"\\title\[[^\]]*\]\{(.*?)\}\s*\n\s*\n", s, re.S)
    if not m:
        m = re.search(r"\\title(?:\[[^\]]*\])?\{((?:[^{}]|\{[^}]*\})*)\}", s, re.S)
    title = " ".join(m.group(1).split()) if m else "Manuscript"
    title = title.replace("~", " ")

    m = re.search(r"\\def\\Authors\{(.*?)\}\s*\n\\def\\Address", s, re.S)
    authors = ""
    if m:
        authors = m.group(1)
        authors = re.sub(r"\\,\$\^\{[^}]*\}\$", "", authors)
        authors = authors.replace("~", " ").replace("\\", "")
        authors = " ".join(authors.split())

    m = re.search(r"\\def\\Address\{(.*?)\}\s*\n", s, re.S)
    address = ""
    if m:
        address = re.sub(r"\$\^\{[^}]*\}\$", "", m.group(1)).strip()

    body_start = s.index("\\begin{document}") + len("\\begin{document}")
    body_end = s.index("\\end{document}")
    body = s[body_start:body_end]
    original_words = prose_words(body)

    # --- abstract becomes an ordinary section -------------------------------
    body = body.replace("\\begin{abstract}\n\\section{}",
                        "\\section*{Abstract}")
    body = body.replace("\\end{abstract}", "")

    # --- drop class-only commands that carry no text ------------------------
    for cmd in (r"\\onecolumn", r"\\linenumbers", r"\\maketitle", r"\\tiny",
                r"\\normalsize", r"\\clearpage", r"\\newpage"):
        body = re.sub(cmd + r"\b", "", body)
    body = re.sub(r"\\firstpage\{[^}]*\}", "", body)
    body = re.sub(r"\\author\[[^\]]*\]\{[^}]*\}", "", body)
    for cmd in ("address", "correspondance", "extraAuth"):
        body = re.sub(r"\\" + cmd + r"\{[^}]*\}", "", body)
    body = re.sub(r"\\title(?:\[[^\]]*\])?\{(?:[^{}]|\{[^}]*\})*\}", "", body, flags=re.S)
    body = expand_macro(body, "keyFont", 1, lambda a: a)
    body = body.replace("\\section{Keywords:}", "\n\n\\textbf{Keywords:} ")

    # --- resolve cross-references from the .aux -----------------------------
    unresolved = set()

    def ref(kind):
        def go(key):
            if key in labels:
                return labels[key]
            unresolved.add(key)
            return "?"
        return go

    body = expand_macro(body, "ref", 1, ref("ref"))
    body = expand_macro(body, "eqref", 1, lambda k: "(%s)" % labels.get(k, "?"))
    body = expand_macro(body, "label", 1, lambda k: "")

    # --- siunitx ------------------------------------------------------------
    body = expand_siunitx(body)

    # --- figures: PDFs do not embed in a .docx, PNGs do ---------------------
    def fig(path):
        png = re.sub(r"\.pdf$", ".png", path.strip())
        cand = os.path.join(os.path.dirname(os.path.abspath(src)), png)
        return png if os.path.isfile(cand) else path

    body = re.sub(r"\\includegraphics(\[[^\]]*\])?\{([^}]*)\}",
                  lambda m: "\\includegraphics%s{%s}" % (m.group(1) or "", fig(m.group(2))),
                  body)

    # --- assemble -----------------------------------------------------------
    head = [
        "\\documentclass{article}",
        "\\usepackage{graphicx}",
        "\\usepackage{booktabs}",
        "\\usepackage{amsmath}",
        "\\begin{document}",
        "",
        "\\section*{%s}" % title,
        "",
    ]
    if authors:
        head += ["%s\\\\" % authors, ""]
    if address:
        head += ["\\emph{%s}" % address, ""]

    out = "\n".join(head) + body + "\n\\end{document}\n"
    open(dst, "w", encoding="utf8").write(out)

    # --- report; nothing here is allowed to fail quietly --------------------
    final_words = prose_words(body)
    kept = 100.0 * final_words / original_words if original_words else 0
    print("      title:   %s..." % title[:58])
    print("      authors: %s" % (authors or "(none found)"))
    print("      refs:    %d resolved from .aux" % len(labels))
    if unresolved:
        print("      WARNING: %d unresolved label(s): %s"
              % (len(unresolved), ", ".join(sorted(unresolved)[:6])))
    print("      prose:   %d -> %d words (%.1f%% kept)"
          % (original_words, final_words, kept))
    if kept < 99.0:
        sys.exit("      ABORT: preprocessing lost prose, not just markup. "
                 "Inspect %s before trusting the .docx." % dst)
    leftovers = re.findall(r"\\(SI|num|SIrange|numlist|si)\b", body)
    if leftovers:
        print("      WARNING: %d unexpanded siunitx macro(s)" % len(leftovers))


if __name__ == "__main__":
    main()
