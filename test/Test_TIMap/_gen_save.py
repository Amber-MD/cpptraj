"""Generate Test_TemplateMatch *.map.save files matching Exec_TIMap printf.

Regenerates expected maps for:
  FLE.lib (rA) onto ERN.lib (dA)
  ERN naorder self-map
  phenol.mol2 onto benzene.mol2
  sec.lib (elmnt 34) onto cys.lib (elmnt 16)
  shipped dat/templatematch/DAA.mol2 identity

Cys/Sec are written as Amber OFF because mol2 atom name SE is sulfur in cpptraj.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "util" / "xna_libalign"))
from io_mol2 import read_mol2
from io_off import read_off, write_off
from match import map_units, output_order, canonical_order
from scaffold import detect_scaffold


def digit_width(n: int) -> int:
    w = 1
    x = abs(n)
    while x >= 10:
        x //= 10
        w += 1
    return w


def write_map(path: Path, tgt_name: str, tpl_name: str, tgt, tpl,
              mapping_dict, order, tgt_kind: str, tpl_kind: str) -> None:
    tgt_names = [a.amber_name() for a in tgt.atoms]
    tpl_names = [a.amber_name() for a in tpl.atoms]
    mapping = [mapping_dict.get(i, -1) for i in range(tgt.n_atom())]
    n_mapped = sum(1 for m in mapping if m >= 0)
    n_ins = sum(1 for m in mapping if m < 0)
    used = [False] * tpl.n_atom()
    for m in mapping:
        if m >= 0:
            used[m] = True
    n_unmap = used.count(False)
    a_width = max(6, digit_width(tgt.n_atom()), digit_width(tpl.n_atom()))
    n_width = 6
    for nm in tgt_names + tpl_names:
        n_width = max(n_width, len(nm))
    lines = [
        f"# timap tgt='{tgt_name}' template='{tpl_name}'",
        f"# kind tgt={tgt_kind} template={tpl_kind}",
        f"# mapped= {n_mapped}  insertion= {n_ins}  unmapped_template= {n_unmap}"
        f"  n_tgt= {tgt.n_atom()}  n_tpl= {tpl.n_atom()}",
        "%-*s %*s %-*s %*s %-*s %s" % (
            a_width, "#Out", a_width, "TgtAt", n_width, "TgtName",
            a_width, "TplAt", n_width, "TplName", "Ins",
        ),
    ]
    for k, oldat in enumerate(order):
        tplidx = mapping[oldat]
        if tplidx < 0:
            lines.append("%*i %*i %-*s %*s %-*s %s" % (
                a_width, k + 1, a_width, oldat + 1, n_width, tgt_names[oldat],
                a_width, "---", n_width, "---", "1",
            ))
        else:
            lines.append("%*i %*i %-*s %*i %-*s %s" % (
                a_width, k + 1, a_width, oldat + 1, n_width, tgt_names[oldat],
                a_width, tplidx + 1, n_width, tpl_names[tplidx], "0",
            ))
    if n_unmap:
        lines.append("# unmapped template atoms:")
        for r, u in enumerate(used):
            if not u:
                lines.append("#   %i %s" % (r + 1, tpl_names[r]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote", path.name, "mapped", n_mapped, "ins", n_ins, "unmap", n_unmap)


def align(tgt, tpl, naorder=False):
    tsc = detect_scaffold(tgt)
    psc = detect_scaffold(tpl)
    mapping = map_units(tgt, tpl, tsc, psc)
    parent = canonical_order(tpl, psc) if naorder else list(range(tpl.n_atom()))
    order = output_order(tgt, tpl, mapping, tsc, psc, parent_order=parent)
    return mapping, order, tsc.kind, psc.kind


def main() -> None:
    fle = read_off(HERE / "FLE.lib")[0]
    ern = read_off(HERE / "ERN.lib")[0]
    m, o, tk, pk = align(fle, ern)
    write_map(HERE / "fle_to_ern.map.save", "FLE[FLE]", "ERN[ERN]", fle, ern, m, o, tk, pk)

    m, o, tk, pk = align(ern, ern, naorder=True)
    # identity mapping for self-naorder
    mapping = {i: i for i in range(ern.n_atom())}
    write_map(HERE / "ern_naorder.map.save", "ERN[ERN]", "ERN[ERN]", ern, ern, mapping, o, tk, pk)

    bnz = read_mol2(HERE / "benzene.mol2")[0]
    phe = read_mol2(HERE / "phenol.mol2")[0]
    m, o, tk, pk = align(phe, bnz)
    write_map(HERE / "phenol_to_benzene.map.save", "phenol", "benzene", phe, bnz, m, o, tk, pk)

    cys = read_mol2(HERE / "cys.mol2")[0]
    sec = read_mol2(HERE / "sec.mol2")[0]
    cys.name, cys.resname = "CYS", "CYS"
    sec.name, sec.resname = "SEC", "SEC"
    write_off(HERE / "cys.lib", [cys])
    write_off(HERE / "sec.lib", [sec])
    # Reload so the expected map matches the OFF units cpptraj will read.
    cys = read_off(HERE / "cys.lib")[0]
    sec = read_off(HERE / "sec.lib")[0]
    m, o, tk, pk = align(sec, cys)
    write_map(HERE / "sec_to_cys.map.save", "SEC[SEC]", "CYS[CYS]", sec, cys, m, o, tk, pk)

    daa = read_mol2(ROOT / "dat" / "templatematch" / "DAA.mol2")[0]
    m, o, tk, pk = align(daa, daa)
    write_map(HERE / "daa_identity.map.save", "DAA", "DAA", daa, daa, m, o, tk, pk)


if __name__ == "__main__":
    main()
