#!/usr/bin/env python3
from pathlib import Path
import json
import os
import sys

ROOT = Path(__file__).resolve().parents[1]
SPEC = json.loads((ROOT/"config"/"Figure1_layout_spec.json").read_text(encoding="utf-8"))
FIG = Path(os.environ.get("VDR_FIGURE_OUTPUT",str(ROOT/"results"/"figures")))
if '--check' in sys.argv:
    assert isinstance(SPEC.get('nodes'),list) and SPEC['nodes']
    print('Figure1: layout input available')
    sys.exit(0)
FIG.mkdir(parents=True, exist_ok=True)

W=SPEC["canvas"]["width_px"]; H=SPEC["canvas"]["height_px"]; MM=SPEC["target_width_mm"]
ST=SPEC["style"]; TY=SPEC["typography"]; FONT=SPEC["font_family"]
N={n["id"]:n for n in SPEC["nodes"]}

def esc(s):
    return str(s).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")

def tx(x,y,s,size,weight=400,fill=None,anchor="start"):
    return f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" fill="{fill or ST["text"]}">{esc(s)}</text>'

def ml(x,y,lines,size,gap=1.22,weight=400,fill=None,anchor="start"):
    p=[f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" fill="{fill or ST["text"]}">']
    for i,s in enumerate(lines):
        dy=0 if i==0 else size*gap
        p.append(f'<tspan x="{x}" dy="{dy}">{esc(s)}</tspan>')
    p.append("</text>")
    return "".join(p)

def box(n, sw=None):
    dash=' stroke-dasharray="8 7"' if n.get("dashed") else ""
    return f'<rect x="{n["x"]}" y="{n["y"]}" width="{n["w"]}" height="{n["h"]}" rx="{ST["corner_radius_px"]}" ry="{ST["corner_radius_px"]}" fill="{n["fill"]}" stroke="{n["stroke"]}" stroke-width="{sw or ST["node_border_px"]}"{dash}/>'

def route(points, arrow=False, dash=None, stroke=None, sw=None):
    pts=" ".join(f"{x},{y}" for x,y in points)
    ma=' marker-end="url(#arrow)"' if arrow else ""
    da=f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke or ST["connector"]}" stroke-width="{sw or ST["connector_width_px"]}" stroke-linecap="round" stroke-linejoin="round"{ma}{da}/>'

head = (
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{MM}mm" height="{MM*H/W}mm" viewBox="0 0 {W} {H}">'
    f'<defs><marker id="arrow" markerWidth="{ST["arrow_length_px"]}" markerHeight="{ST["arrow_width_px"]}" '
    f'refX="{ST["arrow_length_px"]-0.5}" refY="{ST["arrow_width_px"]/2}" orient="auto" markerUnits="userSpaceOnUse">'
    f'<path d="M0,0 L0,{ST["arrow_width_px"]} L{ST["arrow_length_px"]},{ST["arrow_width_px"]/2} z" fill="{ST["connector"]}"/>'
    '</marker></defs>'
    f'<rect x="0" y="0" width="{W}" height="{H}" fill="{SPEC["canvas"]["background"]}"/>'
)
svg=[head]

for s in SPEC["sections"]:
    svg.append(tx(s["x"],s["y"],s["number"],TY["section_number_px"],700,s["number_color"]))
    svg.append(tx(s["x"]+37,s["y"],s["title"],TY["section_title_px"],700))

a=N["experimental_perturbation"]; b=N["rna_candidates"]; c=N["locked_vdr_programs"]; e=N["published_occupancy"]
svg.append(route([(a["x"]+a["w"],a["y"]+a["h"]/2),(b["x"],b["y"]+b["h"]/2)],True))
svg.append(route([(b["x"]+b["w"],b["y"]+b["h"]/2),(c["x"],c["y"]+c["h"]/2)],True))
svg.append(route([(e["x"]+e["w"]/2,e["y"]+e["h"]),(e["x"]+e["w"]/2,c["y"])],True))

p=N["patient_bulk"]; f=N["fertile_luteal"]; r=N["rif_tests"]
centers=[p["x"]+p["w"]/2,f["x"]+f["w"]/2,r["x"]+r["w"]/2]
dy=SPEC["routing"]["stage2_distributor_y"]
svg.append(route([(centers[2],c["y"]+c["h"]),(centers[2],dy)]))
svg.append(route([(centers[0],dy),(centers[2],dy)]))
for cx,n in zip(centers,[p,f,r]):
    svg.append(route([(cx,dy),(cx,n["y"])],True))

cy=SPEC["routing"]["convergence_y"]
for cx,n in zip(centers,[p,f,r]):
    svg.append(route([(cx,n["y"]+n["h"]),(cx,cy)]))
svg.append(route([(centers[0],cy),(centers[2],cy)]))
integ=N["integrated"]; ix=integ["x"]+integ["w"]/2
svg.append(route([(ix,cy),(ix,integ["y"])],True))
bd=N["boundary"]
svg.append(route([(ix,integ["y"]+integ["h"]),(ix,bd["y"])],False,"7 6",ST["boundary_border"],1.8))

order=["experimental_perturbation","rna_candidates","locked_vdr_programs","published_occupancy",
       "patient_bulk","fertile_luteal","rif_tests","integrated","boundary"]
for nid in order:
    n=N[nid]
    sw=3.0 if nid=="integrated" else (ST["boundary_width_px"] if nid=="boundary" else None)
    svg.append(box(n,sw))
    cls=n["node_class"]
    if cls=="external_capsule":
        svg.append(tx(n["x"]+n["w"]/2,n["y"]+37,n["title"],26,700,n["stroke"],"middle"))
        continue
    x=n["x"]+n.get("padding_left_px",34)
    if nid=="integrated":
        svg.append(tx(x,n["y"]+40,n["title"],TY["integration_title_px"],700,n["stroke"]))
        svg.append(tx(x,n["y"]+76,n["body"][0],TY["integration_body_px"],400))
        svg.append(tx(x,n["y"]+108,n["takeaway"],TY["integration_body_px"],700))
    elif nid=="boundary":
        svg.append(tx(x,n["y"]+36,n["title"],TY["boundary_title_px"],700,n["stroke"]))
        svg.append(ml(x,n["y"]+70,n["body"],TY["boundary_body_px"],1.22))
    else:
        svg.append(tx(x,n["y"]+43,n["title"],TY["node_title_px"],700))
        svg.append(ml(x,n["y"]+83,n["body"],TY["node_body_px"],1.25))
        if n.get("takeaway"):
            svg.append(tx(x,n["y"]+n["h"]-30,n["takeaway"],TY["takeaway_px"],700))

svg.append("</svg>")
s="".join(svg)

svg_p=FIG/"Figure1.svg"
if svg_p.exists():
    raise FileExistsError(svg_p)
svg_p.write_text(s,encoding="utf-8")

print("OK")
