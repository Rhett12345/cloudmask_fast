#!/usr/bin/env python3
"""Convert Fortran threshold file to structured YAML."""
import yaml, re, os, sys
from collections import OrderedDict

THR_FILE = os.path.join(os.path.dirname(__file__), "..", "coeff", "fylat_thresholds.mersi.ii3d.v8")
OUT_FILE = os.path.join(os.path.dirname(__file__), "..", "coeff", "thresholds_mersi_ii.yaml")

def parse_thr(filepath):
    thresholds = OrderedDict()
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('!') or line.startswith('rcs_id') or line.startswith('thresholds_file_ver'):
                continue
            if ':' in line:
                name, _, vals = line.partition(':')
                name = name.strip()
                vals = re.sub(r'!.*$', '', vals)
                values = []
                for t in vals.replace(',', ' ').split():
                    try: values.append(float(t))
                    except ValueError: pass
                if values and name:
                    thresholds[name] = values
    return thresholds

# Scene definitions with known threshold keys
scene_defs = OrderedDict([
    ("ocean_day", {
        "desc": "Daytime ocean (non-polar)",
        "keys": ["dobt11","do11_12hi","do11_4lo","doref2","doref3","dovrathi","dovratlo","dotci"],
    }),
    ("ocean_nite", {
        "desc": "Nighttime ocean (non-polar)",
        "keys": ["nobt11","no11_12hi","no11_4lo","no86_73","no_11var"],
    }),
    ("land_day", {
        "desc": "Daytime land (non-polar)",
        "keys": ["dl11_12hi","dl11_4lo","dlref1","dlref3","dlvrat","dltci"],
    }),
    ("land_nite", {
        "desc": "Nighttime land (non-polar)",
        "keys": ["nl4_12hi","nl7_11s","nl11_12hi","nl_11_4l","nl_11_4h","nl_11_4m","bt_diff_bounds"],
    }),
    ("day_snow", {
        "desc": "Daytime snow/ice",
        "keys": ["ds11_12hi","ds11_12adj","ds4_11","ds4_11hel","dsref3","dstci"],
    }),
    ("nite_snow", {
        "desc": "Nighttime snow/ice",
        "keys": ["ns11_12hi","ns11_12adj","ns11_4lo","ns4_12hi"],
    }),
    ("day_desert", {
        "desc": "Daytime desert",
        "keys": ["lds11_12hi","lds11_4hi","lds11_4lo","ldsref2","ldsref3","ldstci"],
    }),
    ("polar_day_land", {
        "desc": "Polar daytime land",
        "keys": ["pdl11_12hi","pdl11_4lo","pdlref1","pdlvrat","pdlref3","pdltci"],
    }),
    ("polar_day_ocean", {
        "desc": "Polar daytime ocean",
        "keys": ["pdobt11","pdo11_12hi","pdo11_4lo","pdoref2","pdoref3","pdovrathi","pdovratlo","pdotci"],
    }),
    ("polar_day_snow", {
        "desc": "Polar daytime snow",
        "keys": ["dps11_12hi","dps11_12adj","dps4_11l","dps4_11h","dps4_11m1","bt_11_bnds3","dpsref3","dpstci"],
    }),
    ("polar_day_desert", {
        "desc": "Polar daytime desert",
        "keys": ["pds11_12hi","pds11_4hi","pds11_4lo","pdsref2","pdsref3","pdstci"],
    }),
    ("polar_nite_land", {
        "desc": "Polar nighttime land",
        "keys": ["pnl11_12hi","pn_11_4l","pn_11_4h","pn_11_4m1","pn_7_11l","pn_7_11h","pn_7_11m1","pn_4_12l","pn_4_12h","pn_4_12m1"],
    }),
    ("polar_nite_ocean", {
        "desc": "Polar nighttime ocean",
        "keys": ["pnobt11","pno11_12hi","pno11_4lo","pno86_73","pno_11var"],
    }),
    ("polar_nite_snow", {
        "desc": "Polar nighttime snow",
        "keys": ["pns11_12hi","pn11_12adj","pn_11_4l","pn_11_4h","pn_11_4m1","pn_7_11l","pn_7_11h","pn_7_11m1","pn_4_12l","pn_4_12h","pn_4_12m1"],
    }),
    ("antarctic_day", {
        "desc": "Antarctic daytime (special)",
        "keys": ["ant4_11l","ant4_11h","ant4_11m1","bt_11_bnds4"],
    }),
    ("shared", {
        "desc": "Shared across multiple scenes (FMFT, sun glint, shadows, land restoral, snow mask)",
        "keys": [
            "pfmft_11maxthre","pfmft_btd_min","pfmft_ocean","pfmft_land","pfmft_cold","pfmft_snow",
            "nfmft_ocean","nfmft_land","nfmft_desert","nfmft_snow","nfmft_maxthre",
            "snglntv","snglntvcl","snglntvch","sg_tbdfl","sg_tbdfh","snglrat",
            "snglnt0","snglnt10","snglnt20","snglnt_bounds",
            "dovar11","novar11",
            "shadnir","shavrat","shad124",
            "nc21","ncrat","ncvrat","ncsig","nc11_12",
            "sm_bt11","sm_ndsi","sm_ref2","sm_ref3","sm85_11","sm37_11","sm37_11hel","sm_mnir",
            "swc_ndvi",
            "ldsbt11","ldsbt11bd","ldsr5_4_thr","ldr5_4_thr","ld20m22","ld22m31","lnbt11",
        ],
    }),
])

def main():
    thr = parse_thr(THR_FILE)

    scenes = OrderedDict()
    for name, defn in scene_defs.items():
        scene_thr = OrderedDict()
        for key in defn["keys"]:
            if key in thr:
                scene_thr[key] = thr[key]
        scenes[name] = OrderedDict({
            "description": defn["desc"],
            "thresholds": scene_thr,
        })

    output = OrderedDict()
    output["metadata"] = OrderedDict([
        ("description", "FY-3D MERSI-II cloud mask thresholds — YAML version"),
        ("source", "coeff/fylat_thresholds.mersi.ii3d.v8"),
        ("version", "1.0"),
        ("note", "Values: [locut, midpt, hicut, power] for conf_test S-curves"),
    ])
    output["scenes"] = OrderedDict()
    for sn, sd in scenes.items():
        if sd["thresholds"]:
            output["scenes"][sn] = sd

    # Convert OrderedDict to plain dict to avoid YAML tag issues
    def _to_plain_dict(obj):
        if isinstance(obj, OrderedDict):
            return {k: _to_plain_dict(v) for k, v in obj.items()}
        elif isinstance(obj, dict):
            return {k: _to_plain_dict(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [_to_plain_dict(v) for v in obj]
        return obj

    plain_output = _to_plain_dict(output)
    with open(OUT_FILE, "w") as f:
        yaml.dump(plain_output, f, default_flow_style=False, allow_unicode=True, sort_keys=False, width=120)

    total = sum(len(s["thresholds"]) for s in scenes.values())
    print(f"Wrote {total} parameters across {len([s for s in scenes.values() if s['thresholds']])} scenes -> {OUT_FILE}")
    for sn, sd in scenes.items():
        n = len(sd["thresholds"])
        if n > 0:
            print(f"  {sn}: {n} params")

if __name__ == "__main__":
    main()
