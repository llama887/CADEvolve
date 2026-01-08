# Export_Test.py
import os
import json
from onshape_client.client import Client
from pathlib import Path
import time
import pprint
from dotenv import load_dotenv

ROOT = Path(__file__).resolve().parents[1]  # repo root
load_dotenv(ROOT / ".env")

BASE_URL = "https://cad.onshape.com"  # change if your Onshape domain is different

def pick(*vals):
    for v in vals:
        if v is not None:
            return v
    return None

def pget(pm: dict, *keys):
    """Return the first parameter-message dict found for any key."""
    for k in keys:
        if k in pm:
            return pm[k]
    return None

def pquantity(pm: dict, *keys):
    """Return (expr, mm) for a quantity parameter from a list of possible ids."""
    q = pget(pm, *keys)
    expr = q.get("expression") if q else None
    return expr, parse_length_mm(expr)

def pbool(pm: dict, *keys):
    b = pget(pm, *keys)
    return b.get("value") if b else None

def penum(pm: dict, *keys):
    e = pget(pm, *keys)
    return e.get("value") if e else None

def pstr(pm: dict, *keys):
    s = pget(pm, *keys)
    return s.get("value") if s else None
import re

def parse_length_mm(expr: str | None) -> float | None:
    """
    Parse Onshape expressions like '25 mm', '0.012 m', '1 in' into mm.
    Returns None if can't parse.
    """
    if not expr:
        return None
    s = expr.strip()

    m = re.match(r"^\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*([a-zA-Z]+)\s*$", s)
    if not m:
        return None

    val = float(m.group(1))
    unit = m.group(2).lower()

    # common length units
    if unit == "mm":
        return val
    if unit == "cm":
        return val * 10.0
    if unit == "m":
        return val * 1000.0
    if unit in ("in", "inch", "inches"):
        return val * 25.4
    if unit in ("ft", "feet"):
        return val * 304.8

    return None


def feature_type(f: dict) -> str:
    return pick(f.get("typeName"), f.get("type_name"), f.get("featureType"), f.get("btType"), f.get("type")) or ""

def feature_id(f: dict) -> str:
    return pick(f.get("featureId"), f.get("feature_id"), f.get("id")) or ""

def feature_name(f: dict) -> str:
    return pick(f.get("name"), f.get("featureName")) or ""

def iter_params(f: dict):
    """Yield normalized parameter dicts from a feature."""
    for p in f.get("parameters", []) or []:
        if isinstance(p, dict):
            yield p

def param_id(p: dict) -> str:
    return pick(p.get("parameterId"), p.get("parameter_id"), p.get("id")) or ""

def param_value(p: dict):
    # some params store numeric in "value", others in nested dicts
    return pick(p.get("value"), p.get("booleanValue"), p.get("stringValue"), p.get("integerValue"))

def param_expr(p: dict):
    return pick(p.get("expression"), p.get("expr"), p.get("valueExpression"))

def param_type(p: dict) -> str:
    return pick(p.get("typeName"), p.get("type_name"), p.get("type")) or ""

def build_param_map(f: dict) -> dict:
    """Map parameterId -> (param dict) for easy lookup."""
    m = {}
    for p in iter_params(f):
        pid = param_id(p)
        if pid:
            m[pid] = p
    return m

def create_step_translation(api_client, did, wvm, wvmid, eid):
    payload = {
        "formatName": "STEP",
        "storeInDocument": False,
        "stepVersion": "AP214",
        "includeExportIds": False,
    }

    resp, status, headers = call_api_raw(
        api_client,
        method="POST",
        resource_path="/api/partstudios/d/{did}/{wvm}/{wvmid}/e/{eid}/translations",
        path_params={"did": did, "wvm": wvm, "wvmid": wvmid, "eid": eid},
        header_params={"Accept": "application/json", "Content-Type": "application/json"},
        body=payload,
    )

    raw = resp.data.decode("utf-8", errors="replace")
    if status >= 400:
        raise RuntimeError(f"Create translation failed HTTP {status}\n{raw}")

    return json.loads(raw)

def wait_for_translation_done(api_client, translation_id, timeout_s=180):
    t0 = time.time()
    while True:
        resp, status, headers = call_api_raw(
            api_client,
            method="GET",
            resource_path="/api/translations/{tid}",
            path_params={"tid": translation_id},
            header_params={"Accept": "application/json"},
            body=None,
        )

        raw = resp.data.decode("utf-8", errors="replace")
        if status >= 400:
            raise RuntimeError(f"Poll translation failed HTTP {status}\n{raw}")

        info = json.loads(raw)
        state = info.get("requestState")
        print(f"[poll] {translation_id} state={state}", flush=True)

        if state == "DONE":
            return info
        if state == "FAILED":
            raise RuntimeError("Translation FAILED:\n" + json.dumps(info, indent=2))

        if time.time() - t0 > timeout_s:
            raise TimeoutError("Timed out.\n" + json.dumps(info, indent=2))

        time.sleep(2)

def download_external_data(api_client, did, external_data_id, out_path):
    resp, status, headers = call_api_raw(
        api_client,
        method="GET",
        resource_path="/api/documents/d/{did}/externaldata/{xid}",
        path_params={"did": did, "xid": external_data_id},
        header_params={"Accept": "application/octet-stream"},
        body=None,
    )

    if status >= 400:
        raw = resp.data.decode("utf-8", errors="replace")
        raise RuntimeError(f"Download failed HTTP {status}\n{raw}")

    with open(out_path, "wb") as f:
        f.write(resp.data)

    print(f"✅ Wrote STEP: {out_path}", flush=True)


def export_step(api_client, did, wvm, wvmid, eid, out_path):
    payload = {
        "formatName": "STEP",
        "storeInDocument": False,
        "stepVersion": "AP214",
        "includeExportIds": False,
    }

    resp, status, headers = api_client.call_api(
        resource_path="/api/partstudios/d/{did}/{wvm}/{wvmid}/e/{eid}/translations",
        method="POST",
        path_params={"did": did, "wvm": wvm, "wvmid": wvmid, "eid": eid},
        query_params={},
        header_params={
            "Accept": "application/octet-stream",
            "Content-Type": "application/json",
        },
        body=payload,
        post_params=[],
        files={},
        response_type=None,
        auth_settings=None,
        _return_http_data_only=False,
        _preload_content=False,
        _check_type=False,
    )

    if status >= 400:
        err = resp.data.decode("utf-8", errors="replace")
        raise RuntimeError(f"STEP export failed HTTP {status}\n{err}")

    with open(out_path, "wb") as f:
        f.write(resp.data)

    print(f"Wrote STEP: {out_path}")


EXPORT_PARAMS = {
    "formatName": "STEP",
    "storeInDocument": False,   # download to client instead
    "stepVersion": "AP214",     # AP203 or AP214
    "includeExportIds": False
}
def extract_extrudes(features: list[dict]) -> list[dict]:
    out = []
    for f in features:
        fm = msg(f)
        meta = get_feature_meta(f)
        if (meta.get("featureType") or "").lower() != "extrude":
            continue

        pm = param_map(fm)

        depth_expr, depth_mm = pquantity(pm, "depth", "endDepth", "distance", "extrudeDepth")
        end_condition = penum(pm, "endBound", "endType", "endCondition")
        opposite = pbool(pm, "oppositeDirection")
        symmetric = pbool(pm, "symmetric") or pbool(pm, "midplane")

        has_second = pbool(pm, "hasSecondDirection")
        second_depth_expr, second_depth_mm = pquantity(pm, "secondDirectionDepth")
        second_end_condition = penum(pm, "secondDirectionBound")

        out.append({
            **meta,
            "depth_expression": depth_expr,
            "depth_mm": depth_mm,
            "end_condition": end_condition,
            "direction": "opposite" if opposite else "normal",
            "symmetric": bool(symmetric),
            "has_second_direction": bool(has_second),
            "second_depth_expression": second_depth_expr,
            "second_depth_mm": second_depth_mm,
            "second_end_condition": second_end_condition,
        })
    return out


def extract_fillets(features: list[dict]) -> list[dict]:
    out = []
    for f in features:
        fm = msg(f)
        meta = get_feature_meta(f)
        if (meta.get("featureType") or "").lower() != "fillet":
            continue

        pm = param_map(fm)

        rad_expr, rad_mm = pquantity(pm, "radius", "filletRadius", "nonCircularRadius")
        fillet_type = penum(pm, "filletType")
        is_variable = pbool(pm, "isVariable")

        out.append({
            **meta,
            "radius_expression": rad_expr,
            "radius_mm": rad_mm,
            "fillet_type": fillet_type,
            "is_variable": bool(is_variable),
        })
    return out


def require_env(name: str) -> str:
    v = os.getenv(name, "").strip()
    if not v:
        raise RuntimeError(
            f"Missing environment variable: {name}\n\n"
            "In Command Prompt (cmd), run:\n"
            "  set ONSHAPE_ACCESS_KEY=YOUR_ACCESS_KEY\n"
            "  set ONSHAPE_SECRET_KEY=YOUR_SECRET_KEY\n"
            "Then run:\n"
            "  python Export_Test.py\n"
        )
    return v


def as_dict(x):
    """Convert OpenAPI model objects (e.g., BTPartMetadataInfo) to plain dict safely."""
    if hasattr(x, "to_dict"):
        try:
            return x.to_dict()
        except Exception:
            pass
    return getattr(x, "_data_store", {}) if hasattr(x, "_data_store") else x


def get_raw_json(api_client, path: str, path_params: dict):
    """
    Fetch raw JSON from Onshape without OpenAPI model deserialization.

    This avoids schema drift issues like:
    - feature_states returned as list instead of dict
    - KeyError: 's' when response_type='str'
    """
    resp, status_code, headers = api_client.call_api(
        resource_path=path,
        method="GET",
        path_params=path_params,
        query_params={},
        header_params={"Accept": "application/json"},
        body=None,
        post_params=[],
        files={},
        response_type=None,
        auth_settings=None,
        _return_http_data_only=False,
        _preload_content=False,
        _check_type=False,
    )

    raw = resp.data.decode("utf-8")
    if status_code >= 400:
        raise RuntimeError(f"HTTP {status_code} when calling {path}: {raw[:500]}")
    return json.loads(raw)


def try_mass_properties(client: Client, DID: str, WID: str, EID: str):
    """
    Different onshape-client builds expose different method names.
    We'll try a few common ones; if none exist, print available methods containing 'mass'.
    """
    api = client.parts_api
    candidates = [
        "get_mass_properties_wmv",
        "get_mass_properties",
        "get_part_mass_properties_wmv",
        "get_part_mass_properties",
        "get_massproperties_wmv",
        "get_massproperties",
    ]

    for name in candidates:
        if hasattr(api, name):
            fn = getattr(api, name)
            try:
                # Most WMV-style endpoints take (did, wvm, wvmid, ...) and sometimes eid
                try:
                    mp = fn(did=DID, wvm="w", wvmid=WID, eid=EID)
                except TypeError:
                    try:
                        mp = fn(did=DID, wvm="w", wvmid=WID)
                    except TypeError:
                        mp = fn(did=DID, wid=WID, eid=EID)
                mpd = as_dict(mp)
                print("\nMass properties (from method:", name + "):")
                if isinstance(mpd, dict):
                    print("Keys:", list(mpd.keys()))
                    print("Mass:", mpd.get("mass"))
                    print("Center of mass:", mpd.get("centerOfMass"))
                else:
                    print(mpd)
                return
            except Exception as e:
                print(f"\nTried {name} but it failed: {type(e).__name__}: {e}")

    # If none worked, show what's available
    mass_methods = [m for m in dir(api) if "mass" in m.lower()]
    print("\nNo mass-properties method succeeded. Methods containing 'mass' on your parts_api:")
    print("\n".join(sorted(mass_methods)))

def msg(x):
    """Return the nested 'message' dict if present, else the dict itself."""
    if isinstance(x, dict) and "message" in x and isinstance(x["message"], dict):
        return x["message"]
    return x if isinstance(x, dict) else {}

def get_feature_meta(f):
    fm = msg(f)
    return {
        "featureId": fm.get("featureId"),
        "name": fm.get("name"),
        "featureType": fm.get("featureType"),   # e.g. 'extrude', 'fillet', 'hole', 'newSketch'
        "typeName": f.get("typeName"),          # e.g. 'BTMFeature', 'BTMSketch'
        "suppressed": fm.get("suppressed"),
    }

def param_map(fm: dict) -> dict:
    """
    Build map: parameterId -> parameter_message_dict
    fm is the feature 'message' dict.
    """
    out = {}
    for p in fm.get("parameters", []) or []:
        pm = msg(p)
        pid = pm.get("parameterId")
        if pid:
            out[pid] = pm
    return out

def extract_quantity(pm: dict):
    """Return (expression, value, units) if available for quantity params."""
    if not pm:
        return (None, None, None)
    return (pm.get("expression"), pm.get("value"), pm.get("units"))

def extract_string(pm: dict):
    return pm.get("value") if pm else None

def extract_bool(pm: dict):
    return pm.get("value") if pm else None

def extract_enum(pm: dict):
    return pm.get("value") if pm else None

def extract_holes(features: list[dict]) -> list[dict]:
    out = []
    for f in features:
        fm = msg(f)
        meta = get_feature_meta(f)
        if (meta.get("featureType") or "").lower() != "hole":
            continue

        pm = param_map(fm)

        # Prefer V2 fields if present, otherwise fallback
        dia_expr, dia_mm = pquantity(pm, "holeDiameterV2", "holeDiameter", "diameter", "holeDiameterV2FitToleranceTable")
        depth_expr, depth_mm = pquantity(pm, "holeDepth", "holeDepthComputed", "holeDepthMultiple", "depth")

        end_style = penum(pm, "endStyleV2", "endStyle")  # often BLIND / THROUGH / UP_TO
        opposite = pbool(pm, "oppositeDirection")

        # Optional: counterbore / countersink
        cbore_d_expr, cbore_d_mm = pquantity(pm, "cBoreDiameter")
        cbore_depth_expr, cbore_depth_mm = pquantity(pm, "cBoreDepth")
        csink_d_expr, csink_d_mm = pquantity(pm, "cSinkDiameter")
        csink_angle = pget(pm, "cSinkAngle").get("expression") if pget(pm, "cSinkAngle") else None

        # Threading (if present)
        thread_standard = pstr(pm, "threadStandard")
        is_tapped_through = pbool(pm, "isTappedThrough")
        tapped_depth_expr, tapped_depth_mm = pquantity(pm, "tappedDepth")

        # Name: your current hole name is '#featureName' because Onshape stores actual in parameter 'featureName'
        actual_name = meta.get("name")
        if actual_name == "#featureName":
            actual_name = pstr(pm, "featureName") or actual_name

        out.append({
            **meta,
            "name": actual_name,
            "diameter_expression": dia_expr,
            "diameter_mm": dia_mm,
            "depth_expression": depth_expr,
            "depth_mm": depth_mm,
            "end_style": end_style,
            "direction": "opposite" if opposite else "normal",
            "counterbore_diameter_mm": cbore_d_mm,
            "counterbore_depth_mm": cbore_depth_mm,
            "countersink_diameter_mm": csink_d_mm,
            "countersink_angle_expression": csink_angle,
            "thread_standard": thread_standard,
            "is_tapped_through": is_tapped_through,
            "tapped_depth_mm": tapped_depth_mm,
        })
    return out

def call_api_raw(api_client, *, method, resource_path, path_params, header_params=None, body=None):
    if header_params is None:
        header_params = {}

    resp, status, headers = api_client.call_api(
        resource_path=resource_path,
        method=method,
        path_params=path_params,
        query_params={},          # <-- CRITICAL: never None
        header_params=header_params,
        body=body,
        post_params=[],           # <-- keep defined
        files={},                 # <-- keep defined
        response_type=None,
        auth_settings=None,
        _return_http_data_only=False,
        _preload_content=False,
        _check_type=False,
    )
    return resp, status, headers

def main():
    client = Client(
        configuration={
            "base_url": BASE_URL,
            "access_key": require_env("ONSHAPE_ACCESS_KEY"),
            "secret_key": require_env("ONSHAPE_SECRET_KEY"),
        }
    )

    # ---- Paste your IDs from the Onshape URL (documents/<did>/w/<wid>/e/<eid>) ----
    DID = "25d67ff80d39c8eea0cc313b"
    WID = "8cf89dd6da3e30940cf16183"
    EID = "98cf0ddaf5244f20468bc887"

    # 1) Feature tree (raw JSON, stable)
    features_resp = get_raw_json(
        client.api_client,
        "/api/partstudios/d/{did}/w/{wid}/e/{eid}/features",
        {"did": DID, "wid": WID, "eid": EID},
    )
    features = features_resp.get("features", [])
    import json
    bodydetails = get_raw_json(
    client.api_client,
    "/api/partstudios/d/{did}/w/{wid}/e/{eid}/bodydetails",
    {"did": DID, "wid": WID, "eid": EID},
    )
    out = {
        "metadata": {
            "feature_count": len(features),
            "note": "All feature JSON dumped from Onshape Part Studio"
        },
        "features": []
    }
    for i, f in enumerate(features):
        fm = msg(f)
        out["features"].append({
            "index": i,
            "featureType": fm.get("featureType"),
            "name": fm.get("name"),
            "suppressed": fm.get("suppressed"),
            "raw": fm
        })

    with open("onshape_features.json", "w", encoding="utf-8") as fp:
        json.dump(out, fp, indent=2)

    print(f"\nWrote {len(features)} features to onshape_features.json")

    with open("onshape_bodydetails.json", "w", encoding="utf-8") as f:
        json.dump(bodydetails, f, indent=2)

    print("Wrote bodydetails to onshape_bodydetails.json")
    extrudes = extract_extrudes(features)
    fillets = extract_fillets(features)
    holes = extract_holes(features)

    print("\n=== EXTRUDES ===")
    for e in extrudes:
        print(e)

    print("\n=== FILLETS ===")
    for fl in fillets:
        print(fl)

    print("\n=== HOLES ===")
    for h in holes:
        print(h)

    print("\n--- DEBUG: features_resp keys ---")
    if isinstance(features_resp, dict):
        print(list(features_resp.keys()))

    print("\n--- DEBUG: first feature raw ---")
    if features:
        pprint.pprint(features[0])
    else:
        print("No features returned")
    print(f"\nFeature count: {len(features)}")
    for f in features[:12]:
        print(
            f"- {f.get('featureId')} | "
            f"type={f.get('typeName')} | "
            f"name={f.get('name')} | "
            f"suppressed={f.get('suppressed')}"
        )

    # 2) Parts list (WMV style in your client)
    if not hasattr(client.parts_api, "get_parts_wmv"):
        raise RuntimeError(
            "Your installed onshape-client doesn't have parts_api.get_parts_wmv.\n"
            "Run: pip show onshape-client\n"
        )

    parts_resp = client.parts_api.get_parts_wmv(did=DID, wvm="w", wvmid=WID)

    # Normalize response to list of parts
    if isinstance(parts_resp, dict):
        parts = parts_resp.get("parts", [])
    else:
        parts = parts_resp

    parts_dicts = [as_dict(p) for p in parts]
    print(f"\nTotal parts returned: {len(parts_dicts)}")

    if parts_dicts:
        # Print one sample's keys to understand structure
        sample = parts_dicts[0]
        if isinstance(sample, dict):
            print("Sample part keys:", list(sample.keys()))
        else:
            print("Sample part (non-dict):", sample)

    # Try to filter parts that belong to this element (field name varies)
    parts_for_eid = []
    for field in ["elementId", "elementID", "eid", "elementIdString"]:
        matches = [p for p in parts_dicts if isinstance(p, dict) and p.get(field) == EID]
        if matches:
            parts_for_eid = matches
            print(f"\nMatched {len(parts_for_eid)} parts using field '{field}'")
            break

    # If filtering didn't work, just print all partIds/names (still useful)
    to_print = parts_for_eid if parts_for_eid else parts_dicts
    # Correct field names are snake_case in your client
    parts_for_eid = [
        p for p in parts_dicts
        if isinstance(p, dict) and p.get("element_id") == EID
    ]      

    to_print = parts_for_eid if parts_for_eid else parts_dicts

    print("\nParts:")
    for p in to_print[:50]:
        print(
            f"- part_id={p.get('part_id')} | "
            f"name={p.get('name')} | "
            f"element_id={p.get('element_id')}"
        )


    # 3) Mass properties (best-effort across client versions)
    if parts_dicts and isinstance(parts_dicts[0], dict):
        partid = parts_dicts[0].get("part_id")

    if partid:
        massprops = client.parts_api.get_mass_properties(
            did=DID,
            wvm="w",
            wvmid=WID,
            eid=EID,
            partid=partid
        )

        mp = as_dict(massprops)
        print("\nMass properties:")
        if isinstance(mp, dict):
            print("Keys:", list(mp.keys()))
            print("Mass:", mp.get("mass"))
            print("Center of mass:", mp.get("center_of_mass"))
        else:
            print(mp)
    else:
        print("\nNo part_id found; cannot query mass properties.")

    extrudes = extract_extrudes(features)
    fillets = extract_fillets(features)
    holes = extract_holes(features)

    print("\n=== EXTRUDES (clean) ===")
    for e in extrudes:
        print(
            e.get("name"),
            "| depth_mm=", e.get("depth_mm"),
            "| end=", e.get("end_condition"),
            "| dir=", e.get("direction"),
            "| symmetric=", e.get("symmetric")
        )

    print("\n=== FILLETS (clean) ===")
    for fl in fillets:
        print(
            fl.get("name"),
            "| radius_mm=", fl.get("radius_mm"),
            "| variable=", fl.get("is_variable"),
            "| fillet_type=", fl.get("fillet_type")
        )

    print("\n=== HOLES (clean) ===")
    for h in holes:
        print(
            h.get("name"),
            "| dia_mm=", h.get("diameter_mm"),
            "| depth_mm=", h.get("depth_mm"),
            "| end=", h.get("end_style"),
            "| cbore_d_mm=", h.get("counterbore_diameter_mm"),
            "| cbore_depth_mm=", h.get("counterbore_depth_mm")
        )
    out_dir = Path(r"C:\Users\daegy\OneDrive\Desktop\Onshape Test")
    WVM = "w"     # or "v" or "m" (match your Onshape URL)
    WVMID = WID   # if WVM="w", WVMID is the workspace id
    # You already have these:
    # DID, WVM ("w"/"v"/"m"), WVMID, EID
    import json
    job = create_step_translation(client.api_client, DID, WVM, WVMID, EID)
    print("Create translation job:\n", json.dumps(job, indent=2), flush=True)

    tid = job.get("id")
    if not tid:
        raise RuntimeError("No translation id returned:\n" + json.dumps(job, indent=2))

    final_info = wait_for_translation_done(client.api_client, tid)

    ext_ids = final_info.get("resultExternalDataIds") or []
    print("resultExternalDataIds:", ext_ids, flush=True)
    if not ext_ids:
        raise RuntimeError("No resultExternalDataIds:\n" + json.dumps(final_info, indent=2))

    download_external_data(client.api_client, DID, ext_ids[0], out_dir / "export.step")

    print("✅ STEP saved:", out_dir / "export.step")

if __name__ == "__main__":
    main()
