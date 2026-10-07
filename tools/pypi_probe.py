import json, ssl, urllib.request

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

PKGS = ["pysdot", "sdot", "geogram", "pygeogram", "POT", "geomloss",
        "ott-jax", "power-diagram", "laguerre", "pyvoro", "scipy",
        "gudhi", "pygalmesh", "torch-ot", "cgal"]

for p in PKGS:
    try:
        with urllib.request.urlopen(
                f"https://pypi.org/pypi/{p}/json", timeout=30, context=CTX) as r:
            j = json.load(r)
        info = j["info"]
        v = info["version"]
        files = j["releases"].get(v, [])
        up = files[0]["upload_time"] if files else "?"
        names = [f["filename"] for f in files]
        win = [n for n in names if "win" in n.lower()]
        lic = (info.get("license") or "N/A").replace("\n", " ")[:40]
        print(f"=== {p}  v{v}   license={lic}")
        print(f"    uploaded={up}  requires_python={info.get('requires_python')}")
        print(f"    homepage={info.get('home_page') or info.get('project_url')}")
        print(f"    files({len(names)}): {names[:6]}")
        print(f"    WINDOWS WHEELS: {win[:6] if win else 'NONE'}")
    except Exception as e:
        print(f"=== {p}  --> ERROR {type(e).__name__}: {e}")
    print()
