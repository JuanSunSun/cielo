import json, sys, numpy as np, xarray as xr, s3fs
from datetime import datetime, timedelta, timezone
sys.path.insert(0, ".")
import procesar as P
LAT, LON = -41.319, -72.985
g = P.GOES()
ahora = datetime.now(timezone.utc)
out = {"ahora": ahora.isoformat()}
acm = g.listar("ABI-L2-ACMF", ahora - timedelta(hours=1), ahora)
ts, f = acm[-1]
out["acm_file"] = f
def win(ds, var, half=3):
    pj, H, lon0 = P.proyeccion(ds)
    xm, ym = pj(LON, LAT)
    ix = int(np.argmin(np.abs(ds.x.values - xm / H))); iy = int(np.argmin(np.abs(ds.y.values - ym / H)))
    sub = ds[var].isel(x=slice(ix-half, ix+half+1), y=slice(iy-half, iy+half+1)).load()
    X, Y = np.meshgrid(ds.x.values[ix-half:ix+half+1]*H, ds.y.values[iy-half:iy+half+1]*H)
    lo, la = pj(X, Y, inverse=True)
    return {"ix": ix, "iy": iy, "vals": np.round(sub.values.astype(float), 2).tolist(),
            "lat_c": float(la[half, half]), "lon_c": float(lo[half, half]),
            "attrs": {k: str(v) for k, v in ds[var].attrs.items()},
            "enc": {k: str(v) for k, v in ds[var].encoding.items() if k in ("dtype","_FillValue","scale_factor","add_offset","_Unsigned")}}
with g.abrir(f) as ds:
    out["acm_vars"] = list(ds.data_vars)
    for v in ("BCM", "ACM", "Cloud_Probabilities", "DQF"):
        if v in ds: out[v] = win(ds, v)
    # region stats
    sy, sx = P.caja_indices(ds, P.REGION, 0)
    b = ds["BCM"].isel(y=sy, x=sx).load().values
    out["bcm_region_unique"] = {str(k): int(n) for k, n in zip(*np.unique(np.nan_to_num(b, nan=-9), return_counts=True))}
    raw = xr.open_dataset(g.fs.open(f, "rb"), engine="h5netcdf", mask_and_scale=False)
    rb = raw["BCM"].isel(y=sy, x=sx).load().values
    out["bcm_raw_unique"] = {str(k): int(n) for k, n in zip(*np.unique(rb, return_counts=True))}
    out["bcm_raw_attrs"] = {k: str(v) for k, v in raw["BCM"].attrs.items()}
for canal in ("C14", "C07"):
    l = g.listar("ABI-L2-CMIPF", ts - timedelta(minutes=12), ts + timedelta(minutes=12), canal + "_")
    ff = P.mas_cercano(l, ts, 12)
    if ff:
        with g.abrir(ff) as ds: out[canal] = win(ds, "CMI")
lh = g.listar("ABI-L2-ACHAF", ts - timedelta(minutes=70), ts + timedelta(minutes=15))
fh = P.mas_cercano(lh, ts, 70)
out["acha_file"] = fh
if fh:
    with g.abrir(fh) as ds: out["HT"] = win(ds, "HT", 2)
json.dump(out, open(sys.argv[1], "w"), indent=1, default=str)
print(json.dumps({k: out[k] for k in out if k in ("acm_file","bcm_region_unique","bcm_raw_unique")}, default=str))
