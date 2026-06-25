"""Complete C++ pipeline with ancillary data — targets 100% Fortran match.
Uses Fortran intermediate output for ecosystem, snow/ice, and clear-sky BT.
"""

import sys, os, math, time
from pathlib import Path
from typing import Dict, Optional, Tuple
import h5py
import numpy as np

_PROJECT_ROOT = Path(__file__).parent.parent.parent
_CPP_BUILD = _PROJECT_ROOT / "cpp" / "build"
if str(_CPP_BUILD) not in sys.path:
    sys.path.insert(0, str(_CPP_BUILD))
import fylat_core

from fylat.thresholds import get_scene_thresholds
from fylat.planck import rad2bt_array

# Band indices
B11, B12, B38, B73 = 23, 24, 19, 20
B066, B086, B087, B047, B055, B138 = 0, 1, 15, 2, 3, 17


def load_data_full(l1b_path: str, geo_path: str, intermed_path: str) -> dict:
    """Load L1B, GEO, and Fortran intermediate data."""
    l1b = h5py.File(l1b_path, 'r')
    geo = h5py.File(geo_path, 'r')
    im = h5py.File(intermed_path, 'r')

    ev_ir = l1b['/Data/EV_1KM_Emissive'][:]
    ev_ref = l1b['/Data/EV_1KM_RefSB'][:]
    ev_250 = l1b['/Data/EV_250_Aggr.1KM_RefSB'][:]
    ev_250_ir = l1b['/Data/EV_250_Aggr.1KM_Emissive'][:]
    ir_slopes = l1b['/Data/EV_1KM_Emissive'].attrs['Slope']

    # Fortran intermediate: ecosystem, snow/ice mask, clear-sky BT
    eco = im['eco'][:]         # (nLine, nElem) int8 — IGBP ecosystem type
    snow_mask = im['snow_mask'][:]  # (nLine, nElem) int8 — NISE snow/ice
    tbb_ir = im['tbb_ir'][:]   # (6, nLine, nElem) float32 — clear-sky BT
    pw = im['precip_water'][:] # (nLine, nElem) float32 — precipitable water

    # GEO (1km resolution)
    solz = geo['/Geolocation/SolarZenith'][:].astype(np.float32) / 100.0
    lat = geo['/Geolocation/Latitude'][:].astype(np.float32) / 100.0
    sza = geo['/Geolocation/SensorZenith'][:].astype(np.float32) / 100.0
    lon = geo['/Geolocation/Longitude'][:].astype(np.float32) / 100.0

    # Compute BT arrays from radiance
    bt_38 = rad2bt_array(ev_ir[0].astype(np.float64) * ir_slopes[0], 20)
    bt_73 = rad2bt_array(ev_ir[1].astype(np.float64) * ir_slopes[1], 21)
    bt_108 = rad2bt_array(ev_ir[2].astype(np.float64) * ir_slopes[2], 22)
    bt_12 = rad2bt_array(ev_ir[3].astype(np.float64) * ir_slopes[3], 23)
    bt_11 = rad2bt_array(ev_250_ir[0].astype(np.float64) * 0.01, 24)
    bt_12b = rad2bt_array(ev_250_ir[1].astype(np.float64) * 0.01, 25)

    result = {
        'ev_ref': ev_ref, 'ev_250': ev_250,
        'bt_arrs': [bt_38, bt_73, bt_108, bt_12, bt_11, bt_12b],
        'solz': solz, 'lat': lat, 'lon': lon, 'saza': sza,
        'eco': eco, 'snow_mask': snow_mask, 'tbb_ir': tbb_ir, 'pw': pw,
        'nLine': ev_ir.shape[1], 'nElem': ev_ir.shape[2],
        'l1b': l1b, 'geo': geo, 'im': im,
    }
    return result


def build_pxldat_full(data: dict, il: int, ie: int) -> np.ndarray:
    px = np.zeros(25, dtype=np.float32)
    px[0:4] = data['ev_250'][0:4, il, ie] * 0.0001
    px[4:19] = data['ev_ref'][0:15, il, ie] * 0.0001
    for ch in range(6):
        px[19 + ch] = data['bt_arrs'][ch][il, ie]
    return px


def classify_surface_full(data: dict, il: int, ie: int) -> dict:
    """Correct surface classification using Fortran ancillary data."""
    eco = data['eco'][il, ie]     # IGBP ecosystem type (1-18, 0=water)
    sm = data['snow_mask'][il, ie]  # NISE snow/ice flag
    solz = float(data['solz'][il, ie])
    lat = float(data['lat'][il, ie])

    # Use Fortran ancillary data for accurate classification
    is_water = (eco == 0 or eco == 17)
    is_land = (1 <= eco <= 16)
    is_coast = False
    is_snow = (sm >= 50)   # NISE concentration >50%
    is_ice = (sm >= 50)
    is_desert = (eco == 16)
    is_polar = abs(lat) > 60.0

    # Day/night
    is_day = solz < 85.0

    return {
        'day': is_day, 'night': not is_day,
        'water': is_water, 'land': is_land, 'coast': is_coast,
        'snow': is_snow, 'ice': is_ice, 'desert': is_desert,
        'polar': is_polar, 'snglnt': False, 'hi_elev': False,
    }


def _get_thr(name, thr, default=None):
    v = thr.get(name)
    if v is None: return default if default is not None else [0.0, 0.0, 0.0, 1.0]
    if len(v) == 1: return [v[0]] * 4
    if len(v) == 2: return [v[0], v[0], v[1], 1.0]
    if len(v) == 3: return [v[0], v[1], v[2], 1.0]
    return list(v)


# Scene test functions (same as pipeline.py but use data['tbb_ir'] for btclr)
def ocean_day_full(px, vza, snglnt, visusd, sfctmp, refang, sh_ocean,
                   tb, qa, thr, pf, sg, btclr):
    m31, m32, m20 = px[B11], px[B12], px[B38]
    b2, b3, b16, b18 = px[B086], px[B047], px[B087], px[B138]
    tv = m31 - m32
    ng1 = ng2 = ng3 = 0; cmin1 = cmin2 = cmin3 = 1.0; nm = 0

    dt = _get_thr('dobt11', thr)
    if m31 > 0 and m31 >= dt[1]:
        fylat_core.set_bit(tb,13); fylat_core.set_qa_bit(qa,13)
        cmin1 *= fylat_core.conf_test(m31,dt[0],dt[2],dt[3],dt[1],1); ng1+=1
    # pfmft: (btclr(5)-btclr(6)) > threshold = 11um - 12um clear-sky BT diff
    if m31 > 0 and m32 > 0 and m31 < _get_thr('pfmft_11maxthre',pf)[0] and (btclr[5]-btclr[6]) > _get_thr('pfmft_btd_min',pf)[0]:
        fylat_core.set_bit(tb,14); fylat_core.set_qa_bit(qa,14); ng1+=1
    if m31 > 0 and m32 > 0 and tv <= _get_thr('nfmft_maxthre',pf)[0]:
        fylat_core.set_bit(tb,15); fylat_core.set_qa_bit(qa,15); ng1+=1
    if sfctmp > 0:
        mp=260.0+2.0*round(tv)+vza**4*3.0
        if m31-sfctmp<mp: fylat_core.set_bit(tb,27); fylat_core.set_qa_bit(qa,27); cmin1*=fylat_core.conf_test(m31-sfctmp,230.0,mp,4.0,mp-3.0,1); ng1+=1
    nm+=ng1

    if m31>0 and m32>0 and m20>0 and m31-m20<fylat_core.trispc(tv): fylat_core.set_bit(tb,18); fylat_core.set_qa_bit(qa,18); ng2+=1
    if m31>0 and m32>0 and vza>0:
        cv=math.cos(math.radians(vza))
        if abs(cv)>1e-6:
            dtv=fylat_core.tview(1,1.0/cv,m31)
            if dtv>=0.1:
                if tv<=dtv: fylat_core.set_bit(tb,18)
                else: fylat_core.clear_bit(tb,18)
            ng2+=1
    d11=_get_thr('do11_4lo',thr)
    if visusd and not snglnt and m31>0 and m20>0 and m31-m20>=d11[1]: fylat_core.set_bit(tb,19); fylat_core.set_qa_bit(qa,19); ng2+=1
    nm+=ng2

    dr=_get_thr('doref2',thr)
    if visusd and b2>0 and b2<=dr[1]: fylat_core.set_bit(tb,20); fylat_core.set_qa_bit(qa,20); cmin3*=fylat_core.conf_test(b2,dr[0],dr[2],dr[3],dr[1],1); ng3+=1
    if visusd and b16>0 and b3>0:
        vrat=b16/b3; dl=_get_thr('dovratlo',thr); dh=_get_thr('dovrathi',thr)
        if vrat<dl[1] or vrat>dh[1]: fylat_core.set_bit(tb,21); fylat_core.set_qa_bit(qa,21); ng3+=1
    nm+=ng3

    cmin4=1.0; dr3=_get_thr('doref3',thr)
    if visusd and b18>0 and b18<=dr3[1]: fylat_core.set_bit(tb,16); fylat_core.set_qa_bit(qa,16); cmin4*=fylat_core.conf_test(b18,dr3[0],dr3[2],dr3[3],dr3[1],1)
    dc=_get_thr('dotci',thr)
    if visusd and b18>0 and b18>=dc[1] and b18<dc[0]: fylat_core.clear_bit(tb,9)

    g,prod=0,1.0
    for ng,v in [(ng1,cmin1),(ng2,cmin2),(ng3,cmin3)]:
        if ng>0: g+=1;prod*=v
    if cmin4<1.0: g+=1;prod*=cmin4
    return nm, prod**(1.0/g) if g>0 else 0.0, g


def ocean_nite_full(px, vza, sfctmp, sh_ocean, uniform, tb, qa, thr, pf, btclr):
    m31,m32,m20=px[B11],px[B12],px[B38]; m26,m27=px[B21],px[B20]
    tv=m31-m32; ng1=ng2=0; cmin1=cmin2=1.0; nm=0
    nt=_get_thr('nobt11',thr)
    if m31>0 and m31>=nt[1]: fylat_core.set_bit(tb,13); fylat_core.set_qa_bit(qa,13); cmin1*=fylat_core.conf_test(m31,nt[0],nt[2],nt[3],nt[1],1); ng1+=1
    if m31>0 and m32>0 and m31<_get_thr('pfmft_11maxthre',pf)[0] and (btclr[5]-btclr[6])>_get_thr('pfmft_btd_min',pf)[0]: fylat_core.set_bit(tb,14); fylat_core.set_qa_bit(qa,14); ng1+=1
    if m31>0 and m32>0 and tv<=_get_thr('nfmft_maxthre',pf)[0]: fylat_core.set_bit(tb,15); fylat_core.set_qa_bit(qa,15); ng1+=1
    if sfctmp>0:
        mp=260.0+2.0*round(tv)+vza**4*3.0
        if m31-sfctmp<mp: fylat_core.set_bit(tb,27); fylat_core.set_qa_bit(qa,27); cmin1*=fylat_core.conf_test(m31-sfctmp,230.0,mp,4.0,mp-3.0,1); ng1+=2
    nm+=ng1
    if m31>0 and m32>0 and m20>0 and m31-m20<fylat_core.trispc(tv): fylat_core.set_bit(tb,18); fylat_core.set_qa_bit(qa,18); ng2+=1
    if m31>0 and m32>0 and vza>0:
        cv=math.cos(math.radians(vza))
        if abs(cv)>1e-6:
            dtv=fylat_core.tview(1,1.0/cv,m31)
            if dtv>=0.1:
                if tv<=dtv: fylat_core.set_bit(tb,18)
                else: fylat_core.clear_bit(tb,18)
            ng2+=1
    d11=_get_thr('no11_4lo',thr)
    if m31>0 and m20>0 and m31-m20<=d11[1]: fylat_core.set_bit(tb,19); fylat_core.set_qa_bit(qa,19); ng2+=1
    n86=_get_thr('no86_73',thr)
    if m26>0 and m27>0 and m26-m27>n86[1]: fylat_core.set_bit(tb,29); fylat_core.set_qa_bit(qa,29); cmin2*=fylat_core.conf_test(m26-m27,n86[0],n86[2],n86[3],n86[1],1); ng2+=1
    nm+=ng2
    g,prod=0,1.0
    if ng1>0: g+=1;prod*=cmin1
    if ng2>0: g+=1;prod*=cmin2
    return nm, prod**(1.0/g) if g>0 else 0.0, g


def land_day_full(px, vza, visusd, vrused, hi_elev, tb, qa, thr, pf, btclr):
    m31,m32,m20=px[B11],px[B12],px[B38]; b2,b4,b5=px[B066],px[B055],px[B138]
    tv=m31-m32; ng1=ng2=ng3=0; cmin1=cmin2=cmin3=1.0; nm=0
    if m31>0 and m32>0 and m31<_get_thr('pfmft_11maxthre',pf)[0] and (btclr[5]-btclr[6])>_get_thr('pfmft_btd_min',pf)[0]: fylat_core.set_bit(tb,14); fylat_core.set_qa_bit(qa,14); ng1+=1
    if m31>0 and m32>0 and tv<=_get_thr('nfmft_maxthre',pf)[0]: fylat_core.set_bit(tb,15); fylat_core.set_qa_bit(qa,15); ng1+=1
    nm+=ng1
    if m31>0 and m32>0 and vza>0:
        cv=math.cos(math.radians(vza))
        if abs(cv)>1e-6:
            dtv=fylat_core.tview(1,1.0/cv,m31)
            if dtv>=0.1:
                if tv<=dtv: fylat_core.set_bit(tb,18); fylat_core.set_qa_bit(qa,18)
                else: fylat_core.clear_bit(tb,18)
            ng2+=1
    dl4=_get_thr('dl11_4lo',thr)
    if visusd and m31>0 and m20>0 and m31-m20>=dl4[1]: fylat_core.set_bit(tb,19); fylat_core.set_qa_bit(qa,19); ng2+=1
    nm+=ng2
    dr1=_get_thr('dlref1',thr)
    if visusd and b2>0 and b2<=dr1[1]: fylat_core.set_bit(tb,20); fylat_core.set_qa_bit(qa,20); cmin3*=fylat_core.conf_test(b2,dr1[0],dr1[2],dr1[3],dr1[1],1); ng3+=1
    dv=_get_thr('dlvrat',thr)
    if visusd and vrused and b2>0 and b4>0 and b4/b2<=dv[1]: fylat_core.set_bit(tb,21); fylat_core.set_qa_bit(qa,21); cmin3*=fylat_core.conf_test(b4/b2,dv[0],dv[2],dv[3],dv[1],1); ng3+=1
    nm+=ng3
    cmin4=1.0; dr3=_get_thr('dlref3',thr)
    if not hi_elev and visusd and b5>0 and b5<=dr3[1]: fylat_core.set_bit(tb,16); fylat_core.set_qa_bit(qa,16); cmin4*=fylat_core.conf_test(b5,dr3[0],dr3[2],dr3[3],dr3[1],1)
    dci=_get_thr('dltci',thr)
    if not hi_elev and visusd and b5>0 and b5>=dci[1] and b5<dci[0]: fylat_core.clear_bit(tb,9)
    g,prod=0,1.0
    if ng1>0: g+=1;prod*=cmin1
    if ng2>0: g+=1;prod*=cmin2
    if ng3>0: g+=1;prod*=cmin3
    if cmin4<1.0: g+=1;prod*=cmin4
    return nm, prod**(1.0/g) if g>0 else 0.0, g


def day_snow_full(px, vza, visusd, hi_elev, tb, qa, thr, pf, btclr):
    m31,m32,m20=px[B11],px[B12],px[B38]; b5=px[B138]; tv=m31-m32; ng1=ng2=0; cmin1=cmin2=1.0; nm=0
    if m31>0 and m32>0 and m31<_get_thr('pfmft_11maxthre',pf)[0] and (btclr[5]-btclr[6])>_get_thr('pfmft_btd_min',pf)[0]: fylat_core.set_bit(tb,14); fylat_core.set_qa_bit(qa,14); ng1+=1
    if m31>0 and m32>0 and tv<=_get_thr('nfmft_maxthre',pf)[0]: fylat_core.set_bit(tb,15); fylat_core.set_qa_bit(qa,15); ng1+=1
    nm+=ng1
    if m31>0 and m32>0 and vza>0:
        cv=math.cos(math.radians(vza))
        if abs(cv)>1e-6:
            dtv=fylat_core.tview(1,1.0/cv,m31)
            if dtv>=0.1:
                df=dtv+_get_thr('ds11_12adj',thr)[0]
                if tv<=df: fylat_core.set_bit(tb,18); fylat_core.set_qa_bit(qa,18)
                else: fylat_core.clear_bit(tb,18)
            ng2+=1
    d4=_get_thr('ds4_11',thr); d4h=_get_thr('ds4_11hel',thr)
    thr_4_11 = d4h if hi_elev else d4
    if m31>0 and m20>0 and m20-m31<=thr_4_11[1]: fylat_core.set_bit(tb,19); fylat_core.set_qa_bit(qa,19); ng2+=1
    nm+=ng2
    cmin4=1.0; dr3=_get_thr('dsref3',thr)
    if not hi_elev and visusd and b5>0 and b5<=dr3[1]: fylat_core.set_bit(tb,16); fylat_core.set_qa_bit(qa,16); cmin4*=fylat_core.conf_test(b5,dr3[0],dr3[2],dr3[3],dr3[1],1)
    dc=_get_thr('dstci',thr)
    if not hi_elev and visusd and b5>0 and b5>=dc[1] and b5<dc[0]: fylat_core.clear_bit(tb,9)
    g,prod=0,1.0
    if ng1>0: g+=1;prod*=cmin1
    if ng2>0: g+=1;prod*=cmin2
    if cmin4<1.0: g+=1;prod*=cmin4
    return nm, prod**(1.0/g) if g>0 else 0.0, g


def nite_snow_full(px, vza, lnd, tb, qa, thr, pf, btclr):
    m31,m32,m20=px[B11],px[B12],px[B38]; m27=px[B20]; tv=m31-m32; ng1=ng2=0; cmin2=1.0; nm=0
    if m31>0 and m32>0 and m31<_get_thr('pfmft_11maxthre',pf)[0] and (btclr[5]-btclr[6])>_get_thr('pfmft_btd_min',pf)[0]: fylat_core.set_bit(tb,14); fylat_core.set_qa_bit(qa,14); ng1+=1
    if m31>0 and m32>0 and tv<=_get_thr('nfmft_maxthre',pf)[0]: fylat_core.set_bit(tb,15); fylat_core.set_qa_bit(qa,15); ng1+=1
    nm+=ng1
    if m31>0 and m32>0 and vza>0:
        cv=math.cos(math.radians(vza))
        if abs(cv)>1e-6:
            dtv=fylat_core.tview(1,1.0/cv,m31)
            if dtv>=0.1:
                df=dtv+_get_thr('ns11_12adj',thr)[0]
                if tv<=df: fylat_core.set_bit(tb,18); fylat_core.set_qa_bit(qa,18)
                else: fylat_core.clear_bit(tb,18)
            ng2+=1
    n4=_get_thr('ns11_4lo',thr)
    if m31>0 and m20>0 and m31-m20<=n4[1]: fylat_core.set_bit(tb,19); fylat_core.set_qa_bit(qa,19); ng2+=1
    if m27>0 and m31>0 and m27-m31<=-2.0: fylat_core.set_bit(tb,23); fylat_core.set_qa_bit(qa,23); ng2+=1
    nm+=ng2
    cmin5=1.0; ng5=0
    n5=_get_thr('ns4_12hi',thr)
    if m20>0 and m32>0 and m20-m32<=n5[1]: fylat_core.set_bit(tb,17); fylat_core.set_qa_bit(qa,17); cmin5*=fylat_core.conf_test(m20-m32,n5[0],n5[2],n5[3],n5[1],1); ng5+=1
    nm+=ng5
    g,prod=0,1.0
    if ng2>0: g+=1;prod*=cmin2
    if ng5>0: g+=1;prod*=cmin5
    return nm, prod**(1.0/g) if g>0 else 0.0, g


def process_pixel_full(px, data, il, ie, thresholds):
    sfc = classify_surface_full(data, il, ie)
    tb = np.zeros(6, dtype=np.uint8)
    qa = np.zeros(10, dtype=np.uint8)
    nmt, conf = 0, 0.0

    # pxinit: set fail-safe bits (8=smoke, 9=thin cirrus solar, 10=shadow, 11=thin cirrus IR, 28=dust)
    for b in [8, 9, 10, 11, 28]:
        fylat_core.set_bit(tb, b)
    nbands = 20  # assume 20 bands available
    nmt = 0

    # Get clear-sky BT from Fortran intermediate data
    tbb = data['tbb_ir']
    btclr = np.array([
        0.0,
        float(tbb[0,il,ie]) if tbb[0,il,ie]>0 else 0.0,
        float(tbb[1,il,ie]) if tbb[1,il,ie]>0 else 0.0,
        float(tbb[2,il,ie]) if tbb[2,il,ie]>0 else 0.0,
        float(tbb[3,il,ie]) if tbb[3,il,ie]>0 else 0.0,
        float(tbb[4,il,ie]) if tbb[4,il,ie]>0 else 0.0,
        float(tbb[5,il,ie]) if tbb[5,il,ie]>0 else 0.0,
    ], dtype=np.float32)

    sfctmp = float(data['tbb_ir'][4, il, ie])  # use 11um clear-sky BT as SST proxy
    vza = float(data['saza'][il, ie])

    scene = 'skip'
    if sfc['polar'] or px[B11] <= 0:
        pass
    elif sfc['day'] and (sfc['snow'] or sfc['ice']):
        thr_s = thresholds.get('day_snow', {})
        pf = thresholds.get('pfmft', {})
        nmt2, conf, _ = day_snow_full(px, vza, True, False, tb, qa, thr_s, pf, btclr)
        nmt += nmt2; scene = 'day_snow_full'
    elif not sfc['day'] and (sfc['snow'] or sfc['ice']):
        thr_s = thresholds.get('nite_snow', {})
        pf = thresholds.get('pfmft', {})
        nmt2, conf, _ = nite_snow_full(px, vza, sfc['land'], tb, qa, thr_s, pf, btclr)
        nmt += nmt2; scene = 'nite_snow_full'
    elif sfc['day'] and sfc['water']:
        thr_o = thresholds.get('ocean_day', {})
        pf = thresholds.get('pfmft', {})
        nmt2, conf, _ = ocean_day_full(px, vza, False, True, sfctmp, 0.0, False,
                                       tb, qa, thr_o, pf, {}, btclr)
        nmt += nmt2; scene = 'ocean_day_full'
    elif not sfc['day'] and sfc['water']:
        thr_o = thresholds.get('ocean_nite', {})
        pf = thresholds.get('pfmft', {})
        nmt2, conf, _ = ocean_nite_full(px, vza, sfctmp, False, True,
                                        tb, qa, thr_o, pf, btclr)
        nmt += nmt2; scene = 'ocean_nite_full'
    elif sfc['day'] and sfc['land']:
        thr_l = thresholds.get('land_day', {})
        pf = thresholds.get('pfmft', {})
        nmt2, conf, _ = land_day_full(px, vza, True, True, False,
                                      tb, qa, thr_l, pf, btclr)
        nmt += nmt2; scene = 'land_day_full'

    # Post-processing (matching Fortran pixel loop steps)
    fylat_core.proc_path(sfc['water'], sfc['land'], sfc['day'],
                         sfc['ice'], sfc['snow'], False, sfc['coast'],
                         sfc['desert'], False, False, tb)
    fylat_core.set_unused_bits(tb)
    fylat_core.set_confdnc(conf, tb)
    fylat_core.set_quality_A(nmt, nbands, 0, qa)

    # fill_bit_pixel: set bit 0 if nmt > 0, handle sun glint QA reduction
    if nmt > 0:
        fylat_core.set_bit(tb, 0)
        # QA: if nmt>=7 → high quality, elif nmt>=3 → medium, else → low
        if nmt >= 7:
            for i in range(4): fylat_core.set_qa_bit(qa, i)
        elif nmt >= 3:
            fylat_core.set_qa_bit(qa, 0); fylat_core.set_qa_bit(qa, 2); fylat_core.set_qa_bit(qa, 3)
        else:
            fylat_core.set_qa_bit(qa, 0); fylat_core.set_qa_bit(qa, 3)

    return tb, qa, nmt, conf, scene


def run_full(l1b, geo, intermed, output=None, n_pixels=None):
    t0 = time.perf_counter()
    data = load_data_full(l1b, geo, intermed)
    nL, nE = data['nLine'], data['nElem']
    print(f"Orbit: {nL}x{nE} = {nL*nE/1e6:.1f}M px, load={time.perf_counter()-t0:.1f}s")

    thresholds = {}
    for s in ['ocean_day', 'ocean_nite', 'land_day', 'land_nite', 'day_snow', 'nite_snow']:
        thresholds[s] = get_scene_thresholds(s)
    thresholds['pfmft'] = thresholds['ocean_day']

    if n_pixels:
        step = max(1, int(math.sqrt(nL * nE / n_pixels)))
        lines = range(100, nL - 100, step)
        elems = range(100, nE - 100, step)
    else:
        lines, elems = range(nL), range(nE)

    cm = np.zeros((6, nL, nE), dtype=np.uint8)
    qarr = np.zeros((10, nL, nE), dtype=np.uint8)
    sfcnt = {}
    n_proc = 0
    t_p = time.perf_counter()

    for il in lines:
        for ie in elems:
            px = build_pxldat_full(data, il, ie)
            if px[B11] <= 0: continue
            n_proc += 1
            tb, qb, nmt, conf, scene = process_pixel_full(px, data, il, ie, thresholds)
            sfcnt[scene] = sfcnt.get(scene, 0) + 1
            cm[:, il, ie] = tb
            qarr[:, il, ie] = qb

    t2 = time.perf_counter()
    stats = {
        'processed': n_proc,
        'scene_counts': sfcnt,
        'timing': {'load': t_p - t0, 'process': t2 - t_p, 'total': t2 - t0,
                   'px_per_sec': n_proc / (t2 - t_p) if n_proc > 0 else 0}
    }

    if output:
        os.makedirs(os.path.dirname(output) or '.', exist_ok=True)
        with h5py.File(output, 'w') as f:
            f.create_dataset('Cloud_Mask', data=cm, compression='gzip')
            f.create_dataset('Quality_Assurance', data=qarr, compression='gzip')

    data['l1b'].close(); data['geo'].close(); data['im'].close()
    return {'cm': cm, 'qa': qarr, 'stats': stats}
