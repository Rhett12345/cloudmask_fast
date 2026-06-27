# FYLAT Cloud Mask — 20220803 Batch Validation Report (v2.4.5)

**Date**: 2026-06-27  
**Data**: FY-3D MERSI-II full-day (63 time slots, 00:45–23:15 UTC, ~5-min cadence)  
**Reference**: Aqua MYD35 L2 (061 collection), ±15 min temporal window  
**Calibration**: Business (onboard) + Recalibration (external monthly) — 126 inversion runs
**Cores**: 8 parallel

## 1. Performance

| Metric | Value |
|--------|-------|
| Time slots | 63 |
| Tasks (2 calibrations) | 126 |
| Success rate | 126/126 |
| Wall clock | ~28 min |
| Per-task latency | 71–286s (median ~80s) |
| Output per slot | ~63 MB CLM HDF5 |
| Total output | ~7.9 GB |

Processing time is NWP-I/O-bound; scene type has minor effect. Sporadic outliers
(1805, 2130 business >280s) likely caused by filesystem contention.

## 2. MYD35 Overlap Coverage

- **49/63 slots** had valid MYD35 overlap (>5% pixel overlap within ±15 min)
- **14 slots failed**: MYD35 scan gaps at 15:00–15:05 UTC, 18:10–18:25 UTC, and
  21:45–23:15 UTC. The overlapping pixels fell below the 100-pixel minimum in
  several twilight scenes.
- Best overlap window: 07:00–09:30 UTC (50–90% overlap, Δt = 0 min)

## 3. Stratified Validation (SZA day/night + latitude + surface type)

### 3.1 By Day/Night (Solar Zenith Angle)

| | Strata | Mean HSS | Median HSS | Worst HSS |
|---|--------|----------|------------|-----------|
| **Day (SZA < 80°)** | 89 | **0.372** | 0.419 | −0.072 |
| **Night (SZA ≥ 90°)** | 45 | **0.272** | 0.267 | −0.070 |

Daytime significantly outperforms nighttime. Reflective solar bands (ch 1–19)
provide strong additional discrimination that IR-only nighttime tests cannot
replicate. Nighttime IR cloud tests (BTD 11–12, 8–11, 11–4) show weak
separability, especially in warm/humid atmospheres.

### 3.2 By Latitude Band

| Band | Strata | Mean HSS |
|------|--------|----------|
| Polar (≥60°) | 9 | **0.392** |
| Mid-lat (30–60°) | 42 | 0.382 |
| Tropical (<30°) | 83 | **0.310** |

Latitude gradient: colder backgrounds and drier atmospheres at high latitudes
produce better IR contrast. Tropical water vapor loading reduces BTD dynamic
range and blurs the cloud/clear boundary.

### 3.3 By Surface Type

| Surface | Strata | Mean HSS |
|---------|--------|----------|
| **Ocean** | 53 | **0.415** |
| Land | 41 | 0.320 |
| **Inland Water** (lakes, coastlines) | 40 | **0.254** |

Ocean is best — uniform emissivity, predictable clear-sky BT. Land has
variable emissivity and hot daytime backgrounds. Inland water is the
worst-case: mixed land/water pixels at coastlines and small lakes produce
unstable emissivity and sub-pixel cloud contamination.

### 3.4 Worst-Performing Strata (HSS < 0.01)

| Time | SZA | Lat | Surface | N | HSS | POD | FAR | Problem |
|------|-----|-----|---------|-----|------|------|------|---------|
| 0530 | Day | MidLat | InlandWater | 7.5k | −0.072 | 42% | 52% | Random skill |
| 1625 | Night | Tropic | Ocean | 2.5k | −0.070 | 80% | 9% | Miss + false |
| 0720 | Day | Tropic | Land | 436k | −0.000 | 100% | 0% | All-cloud (no skill) |
| 0830 | Night | MidLat | Ocean | 2.2k | 0.000 | 85% | 0% | 14.6% cloud miss |
| 1955 | Night | Tropic | Land | 53k | 0.000 | 100% | 0% | All-cloud |
| 0205 | Day | MidLat | Ocean | 19k | 0.000 | 100% | 0% | All-cloud |
| 0830 | Night | Polar | Land | 1.05M | 0.003 | 67% | 4% | 28.7% cloud miss |

**Common failure modes**:
1. **HSS ≈ 0 (all-cloud)**: Algorithm classifies everything as cloudy.
   IR thresholds too aggressive for warm backgrounds. Affects tropical
   land/inland-water at both day and night.
2. **Negative HSS**: Worse than random. POD ~40–80%, FAR ~10–50%.
   Typically inland water pixels where emissivity assumptions break.
3. **Large cloud miss (bias < −15%)**: 0830 orbit misses 15–29% of MYD35
   clouds across multiple surface types. Possible NWP degradation or
   scan-geometry effect specific to this orbit.

### 3.5 Best-Performing Strata (HSS > 0.7)

| Time | SZA | Lat | Surface | N | HSS | POD | FAR |
|------|-----|-----|---------|-----|------|------|------|
| 0740 | Day | Polar | Ocean | 1.67M | **0.822** | 98% | 3% |
| 0920 | Day | Polar | Ocean | 1.93M | **0.807** | 99% | 2% |
| 0410 | Day | MidLat | Ocean | 501k | 0.716 | 98% | 2% |
| 0045 | Day | MidLat | Ocean | 2.05M | 0.710 | 86% | 10% |

Ocean + cold background + daytime multi-spectral = optimal conditions.

## 4. BTD Threshold Analysis

Key BTD separability (cloudy vs clear means, in σ units):

| Scene | BTD 11–12 | BTD 8–11 | BT11 raw |
|-------|-----------|----------|----------|
| Cold (<230K) Ocean | 0.4σ | 0.5σ | 0.8σ |
| Warm (>270K) Land | **1.1σ** | 1.3σ | **1.2σ** |
| Warm (>270K) Ocean | 0.5σ | 1.6σ | 0.1σ |
| Mod (250–270K) Ocean | **1.8σ** | 1.6σ | 1.1σ |

- **Best IR discriminant**: BTD 11–12 for warm/moderate scenes (1.1–1.8σ)
- **Worst IR discriminant**: All BTDs for cold scenes (<230K, all <0.8σ).
  Cold cloud tops have BT values too close to clear-sky cold surface BT.
- Current thresholds differ significantly from MYD35-optimal values,
  especially for BTD 11–12 (current ~3K vs optimal ~30–50K in warm scenes).

## 5. Calibration Comparison

Business (onboard) vs Recalibration (external monthly) cloud masks are
**nearly identical** (ΔHSS < 0.01 across all strata). The recalibration
only affects solar reflective bands (ch 1–19); cloud detection at night
uses IR channels exclusively. Even in daytime strata, the cloud mask
decision tree relies primarily on IR tests and BTD ratios that are
insensitive to absolute VIS calibration.

## 6. Key Issues (Priority-Ordered)

1. **[P0] Tropical daytime inland water over-detection**
   - FAR reaches 52%, HSS negative in multiple time slots
   - Root cause: coastline/lake mixed-pixel emissivity errors in RTM
   - Suggested: relax cloud thresholds near coastlines, add inland water mask check

2. **[P0] Nighttime all-cloud bias in warm scenes**
   - Multiple tropical/mid-lat night slots show HSS ≈ 0 (100% cloud, no skill)
   - Root cause: IR BTD thresholds too aggressive for warm (BT > 270K) scenes
   - Suggested: raise BT11 and BTD 11–12 thresholds for warm-background night

3. **[P1] 0830 orbit systematic cloud miss**
   - Consistent 15–29% cloud miss across land/ocean/polar
   - Possible NWP degradation: 0830 is +2.5h from the 06Z NWP cycle
   - Suggested: verify NWP interpolation quality at 0830, or use 09Z cycle

4. **[P2] Cold cloud (<230K) poor IR separability**
   - All BTDs show <0.8σ separation below 230K
   - Fundamental physical limitation: cold surface ≈ cold cloud top BT
   - Suggested: rely more on spatial texture tests for polar/cold scenes

5. **[P3] NWP single-cycle limitation**
   - All 63 time slots use the same NWP pair (t06z f018/f021)
   - Slots far from 06Z (e.g., 18Z–23Z) use 12+ hour old forecasts
   - Suggested: cycle-aware NWP selection per time slot

## 7. Tools Added

- `diagnostics/batch_myd35_overlap.py` — per-slot MYD35 matching + overall stats
- `diagnostics/batch_geo_analysis.py` — SZA day/night + lat-band + surface stratified validation
- `batch_run.py` — multi-date batch inversion runner (reusable for other dates)
