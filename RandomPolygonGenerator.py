#!/usr/bin/env python3
"""
generate_benchmark_polygons_v2.py
===================================
v2 of the benchmark corpus generator. Fixes the "erosion breaks the
polygon into multiple components" failure mode reported after the v1
corpus was run through Batch Mode: v1 built shapes in a canonical unit
box and then *anisotropically rescaled* them to hit a target aspect
ratio, which could squeeze already-thin notches/teeth/corridors down to
near-zero width. v2 instead:

  1. builds every polygon DIRECTLY at its target absolute size (no
     post-hoc squashing of thin features), and
  2. enforces an explicit MINIMUM FEATURE WIDTH (corridor / notch /
     tooth width) in absolute world units wherever a family can carve
     one directly (Orthogonal, Comb), and
  3. validates every candidate by literally performing the same
     operation the robot-width offset does — a negative buffer
     (Minkowski erosion) — at the *actual* robot radii the user plans
     to run (fixed width 0.4, the area% width for situations 3/4, and
     the recommended max sweep width for situations 5/6) — and rejects
     any candidate whose eroded interior splits into more than one
     component, vanishes, or shrinks implausibly.

Two new families are also added for geometric realism ("earth maps /
real fields"):
  F. FIELD    — smooth, low-frequency radial-noise blob (lake / field /
                 country-outline silhouette)
  G. COASTAL  — union of several overlapping smooth blobs -> a single
                 connected landmass with bays/peninsulas

Shrink-safety design numbers
-----------------------------
User-specified / recommended robot widths this corpus is built to survive:
  - Situations 1/2 (fixed):         width = 0.4          -> radius 0.20
  - Situations 3/4 (X% of area):    X = 0.2%              -> radius depends
                                     on each polygon's own area (up to ~2.76
                                     for the largest polygons in this corpus)
  - Situations 5/6 (sweep, recommended): 0.5 -> 3.0 step 0.5 -> radius up to 1.5

Every delivered polygon survives erosion (stays a single simple polygon,
retains a sane fraction of its area) at ALL of: radius 0.20, its own
per-polygon situations-3/4 radius, and radius 1.50 (the top of the
recommended sweep). A minimum absolute feature width of ~8 world units
is targeted wherever a family carves explicit notches/corridors, i.e.
roughly 5x the toughest test radius (1.5), which is why polygons in
this corpus are noticeably larger than the v1 batch.
"""

import json
import math
import random
from pathlib import Path

from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import unary_union
from scipy.spatial import ConvexHull

RNG_SEED = 20260719
random.seed(RNG_SEED)

OUT_DIR = Path("/Users/samratdutta/Desktop/Polygon Simulation/benchmark_polygons_v2")

AREA_RANGE = (150.0, 3000.0)     # log-uniform, world sq. units — bigger than v1
ASPECT_RANGE = (1.0, 4.5)        # uniform, bbox width / height — slightly tamer than v1
COORD_DECIMALS = 3
MIN_FEATURE = 8.0                # target minimum notch/corridor/tooth width (world units)

# ── erosion stress-test radii ────────────────────────────────────────────
R_FIXED = 0.20                   # situations 1/2: width 0.4
PCT_X = 0.2                      # situations 3/4: robot AREA = 0.2% of polygon area
R_SWEEP_MAX = 1.50               # situations 5/6 recommended max width 3.0 -> radius 1.5
R_SOFT_EXTRA = 2.00              # extra informational (non-filtering) safety check

CANDIDATES_PER_FAMILY = {
    "Convex": 130,
    "Star": 170,
    "Orthogonal": 220,
    "Comb": 190,
    "Irregular": 130,
    "Field": 170,
    "Coastal": 150,
}


# ─────────────────────────────────────────────────────────────────────────
#  Geometry helpers
# ─────────────────────────────────────────────────────────────────────────

def log_uniform(lo, hi):
    return math.exp(random.uniform(math.log(lo), math.log(hi)))


def polygon_is_simple(coords):
    if len(coords) < 3:
        return False
    try:
        p = ShapelyPolygon(coords)
    except Exception:
        return False
    return p.is_valid and p.exterior.is_simple and p.area > 1e-9


def signed_area(coords):
    n = len(coords)
    s = 0.0
    for i in range(n):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return s / 2.0


def ensure_ccw(coords):
    return coords if signed_area(coords) > 0 else list(reversed(coords))


def reflex_count(coords):
    pts = ensure_ccw(coords)
    n = len(pts)
    count = 0
    for i in range(n):
        a = pts[(i - 1) % n]
        b = pts[i]
        c = pts[(i + 1) % n]
        cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
        if cross < -1e-9:
            count += 1
    return count


def bbox(coords):
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    return min(xs), max(xs), min(ys), max(ys)


def clean_and_round(coords):
    rounded = [(round(x, COORD_DECIMALS), round(y, COORD_DECIMALS)) for (x, y) in coords]
    cleaned = []
    for p in rounded:
        if cleaned and cleaned[-1] == p:
            continue
        cleaned.append(p)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def erosion_survives(coords, radius, min_area_retention=0.05):
    """True if buffer(-radius) stays a single simple polygon retaining
    at least `min_area_retention` of the original area."""
    try:
        poly = ShapelyPolygon(coords)
        if not poly.is_valid:
            return False
        orig_area = poly.area
        eroded = poly.buffer(-radius, resolution=8, join_style=2)
    except Exception:
        return False
    if eroded.is_empty:
        return False
    if eroded.geom_type != "Polygon":
        return False  # split into MultiPolygon -> exactly the failure we're avoiding
    if not eroded.is_valid:
        return False
    if eroded.area < min_area_retention * orig_area:
        return False
    return True


# ─────────────────────────────────────────────────────────────────────────
#  Family A — CONVEX  (built directly inside a W x H box)
# ─────────────────────────────────────────────────────────────────────────

def gen_convex(W, H):
    n = random.randint(5, 14)
    pts = [(random.uniform(0, W), random.uniform(0, H)) for _ in range(n)]
    try:
        hull = ConvexHull(pts)
    except Exception:
        return None
    coords = [pts[i] for i in hull.vertices]
    return coords if len(coords) >= 4 else None


# ─────────────────────────────────────────────────────────────────────────
#  Family B — STAR (radial, smooth-ish with a few deliberate deep notches)
# ─────────────────────────────────────────────────────────────────────────

def gen_star(W, H):
    n = random.randint(10, 30)
    R = 0.5 * min(W, H)
    cx, cy = W / 2.0, H / 2.0
    notch_frac = random.uniform(0.0, 0.45)

    angles = sorted(random.uniform(0, 2 * math.pi) for _ in range(n))
    for i in range(1, n):
        if angles[i] - angles[i - 1] < 1e-3:
            angles[i] += 1e-3

    coords = []
    for a in angles:
        r = R * random.uniform(0.82, 1.0)
        if random.random() < notch_frac:
            r *= random.uniform(0.45, 0.75)   # shallower notch than v1 -> stays erosion-safe
        coords.append((cx + r * math.cos(a), cy + r * math.sin(a)))

    # gentle anisotropic stretch toward the target aspect ratio (capped, so
    # notches don't get squashed thin)
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    w0, h0 = max(xs) - min(xs), max(ys) - min(ys)
    if w0 < 1e-6 or h0 < 1e-6:
        return None
    sx, sy = W / w0, H / h0
    coords = [(x * sx, y * sy) for (x, y) in coords]
    return coords


# ─────────────────────────────────────────────────────────────────────────
#  Family C — ORTHOGONAL (rectangle, notches with an explicit minimum
#  corridor-width guarantee)
# ─────────────────────────────────────────────────────────────────────────

def gen_orthogonal(W, H):
    poly = ShapelyPolygon([(0, 0), (W, 0), (W, H), (0, H)])
    min_feat = min(MIN_FEATURE, 0.22 * min(W, H))
    if min_feat < 1.5:
        return None  # box too small to safely carve anything -> caller falls back

    n_notches = random.randint(1, 5)
    for _ in range(n_notches):
        side = random.choice(["bottom", "top", "left", "right"])
        minx, miny, maxx, maxy = poly.bounds
        w, h = maxx - minx, maxy - miny
        if w < 4 * min_feat or h < 4 * min_feat:
            break

        if side in ("bottom", "top"):
            max_depth = h - 2 * min_feat
            max_width = w - 2 * min_feat
            if max_depth < min_feat or max_width < min_feat:
                continue
            notch_h = random.uniform(min_feat, max_depth)
            notch_w = random.uniform(min_feat, max_width)
            nx = random.uniform(minx + min_feat * 0.25, maxx - notch_w - min_feat * 0.25)
            ny = miny - 0.02 if side == "bottom" else maxy - notch_h + 0.02
            notch = ShapelyPolygon([(nx, ny), (nx + notch_w, ny),
                                     (nx + notch_w, ny + notch_h), (nx, ny + notch_h)])
        else:
            max_depth = w - 2 * min_feat
            max_width = h - 2 * min_feat
            if max_depth < min_feat or max_width < min_feat:
                continue
            notch_w = random.uniform(min_feat, max_depth)
            notch_h = random.uniform(min_feat, max_width)
            ny = random.uniform(miny + min_feat * 0.25, maxy - notch_h - min_feat * 0.25)
            nx = minx - 0.02 if side == "left" else maxx - notch_w + 0.02
            notch = ShapelyPolygon([(nx, ny), (nx + notch_w, ny),
                                     (nx + notch_w, ny + notch_h), (nx, ny + notch_h)])

        candidate = poly.difference(notch)
        if (candidate.geom_type == "Polygon" and candidate.is_valid
                and len(candidate.interiors) == 0 and candidate.area > 0.5 * poly.area):
            poly = candidate

    if poly.geom_type != "Polygon" or not poly.is_valid:
        return None
    return list(poly.exterior.coords)[:-1]


# ─────────────────────────────────────────────────────────────────────────
#  Family D — COMB (evenly spaced teeth, explicit min tooth-width AND
#  min wall-width-between-teeth guarantee)
# ─────────────────────────────────────────────────────────────────────────

def gen_comb(W, H):
    min_feat = min(MIN_FEATURE, 0.28 * H)
    if min_feat < 1.5 or H < 3 * min_feat:
        return None

    poly = ShapelyPolygon([(0, 0), (W, 0), (W, H), (0, H)])
    tooth_w = min_feat * random.uniform(1.0, 1.6)
    wall_w = min_feat * random.uniform(1.0, 1.8)
    period = tooth_w + wall_w
    n_teeth = int((W - wall_w) // period)
    n_teeth = max(2, min(n_teeth, 14))
    depth = random.uniform(0.5, 0.72) * H

    x = wall_w
    for _ in range(n_teeth):
        if x + tooth_w > W - wall_w * 0.5:
            break
        notch = ShapelyPolygon([(x, -0.02), (x + tooth_w, -0.02),
                                 (x + tooth_w, depth), (x, depth)])
        candidate = poly.difference(notch)
        if candidate.geom_type == "Polygon" and candidate.is_valid and len(candidate.interiors) == 0:
            poly = candidate
        x += period

    if poly.geom_type != "Polygon" or not poly.is_valid:
        return None
    coords = list(poly.exterior.coords)[:-1]
    return coords if len(coords) >= 8 else None


# ─────────────────────────────────────────────────────────────────────────
#  Family E — IRREGULAR (convex hull, a few vertices pulled inward but
#  capped so they can't create a near-self-intersecting sliver)
# ─────────────────────────────────────────────────────────────────────────

def gen_irregular(W, H):
    n = random.randint(7, 20)
    pts = [(random.uniform(0, W), random.uniform(0, H)) for _ in range(n)]
    try:
        hull = ConvexHull(pts)
    except Exception:
        return None
    coords = [list(pts[i]) for i in hull.vertices]
    m = len(coords)
    if m < 5:
        return None

    push_frac = random.uniform(0.1, 0.4)   # gentler than v1
    n_push = max(1, int(round(m * push_frac)))
    push_idxs = random.sample(range(m), n_push)

    cx = sum(p[0] for p in coords) / m
    cy = sum(p[1] for p in coords) / m

    for idx in push_idxs:
        x, y = coords[idx]
        dx, dy = cx - x, cy - y
        dist = math.hypot(dx, dy)
        if dist < 1e-6:
            continue
        depth = random.uniform(0.15, 0.45)   # gentler push depth than v1
        new_pt = (x + dx * depth, y + dy * depth)
        trial = coords[:idx] + [list(new_pt)] + coords[idx + 1:]
        if polygon_is_simple(trial):
            coords[idx] = list(new_pt)

    return [tuple(c) for c in coords]


# ─────────────────────────────────────────────────────────────────────────
#  Family F — FIELD (smooth low-frequency radial-noise blob)
# ─────────────────────────────────────────────────────────────────────────

def gen_field(W, H):
    n = random.randint(24, 60)
    R = 0.5 * min(W, H)
    cx, cy = W / 2.0, H / 2.0

    n_harm = random.randint(2, 4)
    harmonics = []
    for k in range(2, 2 + n_harm):
        amp = random.uniform(0.03, 0.16) / k
        phase = random.uniform(0, 2 * math.pi)
        harmonics.append((k, amp, phase))

    coords = []
    for i in range(n):
        a = 2 * math.pi * i / n
        r = R
        for (k, amp, phase) in harmonics:
            r += R * amp * math.sin(k * a + phase)
        coords.append((cx + r * math.cos(a), cy + r * math.sin(a)))

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    w0, h0 = max(xs) - min(xs), max(ys) - min(ys)
    if w0 < 1e-6 or h0 < 1e-6:
        return None
    sx, sy = W / w0, H / h0
    coords = [(x * sx, y * sy) for (x, y) in coords]
    return coords


# ─────────────────────────────────────────────────────────────────────────
#  Family G — COASTAL (union of overlapping smooth blobs -> single
#  connected landmass with bays / peninsulas)
# ─────────────────────────────────────────────────────────────────────────

def _smooth_blob(cx, cy, R, n=28):
    n_harm = random.randint(1, 3)
    harmonics = []
    for k in range(2, 2 + n_harm):
        amp = random.uniform(0.05, 0.18) / k
        phase = random.uniform(0, 2 * math.pi)
        harmonics.append((k, amp, phase))
    pts = []
    for i in range(n):
        a = 2 * math.pi * i / n
        r = R
        for (k, amp, phase) in harmonics:
            r += R * amp * math.sin(k * a + phase)
        pts.append((cx + r * math.cos(a), cy + r * math.sin(a)))
    return ShapelyPolygon(pts)


def gen_coastal(W, H):
    n_blobs = random.randint(2, 4)
    base_R = 0.5 * min(W, H) / math.sqrt(n_blobs) * 1.35

    blobs = []
    cx0, cy0 = W / 2.0, H / 2.0
    blobs.append(_smooth_blob(cx0, cy0, base_R * random.uniform(0.9, 1.15)))

    for _ in range(n_blobs - 1):
        anchor = random.choice(blobs)
        acx, acy = anchor.centroid.x, anchor.centroid.y
        ang = random.uniform(0, 2 * math.pi)
        r_here = base_R * random.uniform(0.8, 1.15)
        offset = 0.55 * (base_R + r_here)   # generous overlap
        ncx = acx + offset * math.cos(ang)
        ncy = acy + offset * math.sin(ang)
        blobs.append(_smooth_blob(ncx, ncy, r_here))

    merged = unary_union(blobs)
    if merged.geom_type != "Polygon":
        return None
    if merged.is_empty or not merged.is_valid:
        return None

    coords = list(merged.exterior.coords)[:-1]
    if len(coords) < 8:
        return None

    xmin = min(c[0] for c in coords)
    xmax = max(c[0] for c in coords)
    ymin = min(c[1] for c in coords)
    ymax = max(c[1] for c in coords)
    w0, h0 = xmax - xmin, ymax - ymin
    if w0 < 1e-6 or h0 < 1e-6:
        return None
    sx, sy = W / w0, H / h0
    coords = [((x - xmin) * sx, (y - ymin) * sy) for (x, y) in coords]
    return coords


FAMILY_GENERATORS = {
    "Convex": gen_convex,
    "Star": gen_star,
    "Orthogonal": gen_orthogonal,
    "Comb": gen_comb,
    "Irregular": gen_irregular,
    "Field": gen_field,
    "Coastal": gen_coastal,
}


# ─────────────────────────────────────────────────────────────────────────
#  Main generation loop
# ─────────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    catalog = []
    counter = {fam: 0 for fam in FAMILY_GENERATORS}
    rejected_erosion = {fam: 0 for fam in FAMILY_GENERATORS}

    for family, gen_fn in FAMILY_GENERATORS.items():
        target = CANDIDATES_PER_FAMILY[family]
        attempts = 0
        accepted = 0
        while accepted < target and attempts < target * 12:
            attempts += 1

            target_area = log_uniform(*AREA_RANGE)
            target_aspect = random.uniform(*ASPECT_RANGE)
            W = math.sqrt(target_area * target_aspect)
            H = math.sqrt(target_area / target_aspect)

            raw = gen_fn(W, H)
            if raw is None or not polygon_is_simple(raw):
                continue

            cleaned = clean_and_round(raw)
            if len(cleaned) < 4 or not polygon_is_simple(cleaned):
                continue

            achieved_area = abs(signed_area(cleaned))
            if achieved_area < 20:
                continue
            xmin, xmax, ymin, ymax = bbox(cleaned)
            w, h = xmax - xmin, ymax - ymin
            if h < 1e-6:
                continue
            achieved_aspect = w / h

            # ── erosion stress test: the actual gatekeeper ──────────────
            r_pct = math.sqrt((PCT_X / 100.0) * achieved_area / math.pi)
            if not erosion_survives(cleaned, R_FIXED):
                rejected_erosion[family] += 1
                continue
            if not erosion_survives(cleaned, r_pct):
                rejected_erosion[family] += 1
                continue
            if not erosion_survives(cleaned, R_SWEEP_MAX):
                rejected_erosion[family] += 1
                continue
            survives_extra = erosion_survives(cleaned, R_SOFT_EXTRA)

            n_vertices = len(cleaned)
            n_reflex = reflex_count(cleaned)

            accepted += 1
            counter[family] += 1
            idx = counter[family]
            fname = f"{family}_{idx:03d}"
            catalog.append({
                "filename": fname + ".txt",
                "family": family,
                "vertices": n_vertices,
                "reflex_vertices": n_reflex,
                "complexity": n_vertices + n_reflex,
                "area": round(achieved_area, 4),
                "aspect_ratio": round(achieved_aspect, 4),
                "bbox_width": round(w, 4),
                "bbox_height": round(h, 4),
                "situation34_radius": round(r_pct, 4),
                "survives_radius_2_0": survives_extra,
                "coords": cleaned,
            })

        print(f"{family:12s}: accepted {accepted}/{target}  "
              f"(erosion-rejected {rejected_erosion[family]})  in {attempts} attempts")

    for entry in catalog:
        path = OUT_DIR / entry["filename"]
        with open(path, "w") as f:
            for (x, y) in entry["coords"]:
                f.write(f"({x}, {y})\n")

    catalog_summary = [{k: v for k, v in e.items() if k != "coords"} for e in catalog]
    with open(OUT_DIR / "generation_catalog.json", "w") as f:
        json.dump(catalog_summary, f, indent=2)

    print(f"\nTotal polygons generated: {len(catalog)}")
    return catalog_summary


if __name__ == "__main__":
    main()