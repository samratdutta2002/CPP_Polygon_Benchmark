# Coverage Path Planning — Randomised Polygon Benchmark

This repository contains the randomised polygon benchmark used to evaluate
generalisation of the rotation-based boustrophedon decomposition coverage
path planning (CPP) algorithm described in:

> [Samrat Dutta, Soumen Roy, Rajat Kumar Pal], "Rotation-Optimised Decomposition of
> Irregular Concave Polygons for Coverage Path Planning".

## Contents

├── polygons/                            943 simple polygon files, "(x, y)" per line
├── results/
│   ├── Situation_1.xlsx                 Fixed width, rotation OFF
│   ├── Situation_2.xlsx                 Fixed width, rotation ON
│   ├── Situation_3.xlsx                 Proportional width, rotation OFF
│   └── Situation_4.xlsx                 Proportional width, rotation ON
├── heatmaps/
│   ├── Situation_1/ ... Situation_4/    coverage heatmap PNG per polygon
├── generate_benchmark_polygons.py       corpus generation script
└── README.md

## Benchmark construction

The corpus contains 943 simple (non-self-intersecting) polygons spanning
seven geometrically distinct shape families, plus a small number of
independently contributed polygons:

- **Convex** — convex hull of random points (zero-reflex control group)
- **Star** — angularly sorted radial vertices with a controlled subset of
  radii reduced to introduce reflex vertices
- **Orthogonal** — axis-aligned rectangle with random rectangular notches
  removed (staircase/L/T shapes)
- **Comb** — rectangle with evenly spaced rectangular teeth removed from
  one long edge (high reflex-vertex density)
- **Irregular** — convex hull with a random subset of vertices displaced
  toward the centroid
- **Field** — smooth, low-frequency radial perturbation of a circle
  (agricultural field/lake boundary)
- **Coastal** — union of several overlapping smooth blobs, producing a
  single connected landmass with bays and peninsulas

Every polygon was validated for geometric simplicity and, critically, for
**shrink-safety**: each candidate was eroded (negative morphological
buffer) at the robot radii actually used in this benchmark, and retained
only if the erosion left a single, valid, simply connected polygon. This
prevents the inward robot-footprint offset from fragmenting the polygon
into multiple disjoint components during path planning. The corpus spans
4–65 vertices, 0–24 reflex vertices, area 17.3–2696.8 units², and aspect
ratio 0.45–5.87. Full generation methodology and parameters are documented
in `generate_benchmark_polygons.py`.

## Experimental configurations

All four `Situation_*.xlsx` files report results for the **same 943
polygons, in the same order**, allowing paired statistical comparison.

| Configuration | Robot width | Rotation |
|---|---|---|
| Situation 1 | Fixed, $w = 0.4$ | OFF |
| Situation 2 | Fixed, $w = 0.4$ | ON |
| Situation 3 | Proportional, robot area = 0.20% of polygon area | OFF |
| Situation 4 | Proportional, robot area = 0.20% of polygon area | ON |

Turn penalty $\alpha = 2$ and overlap ratio $\lambda = 20\%$ were held
constant across all four configurations.

## File format

Each polygon `.txt` file lists vertices in order, one per line:

(x1, y1)
(x2, y2)
...

## Columns in `Situation_*.xlsx`

| Column | Description |
|---|---|
| Idx | Row index |
| Polygon Name | Source `.txt` filename (without extension) |
| Vertices (m) | Vertex count |
| Reflex Vertices (n) | Reflex vertex count |
| Complexity (m+n) | m + n |
| Area | Polygon area |
| Robot Width | Robot diameter used for this run |
| Aspect Ratio | Bounding-box width ÷ height |
| Start Point | Robot start coordinate |
| Coverage (%) | Percentage of polygon area covered |
| Redundancy (%) | Percentage of swept area that was re-covered |
| Path Length | Total path length |
| Coverage/Unit (Area/Length) | Coverage efficiency |
| Turns | Total turn count |
| Rotation Angle (deg) | Optimal rotation angle applied (0 if rotation OFF or not beneficial) |
| Norm. Path Length (PL/sqrt(A)) | Path length normalised by √Area |
| Norm. Turns (T/sqrt(A)) | Turn count normalised by √Area |
| Turn Density (T/PL) | Turns per unit path length |
| Time(sec) | Algorithm runtime |

## Licence

This dataset is released under the [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) licence. If you use this benchmark, please cite the paper above.

## Contact

Questions about this dataset can be directed to [samrat2002dutta@gmail.com / samratdutta2002].
