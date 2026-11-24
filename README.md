# VascuSynth

Generates synthetic 3D vascular trees together with their ground truth: a volumetric image
of the vessels, the exact centreline graph (node positions, parent→child edges, radius and
flow per segment) as GraphML, and a voxel skeleton as NIfTI. The GraphML and NIfTI files
follow the conventions of SkelHub's evaluator, so they can be used directly as references
for skeletonisation methods.

For more information about VascuSynth, please visit: https://vascusynth.cs.sfu.ca/

There is no GUI. VascuSynth is a command-line program driven by small text files: a
parameter file per tree, a supply map, an oxygenation (demand) map, and optional noise files.

## Build with conda (no root, no ITK source build)

Tested with ITK 5.4 and CMake 4.5 from conda-forge on Rocky Linux 9.

```bash
conda create -p <env_path> -c conda-forge --override-channels libitk-devel cmake make cxx-compiler
conda activate <env_path>
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
cmake --build build -j4
cmake --install build --prefix $CONDA_PREFIX   # puts VascuSynth on the env's PATH
```

Build out of source (`-B build`) so the committed `Makefile` is not overwritten.

After changing the code, rerun the `cmake` lines above with that env active.

## Quick start

Run from an empty output folder. Paths inside the parameter file (the supply and
oxygenation maps) are resolved relative to the folder you run from, and each image name
becomes a sub-folder there.

```bash
conda activate vascusynth
REPO=/xx/yy/VascuSynth

mkdir -p my_trees && cd my_trees
cp $REPO/{p1.txt,testS.txt,testOx.txt,testNoise.txt} .
echo p1.txt  > params.txt      # one parameter file per line
echo tree001 > names.txt       # one output folder name per line (same order)

VascuSynth params.txt names.txt 0.04 testNoise.txt
```

This builds one 200-terminal tree on a 100³ map in about 6 s, plus one noisy copy of the image.

### Command line

```
VascuSynth [--original-export] <paramList> <nameList> <imageVoxelWidth> [noiseFile ...]
```

| Argument | Meaning |
|---|---|
| `--original-export` | Optional, accepted anywhere. Write only the original VascuSynth outputs (PNG slices and the GXL file `tree_structure.xml`) instead of the default ground-truth export. See [Outputs](#outputs). |
| `paramList` | Text file listing parameter files, one per line. One tree is built per line. |
| `nameList` | Text file listing output folder names, one per line, matched to `paramList` by line order. Must have at least as many lines as `paramList`. Use plain folder names; parent folders are not created. |
| `imageVoxelWidth` | Voxel size of the output image in mm (the same unit as `VOXEL_WIDTH`). Equal to `VOXEL_WIDTH` gives one image voxel per map voxel; twice `VOXEL_WIDTH` halves the resolution, half doubles it. |
| `noiseFile ...` | Optional. Each noise file produces one extra noisy copy of every image (`noise_image_0`, `noise_image_1`, …). |

Blank lines are ignored, and Windows (CRLF) line endings are accepted. The exit status is
1 on any error, so batch scripts can detect failures.

### Generating many trees

Vary the seed (and any other parameter) by writing one parameter file per tree:

```bash
for i in $(seq -w 1 20); do
  sed "s/^RANDOM_SEED: .*/RANDOM_SEED: $i/" p1.txt > p_$i.txt
  echo p_$i.txt >> params.txt
  echo tree_$i  >> names.txt
done
VascuSynth params.txt names.txt 0.04
```

Each run is single-threaded, so for large batches split the lists and run several processes
(or a Slurm array) in parallel.

## Outputs

For each name in `nameList`, the default export writes:

```
<name>/
├── centerline_graph.graphml     # exact centreline graph (GraphML, undirected)
├── skeleton.nii.gz              # voxel skeleton, uint8 0/1
├── vessel_mask.nii.gz           # binary vessel mask, uint8 0/1
├── vessel_image.nii.gz          # partial-volume vessel image, uint8 0-255
├── original_image/imageNNN.png  # same image as 8-bit PNG slices, one per z index
├── noise_image_K.nii.gz         # noisy copy K (only if noise files were given)
└── noise_image_K/imageNNN.png   # same noisy copy as PNG slices
```

With `--original-export` you get exactly the original VascuSynth outputs instead:
`original_image/`, `noise_image_K/` and `tree_structure.xml` (GXL), with no GraphML or NIfTI.
Both modes draw the same tree and identical PNG slices for the same parameter file.

### Coordinates and units

Everything in the default export uses one frame:

- **Units are mm.** `VOXEL_WIDTH` (the map voxel size, $w_\text{map}$) and `imageVoxelWidth`
  ($w_\text{img}$) are read as mm.
- **Voxel index.** `voxel_pos` and `centerline_voxel_points` in the GraphML are $(i, j, k)$
  indices into the NIfTI arrays, as loaded by `nibabel` (`data[i, j, k]`). Voxel centres sit at
  integer indices. A tree position $p$ in map voxels maps to the voxel index
  $p \cdot w_\text{map} / w_\text{img}$.
- **World coordinates.** Every NIfTI has the affine $\mathrm{diag}(w_\text{img}, w_\text{img}, w_\text{img}, 1)$
  with origin 0 and unit mm. The GraphML `X`, `Y`, `Z` and `centerline_world_points` equal the
  affine applied to the voxel index, i.e. index $\times\, w_\text{img}$.
- **PNG slices.** Slice `imageNNN.png` is $k = \text{NNN}$. Its rows are $j$ and its columns
  are $i$, so `png.T == data[:, :, k]`.
- The image has $\lfloor D \cdot w_\text{map} / w_\text{img} \rfloor + 1$ voxels along an axis
  whose map size is $D$ (for example, 101 for a 100-voxel map at equal widths).

### `centerline_graph.graphml`

An undirected GraphML graph that matches SkelHub's GraphML evaluator. It can be passed straight
to `skelhub evaluate --ref centerline_graph.graphml --coordinate-image skeleton.nii.gz ...`.

| Level | Attribute | Meaning |
|---|---|---|
| graph | `generator`, `param_file` | `VascuSynth` and the parameter file used. |
| graph | `spatial_unit`, `map_voxel_width_mm`, `image_voxel_width_mm`, `image_shape` | `mm`, $w_\text{map}$, $w_\text{img}$ and the NIfTI shape. |
| graph | `num_terminals` | Number of terminals $N$. |
| node | `id` | `n<index>`, the same ids as the GXL file. |
| node | `voxel_pos` | JSON `[i, j, k]` voxel index (floats, exact). |
| node | `X`, `Y`, `Z` | World position in mm. |
| node | `map_pos` | Original VascuSynth position in map voxels. |
| node | `kind` | `root`, `bifurcation` or `terminal`. |
| edge | `parent`, `child` | Node ids. Edges are also written `source` = parent, `target` = child. |
| edge | `radius_mm`, `radius_voxels` | Vessel radius, constant along the segment. |
| edge | `flow` | Flow through the segment, in the unit of `PERF_FLOW`. |
| edge | `length_mm` | Segment length. |
| edge | `centerline_voxel_points` | JSON polyline `[[i,j,k],[i,j,k]]`: the exact straight centreline from parent to child. |
| edge | `centerline_world_points` | The same polyline in mm. |
| edge | `centerline_voxels`, `num_centerline_voxels` | This segment's voxel path through the cleaned `skeleton.nii.gz`: the shortest 26-connected path from the parent's node voxel to the child's. |

A tree with `NUM_NODES` = $N$ has $2N$ nodes (root, $N$ terminals, $N-1$ bifurcations) and
$2N-1$ edges. Vessels are straight cylinders between connected nodes, so the polylines are the
exact centreline of the rendered vessels.

### NIfTI volumes

- `vessel_image.nii.gz`: each voxel is sampled at 8 sub-points; the value is
  $32 \times (\text{sub-points inside a vessel}) - 1$, clipped to 0–255. Same data as the
  clean PNG slices.
- `vessel_mask.nii.gz`: 1 where any of the 8 sub-points is inside a vessel
  (`vessel_image > 0`), plus every skeleton voxel. The mask therefore always contains 100% of
  the skeleton, which skeletonisation methods and SkelHub's `--foreground` option expect.
- `skeleton.nii.gz`: a thin, loop-free voxel tree (see below).
- `noise_image_K.nii.gz`: noisy copy $K$, same data as its PNG slices.

### How the skeleton is made

1. Each segment is drawn as a 26-connected line of voxels between its rounded end points, using
   the same rule as SkelHub's GraphML writer.
2. Where lines meet or pass close to each other, the drawn voxels form clumps and small false
   loops. They are cleaned with 3D digital topology (foreground 26-connected, background
   6-connected, the convention SkelHub's evaluator uses):
   - **Thinning.** Voxels whose removal cannot change the topology (simple points) are removed,
     farthest from the true centreline first.
   - **Cutting loops.** A loop is broken by removing a voxel whose neighbours split locally but
     stay joined through the rest of the skeleton. Spots where two unrelated branches touch are
     cut first, then the voxel farthest from the centreline.
   - **Filling tiny rings.** A ring of a few voxels that each carry a branch cannot be cut
     without detaching something, so the background voxel it encloses is filled instead, and
     thinning then pares the clump back.
   - Root and terminal voxels are always kept. Bifurcation voxels are kept unless a loop can
     only be broken through one; the voxel junction then sits one voxel away, and the edges'
     `centerline_voxels` start at that neighbouring voxel.
3. The run log reports the drawn and final voxel counts and the number of loops removed.

Results on test trees (sample parameters with different seeds, sizes and resolutions):

| Tree | Loops left | Components | Voxel end points / graph end points | Max distance from true centreline |
|---|---|---|---|---|
| 200 terminals, 3 seeds, 0.04 mm | 0 | 1 | 191–193 / 201 | 1.06–1.26 voxels |
| 500 terminals, 0.04 mm | 0 | 1 | 468 / 501 | 1.29 voxels |
| 1000 terminals, 0.04 mm | 0 | 1 | 890 / 1001 | 1.29 voxels |
| 200 terminals, 0.02 mm | 0 | 1 | 198 / 201 | 1.18 voxels |
| 200 terminals, 0.08 mm | 0 | 1 | 171 / 201 | 1.16 voxels |

The missing voxel end points are terminals lying within about a voxel of another skeleton line,
mostly on very short terminal segments (1–4 voxels long). The terminal voxel touches that line,
so it cannot be a one-neighbour end point without breaking connectivity. This is a resolution
limit: a finer `imageVoxelWidth` recovers most of them. The GraphML always has the exact
topology, so use it for topology and end-point counts and the NIfTI skeleton for voxel-level
comparisons. The edges' `centerline_voxels` together cover the skeleton, except occasionally a
single corner voxel at a junction.

### Original GXL file (`--original-export`)

`tree_structure.xml` holds nodes with `nodeType` and `position` (map voxels) and parent → child
edges with `radius` (mm) and `flow`. Attribute names have a leading space (`name=" radius"`),
so strict GXL readers may reject the file; a plain XML parser works.

## Parameter file

One `KEY: value` per line. There must be exactly one space after the colon. Every key is
required except `RANDOM_SEED`. Values from the shipped `p1.txt` are shown.

| Key | `p1.txt` | What it controls |
|---|---|---|
| `RANDOM_SEED` | `1` | Seed for where terminals are placed. The same seed and parameters give an identical tree and clean image. If omitted, a time-based seed is used. |
| `SUPPLY_MAP` | `testS.txt` | Supply map file (see below). |
| `OXYGENATION_MAP` | `testOx.txt` | Oxygenation (demand) map file (see below). It also sets the map size. |
| `PERF_POINT` | `0 50 50` | Position of the root (inlet), in map voxel coordinates `x y z`. |
| `PERF_PRESSURE` | `133000` | Pressure at the inlet. |
| `TERM_PRESSURE` | `83000` | Pressure at every terminal. Only the drop $\Delta P = P_\text{perf} - P_\text{term}$ matters: it sets the root radius (below). |
| `PERF_FLOW` | `8.33` | Total inflow $Q_\text{perf}$. Each terminal receives $Q_\text{perf} / N$. |
| `RHO` | `0.036` | Blood viscosity $\rho$, used in the Poiseuille resistance of each segment, $8\rho\ell/\pi$. |
| `GAMMA` | `3` | Bifurcation exponent $\gamma$ (Murray's law): $r_\text{parent}^\gamma = r_\text{left}^\gamma + r_\text{right}^\gamma$. Lower values make daughters thinner relative to the parent. |
| `LAMBDA` | `2` | Exponent $\lambda$ on radius in the cost $\sum_i \ell_i^{\mu} r_i^{\lambda}$ that is minimised when choosing where a new terminal joins the tree. $\lambda = 2$, $\mu = 1$ minimises total vessel volume. |
| `MU` | `1` | Exponent $\mu$ on segment length in the same cost. |
| `MIN_DISTANCE` | `1` | A new terminal is rejected if it is closer than this to any existing segment (map voxels). Larger values give sparser trees. |
| `NUM_NODES` | `200` | Number of terminals $N$. This is the main knob for tree size and run time. |
| `VOXEL_WIDTH` | `0.04` | Physical size of one map voxel in mm, $w_\text{map}$. It sets the unit of all lengths and radii. (The original VascuSynth samples were in CGS, i.e. cm; this fork reads it as mm.) |
| `CLOSEST_NEIGHBOURS` | `5` | How many nearby segments are tried as connection points for each new terminal. More is slower but finds lower-cost connections. |

Radii are computed after the tree is grown. The root segment gets
$r_\text{root} = \left( Q\,\tilde R / \Delta P \right)^{1/4}$, where $\tilde R$ is the
reduced resistance of the whole tree, and each bifurcation splits radius by $\gamma$ and the
flow and resistance of the two subtrees.

## Oxygenation (demand) map

Sets where terminals are wanted and where vessels may go.

```
X Y Z                 # map size in voxels
x0 y0 z0 x1 y1 z1     # a box, upper bounds exclusive
value                 # demand inside that box
...                   # more box/value pairs; later boxes overwrite earlier ones
```

Each new terminal is drawn at random with probability proportional to the current demand.
Voxels with demand 0 never get terminals, and connections are not allowed to pass through
them, so zero boxes act as obstacles. The shipped `testOx.txt` fills the 100³ volume with 1
and then cuts twelve zero-valued bars through it.

## Supply map

Sets how much demand a terminal removes around itself after it is placed. This spreads
terminals out.

```
X Y Z K               # map size and number of coefficients K
x0 y0 z0 x1 y1 z1     # a box, upper bounds exclusive
a0 a1 ... a(K-1)      # K coefficients for that box
...
```

After a terminal is placed, the demand at every voxel at distance $d$ (map voxels) is
multiplied by

$$
f(d) =
\begin{cases}
1 & d \ge a_{K-1} \\
\max\!\left(0,\; 1 - \dfrac{1}{\sum_{i=0}^{K-1} a_i d^{i}}\right) & \text{otherwise}
\end{cases}
$$

The last coefficient does double duty: it is the cut-off distance and also the highest-order
coefficient. The shipped `testS.txt` uses `0.65 0.34 0.01 7`, so demand is reduced within 7
voxels of each terminal.

## Noise files

One noise operation per line, applied in file order to a copy of the clean image:

```
GAUSSIAN: mean sigma                       # additive Gaussian noise
UNIFORM: low high                          # additive uniform noise
SALTPEPPER: saltValue pSalt pepperValue pPepper
SHADOW: count                              # 'count' dark spots, see below
```

Each shadow picks a random segment and darkens a sphere centred on its midpoint, with radius
equal to the segment length. Intensity is scaled by a random factor in $[0, 1)$ at the
centre, rising linearly back to 1 at the sphere's edge.

The noise is not seeded, so noisy images differ between runs even with the same
`RANDOM_SEED`. The shipped `testNoise.txt` is `GAUSSIAN: 0 20`.

## Run time

Measured with 100³ map, `imageVoxelWidth` = `VOXEL_WIDTH` = 0.04,
no noise:

| `NUM_NODES` | Time | Peak memory |
|---|---|---|
| 200 | ~6 s | — |
| 500 | ~25 s | 135 MB |
| 1000 | ~98 s | 190 MB |

Time grows roughly with $N^2$. A finer `imageVoxelWidth` makes the drawing step slower too,
since every image voxel is tested against every segment.

Terminals are added until `NUM_NODES` have been placed, with no attempt limit. From reading
the code, settings that make placement impossible (for example too many terminals for the
available demand, or a very large `MIN_DISTANCE`) can make the program run forever. Wrap
large or unusual runs in `timeout`.

## Original build instructions (ITK 4.13, Ubuntu)

Follow the steps below meticulously to compile the code base on Ubuntu:

0. `git clone https://github.com/sfu-mial/VascuSynth`
1. `sudo apt-get install libinsighttoolkit4-dev`
2. `git clone -b v4.13.3 https://github.com/InsightSoftwareConsortium/ITK ITK-4.13.3`
3. `cd ITK-4.13.3`
4. `cmake -B build `
5. `cd build && make -j32 && cd ..`
6. `cmake . # to build the vascusynth bianry`
7. `make -j32`


**NOTE**: Dependencies for Windows 10: CMake 3.29.1 + VS Studio 2022 community + MSVC v140 C++ build tools (v14.00) + ITK 5.1.0

To run the built binary file:

`./VascuSynth paramFiles.txt imageNames.txt 0.04`
