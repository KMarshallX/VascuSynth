# VascuSynth

Generates synthetic 3D vascular trees together with their ground truth: a volumetric image
of the vessels and a graph of the tree (node positions, parent→child edges, radius and flow
per segment).

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
VascuSynth <paramList> <nameList> <imageVoxelWidth> [noiseFile ...]
```

| Argument | Meaning |
|---|---|
| `paramList` | Text file listing parameter files, one per line. One tree is built per line. |
| `nameList` | Text file listing output folder names, one per line, matched to `paramList` by line order. Must have at least as many lines as `paramList`. Use plain folder names; parent folders are not created. |
| `imageVoxelWidth` | Voxel size of the output image, in the same physical unit as `VOXEL_WIDTH` in the parameter file. Equal to `VOXEL_WIDTH` gives one image voxel per map voxel; twice `VOXEL_WIDTH` halves the resolution. |
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

For each name in `nameList`:

```
<name>/
├── tree_structure.xml          # GXL graph of the tree
├── original_image/imageNNN.png # clean volume, one 8-bit PNG slice per z index
└── noise_image_K/imageNNN.png  # only if noise files were given
```

**Graph (`tree_structure.xml`).** Nodes carry `nodeType` (root, bifurcation, terminal) and
`position`. Edges go parent → child and carry `radius` and `flow`. Vessels are straight
cylinders between connected nodes, so the polyline through the nodes is the exact
centreline of the rendered vessels. A tree with `NUM_NODES` = $N$ has $2N$ nodes (root,
$N$ terminals, $N-1$ bifurcations) and $2N-1$ edges. Attribute names have a leading space
(`name=" radius"`), so strict GXL readers may reject the file; a plain XML parser works.

**Units and coordinates.**

- `position` is in oxygenation-map voxel coordinates. With $w_\text{map}$ = `VOXEL_WIDTH` and
  $w_\text{img}$ = `imageVoxelWidth`, the image voxel index is
  $\text{pos} \cdot w_\text{map} / w_\text{img}$.
- `radius` is in physical units (the unit of `VOXEL_WIDTH`), so the radius in image voxels is
  $r / w_\text{img}$.
- The image has $\lfloor D \cdot w_\text{map} / w_\text{img} \rfloor + 1$ voxels along an axis
  whose map size is $D$ (for example, 101 for a 100-voxel map at equal widths).
- Image intensity is a partial-volume value: each voxel is sampled at 8 sub-points and the
  value is about $255 \times$ (fraction of sub-points inside a vessel). Threshold it for a
  binary mask.
- The axis order between the PNG slices and the graph positions has not been checked here;
  check it on one tree before overlaying.

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
| `VOXEL_WIDTH` | `0.04` | Physical size of one map voxel, $w_\text{map}$. It sets the unit of all lengths and radii. |
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
