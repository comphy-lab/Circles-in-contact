# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Initial condition generator for drop/bubble coalescence simulations using Basilisk CFD. Creates smooth interface geometries for unequal-sized drops/bubbles with a fillet transition at the contact point.

## Build Commands

Compile the Basilisk C code:
```bash
qcc -O2 -Wall InitialCondition.c -o InitialCondition -lm
```

Run the compiled executable (L0 is domain size):
```bash
./InitialCondition 8
```

## Architecture

### Workflow
1. **Python (Jupyter)**: `InitialCondition.ipynb` generates interface coordinates using `GetCircles(delta, Rr)`
2. **Data file**: Coordinates saved to `f_Testing.dat` (space-separated x,y pairs)
3. **Basilisk C**: `InitialCondition.c` reads the data file and creates volume fraction field using the `distance.h` function

### Key Parameters
- `delta`: Neck radius / minimum length scale (separation = 2*delta)
- `Rr`: Radius ratio (large drop / small drop)

### Geometry Construction
Three connected circular arcs form the interface:
1. **Circle 1** (orange): Main drop, radius 1, centered at (-(1+delta), 0)
2. **Fillet circle** (red): Smooth transition at contact, computed radius Rf
3. **Circle 2** (blue): Second drop, radius Rr, positioned for tangent continuity

See `Initial condition unequal sized drops_bubbles.pdf` for schematic.

### Output Files
- `f_Testing.dat`: Interface coordinates for Basilisk input
- `f_facets.dat`: Facets output from Basilisk
- `Init_f`: Basilisk dump file with volume fraction field
- `DataFiles_*` / `ImageFiles_*`: Batch outputs for different delta values

## Python Dependencies

```python
numpy, matplotlib, pandas
```

Requires LaTeX for plot labels (`plt.rc('text', usetex=True)`).

## Basilisk Headers

```c
#include "axi.h"           // Axisymmetric coordinates
#include "fractions.h"     // Volume fraction computation
#include "distance.h"      // Signed distance function
#include "navier-stokes/centered.h"
```
