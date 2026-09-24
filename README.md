# slCBM — smooth Lagrangian Crack-Band Model

This repository contains numerical models and implementation files associated with the smooth Lagrangian Crack-Band Model (slCBM) and the paper:

> Houlin Xu, Anh Tay Nguyen, and Zdeněk P. Bažant,  
> **“Sprain energy consequences for damage localization and fracture mechanics.”**  
> *Proceedings of the National Academy of Sciences* 121(40), e2410668121 (2024).  
> DOI: https://doi.org/10.1073/pnas.2410668121

The repository provides Abaqus user-element implementations, input data, post-processing scripts, and COMSOL models used to reproduce numerical examples corresponding to Figs. 2–4 of the paper.

## Background

The slCBM formulation is designed to prevent spurious localization of fracture damage into an arbitrarily narrow zone. Its localization limiter is based on the second gradient of the displacement field, called **sprain**. In the finite-element implementation, the displacement field and its gradient are represented as independent fields and constrained using Lagrange multipliers.

For the concrete simulations in this repository, the slCBM formulation is coupled with the **microplane model M7**.

The repository also contains comparisons between the full slCBM formulation and a strain-gradient formulation in which material-rotation gradients are suppressed.

## Repository structure

```text
slCBM/
├── Fig 2/
│   ├── small C/
│   ├── middle C/
│   └── large C/
│
├── Fig. 3/
│   ├── small/
│   │   ├── 0/
│   │   ├── 5/
│   │   └── 10/
│   ├── middle/
│   │   ├── 0/
│   │   ├── 5/
│   │   └── 10/
│   └── large/
│       ├── 0/
│       ├── 5/
│       └── 10/
│
└── Fig 4/
    ├── Fig 4B/
    ├── Fig 4C/
    └── Fig 4DE/
```

### Fig. 2

The `Fig 2` directory contains Abaqus models used to study the influence of the sprain-damage threshold on localization near a notch.

The three cases correspond to different values used in the fifth and sixth UEL properties:

| Case | Property 5 | Property 6 |
|---|---:|---:|
| `small C` | `1e-6` | `2e-6` |
| `middle C` | `1e-5` | `2e-5` |
| `large C` | `1.5e-3` | `3e-3` |

In `SprainE.for`, these values are read as `gamma0` and `gammac`. The first value corresponds to the threshold cases discussed in Fig. 2 of the paper.

Representative files include:

- `sample.inp` — Abaqus input file.
- `UEL_new1.for` — 2D Abaqus user-element implementation.
- `SprainE.for` — sprain/spress contribution, damage evolution, and corresponding tangent terms.
- `M7fMATERIAL.for` — microplane model M7 constitutive implementation.
- `geninp.py` — helper script that generates element connectivity files from `ori.txt`.
- `nodes.txt` — nodal data.
- `ori.txt` — original connectivity/data source used by `geninp.py`.
- `real.txt` — connectivity for the user elements.
- `fake.txt` — duplicated conventional elements used for visualization/output.
- `outputRF.PY` — Abaqus/CAE post-processing script for extracting reaction-force data.

### Fig. 3

The `Fig. 3` directory contains the 3D **gap-test** simulations for three geometrically scaled specimen sizes:

- `small`
- `middle`
- `large`

Within each specimen-size directory, the `0`, `5`, and `10` folders correspond to different imposed surface-pressure levels in the first analysis step:

```text
0  -> Surf-1, P, 0
5  -> Surf-1, P, 5
10 -> Surf-1, P, 10
```

These cases are used to study the effect of crack-parallel compression on fracture behavior.

Representative files include:

- `sample.inp` — 3D Abaqus input file.
- `UEL_new2.for` — 3D slCBM Abaqus user element.
- `SprainE.for` — sprain/spress formulation and damage evolution.
- `M7fMATERIAL.for` — microplane model M7.
- `elastic_pads.txt` — elastic/plastic pad elements used in the gap-test model.
- `nodes.txt`, `ori.txt`, `real.txt`, and `fake.txt` — model geometry/connectivity files.
- size-specific Abaqus input files such as `Small_elastic.inp`, `real_middle.inp`, or `Middle_dense.inp`.

### Fig. 4

The `Fig 4` directory contains COMSOL Multiphysics models used to compare the full slCBM formulation with a strain-gradient formulation.

Each case includes two `.mph` files:

- `shear weak dynamic.mph`
- `shear weak dynamic strain gradient.mph`

These models are associated with the Mode-II/shear comparisons discussed in Fig. 4 of the paper.

## Requirements

### Abaqus models

The Abaqus input files in this repository were generated with **Abaqus/CAE 2020**, as recorded in the `.inp` file headers.

To compile and run the user-element implementation, you need:

- Abaqus/Standard.
- A Fortran compiler supported by your Abaqus installation.
- Python for the optional connectivity-generation script.
- Sufficient memory for the 3D Fig. 3 models.

Abaqus/compiler compatibility depends on the operating system and Abaqus release. Newer Abaqus versions may work, but they have not been verified in this repository.

### COMSOL models

For Fig. 4, use COMSOL Multiphysics to open the supplied `.mph` files.

## Running an Abaqus case

Clone the repository:

```bash
git clone https://github.com/Houlin2018/slCBM.git
cd slCBM
```

Then enter one case directory. For example:

```bash
cd "Fig 2/small C"
```

The supplied `sample.inp` files use relative paths such as `.\nodes.txt`, so the Abaqus job should be launched from the corresponding case directory.

### 1. Optional: regenerate connectivity files

The repository already contains `real.txt` and `fake.txt`. Regeneration is only needed if `ori.txt` has been modified.

```bash
python geninp.py
```

`geninp.py` reads `ori.txt` and writes `real.txt` and `fake.txt`.

### 2. Prepare the Fortran user subroutines

For Fig. 2, the primary sources are:

```text
UEL_new1.for
SprainE.for
M7fMATERIAL.for
```

For Fig. 3, use:

```text
UEL_new2.for
SprainE.for
M7fMATERIAL.for
```

One convenient approach is to concatenate the required source files into one Abaqus user-subroutine file.

On Linux/macOS:

```bash
cat UEL_new1.for SprainE.for M7fMATERIAL.for > user_subroutines.for
```

For a Fig. 3 case, replace `UEL_new1.for` with `UEL_new2.for`.

On Windows Command Prompt:

```bat
copy /b UEL_new1.for+SprainE.for+M7fMATERIAL.for user_subroutines.for
```

### 3. Run Abaqus

Example:

```bash
abaqus job=sample input=sample.inp user=user_subroutines.for interactive
```

Exact command-line syntax can vary with the local Abaqus installation.

## User-element formulation

The Abaqus input files define custom user elements with additional nodal degrees of freedom beyond the conventional displacement field.

For the 2D Fig. 2 implementation:

```text
4-node user element
2 spatial dimensions
40 element DOFs
4 integration points
```

For the 3D Fig. 3 implementation:

```text
8-node user element
3 spatial dimensions
88 element DOFs
8 integration points
```

The implementation contains three coupled fields:

1. displacement,
2. displacement-gradient variables,
3. Lagrange-multiplier variables.

The sprain contribution is evaluated in `SprainE.for`. The code forms the sprain/spress contribution to the residual and tangent stiffness and tracks the associated sprain-energy term.

## Material model

`M7fMATERIAL.for` contains an implicit implementation of the **M7 microplane model for concrete**.

The slCBM user element calls this constitutive routine to obtain the stress response and material tangent while the sprain formulation provides the higher-order localization control.

## Post-processing

For Fig. 2, `outputRF.PY` is an Abaqus/CAE script for extracting reaction forces from the output database.

The current script contains a hard-coded local `.odb` path from the original simulation environment. Before using it, replace that path with the location of your own output database.

The script can be run from Abaqus/CAE or adapted for Abaqus Python.

## Notes on reproducing the paper figures

- **Fig. 2:** compare the `small C`, `middle C`, and `large C` cases to examine the influence of the sprain threshold on strain localization.
- **Fig. 3:** run the three specimen sizes and the `0`, `5`, and `10` pressure cases to reproduce the numerical gap-test trends.
- **Fig. 4:** open the paired COMSOL models to compare slCBM with the strain-gradient formulation.

The paper reports that the slCBM and strain-gradient formulations are similar for symmetric Mode-I fracture, while important differences appear for shear and mixed-mode fracture because the full slCBM includes material-rotation gradients.

## Citation

If you use this repository in academic work, please cite:

```bibtex
@article{Xu2024Sprain,
  author  = {Xu, Houlin and Nguyen, Anh Tay and Bažant, Zdeněk P.},
  title   = {Sprain energy consequences for damage localization and fracture mechanics},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {121},
  number  = {40},
  pages   = {e2410668121},
  year    = {2024},
  doi     = {10.1073/pnas.2410668121}
}
```

The repository itself is cited in the paper as:

> H. Xu, *Smooth Crack Band Model*, GitHub, deposited 6 September 2024.

## Repository

https://github.com/Houlin2018/slCBM

## License

No software license file is currently included in this repository. Until a license is added, reuse and redistribution of the repository contents are not automatically granted by an open-source license.

The publication has its own publication license; that license should not be assumed to apply automatically to the source files in this repository.
