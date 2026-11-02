# Manuscript figure scripts

Reproducible numerical figures for the Perspective on optimal control in soil and plant models. These research scripts are separate from the installable `src/eco_opt` library. No LaTeX manuscript source is stored here.

## Run

From the repository root, install the library and plotting dependencies:

```bash
python -m pip install -e ".[examples]"
python scripts/manuscript/fig_box1_microbial.py
python scripts/manuscript/fig_box2_lignin.py
python scripts/manuscript/fig_box3_root.py
```

The scripts solve the published model formulations through the public `eco_opt` API. The microbial figure compares the numerical solution against the published analytical solution. They require an eco_opt version providing costate estimates (`sol.costate`).

PDFs and PNG previews are written to `scripts/manuscript/figures/` by default. To write directly to a local manuscript's figures directory, set an explicit output location, for example:

```bash
ECO_OPT_FIGURE_DIR=/path/to/manuscript/figures python scripts/manuscript/fig_box1_microbial.py
```

The lignin and plant–soil scripts optionally accept `--plot-data /path/to/cache.npz`. An existing cache is used for plotting; if absent, the script solves the model and saves the trajectories there. Keep generated figures and numerical caches out of Git.

## Contents

- `fig_box1_microbial.py` and `models/microbial.py`: microbial uptake, numerical/analytical comparison and shadow-price interpretation (Manzoni et al., 2023).
- `fig_box2_lignin.py` and `models/lignin.py`: optimal ligninolytic activity (Chakrawal et al., 2024).
- `fig_box3_root.py` and `models/plant_soil.py`: plant–soil carbon and nitrogen allocation (Chakrawal et al., 2026).
- `style.py`: shared plotting, schematic and annotation-overlap helpers.

Publication references and model provenance are retained in the model docstrings. Full solves can take substantially longer than plotting cached trajectories.
