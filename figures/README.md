# Figures

Generate the six figures from the repository root:

```text
python figures/make_figure1_publication.py
python figures/make_figures2_6_publication.py
```

Requires Python 3 with numpy, pandas, matplotlib, and Pillow. Figure 1 is a
schematic; Figures 2–6 read result tables in `analysis/results`.

Each figure has a vector PDF, SVG, 300 dpi PNG, and 600 dpi LZW TIFF in this
directory. Plotting tables are saved under `source_data/`.
See `figure_file_map.md` for the result sources and figure subjects.
