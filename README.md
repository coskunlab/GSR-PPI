## GSR-PPI

This repository contains code for the manuscript "Graph-based spatial proximity of super-resolved protein-protein interactions predict cancer drug responses in single cells". Codes are run under the specified Anaconda environment.

To set up environments, run the following command: `conda env create -f environ.yml`

The data can be found at: https://figshare.com/projects/Signaling_Project_PLA/195958 


## Citation

Please cite: Zhang, N. et al. Graph-Based Spatial Proximity of Super-Resolved Protein–Protein Interactions Predicts Cancer Drug Responses in Single Cells. Cel. Mol. Bioeng. https://doi.org/10.1007/s12195-024-00822-1 (2024) doi:10.1007/s12195-024-00822-1.

## Project Overview

This is a **PLA (Proximity Ligation Assay) super-resolution microscopy analysis** repository for analyzing protein-protein interactions (PPIs) in HCC827 lung cancer cells. The project compares PPI detection across different microscopy modalities (widefield, localized widefield, confocal/SRRF) and drug treatment conditions.

## Environment Setup

**Python Version:** 3.12.7

**Environment Management:** This project uses Conda. Create the environment from the provided file:

```bash
conda env create -f environment.yml
conda activate base
```

The environment includes key dependencies:
- **Image Analysis:** napari, scikit-image, opencv-python, nd2, cellpose
- **Data Processing:** pandas, numpy, dask, h5py
- **Deep Learning:** PyTorch (2.6.0+cu118), JAX, torchvision
- **Visualization:** matplotlib, seaborn, plotly
- **ML/Analysis:** scikit-learn, scanpy, anndata, umap-learn
- **Graph Neural Networks:** PyTorch Geometric libraries

## Repository Structure

```
.
├── data/                    # Experimental data (local relative paths)
│   ├── HCC827 P6 plate 003/ # Main dataset with microscopy comparisons
│   ├── HCC827 P7 plate 002/ # Treatment comparison dataset
│   └── Human Lung */        # Human tissue samples
├── notebooks/               # Jupyter analysis notebooks (numbered workflow)
│   ├── 01_count_dots_per_cell_untreated_vs_treated.ipynb
│   ├── 02_plot_GNN_metrics.ipynb
│   └── 03_plot_HCC827_cells_PPIs_across_microscopes.ipynb
├── figures/                 # Output plots and visualizations
├── environment.yml          # Conda environment specification
└── update_paths.py         # Utility script for data path management
```

## Key Concepts and Architecture

### Data Organization

**Hierarchical Structure:**
- **Plate-level:** Experimental plates (e.g., "HCC827 P6 plate 003")
- **Cycle-level:** Imaging cycles with modality info (e.g., "15Jul2024 cycle 1 Nicky widefield")
- **Well-level:** Individual wells with treatment conditions tracked in Excel metadata
- **FOV-level:** Field of view (FOV) images
- **Single-cell:** Pickle files (`.pkl`) containing per-cell quantification data

**Data Formats:**
- **PKL files:** Serialized pandas DataFrames with single-cell measurements
- **Excel metadata:** Well information, drug treatments, microscope parameters
- **GNN saved models:** PyTorch model checkpoints and training logs

### Analysis Workflow (Notebooks)

**01_count_dots_per_cell_untreated_vs_treated.ipynb**
- Loads single-cell PKL files from multiple datasets
- Aggregates PPI dot counts per cell across FOVs
- Compares untreated vs treated conditions
- Statistical testing with Mann-Whitney U test
- Generates comparative barplots for widefield vs confocal

**02_plot_GNN_metrics.ipynb**
- Analyzes Graph Neural Network model performance
- Loads model training metrics and predictions

**03_plot_HCC827_cells_PPIs_across_microscopes.ipynb**
- Cross-microscope comparison (Widefield, Localized/SRW, Confocal/SRRF)
- Reads Excel configuration files for path mapping
- Computes PPI counts and cell confluency metrics
- Parallel processing with joblib for performance

### Data Processing Patterns

**Common Analysis Flow:**
1. Load PKL files with `pd.read_pickle(path)`
2. Filter/validate data (check for empty dataframes)
3. Map well labels to treatment conditions using Excel metadata
4. Group by FOV and CellLabel, aggregate with `.sum()` or custom aggregation
5. Parallel processing with `joblib.Parallel` for multiple files
6. Statistical comparisons and visualization

**Standard Columns:**
- **ID columns:** Y, X, Z, FOV, MaskCytoLabel, MaskNucLabel, CellLabel, Cycle, CellRegion
- **Marker columns:** DAPI, Phalloidin, PPI pairs (e.g., "FGFR1_PIK3R1", "Combo5")
- **Metadata:** Treatment, Drug, WellLabel, Microscope

### Path Management

**CRITICAL:** All data paths in notebooks use **relative paths from the notebook directory**:
```python
Path(r"../data/HCC827 P6 plate 003/...")
```

The `update_paths.py` utility script handles:
- Converting absolute paths to relative paths in Excel metadata
- Copying data from source locations to local `data/` folder
- Preserving nested directory structure
- Run with `--execute` flag to perform file operations

**Excel Path Format:**
- Original: `Y:\coskun-lab\Nicky\47 PLA super resolution\Data\...`
- Updated: `data/HCC827 P6 plate 003/...` (relative to Github root)

### Visualization Standards

**Figure Saving:**
All notebooks use a common pattern with timestamp-based filenames:
```python
def saveFigLabelTime(fig):
    now = datetime.now().strftime('%d%b%Y_%H%M%S')
    fileOut = os.path.join(screenshotSavePath, now + '.png')
    fig.savefig(fileOut, dpi='figure', bbox_inches='tight', pad_inches=0)
    time.sleep(2)  # Prevent filename collisions
```

**Standard Plotting:**
- DPI: 300 for publication quality
- Style: `sns.set(font_scale=2)`, `sns.set_style('whitegrid')`
- Statistical annotations: Use `statannotations` package with Mann-Whitney tests
- Output directory: `../figures`

### Performance Considerations

**Parallel Processing:**
- Use `joblib.Parallel(n_jobs=-1, prefer='threads')` for I/O-bound tasks
- Typical pattern: load and process multiple PKL files concurrently
- Set `verbose=10` for progress monitoring

**Large Dataset Handling:**
- GNN model directories contain 5000+ files each
- Use Dask for larger-than-memory operations
- PKL files are per-well, enabling modular processing

## Running the Analysis

**Typical workflow:**

1. Ensure data is in `data/` with relative paths set correctly
2. Run notebooks in order (01 → 02 → 03)
3. Check `figures/` for timestamped output plots

**Interactive Analysis:**
- Notebooks are designed for Jupyter with `%gui qt5` for napari visualization
- Use `napari` for interactive cell segmentation and dot visualization
- Cellpose is used for automated cell segmentation

## Important Notes

- **Excel metadata files** are critical for mapping well labels to treatments
- **FOV numbering** may not be sequential; use FOV column from dataframes
- **Treatment conditions** vary by experiment; always verify against well info Excel files
- **Microscope parameters** (DimY, DimX) differ between modalities and are stored in Excel metadata
- **PPI naming:** Single pairs use underscore (e.g., "FGFR1_PIK3R1"), multiplex uses "Combo" prefix
