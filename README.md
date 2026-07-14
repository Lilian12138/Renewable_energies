# Renewable_energies

A collection of processing scripts for China wind/solar installable-capacity potential
assessment, gridded allocation of province-level planning targets, and economic /
emission-reduction benefit accounting.

> All scripts in this repository were restructured from the original Jupyter Notebooks /
> earlier scripts. **The calculation logic is identical to the original code** — only
> structural reorganization was done (splitting into functions, centralizing path/parameter
> configuration, adding comments); no calculation results were changed.
> The original files are kept in [`original_notebooks/`](original_notebooks) for traceability.

## Pipeline overview

```
                         ┌────────────────────────────┐
                         │  mapping_solar_panel_qa.py  │  Solar-panel polygon QA basemaps
                         └────────────────────────────┘
                                (standalone QA tool, not part of the downstream pipeline)

  allocate_wind_capacity.py  ──┐
                                ├──►  economy_carbon_analysis.py  ──►  Province-level cost/emission/profit tables
  allocate_solar_capacity.py ──┘

  plan_install_allocation.py   (a separate grid-allocation implementation, see below;
                                 not currently wired into the economics script)
```

## Script descriptions

### 1. `mapping_solar_panel_qa.py`
For selected (`Select == 1`) PV (photovoltaic) polygons, generates a per-province QA image
with a Tianditu satellite basemap overlay, used to manually spot-check whether polygon
boundaries look reasonable. Provinces with fewer than 2 polygons get only an overview map;
provinces with 2 or more additionally get 2 randomly sampled polygons rendered as zoomed-in
insets, combined into one image.

By default it runs for a single province set via `TARGET_PROVINCE` (matching the original
code's behavior of manual, one-province-at-a-time review). To batch-process every province
that still lacks a QA image, see the commented-out block at the end of the script, which
loops over the return value of `list_pending_provinces()`.

Dependencies: `cartopy` (downloads Tianditu WMTS tiles, requires network access), `geopandas`,
`matplotlib`.

### 2. `allocate_wind_capacity.py` / `allocate_solar_capacity.py`
Allocates each province's planned wind / solar capacity targets (2025/2030/2035/2040/2050/
2060, in units of 10,000 kW / "wan-kW") to the 10km grid, year by year and independently,
under three development scenarios: **Low / Median / High**.

**Allocation algorithm (`plan_install`)**:
1. Sort all grid cells for the target province in descending priority order:
   - Wind: `[already installed in 2022 (s2022), score (Score), potential cap (KW2)]`
   - Solar: `[already installed in 2022 (s2022), solar score (score_sola), 2022 generation
     hours (hour2022), potential cap (KW2)]`
2. That year's target (converted from wan-kW to kW) has the 2022 baseline capacity
   subtracted from it to get the amount still to allocate.
3. Fill grid cells in sorted order, capping each cell at its potential ceiling `KW2`;
   stop once the year's target is met or grid potential is exhausted.

**Important note**: this version computes each year's incremental target **independently**;
there is no cascading relationship where one year's result is carried forward as the floor
for the next — this differs from the newer algorithm in `plan_install_allocation.py`
(see the comparison below).

Output: a grid shapefile per scenario, plus a province-level summary check table (Excel).

> The solar script preserves two historical quirks from the original code (neither was
> modified, both are called out in comments):
> - The planning-table file path used at allocation time differs from the one used at QA
>   time (the two files are presumably meant to hold the same content, but live in
>   different locations).
> - The "low" scenario's QA summary actually reads `Solar_high_20240226.shp` (from a
>   path with an extra space) — i.e. the table titled "low" actually contains a summary
>   of the "high" scenario's data.

### 3. `plan_install_allocation.py`
A separate grid-allocation implementation that uses ArcPy to read score feature classes
and grid current/potential data from an ArcGIS Pro project (`.gdb`), rather than shapefiles.
Main differences from `allocate_wind_capacity.py` / `allocate_solar_capacity.py`:

| | legacy (`allocate_*_capacity.py`) | this script |
|---|---|---|
| Scenarios | Low / Median / High (3 development scenarios) | "province" plan / "low_carbon" plan (2 scenarios) |
| Cross-year relationship | Each year's target computed independently | **Cascading**: each year's allocation uses the previous year's result as a floor and only allocates the new increment |
| Data source | Shapefile | ArcGIS `.gdb` feature classes (requires an ArcPy / ArcGIS Pro environment) |
| Wind handling | No distinct score field for onshore vs. offshore | Splits onshore/offshore wind by `Shengcode == 100`, using the matching score field for each |

Requires ArcPy (i.e. must run inside a Python environment with ArcGIS Pro installed).
Outputs a grid-level shapefile (with `kw2025` ~ `kw2060` columns for both scenarios) and a
province-level summary CSV.

**This script's output is not currently wired into `economy_carbon_analysis.py`**
(which instead reads the "Low/Median/High" scenario results produced by
`allocate_wind_capacity.py` / `allocate_solar_capacity.py`).

### 4. `economy_carbon_analysis.py`
Using the wind / solar grid capacity allocation results (`KW2022` ~ `KW2060`), combined
with each province's unit installed cost, grid emission factor, and carbon price
parameters, computes per grid cell:

- **Installed Cost** = grid capacity × province's unit installed cost (price per watt)
- **Emission Reduction** = grid capacity × province's generation hours × grid emission factor
- **Profit** = Emission Reduction / 1000 × carbon price − Installed Cost

If a province is missing its unit installed cost / grid emission factor parameter, the
value is filled in using the average for its region.

The four computation runs pair installed-capacity scenarios with cost/carbon-price
scenarios as follows (carried over unchanged from the original code):

| Installed-capacity grid | Cost/carbon-price scenario |
|---|---|
| Wind - High-speed scenario (`wind_high`) | High-speed |
| Wind - Policy scenario (`wind_median`) | Policy |
| Solar - Policy scenario (`Solar_median`) | Policy |
| Solar - High-speed scenario (`Solar_high`) | High-speed |

Output:
- A grid-level result shapefile per year/scenario (with `C`/`E`/`P` columns);
- Parameter-check tables `feng.xlsx` / `guang.xlsx` (per-province wind/solar parameter
  averages, for checking whether the fill-in looks reasonable);
- Province-level summary workbook `province_results_*.xlsx` (one sheet per result shp);
- Sensitivity tables `province_solar_median_kw2.xlsx` / `province_solar_high_kw2.xlsx`
  (per-province, per-year grid-cell counts and totals).

## Directory structure

```
processing/code_for_github/
├── README.md
├── mapping_solar_panel_qa.py       # Solar-panel polygon QA basemaps
├── allocate_wind_capacity.py       # Wind grid allocation (Low/Median/High scenarios, independent years)
├── allocate_solar_capacity.py      # Solar grid allocation (Low/Median/High scenarios, independent years)
├── plan_install_allocation.py      # Wind/solar grid allocation (province/low-carbon scenarios, cascading years, ArcPy)
├── economy_carbon_analysis.py      # Installed cost / emission-reduction / profit accounting
└── original_notebooks/             # Original notebooks/scripts (unmodified, kept for traceability)
    ├── mappingSolarPanel.py
    ├── ProcessingWind.ipynb
    ├── SolarProcessing.ipynb
    └── EconomyWindSolar.ipynb
```

## Environment

- Python 3.x
- `pandas`, `numpy`, `geopandas`, `openpyxl` (reading/writing `.xlsx`)
- `matplotlib`, `cartopy` (only needed for `mapping_solar_panel_qa.py`)
- `arcpy` (only needed for `plan_install_allocation.py`, requires an ArcGIS Pro environment)

## Usage notes

- Every script's input/output paths are centralized at the top of the file (hardcoded local
  paths that were scattered throughout the original code have been consolidated into a
  "Paths and parameters" section at the start of each file); edit these to match your actual
  data locations before running.
- Each script can be run independently via `python <script>.py`. The recommended order is
  `allocate_wind_capacity.py` / `allocate_solar_capacity.py` → `economy_carbon_analysis.py`,
  so that the intermediate result shapefiles the downstream script depends on already exist.
- This reorganization pass **did not actually run** any of the scripts to verify results —
  it is a static restructuring of the code only. If any paths have since changed, please
  verify them before running.
