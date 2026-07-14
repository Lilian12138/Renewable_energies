"""
Wind/Solar Capacity Economics and Emission-Reduction Analysis
================================================================

Using the wind and solar grid capacity allocation results (KW2022 ~ KW2060,
produced by allocate_wind_capacity.py / allocate_solar_capacity.py), combined
with each province's unit cost, grid emission factor, and carbon price
parameters, computes per grid cell:

    Installed Cost      Cost     = grid capacity * province's unit installed cost (price per watt)
    Emission Reduction  Emission = grid capacity * province's generation hours * grid emission factor
    Profit               Profit  = Emission / 1000 * carbon price - Cost

If a province is missing its unit installed cost / grid emission factor
parameter, the value is filled in using the average for its region (if the
entire region is also missing the value, that field stays blank for the
province).

The four grid-level computation runs correspond to:
    Wind - High-speed scenario (wind_high)     x  Cost scenario "高速" (High-speed)
    Wind - Policy scenario (wind_median)       x  Cost scenario "政策" (Policy)
    Solar - Policy scenario (Solar_median)     x  Cost scenario "政策" (Policy)
    Solar - High-speed scenario (Solar_high)   x  Cost scenario "高速" (High-speed)
(This pairing between installed-capacity scenario and cost/carbon-price
scenario is carried over unchanged from the original code.)

Ported from: EconomyWindSolar.ipynb (restructured for readability; the
original calculation logic and results are unchanged).

Additional notes (carried over from the original notebook's analysis notes,
for reference only, no effect on the calculation logic): the intent was to
gauge what order of magnitude the ratio of total wind/solar capacity could
reach compared with the Policy- and High-speed-scenario capacities, with an
expected range of 20%-40%; when aggregating, negative values take the larger
of the two, positive values are summed, and existing (baseline) capacity is
kept as-is.
"""

import glob
import os

import geopandas as gpd
import pandas as pd

# ============================================================
# 1. Paths
# ============================================================
PARAM_TABLE = r"C:\Users\ruome\Desktop\New folder\Renewable\ProcessingFiles\table\电网排放因子和成本-碳价参数0110V3.xlsx"

RESULT2_DIR = r"C:\Users\ruome\Desktop\New folder\Renewable\论文编辑\数据输出\结果二"
RESULT3_DIR = r"C:\Users\ruome\Desktop\New folder\Renewable\论文编辑\数据输出\结果三-全生命周期折减"
TABLE_DIR = r"C:\Users\ruome\Desktop\New folder\Renewable\ProcessingFiles\table"
SENSITIVITY_DIR = r"C:\Users\ruome\Desktop\New folder\Renewable\ProcessingFiles\table\敏感性统计表格"

FUTURE_YEARS = ["2025", "2030", "2035", "2040", "2050", "2060"]
ALL_YEARS = ["2022"] + FUTURE_YEARS

# The four grid-level computation runs: (input grid shp, resource type (风=wind/光=solar),
# cost/carbon-price scenario, output shp)
RUNS = [
    {
        "shp_path": os.path.join(RESULT2_DIR, "wind_high_20240320.shp"),
        "sel_str": "风",
        "scenario": "高速",
        "output": os.path.join(RESULT3_DIR, "wind_high_results_20250110gj.shp"),
    },
    {
        "shp_path": os.path.join(RESULT2_DIR, "wind_median_20240320.shp"),
        "sel_str": "风",
        "scenario": "政策",
        "output": os.path.join(RESULT3_DIR, "wind_median_results_20250110gj.shp"),
    },
    {
        "shp_path": os.path.join(RESULT2_DIR, "Solar_median_20240226.shp"),
        "sel_str": "光",
        "scenario": "政策",
        "output": os.path.join(RESULT3_DIR, "Solar_median_results_20250110gj.shp"),
    },
    {
        "shp_path": os.path.join(RESULT2_DIR, "Solar_high_20240226.shp"),
        "sel_str": "光",
        "scenario": "高速",
        "output": os.path.join(RESULT3_DIR, "Solar_high_results_20250110gj.shp"),
    },
]

OUTPUT_COLS = (
    ["NID10", "Shengcode", "NID50", "NID10a", "KW2"]
    + [f"KW{y}" for y in ALL_YEARS]
    + [f"C{y}" for y in ALL_YEARS]
    + [f"E{y}" for y in ALL_YEARS]
    + [f"P{y}" for y in ALL_YEARS]
    + ["geometry"]
)


# ============================================================
# 2. Parameter table loading
# ============================================================
def load_parameter_tables():
    """Load the per-province generation-hours table, the cost/carbon-price
    parameter table, and the province-to-region mapping table."""
    df_hours = pd.read_excel(PARAM_TABLE, sheet_name="2022年分省风光小时数")
    df_hours = df_hours.iloc[4:, -4:]
    df_hours.columns = ["Shengcode", "省份", "风", "光"]
    df_hours.fillna(0, inplace=True)
    df_hours = df_hours.astype({"Shengcode": int, "风": float, "光": float})

    df_parameter_summary = pd.read_excel(PARAM_TABLE, sheet_name="成本参数及碳价参数")
    col = pd.read_excel(PARAM_TABLE, sheet_name="Columns")
    df_parameter_summary.columns = col.iloc[:, 1].to_list()
    df_parameter_summary = df_parameter_summary.iloc[2:]

    df_shengcode = pd.read_excel(PARAM_TABLE, sheet_name="Shengcode")

    return df_hours, df_parameter_summary, df_shengcode


# ============================================================
# 3. Core calculation functions
# ============================================================
def find_str(search_list, olist):
    """Return the first column name that contains every term in search_list."""
    return [col for col in olist if all(term in col for term in search_list)][0]


def _fill_missing_by_region(df_mean, df_mean_by_region, cols, df_shengcode):
    """Fill missing province-level parameter values using the average for their region."""
    df_mean = df_mean.merge(df_shengcode.loc[:, ["Shengcode", "大区"]], how="right", on="Shengcode")
    for col in cols:
        for region in df_mean_by_region.index:
            mean_value = df_mean_by_region.loc[region, col]
            mask = df_mean["大区"] == region
            df_mean.loc[mask, col] = df_mean.loc[mask, col].fillna(mean_value)
    return df_mean


def cal3(df, sel_str, df_shengcode):
    """
    Filter the parameter table by resource type (wind/solar) and compute
    per-province averages (missing values filled with the region average).
    Used to export parameter-check tables (feng.xlsx / guang.xlsx).
    """
    df_value_sel = df.loc[df["风/光"].str.contains(sel_str, case=False, na=False)].copy()
    cols = df_value_sel.columns[12:]
    df_value_sel[cols] = df_value_sel[cols].apply(pd.to_numeric, errors="coerce")

    df_value_sel_mean = df_value_sel.groupby("Shengcode")[cols].mean().reset_index()
    df_value_sel_mean_daqu = df_value_sel.groupby("大区")[cols].mean()

    return _fill_missing_by_region(df_value_sel_mean, df_value_sel_mean_daqu, cols, df_shengcode)


def cal2(df, df_hours, shp_path, sel_str, scenario, df_shengcode):
    """
    For a single grid shapefile, compute the 2022 baseline and each future
    year's (under the given cost/carbon-price scenario) Cost, Emission
    (reduction), and Profit.
    """
    df_value_sel = df.loc[df["风/光"].str.contains(sel_str, case=False, na=False)].copy()
    cols = df_value_sel.columns[12:]
    df_value_sel[cols] = df_value_sel[cols].apply(pd.to_numeric, errors="coerce")

    df_value_sel_mean = df_value_sel.groupby("Shengcode")[cols].mean().reset_index()
    df_value_sel_mean_daqu = df_value_sel.groupby("大区")[cols].mean()
    df_value_sel_mean = _fill_missing_by_region(df_value_sel_mean, df_value_sel_mean_daqu, cols, df_shengcode)

    shp_df = gpd.read_file(shp_path)
    shp_crs = shp_df.crs
    shp_df["Shengcode"] = shp_df["Shengcode"].astype("int")
    shp_df = shp_df.merge(df_value_sel_mean, how="left", on="Shengcode")
    shp_df = shp_df.merge(df_hours.loc[:, ["Shengcode", sel_str]], how="left", on="Shengcode")

    # 2022 baseline
    olist = df.columns
    shp_df["C2022"] = shp_df["KW2022"] * shp_df[find_str(["2022", "单位瓦价格"], olist)] * 1000
    shp_df["E2022"] = shp_df["KW2022"] * shp_df[find_str(["2022", "电力排放因子"], olist)] * shp_df[sel_str]
    shp_df["P2022"] = (shp_df["E2022"] / 1000) * shp_df[find_str(["2022", "碳价"], olist)] - shp_df["C2022"]

    # Future years (under the given scenario)
    for year in FUTURE_YEARS:
        shp_df[f"C{year}"] = shp_df[f"KW{year}"] * shp_df[find_str([scenario, year, "单位瓦价格"], olist)] * 1000
        shp_df[f"E{year}"] = shp_df[f"KW{year}"] * shp_df[sel_str] * shp_df[find_str([scenario, year, "电网因子"], olist)]
        shp_df[f"P{year}"] = (shp_df[f"E{year}"] / 1000) * shp_df[find_str([scenario, year, "碳价"], olist)] - shp_df[f"C{year}"]

    return shp_df.loc[:, OUTPUT_COLS].set_crs(shp_crs)


# ============================================================
# 4. Parameter-check table export
# ============================================================
def export_parameter_check_tables(df_parameter_summary, df_shengcode):
    """Export per-province average wind/solar parameter tables, for manually
    checking whether the region-average fill-in looks reasonable."""
    cal3(df_parameter_summary, "风", df_shengcode).to_excel(os.path.join(TABLE_DIR, "feng.xlsx"))
    cal3(df_parameter_summary, "光", df_shengcode).to_excel(os.path.join(TABLE_DIR, "guang.xlsx"))


# ============================================================
# 5. Per-grid Cost / Emission / Profit calculation
# ============================================================
def run_grid_calculations(df_parameter_summary, df_hours, df_shengcode):
    for run in RUNS:
        shp_df = cal2(
            df_parameter_summary, df_hours, run["shp_path"], run["sel_str"], run["scenario"], df_shengcode
        )
        shp_df.to_file(run["output"])
        print(f"Saved: {run['output']}")


# ============================================================
# 6. Summary statistics
# ============================================================
def build_province_summary_workbook():
    """Walk every result shp in the RESULT3_DIR folder, sum each year's C/E/P
    totals by province, and combine into a single Excel workbook (one sheet
    per shp)."""
    path_lst = glob.glob(os.path.join(RESULT3_DIR, "*.shp"))
    sum_cols = [f"C{y}" for y in ALL_YEARS] + [f"E{y}" for y in ALL_YEARS] + [f"P{y}" for y in ALL_YEARS]

    name_lst, df_lst = [], []
    for path in path_lst:
        name = os.path.splitext(os.path.basename(path))[0]
        name_lst.append(name)
        df = gpd.read_file(path)
        df_province = df.groupby("Shengcode")[sum_cols].sum().reset_index()
        df_lst.append(df_province)

    output_path = os.path.join(TABLE_DIR, "province_results_20250625.xlsx")
    with pd.ExcelWriter(output_path) as writer:
        for name, df_province in zip(name_lst, df_lst):
            df_province.to_excel(writer, sheet_name=name, index=False)
    print(f"Saved: {output_path}")


def build_sensitivity_table(shp_path, output_path, df_shengcode):
    """For a single result shp, count (grid cells with capacity > 0) and sum
    (each metric) per year, by province - used for sensitivity analysis."""
    shp_df = gpd.read_file(shp_path)

    cols_to_count = [f"KW{y}" for y in ALL_YEARS]
    count_df = shp_df.groupby("Shengcode")[cols_to_count].apply(lambda x: (x > 0).sum()).reset_index()
    count_df.columns = ["Shengcode"] + [f"{col}_count" for col in cols_to_count]

    cols_to_sum = (
        [f"KW{y}" for y in ALL_YEARS]
        + [f"C{y}" for y in ALL_YEARS]
        + [f"E{y}" for y in ALL_YEARS]
        + [f"P{y}" for y in ALL_YEARS]
    )
    sum_df = shp_df.groupby("Shengcode")[cols_to_sum].sum().reset_index()
    sum_df.columns = ["Shengcode"] + [f"{col}_sum" for col in cols_to_sum]

    summary_df = pd.merge(count_df, sum_df, on="Shengcode")
    summary_df = pd.merge(df_shengcode, summary_df, on="Shengcode")
    summary_df.to_excel(output_path)
    print(f"Saved: {output_path}")


# ============================================================
# 7. Main entry point
# ============================================================
if __name__ == "__main__":
    df_hours, df_parameter_summary, df_shengcode = load_parameter_tables()

    export_parameter_check_tables(df_parameter_summary, df_shengcode)
    run_grid_calculations(df_parameter_summary, df_hours, df_shengcode)
    build_province_summary_workbook()

    build_sensitivity_table(
        os.path.join(RESULT3_DIR, "Solar_median_results_20250110gj.shp"),
        os.path.join(SENSITIVITY_DIR, "province_solar_median_kw2.xlsx"),
        df_shengcode,
    )
    build_sensitivity_table(
        os.path.join(RESULT3_DIR, "Solar_high_results_20250110gj.shp"),
        os.path.join(SENSITIVITY_DIR, "province_solar_high_kw2.xlsx"),
        df_shengcode,
    )
