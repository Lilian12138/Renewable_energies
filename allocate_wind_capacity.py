"""
Wind Grid Capacity Allocation (per-scenario, independent-year version)
========================================================================

Allocates each province's planned wind capacity targets (2025/2030/2035/2040/
2050/2060, in units of 10,000 kW / "wan-kW") to the 10km wind grid, year by
year and independently, under three development scenarios: Low / Median / High.

Allocation algorithm (plan_install):
    For a given province's grid cells, sort descending by
    [already installed in 2022 (s2022), score (Score), potential cap (KW2)].
    Subtract the 2022 baseline capacity from that year's target, then fill the
    remaining target into cells in sorted order, capping each cell at its
    potential ceiling KW2.

Note: this script computes each year's incremental target independently -
there is no cascading relationship where one year's allocation result is used
as the floor for the next year (unlike the newer algorithm in
plan_install_allocation.py; see that file's docstring for the comparison).
This version is kept for reference and to reproduce historical results.

Ported from: ProcessingWind.ipynb (restructured for readability; the original
calculation logic and results are unchanged).
"""

import os

import geopandas as gpd
import pandas as pd

# ============================================================
# 1. Paths and parameters
# ============================================================
WIND_GRID_SHP = r"D:\GoldWindProject2023\RENEWABLE_RESOURCE\Processing\MapsData\WIND_GRID10KM_20230914.shp"

PLAN_TABLE = r"D:\GoldWindProject2023\RENEWABLE_RESOURCE\00 原始数据\中国海上风光装机容量_20230319.xlsx"
PLAN_SHEET = "风电分省装机容量"

OUTPUT_DIR = r"D:\GoldWindProject2023\RENEWABLE_RESOURCE\Processing\MapsData"
TABLE_OUTPUT_DIR = r"D:\GoldWindProject2023\RENEWABLE_RESOURCE\Processing\TableData"

# Planning-table column names for each of the three development scenarios
SCENARIO_PLAN_COLUMNS = {
    "low": ["2025年_低发展情景", "2030年_低发展情景", "2035年_低发展情景",
            "2040年_低发展情景", "2050年_低发展情景", "2060年_低发展情景"],
    "median": ["2025年_中", "2030年_中", "2035年_中", "2040年_中", "2050年_中", "2060年_中"],
    "high": ["2025年_高发展情景", "2030年_高发展情景", "2035年_高发展情景",
             "2040年_高发展情景", "2050年_高发展情景", "2060年_高发展情景"],
}

SCENARIO_OUTPUT_SHP = {
    "low": os.path.join(OUTPUT_DIR, "wind_low_20240226.shp"),
    "median": os.path.join(OUTPUT_DIR, "wind_median_20240226.shp"),
    "high": os.path.join(OUTPUT_DIR, "wind_high_20240226.shp"),
}


# ============================================================
# 2. Core allocation function
# ============================================================
def plan_install(df, shengcode, kw_year, year_col):
    """Greedily allocate one province's, one year's planned target (in wan-kW) to grid cells."""
    # Restrict to a single province
    df = df.loc[df["Shengcode"] == shengcode]
    kw_year = kw_year * 10000
    df_sorted = df.sort_values(by=["s2022", "Score", "KW2"], ascending=False)

    # The remaining amount to allocate updates dynamically as cells are filled; the
    # starting amount available for allocation must also account for what was already
    # installed in 2022.
    kwyear_updated = kw_year - df_sorted["KW2022"].sum()
    print(kwyear_updated)

    # Before allocating, seed the year column with whatever was already installed in 2022
    df_sorted[year_col] = df_sorted["KW2022"]

    for index, cell in df_sorted.iterrows():
        if kwyear_updated > 0:  # allocate while there is still a remaining amount
            # How much room is left in this cell?
            cell_room = max(cell["KW2"] - cell[year_col], 0)
            if kwyear_updated > cell_room:
                # The remainder is enough to fill this cell completely
                cell[year_col] = cell[year_col] + cell_room
                kwyear_updated = kwyear_updated - cell_room
            else:
                # The remainder doesn't fill the cell's max capacity (KW2); put all of it in
                cell[year_col] = cell[year_col] + kwyear_updated
                kwyear_updated = kwyear_updated - kwyear_updated
            df_sorted.at[index, year_col] = cell[year_col]
            print(f"Cell: {index}, remaining to allocate: {kwyear_updated}, this cell's room: {cell_room}")
        else:
            break
    return df_sorted


def allocate_scenario(scenario):
    """For one development scenario (low/median/high), allocate the province-level plan
    to the wind grid year by year and write the result to a shapefile."""
    print(f"\n{'=' * 60}\nWind allocation scenario: {scenario}\n{'=' * 60}")

    points = gpd.read_file(WIND_GRID_SHP)

    # Drop unneeded columns.
    # (The original "low" scenario cell also had a commented-out line here -
    #  points.drop(columns=['hours','KWH2_2022','KW2_Unif_K'], inplace=True) -
    #  which was never actually executed, so it is not carried over.)
    points.drop(columns=["kw2025", "kw2030", "kw2035"], inplace=True)

    # Exclude grid cells that already exceed their 2022 potential; clamp negatives to 0
    points["exclude_2022"] = points["KW2"] - points["KW2022"]
    points.loc[points["exclude_2022"] < 0, ["exclude_2022"]] = 0

    # Open the planning-target table
    plan_df = pd.read_excel(PLAN_TABLE, sheet_name=PLAN_SHEET)
    shengcode_list = plan_df["Shengcode"].drop_duplicates().dropna().tolist()

    # Set the 2022 baseline flag
    points.loc[(points["KW2022"] > 0) & (points["exclude_2022"] > 0), ["s2022"]] = 1
    points["s2022"].fillna(0, inplace=True)

    points_copy = points.copy()
    for plan in SCENARIO_PLAN_COLUMNS[scenario]:
        df_list = []
        for shengcode in shengcode_list:
            kw_year = plan_df.loc[plan_df["Shengcode"] == shengcode, plan].values[0]
            year_col = "KW" + plan.split("年")[0]
            df = points_copy.loc[points_copy["Shengcode"] == shengcode]
            df_list.append(plan_install(df, shengcode, kw_year, year_col))
        df_kw = pd.concat(df_list)
        points_copy = points_copy.merge(df_kw.loc[:, ["NID10a", year_col]], how="left", on="NID10a")

    output_path = SCENARIO_OUTPUT_SHP[scenario]
    points_copy.to_file(output_path)
    print(f"Saved: {output_path}")
    return points_copy


# ============================================================
# 3. Province-level summary check (QA)
# ============================================================
def summarize_scenario(output_shp_path, table_output_path):
    """Summarize each year's total capacity by province, to check the allocation
    result against the planning targets."""
    df = gpd.read_file(output_shp_path)
    year_cols = ["KW2", "KW2022", "KW2025", "KW2030", "KW2035", "KW2040", "KW2050", "KW2060"]
    groupdf = df.loc[:, ["Shengcode"] + year_cols].groupby("Shengcode").agg("sum").reset_index()

    plan_df = pd.read_excel(PLAN_TABLE, sheet_name=PLAN_SHEET)
    groupdf = groupdf.merge(plan_df.loc[:, ["地区", "Shengcode"]], on="Shengcode", how="left")
    groupdf.to_excel(table_output_path)
    print(f"Province summary saved: {table_output_path}")
    return groupdf


# ============================================================
# 4. Main entry point
# ============================================================
if __name__ == "__main__":
    allocate_scenario("low")
    allocate_scenario("median")
    allocate_scenario("high")

    # In the original code, the "low" scenario's province summary check was commented out
    # (never executed); the same default behavior is preserved here.
    # summarize_scenario(SCENARIO_OUTPUT_SHP["low"], os.path.join(TABLE_OUTPUT_DIR, "groupdf_low.xlsx"))
    summarize_scenario(SCENARIO_OUTPUT_SHP["median"], os.path.join(TABLE_OUTPUT_DIR, "groupdf_median.xlsx"))
    summarize_scenario(SCENARIO_OUTPUT_SHP["high"], os.path.join(TABLE_OUTPUT_DIR, "groupdf_high.xlsx"))
