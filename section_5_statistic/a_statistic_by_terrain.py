"""
按地形类型 + 省份分类统计风电/光伏装机容量
逻辑：取网格中心点，提取地形栅格(Reclass_geomor.tif)的值 1/2/3/4，
      分别对应 Plain, Hills, Mountainous, Complex terrain，
      再按 (省份 Shengcode, 地形类型) 对各年份/情景装机容量求和。
"""

from pathlib import Path
import os
import arcpy
import pandas as pd
import os 

base_folder = Path(__file__).resolve().parents[3]

wind_terrain_shp = os.path.join(base_folder, r"processing\gisfiles\GridValidArea\grid10km_wind_planned.shp")
solar_terrain_shp = os.path.join(base_folder, r"processing\gisfiles\GridValidArea\grid10km_solar_planned.shp")
terrain_type = os.path.join(base_folder, r"processing\gisfiles\terrain_type\Reclass_geom1.tif")


# 按照网格中心点提取地形类型 1，2，3，4 分别代表Plain, Hills, Mountainous, Complex terrain
NID10_INT = "NID10_INT"
SHENGCODE_FIELD = "Shengcode"
TERRAIN_FIELD = "terrain"
TERRAIN_MAP = {1: "Plain", 2: "Hills", 3: "Mountainous", 4: "Complex terrain"}
TERRAIN_ORDER = ["Plain", "Hills", "Mountainous", "Complex terrain"]

# 省的id是Shengcode
PROVINCE_MAP = {
    100: "Offshore", 65: "Xinjiang", 15: "Inner Mongolia", 23: "Heilongjiang",
    62: "Gansu", 63: "Qinghai", 22: "Jilin", 13: "Hebei", 37: "Shandong",
    41: "Henan", 21: "Liaoning", 45: "Guangxi", 53: "Yunnan", 34: "Anhui",
    32: "Jiangsu", 14: "Shanxi", 61: "Shaanxi", 43: "Hunan", 42: "Hubei",
    44: "Guangdong", 52: "Guizhou", 36: "Jiangxi", 64: "Ningxia", 46: "Hainan",
    33: "Zhejiang", 51: "Sichuan", 35: "Fujian", 54: "Xizang", 12: "Tianjin",
    50: "Chongqing", 31: "Shanghai", 11: "Beijing",
}
PROVINCE_ORDER = list(PROVINCE_MAP.values())

scratch_gdb = os.path.join(base_folder, r"processing\arcprojects\MyProject1\scratch.gdb")
terrain_filled_raster = os.path.join(scratch_gdb, "terrain_filled")
OFFSHORE_SHENGCODE = 100
TERRAIN_FIELD_FILLED = "terrain_fill"

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")


def read_fc_to_df(fc_path, fields):
    """用 arcpy 将 feature class / shapefile 读为 DataFrame"""
    all_fields = [f.name for f in arcpy.ListFields(fc_path)]
    use_fields = [f for f in fields if f in all_fields]
    data = [row for row in arcpy.da.SearchCursor(fc_path, use_fields)]
    return pd.DataFrame(data, columns=use_fields)


def build_filled_terrain_raster():
    """用欧氏分配(最近邻)填补地形栅格的 NoData 区域，避免网格中心点落在 NoData 处而提取不到地形值"""
    if not arcpy.Exists(terrain_filled_raster):
        filled = arcpy.sa.EucAllocation(terrain_type)
        filled.save(terrain_filled_raster)
    return terrain_filled_raster


def extract_terrain_by_centroid(grid_shp, out_name):
    """取网格中心点，提取地形栅格值；若中心点落在 NoData 处，除 Shengcode=100(海上)外，
    其余省份改用最近有效栅格的地形值"""
    centroids = os.path.join(scratch_gdb, out_name)
    arcpy.management.FeatureToPoint(grid_shp, centroids, "CENTROID")

    filled_raster = build_filled_terrain_raster()
    arcpy.sa.ExtractMultiValuesToPoints(
        centroids, [[terrain_type, TERRAIN_FIELD], [filled_raster, TERRAIN_FIELD_FILLED]]
    )

    df_terrain = read_fc_to_df(centroids, [NID10_INT, SHENGCODE_FIELD, TERRAIN_FIELD, TERRAIN_FIELD_FILLED])
    shengcode_num = pd.to_numeric(df_terrain[SHENGCODE_FIELD], errors="coerce")

    terrain_code = df_terrain[TERRAIN_FIELD]
    need_fill = terrain_code.isna() & (shengcode_num != OFFSHORE_SHENGCODE)
    terrain_code = terrain_code.where(~need_fill, df_terrain[TERRAIN_FIELD_FILLED])

    df_terrain["terrain_type"] = terrain_code.map(TERRAIN_MAP)
    return df_terrain[[NID10_INT, "terrain_type"]]


# 分地形统计装机量
wind_province_columns = ['kw2025'] + [f'prov_{i}' for i in [2030,2035,2040,2050,2060]]
wind_low_carbon_columns = ['kw2025'] + [f'lc_{i}' for i in [2030,2035,2040,2050,2060]]
solar_province_columns = ['kw2025'] + [f'prov_{i}' for i in [2030,2035,2040,2050,2060]]
solar_low_carbon_columns = ['kw2025'] + [f'lc_{i}' for i in [2030,2035,2040,2050,2060]]


# 输出tables
output_folder = os.path.join(base_folder, r"processing\tables")
output_excel_path = os.path.join(output_folder, "statistic_by_province_terrain.xlsx")


YEARS = [2025, 2030, 2035, 2040, 2050, 2060]


def statistic_by_terrain(grid_shp, energy_type):
    """提取网格中心点地形类型，并按 (Shengcode, 省份, 地形) 统计 province / low_carbon 两套情景装机量"""
    print(f"\n{'='*60}")
    print(f"  分省分地形统计装机量: {energy_type.upper()}")
    print(f"{'='*60}")

    province_columns = wind_province_columns if energy_type == "wind" else solar_province_columns
    low_carbon_columns = wind_low_carbon_columns if energy_type == "wind" else solar_low_carbon_columns
    all_columns = sorted(set(province_columns) | set(low_carbon_columns), key=(province_columns + low_carbon_columns).index)

    df_terrain = extract_terrain_by_centroid(grid_shp, f"{energy_type}_centroids")
    df_attr = read_fc_to_df(grid_shp, [NID10_INT, SHENGCODE_FIELD] + all_columns)

    df = df_attr.merge(df_terrain, on=NID10_INT, how="left")
    df["terrain_type"] = df["terrain_type"].fillna("Unknown")
    df["shengcode"] = pd.to_numeric(df[SHENGCODE_FIELD], errors="coerce").astype("Int64")
    df["province"] = df["shengcode"].map(lambda c: PROVINCE_MAP.get(int(c), "Unknown") if pd.notna(c) else "Unknown")

    province_cat = pd.Categorical(df["province"], categories=PROVINCE_ORDER + ["Unknown"], ordered=True)
    terrain_cat = pd.Categorical(df["terrain_type"], categories=TERRAIN_ORDER + ["Unknown"], ordered=True)
    df["province"] = province_cat
    df["terrain_type"] = terrain_cat

    summary = df.groupby(["province", "shengcode", "terrain_type"], observed=True)[all_columns].sum()
    summary = summary.sort_index().reset_index()

    out = summary[["shengcode", "province", "terrain_type"]].rename(columns={"shengcode": "Shengcode", "province": "Province"})
    for year in YEARS:
        prov_col = "kw2025" if year == 2025 else f"prov_{year}"
        out[f"{year}_province"] = summary[prov_col]
    for year in YEARS:
        lc_col = "kw2025" if year == 2025 else f"lc_{year}"
        out[f"{year}_low_carbon"] = summary[lc_col]

    print(out.to_string(index=False))
    return df, out


if __name__ == "__main__":
    wind_df, wind_summary = statistic_by_terrain(wind_terrain_shp, "wind")
    solar_df, solar_summary = statistic_by_terrain(solar_terrain_shp, "solar")

    with pd.ExcelWriter(output_excel_path, engine="openpyxl") as writer:
        wind_summary.to_excel(writer, sheet_name="wind", index=False)
        solar_summary.to_excel(writer, sheet_name="solar", index=False)

    print(f"\n✅ 分地形统计完成！已导出 -> {output_excel_path}")