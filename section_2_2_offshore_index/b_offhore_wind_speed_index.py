# -*- coding: utf-8 -*-
"""
计算 100m 平均风速适宜性指数 (off_idx_wind)
流程:
  0. 清洗风速栅格 (剔除 NoData 填充值 / 非法值)
  1. 分区统计: 每个 10km 格网的平均风速
  2. 质心提取: 作为无分区统计结果格网的兜底
  3. 阈值筛选 (>= 6 m/s) + 分位数截断归一化到 0.2-1
  4. 写回原始要素类 CL_WGS84
"""

from pathlib import Path
import os

base_folder = Path(__file__).resolve().parents[3]
os.environ['GDAL_DATA'] = r'D:\installs\ArcGIS\Pro\Resources\pedata\gdaldata'

import arcpy
import numpy as np

# ----------------------------- 参数 -----------------------------
wind_power_path = os.path.join(base_folder, r"processing\gisfiles\windspeed100m\merged.tif")
grid_10km = os.path.join(base_folder, r"processing\arcprojects\MyProject1\MyProject1.gdb\CL_WGS84")

scratch_gdb = os.path.join(base_folder, r"processing\arcprojects\MyProject1\scratch.gdb")

GRID_ID_FIELD = "NID10_INT"
GRID_FILTER = "Shengcode = 100"

WIND_THRESHOLD = 6.0        # 达标下限 (m/s)
WIND_VALID_MIN = 0.0        # 物理合法区间下限
WIND_VALID_MAX = 60.0       # 物理合法区间上限 (超出视为填充值/异常)
PCT_LOW, PCT_HIGH = 2, 98   # 归一化分位数截断,防止离群点压扁主体分布

OUT_FIELD = "off_idx_wind"

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
arcpy.env.workspace = scratch_gdb


def add_field_safe(fc, field_name, field_type):
    existing = [f.name for f in arcpy.ListFields(fc)]
    if field_name not in existing:
        arcpy.management.AddField(fc, field_name, field_type)


def delete_if_exists(path):
    if arcpy.Exists(path):
        arcpy.management.Delete(path)


# ------------------- Step 0: 栅格诊断与清洗 -------------------
ras = arcpy.Raster(wind_power_path)
print("=== 原始栅格诊断 ===")
print(f"  noDataValue : {ras.noDataValue}")
print(f"  min / max   : {ras.minimum} / {ras.maximum}")

wind_clean = os.path.join(scratch_gdb, "wind_clean")

if ras.minimum is not None and ras.minimum >= WIND_VALID_MIN \
        and ras.maximum is not None and ras.maximum <= WIND_VALID_MAX:
    print("  栅格取值在合法区间内, 跳过清洗")
    wind_raster = wind_power_path
else:
    print(f"  检测到非法值, 将 (<{WIND_VALID_MIN} 或 >{WIND_VALID_MAX}) 置为 NoData ...")
    delete_if_exists(wind_clean)
    cleaned = arcpy.sa.SetNull((ras < WIND_VALID_MIN) | (ras > WIND_VALID_MAX), ras)
    cleaned.save(wind_clean)
    wind_raster = wind_clean
    ras2 = arcpy.Raster(wind_clean)
    print(f"  清洗后 min / max: {ras2.minimum} / {ras2.maximum}")

# ------------------- 落盘筛选后的格网 -------------------
grid_layer = "grid_filtered_a"
arcpy.MakeFeatureLayer_management(grid_10km, grid_layer, GRID_FILTER)

grid_fc = os.path.join(scratch_gdb, "grid_filtered_a")
delete_if_exists(grid_fc)
arcpy.management.CopyFeatures(grid_layer, grid_fc)

total_cells = int(arcpy.management.GetCount(grid_fc)[0])
print(f"\n筛选后格网数: {total_cells}")

# ------------------- Step 1: 分区统计 (MEAN) -------------------
wind_table = os.path.join(scratch_gdb, "wind_zonal")
delete_if_exists(wind_table)
arcpy.sa.ZonalStatisticsAsTable(grid_fc, GRID_ID_FIELD, wind_raster,
                                wind_table, "DATA", "MEAN")

wind_zonal_dict = {}
for fid, mean_v in arcpy.da.SearchCursor(wind_table, [GRID_ID_FIELD, "MEAN"]):
    # 二次保险: 即使清洗后仍出现非法均值, 一律丢弃走兜底
    if mean_v is not None and WIND_VALID_MIN <= mean_v <= WIND_VALID_MAX:
        wind_zonal_dict[fid] = mean_v

print(f"分区统计有效格网数: {len(wind_zonal_dict)}")

# ------------------- Step 2: 质心提取兜底 -------------------
centroids = os.path.join(scratch_gdb, "grid_centroids_a")
delete_if_exists(centroids)
arcpy.management.FeatureToPoint(grid_fc, centroids, "CENTROID")
arcpy.sa.ExtractMultiValuesToPoints(centroids, [[wind_raster, "wind_pt"]])

wind_pt_dict = {}
with arcpy.da.SearchCursor(centroids, [GRID_ID_FIELD, "wind_pt"]) as cur:
    for fid, v in cur:
        if v is not None and WIND_VALID_MIN <= v <= WIND_VALID_MAX:
            wind_pt_dict[fid] = v

# ------------------- 合并 + 阈值筛选 -------------------
wind_values = {}        # 达标格网: fid -> 风速
below_threshold = 0     # 有数据但不达标
no_data_cells = 0       # 完全无数据

with arcpy.da.SearchCursor(grid_fc, [GRID_ID_FIELD]) as cur:
    for (fid,) in cur:
        val = wind_zonal_dict.get(fid)
        if val is None:
            val = wind_pt_dict.get(fid)
        if val is None:
            no_data_cells += 1
        elif val >= WIND_THRESHOLD:
            wind_values[fid] = val
        else:
            below_threshold += 1

print(f"\n=== 阈值筛选结果 (>= {WIND_THRESHOLD} m/s) ===")
print(f"  达标格网   : {len(wind_values)}")
print(f"  不达标格网 : {below_threshold}")
print(f"  无数据格网 : {no_data_cells}")

if not wind_values:
    raise RuntimeError("没有任何格网达到风速阈值, 请检查栅格单位/阈值设置!")

# ------------------- Step 3: 分位数截断归一化 (0.2-1) -------------------
valid_arr = np.array(list(wind_values.values()), dtype=float)
print(f"\n=== 达标格网风速分布 (m/s) ===")
print(f"  min={valid_arr.min():.2f}  max={valid_arr.max():.2f}")
print(f"  P5={np.percentile(valid_arr, 5):.2f}  "
      f"P50={np.percentile(valid_arr, 50):.2f}  "
      f"P95={np.percentile(valid_arr, 95):.2f}")

low = float(np.percentile(valid_arr, PCT_LOW))
high = float(np.percentile(valid_arr, PCT_HIGH))
val_range = high - low
print(f"  归一化区间 (P{PCT_LOW}-P{PCT_HIGH}): [{low:.2f}, {high:.2f}]")

add_field_safe(grid_fc, OUT_FIELD, "DOUBLE")
with arcpy.da.UpdateCursor(grid_fc, [GRID_ID_FIELD, OUT_FIELD]) as cur:
    for row in cur:
        val = wind_values.get(row[0])
        if val is None:
            row[1] = 0.0                    # 不达标或无数据
        elif val_range > 0:
            v = min(max(val, low), high)    # 截断到分位区间
            row[1] = 0.2 + (v - low) / val_range * 0.8
        else:
            row[1] = 1.0                    # 所有达标格网风速相同
        cur.updateRow(row)

# ------------------- Step 4: 写回原始要素类 -------------------
add_field_safe(grid_10km, OUT_FIELD, "DOUBLE")
join_dict = {r[0]: r[1] for r in
             arcpy.da.SearchCursor(grid_fc, [GRID_ID_FIELD, OUT_FIELD])}

with arcpy.da.UpdateCursor(grid_10km, [GRID_ID_FIELD, OUT_FIELD], GRID_FILTER) as cur:
    for row in cur:
        row[1] = join_dict.get(row[0], 0.0)
        cur.updateRow(row)

# ------------------- 结果分布检查 -------------------
idx_arr = np.array([v for v in join_dict.values()], dtype=float)
nonzero = idx_arr[idx_arr > 0]
print(f"\n=== off_idx_wind 结果分布 ===")
print(f"  为 0 (不达标/无数据): {int((idx_arr == 0).sum())}")
if nonzero.size:
    print(f"  非 0 部分: min={nonzero.min():.3f}  "
          f"P50={np.percentile(nonzero, 50):.3f}  max={nonzero.max():.3f}")

print(f"\nDone: {OUT_FIELD} (0.2-1, 风速 >= {WIND_THRESHOLD} m/s, "
      f"P{PCT_LOW}-P{PCT_HIGH} 截断归一化) 已写回 CL_WGS84")