from pathlib import Path
base_folder = Path(__file__).resolve().parents[3]

import arcpy
import os

wind_turbines_path = os.path.join(base_folder, r"processing\gisfiles\wind_solar_distribution_202605\windturbines.shp")
solar_panel_path = os.path.join(base_folder, r"processing\gisfiles\wind_solar_distribution_202605\solar_panel.shp")
grid_10km = os.path.join(base_folder, r"processing\gisfiles\limited_factors\RedlandGrid10km\CL_WGS84.shp")
scratch_gdb = os.path.join(base_folder, r"processing\arcprojects\MyProject1\scratch.gdb")
installation_gdb = os.path.join(base_folder, r"processing\arcprojects\MyProject1\installation.gdb")

arcpy.env.overwriteOutput = True

# ---------- a. 统计每个网格内风机数量 ----------

wind_output = os.path.join(installation_gdb, "grid_wind_turbine_count")

arcpy.analysis.SpatialJoin(
    target_features=grid_10km,
    join_features=wind_turbines_path,
    out_feature_class=wind_output,
    join_operation="JOIN_ONE_TO_ONE",
    join_type="KEEP_ALL",
    match_option="CONTAINS"
)
# SpatialJoin 默认生成 Join_Count 字段，即每个网格包含的风机数量
# 重命名字段便于识别
arcpy.management.AlterField(wind_output, "Join_Count", "Wind_Turbine_Count", "Wind_Turbine_Count")

print("风机数量统计完成 ->", wind_output)

# ---------- b. 统计每个网格内光伏板面积(m²) ----------

# 先为 solar_panel 计算面积字段（投影到等面积坐标系计算 m²）
# 添加面积字段
arcpy.management.AddField(solar_panel_path, "Area_m2", "DOUBLE")

# 用 WGS84 数据时，geometry area 默认是度，需要指定单位
# 使用 CalculateGeometryAttributes 以 m² 为单位计算
arcpy.management.CalculateGeometryAttributes(
    solar_panel_path,
    [["Area_m2", "AREA_GEODESIC"]],
    area_unit="SQUARE_METERS"
)

solar_output = os.path.join(installation_gdb, "grid_solar_panel_area")

# 空间连接，对 Area_m2 求和
field_mappings = arcpy.FieldMappings()
field_mappings.addTable(grid_10km)
field_mappings.addTable(solar_panel_path)

# 找到 Area_m2 字段，设置合并规则为 Sum
area_idx = field_mappings.findFieldMapIndex("Area_m2")
area_fm = field_mappings.getFieldMap(area_idx)
area_fm.mergeRule = "Sum"
field_mappings.replaceFieldMap(area_idx, area_fm)

arcpy.analysis.SpatialJoin(
    target_features=grid_10km,
    join_features=solar_panel_path,
    out_feature_class=solar_output,
    join_operation="JOIN_ONE_TO_ONE",
    join_type="KEEP_ALL",
    field_mapping=field_mappings,
    match_option="CONTAINS"
)

# 重命名
arcpy.management.AlterField(solar_output, "Area_m2", "Solar_Area_m2", "Solar_Area_m2")

print("光伏面积统计完成 ->", solar_output)