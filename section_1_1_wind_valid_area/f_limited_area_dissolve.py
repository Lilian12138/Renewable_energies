from pathlib import Path
base_folder = Path(__file__).resolve().parents[3]

import os
import time
import arcpy

t0 = time.time()
def elapsed():
    return f"{time.time() - t0:.1f}s"

# ---------------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------------
path_list = [
    r"processing\gisfiles\limited_factors\water\water.shp",
    r"processing\gisfiles\limited_factors\nature_area\natureArea.shp",
    r"processing\gisfiles\scratch\limitedAreasShp\windTurbine_points_buffer500m_singparts.shp",
    r"processing\gisfiles\wind_solar_distribution_202605\solar_panel.shp",
    r"processing\gisfiles\scratch\limitedAreasShp\dem_above_4000m.shp",
    r"processing\gisfiles\scratch\limitedAreasShp\dissolve_roads_railways_buffer1000m.shp",
    r"processing\gisfiles\scratch\limitedAreasShp\windspeed_below_4p5.shp",
    r"processing\gisfiles\scratch\limitedAreasShp\slope_above_30.shp"
]

# ---------------------------------------------------------------------------
# Output files
# ---------------------------------------------------------------------------
output_folder = os.path.join(base_folder, r"processing\gisfiles\scratch\limitedAreasShp")
output_path = os.path.join(output_folder, "wind_limited_areas_dissolved.shp")
os.makedirs(output_folder, exist_ok=True)

arcpy.env.overwriteOutput = True

# ---------------------------------------------------------------------------
# Check input files, use the first valid file's coordinate system as reference
# ---------------------------------------------------------------------------
valid_inputs = []
for p in path_list:
    full_path = os.path.join(base_folder, p)
    if arcpy.Exists(full_path):
        valid_inputs.append(full_path)
        print(f"[找到] {os.path.basename(full_path)}")
    else:
        print(f"[跳过] 文件不存在: {full_path}")

if not valid_inputs:
    print("未找到任何输入文件，退出。")
else:
    arcpy.env.outputCoordinateSystem = arcpy.Describe(valid_inputs[0]).spatialReference

    # Step 1: Merge all layers into memory (avoid writing temporary disk files)
    merged = r"in_memory\temp_merged"
    print(f"\n合并 {len(valid_inputs)} 个图层...  {elapsed()}")
    arcpy.management.Merge(valid_inputs, merged)
    print(f"合并完成  {elapsed()}")

    # Step 2: Dissolve into a single MultiPolygon feature
    print(f"开始溶解...  {elapsed()}")
    arcpy.management.Dissolve(merged, output_path, multi_part="MULTI_PART")
    print(f"溶解完成  {elapsed()}")

    arcpy.management.Delete(merged)
    print(f"\n输出文件: {output_path}  总耗时 {elapsed()}")
