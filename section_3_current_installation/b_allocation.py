from pathlib import Path
base_folder = Path(__file__).resolve().parents[3]

import arcpy
import os
import pandas as pd
import numpy as np

wind_turbines_grid = os.path.join(base_folder, r"processing\arcprojects\MyProject1\installation.gdb\grid_wind_turbine_count")
solar_panel_grid = os.path.join(base_folder, r"processing\arcprojects\MyProject1\installation.gdb\grid_solar_panel_area")
installation_table = os.path.join(base_folder, r"processing\tables\wind and solar development - EN.xlsx")

scratch_gdb = os.path.join(base_folder, r"processing\arcprojects\MyProject1\scratch.gdb")
installation_gdb = os.path.join(base_folder, r"processing\arcprojects\MyProject1\installation.gdb")

arcpy.env.overwriteOutput = True

# ---- 读取省级装机表 ----
installation_df = pd.read_excel(installation_table, sheet_name="2025", header=1)

pv_capacity_col = 'PV installed capacity(10MW)'
pv_generation_col = 'PV power generation(10^8kwh)'
onshore_wind_capacity_col = 'Onshore wind installed capacity (10 MW)'
wind_full_load_hours = 'Wind full load hours'
solar_full_load_hours = 'PV full load hours'
offshore_wind_capacity_col = 'Offshore wind installed capacity (10 MW)'

shengcode_col = 'Shengcode'
shengname_col = 'Shengname_cn'

# 单位换算常量
UNIT_10MW_TO_KW = 1e4          # 10 MW = 10,000 kW
UNIT_10E8KWH_TO_KWH = 1e8     # 10^8 kWh = 100,000,000 kWh

# 过滤有效省份（排除 TOTAL 行和无效代码）
province_df = installation_df[
    (installation_df[shengcode_col] > 0) &
    (installation_df[shengcode_col] != 100)
].copy()
province_df[shengcode_col] = province_df[shengcode_col].astype(int)

# ============================================================
# 1. 陆上风电：按省内网格风机数量占比分配装机，按满发小时算发电量
# ============================================================

# --- 读取风电网格属性表 ---
wind_fields = [f.name for f in arcpy.ListFields(wind_turbines_grid)]
wind_arr = arcpy.da.TableToNumPyArray(wind_turbines_grid, ['NID10', 'Shengcode', 'Wind_Turbine_Count'], skip_nulls=False)
wind_gdf = pd.DataFrame(wind_arr, columns=['NID10', 'Shengcode', 'Wind_Turbine_Count'])
wind_gdf['NID10'] = wind_gdf['NID10'].astype(int)
wind_gdf['Shengcode'] = wind_gdf['Shengcode'].astype(int)
wind_gdf['Wind_Turbine_Count'] = wind_gdf['Wind_Turbine_Count'].fillna(0).astype(float)

# 分离陆上 / 海上
onshore_wind = wind_gdf[wind_gdf['Shengcode'] != 100].copy()
offshore_wind = wind_gdf[wind_gdf['Shengcode'] == 100].copy()

# 省内风机数合计 & 占比
onshore_prov_sum = onshore_wind.groupby('Shengcode')['Wind_Turbine_Count'].transform('sum')
onshore_wind['prov_ratio'] = np.where(onshore_prov_sum > 0,
                                       onshore_wind['Wind_Turbine_Count'] / onshore_prov_sum, 0)

# 合并省级装机数据
prov_wind = province_df[[shengcode_col, onshore_wind_capacity_col, wind_full_load_hours]].copy()
prov_wind.rename(columns={shengcode_col: 'Shengcode'}, inplace=True)
prov_wind['Shengcode'] = prov_wind['Shengcode'].astype(int)

onshore_wind = onshore_wind.merge(prov_wind, on='Shengcode', how='left')

# 装机 kW = 省装机(10MW) × 占比 × 10000
onshore_wind['kw2025'] = (onshore_wind[onshore_wind_capacity_col].fillna(0)
                          * onshore_wind['prov_ratio']
                          * UNIT_10MW_TO_KW)

# 发电量 kWh = 装机 kW × 满发小时
onshore_wind['kwh2025'] = (onshore_wind['kw2025']
                           * onshore_wind[wind_full_load_hours].fillna(0))

# ============================================================
# 2. 海上风电：Shengcode=100 的网格按风机数量占比分配总量
# ============================================================

total_row = installation_df[installation_df['Shengname_cn'] == 'TOTAL'].iloc[0]
offshore_capacity_total = total_row[offshore_wind_capacity_col]
offshore_hours = installation_df[wind_full_load_hours].mean()

offshore_count_sum = offshore_wind['Wind_Turbine_Count'].sum()
offshore_wind['prov_ratio'] = np.where(offshore_count_sum > 0,
                                        offshore_wind['Wind_Turbine_Count'] / offshore_count_sum, 0)

offshore_wind['kw2025'] = (offshore_capacity_total
                           * offshore_wind['prov_ratio']
                           * UNIT_10MW_TO_KW)
offshore_wind['kwh2025'] = offshore_wind['kw2025'] * offshore_hours

# ============================================================
# 合并陆上+海上，写回风电网格
# ============================================================

wind_result = pd.concat([onshore_wind, offshore_wind], ignore_index=True)
wind_result = wind_result[['NID10', 'kw2025', 'kwh2025']].sort_values('NID10')
wind_result['NID10'] = wind_result['NID10'].astype(int)

for col in ['kw2025', 'kwh2025']:
    if col not in wind_fields:
        arcpy.management.AddField(wind_turbines_grid, col, 'DOUBLE')

with arcpy.da.UpdateCursor(wind_turbines_grid, ['NID10', 'kw2025', 'kwh2025']) as cur:
    for row in cur:
        nid = int(row[0])
        match = wind_result.loc[wind_result['NID10'] == nid]
        if not match.empty:
            row[1] = float(match['kw2025'].iloc[0])
            row[2] = float(match['kwh2025'].iloc[0])
        else:
            row[1], row[2] = 0.0, 0.0
        cur.updateRow(row)

print("风电网格写入完成（陆上+海上）")

# ============================================================
# 3. 光伏：按省内网格面积占比分配装机，按满发小时算发电量
# ============================================================

solar_fields = [f.name for f in arcpy.ListFields(solar_panel_grid)]
solar_arr = arcpy.da.TableToNumPyArray(solar_panel_grid, ['NID10', 'Shengcode', 'Solar_Area_m2'], skip_nulls=False)
solar_gdf = pd.DataFrame(solar_arr, columns=['NID10', 'Shengcode', 'Solar_Area_m2'])
solar_gdf['NID10'] = solar_gdf['NID10'].astype(int)
solar_gdf['Shengcode'] = solar_gdf['Shengcode'].astype(int)
solar_gdf['Solar_Area_m2'] = solar_gdf['Solar_Area_m2'].fillna(0).astype(float)

# 只处理陆上
solar_gdf = solar_gdf[solar_gdf['Shengcode'] != 100].copy()

# 省内面积合计 & 占比
prov_area_sum = solar_gdf.groupby('Shengcode')['Solar_Area_m2'].transform('sum')
solar_gdf['prov_ratio'] = np.where(prov_area_sum > 0,
                                    solar_gdf['Solar_Area_m2'] / prov_area_sum, 0)

# 合并省级光伏数据
prov_solar = province_df[[shengcode_col, pv_capacity_col, solar_full_load_hours]].copy()
prov_solar.rename(columns={shengcode_col: 'Shengcode'}, inplace=True)
prov_solar['Shengcode'] = prov_solar['Shengcode'].astype(int)

solar_gdf = solar_gdf.merge(prov_solar, on='Shengcode', how='left')

# 装机 kW
solar_gdf['kw2025'] = (solar_gdf[pv_capacity_col].fillna(0)
                        * solar_gdf['prov_ratio']
                        * UNIT_10MW_TO_KW)

# 发电量 kWh = 装机 × 满发小时
solar_gdf['kwh2025'] = (solar_gdf['kw2025']
                        * solar_gdf[solar_full_load_hours].fillna(0))

# 写回光伏网格
for col in ['kw2025', 'kwh2025']:
    if col not in solar_fields:
        arcpy.management.AddField(solar_panel_grid, col, 'DOUBLE')

solar_write = solar_gdf[['NID10', 'kw2025', 'kwh2025']].sort_values('NID10')
solar_write['NID10'] = solar_write['NID10'].astype(int)

with arcpy.da.UpdateCursor(solar_panel_grid, ['NID10', 'kw2025', 'kwh2025']) as cur:
    for row in cur:
        nid = int(row[0])
        match = solar_write.loc[solar_write['NID10'] == nid]
        if not match.empty:
            row[1] = float(match['kw2025'].iloc[0])
            row[2] = float(match['kwh2025'].iloc[0])
        else:
            row[1], row[2] = 0.0, 0.0
        cur.updateRow(row)

print("光伏网格写入完成")