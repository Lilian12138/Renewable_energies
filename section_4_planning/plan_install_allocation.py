"""
分省分年份网格装机容量分配脚本
逻辑：按 2030 → 2035 → 2040 → 2050 → 2060 顺序，每个年份以上一年为保底，增量按评分排序贪心填充。
输出列：NID50, NID10_INT, Shengcode, kw2025, kw2030, kw2035, kw2040, kw2050, kw2060
"""

from pathlib import Path
import os
import pandas as pd
import numpy as np

# ============================================================
# 1. 路径配置（与你原始代码一致）
# ============================================================
base_folder = Path(__file__).resolve().parents[3]

# 评分网格 (feature class in GDB)
score_path = os.path.join(
    base_folder,
    r"processing\arcprojects\MyProject1\MyProject1.gdb\CL_WGS84"
)
onshore_wind_score = "score_wind_onshore"
offshore_wind_score = "score_wind_offshore"
onshore_solar_score = "score_pv_onshore"

# 省级规划表
planning_path = os.path.join(
    base_folder,
    r"processing\tables\Planned installed capacity.xlsx"
)
wind_sheetname = "Wind"
solar_sheetname = "Solar"

# 网格装机潜力 shapefile（最大可装机容量 KW2）
potential_wind_shp = os.path.join(
    base_folder,
    r"processing\gisfiles\GridValidArea\grid10km_wind_valid_area_statistic.shp"
)
potential_cap_wind = "cap_kw_new"  # 潜力字段 → 重命名为 KW2

potential_solar_shp = os.path.join(
    base_folder,
    r"processing\gisfiles\GridValidArea\grid10km_solar_valid_area_statistic.shp"
)
potential_cap_solar = "cap_kw"  # 潜力字段 → 重命名为 KW2

# 网格装机现状 (feature class in GDB)，字段名均为 kw2025
current_wind_fc = os.path.join(
    base_folder,
    r"processing\arcprojects\MyProject1\installation.gdb\grid_wind_turbine_count"
)
current_solar_fc = os.path.join(
    base_folder,
    r"processing\arcprojects\MyProject1\installation.gdb\grid_solar_panel_area"
)
current_cap_field = "kw2025"  # 两个现状 fc 中的装机字段名

NID10_INT = "NID10_INT"

# 年份配置
YEARS = [2030, 2035, 2040, 2050, 2060]

# 两套情景的规划表列名
PROVINCE_PLAN_COLS = [f"{y}_province" for y in YEARS]
LOWCARBON_PLAN_COLS = [f"{y}_low_carbon" for y in YEARS]

# 两套情景的输出列名（网格级）
PROV_YEAR_COLS = [f"prov_{y}" for y in YEARS]       # province 情景
LC_YEAR_COLS   = [f"lc_{y}" for y in YEARS]         # low_carbon 情景
ALL_YEAR_COLS  = ["kw2025"] + PROV_YEAR_COLS + LC_YEAR_COLS

# 最终保留的输出列
OUTPUT_COLS = ["NID50", NID10_INT, "Shengcode"] + ALL_YEAR_COLS


# ============================================================
# 2. 数据读取（使用 arcpy 读 shapefile / feature class）
# ============================================================
def read_feature_to_df(fc_path, fields):
    """用 arcpy 将 feature class / shapefile 读为 DataFrame"""
    import arcpy
    all_fields = [f.name for f in arcpy.ListFields(fc_path)]
    use_fields = [f for f in fields if f in all_fields]
    data = [row for row in arcpy.da.SearchCursor(fc_path, use_fields)]
    return pd.DataFrame(data, columns=use_fields)


def ensure_nid10_int(fc_path, df):
    """
    确保 DataFrame 中有 NID10_INT 列（long 整型）。
    如果原数据只有 NID10（文本），则转换并重命名。
    """
    import arcpy
    all_fields = [f.name for f in arcpy.ListFields(fc_path)]

    if NID10_INT in df.columns:
        # 已有 NID10_INT，确保整型
        df[NID10_INT] = df[NID10_INT].astype(int)
    elif "NID10" in df.columns:
        # 只有 NID10，转为整型并重命名
        df[NID10_INT] = df["NID10"].astype(int)
        df = df.drop(columns=["NID10"])
    return df


def load_grid_data(potential_shp, potential_field, current_fc):
    """
    读取网格数据，合并潜力和现状两个数据源:
    返回含以下列的 DataFrame:
      - NID10_INT, NID50, Shengcode, KW2, kw2025
    """
    # 读取潜力数据（含 NID50），优先读 NID10_INT，回退到 NID10
    import arcpy
    pot_all_fields = [f.name for f in arcpy.ListFields(potential_shp)]
    nid_field_pot = NID10_INT if NID10_INT in pot_all_fields else "NID10"
    pot_fields = [nid_field_pot, "NID50", "Shengcode", potential_field]
    df_pot = read_feature_to_df(potential_shp, pot_fields)
    df_pot = ensure_nid10_int(potential_shp, df_pot)
    df_pot = df_pot.rename(columns={potential_field: "KW2"})
    df_pot["KW2"] = df_pot["KW2"].fillna(0)

    # 读取现状装机数据，优先读 NID10_INT，回退到 NID10
    cur_all_fields = [f.name for f in arcpy.ListFields(current_fc)]
    nid_field_cur = NID10_INT if NID10_INT in cur_all_fields else "NID10"
    cur_fields = [nid_field_cur, current_cap_field]
    df_cur = read_feature_to_df(current_fc, cur_fields)
    df_cur = ensure_nid10_int(current_fc, df_cur)
    df_cur = df_cur.rename(columns={current_cap_field: "kw2025"})
    df_cur["kw2025"] = df_cur["kw2025"].fillna(0)

    # 去重检查：如有重复 NID10_INT，按合计处理
    pot_dup = df_pot[NID10_INT].duplicated().sum()
    cur_dup = df_cur[NID10_INT].duplicated().sum()
    if pot_dup > 0:
        print(f"  ⚠ 潜力数据 NID10_INT 有 {pot_dup} 条重复，按合计处理")
        df_pot = df_pot.groupby([NID10_INT, "NID50", "Shengcode"], as_index=False)["KW2"].sum()
    if cur_dup > 0:
        print(f"  ⚠ 现状数据 NID10_INT 有 {cur_dup} 条重复，按合计处理")
        df_cur = df_cur.groupby(NID10_INT, as_index=False)["kw2025"].sum()

    # 合并：以潜力表为主，左连接现状
    df = df_pot.merge(df_cur, on=NID10_INT, how="left")
    df["kw2025"] = df["kw2025"].fillna(0)

    print(f"  潜力网格数: {len(df_pot)}, 现状网格数: {len(df_cur)}, 合并后: {len(df)}")

    return df


def load_score_data(score_fc, score_field):
    """从评分要素类读取评分，用于 join"""
    import arcpy
    all_fields = [f.name for f in arcpy.ListFields(score_fc)]
    nid_field = NID10_INT if NID10_INT in all_fields else "NID10"
    fields = [nid_field, score_field]
    df = read_feature_to_df(score_fc, fields)
    df = ensure_nid10_int(score_fc, df)
    df = df.rename(columns={score_field: "Score"})
    # 去重：如有重复取均值
    if df[NID10_INT].duplicated().any():
        df = df.groupby(NID10_INT, as_index=False)["Score"].mean()
    return df


# ============================================================
# 3. 核心分配函数（多年份递增）
# ============================================================
def plan_install_multiyear(df, shengcode, kw_targets, year_cols):
    """
    对单个省份，按年份顺序分配装机到网格。

    参数:
        df         : 该能源类型的全部网格 DataFrame
        shengcode  : 省代码
        kw_targets : list，各年份省级规划目标，单位 **万kW**
        year_cols  : list，对应输出列名 ['kw2030','kw2035','kw2040','kw2060']

    返回:
        该省分配完成的 DataFrame（含所有年份列）
    """
    df_prov = df.loc[df["Shengcode"] == shengcode].copy()

    if df_prov.empty:
        print(f"  ⚠ 省代码 {shengcode} 无对应网格，跳过")
        return df_prov

    # 按现状装机、评分、最大容量降序排列（优先填高分格子）
    df_prov = df_prov.sort_values(
        by=["kw2025", "Score", "KW2"], ascending=False
    ).reset_index(drop=True)

    prev_col = "kw2025"  # 起点列：现状装机

    for year_col, kw_target_wan in zip(year_cols, kw_targets):
        kw_target = kw_target_wan * 10000  # 万kW → kW

        # 以上一年结果为保底
        df_prov[year_col] = df_prov[prev_col].copy()

        already = df_prov[year_col].sum()
        remain = kw_target - already

        if remain <= 0:
            print(f"  {year_col}: 目标 {kw_target_wan:.1f}万kW 已由上一年满足，无需新增")
            prev_col = year_col
            continue

        allocated = 0
        for idx in df_prov.index:
            if remain <= 0:
                break
            cell_room = max(df_prov.at[idx, "KW2"] - df_prov.at[idx, year_col], 0)
            alloc = min(remain, cell_room)
            df_prov.at[idx, year_col] += alloc
            remain -= alloc
            allocated += alloc

        print(
            f"  {year_col}: 目标 {kw_target_wan:.1f}万kW, "
            f"新增 {allocated/10000:.1f}万kW, "
            f"未分配 {max(remain,0)/10000:.1f}万kW"
        )
        prev_col = year_col

    return df_prov


# ============================================================
# 4. 批量运行：对所有省份 × 所有年份
# ============================================================
def allocate_by_plan(df, plan_df, score_field, label):
    """按省份和评分字段分配，返回结果列表。"""
    print(f"读取{label}评分数据: {score_field}...")
    score_df = load_score_data(score_path, score_field)
    df = df.merge(score_df, on=NID10_INT, how="left")
    df["Score"] = df["Score"].fillna(0)

    results = []
    for _, row in plan_df.iterrows():
        shengcode = row["Shengcode"]

        prov_targets = [row[col] for col in PROVINCE_PLAN_COLS]
        print(f"\n省份: {shengcode} | {label} Province目标(万kW): {prov_targets}")
        df_prov = plan_install_multiyear(df, shengcode, prov_targets, PROV_YEAR_COLS)

        lc_targets = [row[col] for col in LOWCARBON_PLAN_COLS]
        print(f"省份: {shengcode} | {label} LowCarbon目标(万kW): {lc_targets}")
        df_lc = plan_install_multiyear(df, shengcode, lc_targets, LC_YEAR_COLS)

        if not df_prov.empty and not df_lc.empty:
            for col in LC_YEAR_COLS:
                df_prov[col] = df_lc[col].values
        results.append(df_prov)

    return results


def run_allocation(energy_type="wind"):
    """
    energy_type: 'wind' 或 'solar'
    """
    print(f"\n{'='*60}")
    print(f"  开始分配: {energy_type.upper()}")
    print(f"{'='*60}")

    # ---- 读取规划表 ----
    sheet = wind_sheetname if energy_type == "wind" else solar_sheetname
    print("读取规划表...")
    plan_df = pd.read_excel(planning_path, sheet_name=sheet)

    # ---- 读取网格潜力 + 现状数据 ----
    if energy_type == "wind":
        pot_shp = potential_wind_shp
        pot_field = potential_cap_wind
        cur_fc = current_wind_fc
    else:
        pot_shp = potential_solar_shp
        pot_field = potential_cap_solar
        cur_fc = current_solar_fc

    print("读取网格潜力 + 现状数据...")
    grid_df = load_grid_data(pot_shp, pot_field, cur_fc)

    results = []
    if energy_type == "wind":
        # 分离陆上/海上网格并分别使用对应评分字段
        onshore_grid = grid_df.loc[grid_df["Shengcode"] != 100].copy()
        offshore_grid = grid_df.loc[grid_df["Shengcode"] == 100].copy()

        onshore_plan = plan_df.loc[plan_df["Shengcode"] != 100].copy()
        offshore_plan = plan_df.loc[plan_df["Shengcode"] == 100].copy()

        if not onshore_plan.empty:
            results.extend(allocate_by_plan(onshore_grid, onshore_plan, onshore_wind_score, "Onshore"))
        else:
            print("  ⚠ Wind 规划表中未找到陆上省份行，跳过陆上分配")

        if not offshore_plan.empty:
            results.extend(allocate_by_plan(offshore_grid, offshore_plan, offshore_wind_score, "Offshore"))
        else:
            print("  ⚠ Wind 规划表中未找到 Shengcode=100 的海上行，跳过海上分配")
    else:
        results.extend(allocate_by_plan(grid_df, plan_df, onshore_solar_score, "Solar"))

    if not results:
        raise ValueError("没有生成任何分配结果，请检查规划表和网格数据")

    final_df = pd.concat(results, ignore_index=True)

    # ---- 只保留指定输出列 ----
    keep_cols = [c for c in OUTPUT_COLS if c in final_df.columns]
    final_df = final_df[keep_cols]

    # ---- 汇总检查 ----
    print(f"\n{'='*60}")
    print("分配结果汇总（万kW）:")
    print(f"{'='*60}")
    summary = final_df.groupby("Shengcode")[ALL_YEAR_COLS].sum() / 10000
    print(summary.to_string())

    total = final_df[ALL_YEAR_COLS].sum() / 10000
    print(f"\n全国合计(万kW):")
    for col in ALL_YEAR_COLS:
        print(f"  {col}: {total[col]:.1f}")

    summary_path = os.path.join(
        base_folder,
        rf"processing\tables\{energy_type}_allocation_summary.csv"
    )
    summary.to_csv(summary_path, encoding="utf-8-sig")
    print(f"分省汇总已导出: {summary_path}")

    return final_df


# ============================================================
# 5. 结果回写到 shapefile（可选）
# ============================================================
def write_results_to_shp(final_df, output_shp, grid_shp):
    """将分配结果 join 回 shapefile 并输出（含 kw2025~kw2060）"""
    import arcpy
    arcpy.env.overwriteOutput = True

    # 先导出为 CSV，再通过 arcpy join
    csv_path = output_shp.replace(".shp", "_result.csv")
    final_df.to_csv(csv_path, index=False)

    # 复制原始 shp
    arcpy.management.CopyFeatures(grid_shp, output_shp)

    # 删除多余字段，只保留 NID50, NID10_INT, Shengcode
    keep_fields = {"NID50", NID10_INT, "Shengcode", "cap_kw", "cap_kw_new",
                   "FID", "Shape", "OBJECTID", "OID"}
    all_fields = arcpy.ListFields(output_shp)
    drop_fields = [f.name for f in all_fields
                   if f.name not in keep_fields and f.type not in ("OID", "Geometry")]
    if drop_fields:
        arcpy.management.DeleteField(output_shp, drop_fields)

    # 添加年份字段并赋值
    for col in ALL_YEAR_COLS:
        arcpy.management.AddField(output_shp, col, "DOUBLE")

    # 用 UpdateCursor 写入（处理 NID10_INT 可能重复的情况）
    # 按 NID10_INT 汇总（如有重复取合计）
    agg_df = final_df.groupby(NID10_INT)[ALL_YEAR_COLS].sum()
    lookup = agg_df.to_dict("index")

    update_fields = [NID10_INT] + ALL_YEAR_COLS

    with arcpy.da.UpdateCursor(output_shp, update_fields) as cursor:
        for row in cursor:
            nid = row[0]
            if nid in lookup:
                for i, col in enumerate(ALL_YEAR_COLS, start=1):
                    row[i] = lookup[nid][col]
                cursor.updateRow(row)

    print(f"结果已写入: {output_shp}")


# ============================================================
# 6. 主入口
# ============================================================
if __name__ == "__main__":
    # ---- 风电分配 ----
    wind_result = run_allocation("wind")

    # ---- 光伏分配 ----
    solar_result = run_allocation("solar")

    # ---- 可选：导出结果 ----
    output_folder = os.path.join(base_folder, r"processing\gisfiles\GridValidArea")

    wind_output = os.path.join(output_folder, "grid10km_wind_planned.shp")
    write_results_to_shp(wind_result, wind_output, potential_wind_shp)

    solar_output = os.path.join(output_folder, "grid10km_solar_planned.shp")
    write_results_to_shp(solar_result, solar_output, potential_solar_shp)

    print("\n✅ 全部分配完成！")