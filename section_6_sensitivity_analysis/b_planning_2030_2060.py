"""
分省分年份网格装机容量分配脚本
逻辑：按 2030 → 2035 → 2040 → 2050 → 2060 顺序，每个年份以上一年为保底，增量按评分排序贪心填充。
支持不同评分情景；风电按陆上/海上拆分，海上目标从规划表中 Shengcode = 100 行读取。
输出：写入 sensitivity.gdb，按情景名称生成要素类。
"""

from pathlib import Path
import os
import pandas as pd

base_folder = Path(__file__).resolve().parents[3]

# 评分网格 (feature class in GDB)
score_path = os.path.join(
    base_folder,
    r"processing\arcprojects\MyProject1\sensitivity.gdb\CL_WGS84"
)
output_gdb = os.path.join(
    base_folder,
    r"processing\arcprojects\MyProject1\sensitivity.gdb"
)

wind_onshore_score_field = [
    'score_wind_onshore_baseline',
    'score_wind_onshore_resource',
    'score_wind_onshore_infra',
    'score_wind_onshore_terrain',
    'score_wind_onshore_equal'
]
# offshore 只使用你配置的 4 个情景字段，不再包含未配置的 terrain 评分。
wind_offshore_score_field = [
    'score_wind_offshore_baseline',
    'score_wind_offshore_resource',
    'score_wind_offshore_infra',
    'score_wind_offshore_equal'
]
solar_score_field = [
    'score_solar_baseline',
    'score_solar_resource',
    'score_solar_infra',
    'score_solar_terrain',
    'score_solar_equal'
]

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
LC_YEAR_COLS = [f"lc_{y}" for y in YEARS]          # low_carbon 情景
ALL_YEAR_COLS = ["kw2025"] + PROV_YEAR_COLS + LC_YEAR_COLS

# 输出要素类中保留的字段（按情景动态拼接）
BASE_OUTPUT_COLS = [NID10_INT, "Shengcode", "KW2", "kw2025"]


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
        df[NID10_INT] = df[NID10_INT].astype(int)
    elif "NID10" in df.columns:
        df[NID10_INT] = df["NID10"].astype(int)
        df = df.drop(columns=["NID10"])
    return df


def load_grid_data(potential_shp, potential_field, current_fc, province_filter=None):
    """
    读取网格数据，合并潜力和现状两个数据源:
    返回含以下列的 DataFrame:
      - NID10_INT, NID50, Shengcode, KW2, kw2025
    """
    import arcpy

    pot_all_fields = [f.name for f in arcpy.ListFields(potential_shp)]
    nid_field_pot = NID10_INT if NID10_INT in pot_all_fields else "NID10"
    pot_fields = [nid_field_pot, "NID50", "Shengcode", potential_field]
    df_pot = read_feature_to_df(potential_shp, pot_fields)
    df_pot = ensure_nid10_int(potential_shp, df_pot)
    df_pot = df_pot.rename(columns={potential_field: "KW2"})
    df_pot["KW2"] = df_pot["KW2"].fillna(0)

    cur_all_fields = [f.name for f in arcpy.ListFields(current_fc)]
    nid_field_cur = NID10_INT if NID10_INT in cur_all_fields else "NID10"
    cur_fields = [nid_field_cur, current_cap_field]
    df_cur = read_feature_to_df(current_fc, cur_fields)
    df_cur = ensure_nid10_int(current_fc, df_cur)
    df_cur = df_cur.rename(columns={current_cap_field: "kw2025"})
    df_cur["kw2025"] = df_cur["kw2025"].fillna(0)

    pot_dup = df_pot[NID10_INT].duplicated().sum()
    cur_dup = df_cur[NID10_INT].duplicated().sum()
    if pot_dup > 0:
        print(f"  ⚠ 潜力数据 NID10_INT 有 {pot_dup} 条重复，按合计处理")
        df_pot = df_pot.groupby([NID10_INT, "NID50", "Shengcode"], as_index=False)["KW2"].sum()
    if cur_dup > 0:
        print(f"  ⚠ 现状数据 NID10_INT 有 {cur_dup} 条重复，按合计处理")
        df_cur = df_cur.groupby(NID10_INT, as_index=False)["kw2025"].sum()

    df = df_pot.merge(df_cur, on=NID10_INT, how="left")
    df["kw2025"] = df["kw2025"].fillna(0)

    if province_filter is not None:
        df = df.loc[province_filter(df["Shengcode"])].copy()

    print(f"  潜力网格数: {len(df_pot)}, 现状网格数: {len(df_cur)}, 合并后: {len(df)}")
    return df


def load_score_data(score_fc, score_field):
    """从评分要素类读取评分，用于 join；若字段不存在则按 0 处理。"""
    import arcpy
    all_fields = [f.name for f in arcpy.ListFields(score_fc)]
    nid_field = NID10_INT if NID10_INT in all_fields else "NID10"

    if score_field not in all_fields:
        print(f"  ⚠ 评分字段 {score_field} 不存在，按 0 处理")
        df = read_feature_to_df(score_fc, [nid_field])
        df = ensure_nid10_int(score_fc, df)
        df["Score"] = 0.0
        return df

    fields = [nid_field, score_field]
    df = read_feature_to_df(score_fc, fields)
    df = ensure_nid10_int(score_fc, df)
    df = df.rename(columns={score_field: "Score"})
    if df[NID10_INT].duplicated().any():
        df = df.groupby(NID10_INT, as_index=False)["Score"].mean()
    return df


# ============================================================
# 3. 核心分配函数（多年份递增）
# ============================================================
def plan_install_multiyear(df, shengcode, kw_targets, year_cols):
    """
    对单个省份，按年份顺序分配装机到网格。
    """
    df_prov = df.loc[df["Shengcode"] == shengcode].copy()

    if df_prov.empty:
        print(f"  ⚠ 省代码 {shengcode} 无对应网格，跳过")
        return df_prov

    df_prov = df_prov.sort_values(
        by=["kw2025", "Score", "KW2"], ascending=False
    ).reset_index(drop=True)

    prev_col = "kw2025"

    for year_col, kw_target_wan in zip(year_cols, kw_targets):
        kw_target = kw_target_wan * 10000
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


def ensure_output_gdb():
    import arcpy
    if not arcpy.Exists(output_gdb):
        arcpy.management.CreateFileGDB(os.path.dirname(output_gdb), os.path.basename(output_gdb))
        print(f"[Created] {output_gdb}")


def prepare_output_feature_class(base_fc, output_fc):
    import arcpy
    arcpy.env.overwriteOutput = True
    if arcpy.Exists(output_fc):
        arcpy.management.Delete(output_fc)
    arcpy.management.CopyFeatures(base_fc, output_fc)
    print(f"[Prepared] {output_fc}")
    return output_fc


def write_results_to_fc(final_df, output_fc, field_cols):
    import arcpy
    arcpy.env.overwriteOutput = True

    if not arcpy.Exists(output_fc):
        raise FileNotFoundError(f"输出要素类不存在: {output_fc}")

    existing_fields = [f.name for f in arcpy.ListFields(output_fc)]
    for col in field_cols:
        if col not in existing_fields:
            if col in [NID10_INT, "Shengcode"]:
                arcpy.management.AddField(output_fc, col, "LONG")
            else:
                arcpy.management.AddField(output_fc, col, "DOUBLE")

    protected_fields = {"OBJECTID", "OID", "FID", "Shape", "SHAPE", "Shape_Length", "Shape_Area"}
    keep_fields = protected_fields | {NID10_INT, "Shengcode", "KW2", "kw2025"} | set(field_cols)

    for field in arcpy.ListFields(output_fc):
        if field.name in keep_fields or field.required:
            continue
        try:
            arcpy.management.DeleteField(output_fc, field.name)
        except Exception as ex:
            print(f"  跳过字段 {field.name}: {ex}")

    lookup = final_df.groupby(NID10_INT)[field_cols].first().to_dict("index")
    update_fields = [NID10_INT] + field_cols

    with arcpy.da.UpdateCursor(output_fc, update_fields) as cursor:
        for row in cursor:
            nid = row[0]
            if nid in lookup:
                values = lookup[nid]
                for i, col in enumerate(field_cols, start=1):
                    row[i] = values.get(col, 0.0)
                cursor.updateRow(row)

    print(f"[Wrote] {output_fc}")


def run_scenario(energy_type, scenario_name, score_field, pot_shp, pot_field, cur_fc, sheet_name, province_filter):
    """按单个评分情景执行一次分配，并输出到 sensitivity.gdb。"""
    print(f"\n{'=' * 60}")
    print(f"开始分配: {energy_type} | 情景: {scenario_name}")
    print(f"评分字段: {score_field}")
    print(f"{'=' * 60}")

    print("读取网格潜力 + 现状数据...")
    grid_df = load_grid_data(pot_shp, pot_field, cur_fc, province_filter=province_filter)

    print("读取评分数据...")
    score_df = load_score_data(score_path, score_field)
    if "Score" in score_df.columns:
        score_df = score_df.rename(columns={"Score": score_field})
    grid_df = grid_df.merge(score_df, on=NID10_INT, how="left")
    if score_field not in grid_df.columns:
        grid_df[score_field] = 0.0
    else:
        grid_df[score_field] = grid_df[score_field].fillna(0)
    grid_df["Score"] = grid_df[score_field]

    print("读取规划表...")
    plan_df = pd.read_excel(planning_path, sheet_name=sheet_name)
    if province_filter is not None:
        plan_df = plan_df.loc[plan_df["Shengcode"].apply(province_filter)].copy()

    if plan_df.empty:
        print("  ⚠ 规划表中没有满足条件的省份，跳过")
        return None

    results = []
    for _, row in plan_df.iterrows():
        shengcode = row["Shengcode"]

        prov_targets = [row[col] for col in PROVINCE_PLAN_COLS]
        print(f"\n省份: {shengcode} | Province目标(万kW): {prov_targets}")
        df_prov = plan_install_multiyear(grid_df, shengcode, prov_targets, PROV_YEAR_COLS)

        lc_targets = [row[col] for col in LOWCARBON_PLAN_COLS]
        print(f"省份: {shengcode} | LowCarbon目标(万kW): {lc_targets}")
        df_lc = plan_install_multiyear(grid_df, shengcode, lc_targets, LC_YEAR_COLS)

        if not df_prov.empty and not df_lc.empty:
            for col in LC_YEAR_COLS:
                df_prov[col] = df_lc[col].values
        results.append(df_prov)

    if not results:
        print("  ⚠ 没有生成任何结果，跳过")
        return None

    final_df = pd.concat(results, ignore_index=True)
    keep_cols = [c for c in BASE_OUTPUT_COLS + [score_field] + PROV_YEAR_COLS + LC_YEAR_COLS if c in final_df.columns]
    final_df = final_df[keep_cols]

    print(f"\n{'=' * 60}")
    print("分配结果汇总（万kW）:")
    print(f"{'=' * 60}")
    summary = final_df.groupby("Shengcode")[ALL_YEAR_COLS].sum() / 10000
    print(summary.to_string())

    total = final_df[ALL_YEAR_COLS].sum() / 10000
    print("\n全国合计(万kW):")
    for col in ALL_YEAR_COLS:
        print(f"  {col}: {total[col]:.1f}")

    output_name = f"{energy_type}_{scenario_name}"
    output_fc = os.path.join(output_gdb, output_name)
    prepare_output_feature_class(score_path, output_fc)
    write_results_to_fc(final_df, output_fc, keep_cols)

    summary_path = os.path.join(
        base_folder,
        rf"processing\tables\{output_name}_allocation_summary.csv"
    )
    summary.to_csv(summary_path, encoding="utf-8-sig")
    print(f"分省汇总已导出: {summary_path}")
    return final_df


# ============================================================
# 4. 批量运行：按情景生成输出要素类
# ============================================================
def main():
    import arcpy
    ensure_output_gdb()
    arcpy.env.overwriteOutput = True

    scenario_configs = [
        ("wind_onshore", wind_onshore_score_field, potential_wind_shp, potential_cap_wind, current_wind_fc, wind_sheetname, lambda x: x != 100),
        ("wind_offshore", wind_offshore_score_field, potential_wind_shp, potential_cap_wind, current_wind_fc, wind_sheetname, lambda x: x == 100),
        ("solar", solar_score_field, potential_solar_shp, potential_cap_solar, current_solar_fc, solar_sheetname, None),
    ]

    for energy_type, score_fields, pot_shp, pot_field, cur_fc, sheet_name, province_filter in scenario_configs:
        for score_field in score_fields:
            scenario_name = score_field.rsplit("_", 1)[-1]
            run_scenario(
                energy_type=energy_type,
                scenario_name=scenario_name,
                score_field=score_field,
                pot_shp=pot_shp,
                pot_field=pot_field,
                cur_fc=cur_fc,
                sheet_name=sheet_name,
                province_filter=province_filter,
            )


if __name__ == "__main__":
    main()
    print("\n✅ 全部分配完成！")