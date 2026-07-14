"""
Solar Panel Polygon QA Basemap Generator
=========================================

For selected PV (photovoltaic) polygons (Select == 1), generate a QA image per
province with a Tianditu satellite basemap overlaid, so the polygon boundaries
can be manually spot-checked:
  - If a province has fewer than 2 selected polygons: produce a single overview map.
  - If a province has >= 2 selected polygons: additionally sample 2 random polygons
    and produce zoomed-in inset maps, combined with the overview map into one image.

By default only the province given in `TARGET_PROVINCE` is processed (matching the
original code's behavior); to batch-process every province still pending review,
see the commented-out block at the bottom of the file, which loops over `left_lst`.

Ported from: mappingSolarPanel.py (restructured for readability; the original
calculation logic and results are unchanged).
"""

import os
from glob import glob

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import geopandas as gpd
import matplotlib.pyplot as plt

# ============================================================
# 1. Paths and parameters
# ============================================================
PV_POLYGON_SHP = r"C:\Users\ruome\Desktop\Renewable\论文编辑\数据输出\结果一\PV_polygonMerge.shp"
OUTPUT_DIR = r"C:\Users\ruome\Desktop\Renewable\ProcessingFiles\images\download\solar"

# Province processed in a single run (matches the original: manual, one province at a time,
# rather than batch-processing every province).
TARGET_PROVINCE = "吉林省"

# Tianditu tile zoom level
MAIN_MAP_DEGREE = 18
SUBFIG_DEGREE = 18


# ============================================================
# 2. Tianditu WMTS basemap
# ============================================================
class TDT_img(cimgt.GoogleWTS):
    """Tianditu satellite imagery WMTS tile source."""

    _TDT_KEYS = [
        "82497ed7105826b2c770bc58ccbf7a02",
        "edf5ed862d6dca477ce3293f0d621827",
        "e44b6797f396ec04c056bca92a9379c1",
        "ec0eb06aec43d59ccb1ce1221bcc75ee",
        "ead6e85c4ec34c1a671e1a403075ea56",
        "ec0eb06aec43d59ccb1ce1221bcc75ee",
        "5d879209289e61417b0a06e98abfbe19",
        "a820f85ff87edf1e824056479e196ecd",
    ]

    def _image_url(self, tile):
        x, y, z = tile
        return "http://t3.tianditu.gov.cn/DataServer?T=img_w&x=%s&y=%s&l=%s&tk=%s" % (
            x, y, z, self._TDT_KEYS[1],
        )


def make_map(fig, projection=ccrs.PlateCarree()):
    return fig.add_subplot(1, 1, 1, projection=projection)


def make_submap(subfig, projection=ccrs.PlateCarree()):
    ax12 = subfig.add_subplot(2, 1, 1, projection=projection)
    ax13 = subfig.add_subplot(2, 1, 2, projection=projection)
    return ax12, ax13


# ============================================================
# 3. Data preparation
# ============================================================
def load_selected_polygons(province):
    """Load the selected (Select==1), de-duplicated PV polygons for one province."""
    select_gdf = gpd.read_file(PV_POLYGON_SHP)
    select_gdf = select_gdf.loc[select_gdf["Select"] == 1]

    gdf = select_gdf.loc[select_gdf["name"] == province]
    gdf = gdf.drop_duplicates(subset=["geometry"])
    return select_gdf, gdf


def list_pending_provinces(select_gdf):
    """Compare the selected-province list against provinces that already have a QA image;
    return the provinces still pending (informational only)."""
    province_lst = select_gdf["name"].drop_duplicates().to_list()
    exit_lst = glob(os.path.join(OUTPUT_DIR, "*.png"))
    province_exit_lst = [os.path.split(i)[-1].split(".")[0] for i in exit_lst]
    return [i for i in province_lst if i not in province_exit_lst]


def compute_extent(gdf, pad_ratio_denominator=5):
    """Compute a square display extent around the polygons' bounding box, padded by
    1 / `pad_ratio_denominator` of the half-width on each side."""
    lon_min, lat_min, lon_max, lat_max = gdf.total_bounds
    d = max(lon_max - lon_min, lat_max - lat_min) / 2
    center_x = lon_min + d
    center_y = lat_min + d
    pad = d / pad_ratio_denominator
    extent = [center_x - d - pad, center_x + d + pad, center_y - d - pad, center_y + d + pad]
    return extent, d


# ============================================================
# 4. Plotting
# ============================================================
def plot_single_polygon_overview(gdf, province, extent):
    """Fewer than 2 polygons: output a single overview map only."""
    print(province, "less than 2")
    request = TDT_img()
    fig = plt.figure(layout="constrained", figsize=(10, 8))
    ax = make_map(fig, projection=request.crs)
    ax.set_extent(extent)
    ax.add_image(request, MAIN_MAP_DEGREE)
    gdf.plot(ax=ax, facecolor="none", edgecolor="red", transform=ccrs.Geodetic())
    plt.savefig(os.path.join(OUTPUT_DIR, f"{province}.png"))


def plot_sampled_polygons_overview(gdf, province, extent, d):
    """2 or more polygons: overview map + zoomed-in insets for 2 randomly sampled polygons."""
    sample_gdf = gdf.sample(n=2).reset_index()
    sample1_gdf = gpd.GeoDataFrame(sample_gdf.loc[0].to_frame().T, crs=gdf.crs)
    sample2_gdf = gpd.GeoDataFrame(sample_gdf.loc[1].to_frame().T, crs=gdf.crs)

    d1 = d / 2

    lon1_min, lat1_min, lon1_max, lat1_max = sample1_gdf.total_bounds
    center1_x = lon1_min + ((lon1_max - lon1_min) / 2)
    center1_y = lat1_min + ((lat1_max - lat1_min) / 2)
    sample1_extent = [center1_x - d1, center1_x + d1, center1_y - (d1 / 2), center1_y + (d1 / 2)]

    lon2_min, lat2_min, lon2_max, lat2_max = sample2_gdf.total_bounds
    center2_x = lon2_min + ((lon2_max - lon2_min) / 2)
    center2_y = lat2_min + ((lat2_max - lat2_min) / 2)
    sample2_extent = [center2_x - d1, center2_x + d1, center2_y - (d1 / 2), center2_y + (d1 / 2)]

    request = TDT_img()
    fig = plt.figure(figsize=(10, 10), layout="tight")
    subfigs = fig.subfigures(1, 2, width_ratios=[1, 1])
    subfig1, subfig2 = subfigs[0], subfigs[1]

    # Overview map: all polygons + the 2 sampled polygons highlighted
    ax11 = make_map(subfig1, projection=ccrs.PlateCarree())
    ax11.set_extent(extent)
    ax11.add_image(request, MAIN_MAP_DEGREE)
    gdf.plot(ax=ax11, facecolor="none", edgecolor="red", transform=ccrs.Geodetic())
    sample1_gdf.plot(ax=ax11, facecolor="none", edgecolor="blue", transform=ccrs.Geodetic())
    sample2_gdf.plot(ax=ax11, facecolor="none", edgecolor="yellow", transform=ccrs.Geodetic())

    # Zoomed-in insets: sample 1, sample 2
    ax12, ax13 = make_submap(subfig2, projection=ccrs.PlateCarree())
    ax12.set_extent(sample1_extent)
    ax12.add_image(request, SUBFIG_DEGREE)
    sample1_gdf.plot(ax=ax12, facecolor="none", edgecolor="blue", transform=ccrs.Geodetic())

    ax13.set_extent(sample2_extent)
    ax13.add_image(request, SUBFIG_DEGREE)
    sample2_gdf.plot(ax=ax13, facecolor="none", edgecolor="yellow", transform=ccrs.Geodetic())

    plt.subplots_adjust(left=0.05, right=2, top=0.95, bottom=0.05)
    fig.set_figwidth(10)
    fig.set_figheight(5)

    plt.savefig(os.path.join(OUTPUT_DIR, f"{province}.png"))
    print("complete:", province)


def generate_qa_map(province):
    select_gdf, gdf = load_selected_polygons(province)

    left_lst = list_pending_provinces(select_gdf)
    print(left_lst)

    print("Start Province:", province)
    extent, d = compute_extent(gdf)

    if len(gdf) < 2:
        plot_single_polygon_overview(gdf, province, extent)
    else:
        plot_sampled_polygons_overview(gdf, province, extent, d)


# ============================================================
# 5. Main entry point
# ============================================================
if __name__ == "__main__":
    generate_qa_map(TARGET_PROVINCE)

    # To batch-process every province that still lacks a QA image, use instead:
    # select_gdf, _ = load_selected_polygons(TARGET_PROVINCE)
    # for province in list_pending_provinces(select_gdf):
    #     generate_qa_map(province)
