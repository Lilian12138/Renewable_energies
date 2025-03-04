import matplotlib.pyplot as plt
import numpy as np
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import geopandas as gpd
from glob import glob
import os

class TDT_img(cimgt.GoogleWTS):
    def _image_url(self, tile):
        x, y, z = tile
        keys = [
            "82497ed7105826b2c770bc58ccbf7a02",
            "edf5ed862d6dca477ce3293f0d621827",
            "e44b6797f396ec04c056bca92a9379c1",
            "ec0eb06aec43d59ccb1ce1221bcc75ee",
            "ead6e85c4ec34c1a671e1a403075ea56",
            "ec0eb06aec43d59ccb1ce1221bcc75ee",
            "5d879209289e61417b0a06e98abfbe19",
            "a820f85ff87edf1e824056479e196ecd"
        ]
        url = 'http://t3.tianditu.gov.cn/DataServer?T=img_w&x=%s&y=%s&l=%s&tk=%s' % (x, y, z, keys[1])
        # "https://t0.tianditu.gov.cn/img_w/wmts?SERVICE=WMTS&REQUEST=GetTile&VERSION=1.0.0&LAYER=img&STYLE=default&TILEMATRIXSET=w&FORMAT=tiles&TILEMATRIX=1&TILEROW=1&TILECOL=1&tk=9fe4dd1ce46bf96194ee075b13474f76&level=2"
        return url

def make_map(fig, projection=ccrs.PlateCarree()):
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    return ax

def make_submap(subfig2, projection=ccrs.PlateCarree()):
    ax12 = subfig2.add_subplot(2, 1, 1, projection=projection)  
    ax13 = subfig2.add_subplot(2, 1, 2, projection=projection)  
    return ax12, ax13

path = r"C:\Users\ruome\Desktop\Renewable\论文编辑\数据输出\结果一\PV_polygonMerge.shp"
select_gdf = gpd.read_file(path)
select_gdf = select_gdf.loc[select_gdf['Select']==1]
province_lst = select_gdf['name'].drop_duplicates().to_list()
exit_lst = glob(r"C:\Users\ruome\Desktop\Renewable\ProcessingFiles\images\download\solar\*.png")
province_exit_lst = [os.path.split(i)[-1].split('.')[0] for i in exit_lst]
left_lst = [i for i in province_lst if i not in province_exit_lst]
print(left_lst)
province = "吉林省"
# for province in province_lst:
print("Start Province:", province)
gdf = select_gdf.loc[select_gdf['name']==province]
# drop duplicates
gdf = gdf.drop_duplicates(subset=['geometry'])

degree = 18
subfig_degree = 18
lon_min, lat_min, lon_max, lat_max = gdf.total_bounds
d = max(lon_max-lon_min, lat_max-lat_min)/2
center_x = lon_min + d
center_y = lat_min + d
extent = [center_x-d-(d/5), center_x+d+(d/5), center_y-d-(d/5), center_y+d+(d/5)]

if len(gdf) < 2:
    print(province, "less than 2")
    # longitude, latitude
    request = TDT_img()
    # Create a main figure
    fig = plt.figure(layout='constrained', figsize=(10, 8))
    ax = make_map(fig, projection=request.crs)
    ax.set_extent(extent)
    ax.add_image(request,18)
    gdf.plot(ax = ax, facecolor='none', edgecolor='red', transform=ccrs.Geodetic())
    # plt.show()
    plt.savefig(r'C:\Users\ruome\Desktop\Renewable\ProcessingFiles\images\download\solar\{}.png'.format(province))
else:
    sample_gdf = gdf.sample(n=2).reset_index()
    sample1_gdf = gpd.GeoDataFrame(sample_gdf.loc[0].to_frame().T, crs=gdf.crs)
    sample2_gdf = gpd.GeoDataFrame(sample_gdf.loc[1].to_frame().T, crs=gdf.crs)

    lon1_min, lat1_min, lon1_max, lat1_max = sample1_gdf.total_bounds
    d1 = d/2
    center1_x = lon1_min + ((lon1_max-lon1_min)/2)
    center1_y = lat1_min + ((lat1_max-lat1_min)/2)
    sample1_extent = [center1_x-d1, center1_x+d1, center1_y-(d1/2), center1_y+(d1/2)]

    lon2_min, lat2_min, lon2_max, lat2_max = sample2_gdf.total_bounds
    d2 = d/2
    center2_x = lon2_min + ((lon2_max-lon2_min)/2)
    center2_y = lat2_min + ((lat2_max-lat2_min)/2)
    sample2_extent =  [center2_x-d1, center2_x+d1, center2_y-(d1/2), center2_y+(d1/2)]

    # Create a main figure
    fig = plt.figure(figsize=(10,10), layout="tight")
    subfigs = fig.subfigures(1, 2, width_ratios=[1,1])
    # Access individual subfigures
    subfig1 = subfigs[0]
    subfig2 = subfigs[1]
    # First map in subfigure 1
    ax11 = make_map(subfig1, projection=ccrs.PlateCarree())
    request = TDT_img()
    ax11.set_extent(extent)
    ax11.add_image(request, degree)
    gdf.plot(ax = ax11, facecolor='none', edgecolor='red', transform=ccrs.Geodetic())
    sample1_gdf.plot(ax = ax11, facecolor='none', edgecolor='blue',transform=ccrs.Geodetic())
    sample2_gdf.plot(ax = ax11, facecolor='none', edgecolor='yellow',transform=ccrs.Geodetic())

    # Second map in subfigure 2
    ax12, ax13 = make_submap(subfig2, projection=ccrs.PlateCarree())
    ax12.set_extent(sample1_extent)
    ax12.add_image(request, subfig_degree)
    sample1_gdf.plot(ax = ax12, facecolor='none', edgecolor='blue', transform=ccrs.Geodetic())

    ax13.set_extent(sample2_extent)
    ax13.add_image(request, subfig_degree)
    sample2_gdf.plot(ax = ax13, facecolor='none', edgecolor='yellow', transform=ccrs.Geodetic())

    # Adjust the layout
    plt.subplots_adjust(left=0.05, right=2, top=0.95, bottom=0.05)

    # Set figure window size and position (optional, works in interactive environments)
    fig.set_figwidth(10)
    fig.set_figheight(5)

    # plt.show()
    plt.savefig(r'C:\Users\ruome\Desktop\Renewable\ProcessingFiles\images\download\solar\{}.png'.format(province))
    print("complete:", province)