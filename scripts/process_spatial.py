
from shapely.geometry import LineString, LinearRing, Point,Polygon, MultiPolygon, MultiPoint
from shapely.ops import unary_union
import geopandas as gpd
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt
import loopflopy
from loopflopy.mesh_routines import resample_linestring, resample_shapely_poly, resample_gdf_poly
from shapely.affinity import translate

def remove_duplicate_points(polygon):
    # Extract unique points using LinearRing
    linear_ring = LinearRing(polygon.exterior.coords)
    unique_coords = list(linear_ring.coords)
    # Reconstruct the Polygon without duplicates
    return Polygon(unique_coords)

def make_shp(spatial, xcoords, ycoords):
    shp = Polygon(list(zip(xcoords, ycoords)))
    shp_gdf = gpd.GeoDataFrame(geometry=[bbox], crs = spatial.epsg)
    shp_gdf.to_file('../data/data_shp/Gingin_shape_file.shp')

def model_boundary(spatial, boundary_buff, simplify_tolerance, node_spacing):

    model_boundary_gdf = gpd.read_file('../data/data_shp/Gingin_transect_1.shp')
    model_boundary_gdf.to_crs(epsg=28350, inplace=True)
    model_boundary_gs = model_boundary_gdf.geometry.simplify(tolerance=simplify_tolerance, preserve_topology=True)
    model_boundary_poly = resample_gdf_poly(model_boundary_gs, node_spacing)

    # Create an inner boundary for meshing
    inner_boundary_poly = model_boundary_poly.buffer(-boundary_buff)
    inner_boundary_gs = gpd.GeoSeries([inner_boundary_poly])
    inner_boundary_poly = resample_gdf_poly(inner_boundary_gs, node_spacing) #  A few less nodes in inside boundary
    
    
    spatial.model_boundary_gdf = model_boundary_gdf
    model_boundary_gdf.to_file('../modelfiles/model_boundary.shp')
    spatial.model_boundary_poly = model_boundary_poly
    spatial.inner_boundary_poly = inner_boundary_poly
    spatial.x0, spatial.y0, spatial.x1, spatial.y1 = model_boundary_poly.bounds

def head_boundary(spatial):    

    # CHD EAST BOUNDARY
    inner_boundary_poly = spatial.model_boundary_poly.buffer(-1)
    coords = list(inner_boundary_poly.exterior.coords)
    new_coords = []
    for coord in coords[:-1]: # Don't include last point
        x,y = coord[0], coord[1]
        if x > spatial.x1 - 100:
            new_coords.append(coord)
    
    chd_east_ls = LineString(new_coords)
    chd_east_gdf = gpd.GeoDataFrame({'geometry': [chd_east_ls]}, crs = spatial.epsg)
    spatial.chd_east_gdf = chd_east_gdf
    spatial.chd_east_ls = chd_east_ls

    # CHD WEST BOUNDARY
    inner_boundary_poly = spatial.model_boundary_poly.buffer(-1)
    coords = list(inner_boundary_poly.exterior.coords)
    new_coords = []
    for coord in coords[:-1]: # Don't include last point
        x,y = coord[0], coord[1]
        if x < 370500:
            if y < (spatial.y1 - 5):
                if y > (spatial.y0 + 5):
                    new_coords.append(coord)
    
    chd_west_ls = LineString(new_coords)
    chd_west_gdf = gpd.GeoDataFrame({'geometry': [chd_west_ls]}, crs = spatial.epsg)
    spatial.chd_west_gdf = chd_west_gdf
    spatial.chd_west_ls = chd_west_ls

def geo_bores(spatial):   
    df = pd.read_excel('../data/data_geology/data_geology.xlsx', sheet_name = 'geo_bores')
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Easting, df.Northing), crs=spatial.epsg)
    gdf = gpd.clip(gdf, spatial.model_boundary_poly).reset_index(drop=True)
    spatial.geobore_gdf = gdf
    spatial.idgeobores = list(gdf.ID)
    spatial.xygeobores = list(zip(gdf.Easting, gdf.Northing))
    spatial.nobs = len(spatial.xygeobores)

def obs_bores(spatial):   
    df = pd.read_excel('../data/data_geology/data_geology.xlsx', sheet_name = 'obs_bores')
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Easting, df.Northing), crs=spatial.epsg)
    gdf = gpd.clip(gdf, spatial.model_boundary_poly).reset_index(drop=True)
    spatial.obsbore_gdf = gdf
    spatial.idobsbores = list(gdf.ID)
    spatial.xyobsbores = list(zip(gdf.Easting, gdf.Northing))
    spatial.nobs = len(spatial.xyobsbores)
    
def pump_bores(spatial):    
    df = pd.read_excel('../data/data_geology/data_geology.xlsx', sheet_name = 'pumping_bores')
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Easting, df.Northing), crs=spatial.epsg)
    gdf = gpd.clip(gdf, spatial.model_boundary_poly).reset_index(drop=True)
    spatial.pumpbore_gdf = gdf
    spatial.idpumpbores = list(gdf.ID)
    spatial.xypumpbores = list(zip(gdf.Easting, gdf.Northing))
    spatial.npump = len(spatial.xypumpbores)

def faults(spatial):  
    # Import fault shape file
    faults_gdf = gpd.read_file('../data/data_shp/Gingin_transect_1_fault.shp')
    print(faults_gdf)
    faults_gdf.to_crs(epsg=28350, inplace=True)
    #faults_gdf = gpd.clip(faults_gdf, spatial.model_boundary_poly).reset_index(drop=True)

    #get nodes of fault line
    fault_nodes = []
    for i in range(len(faults_gdf.geometry[0].coords)):
        fault_nodes.append(faults_gdf.geometry[0].coords[i])       

    r = 200 # distance between points
    ls = faults_gdf.geometry[0]    
    ls_resample = resample_linestring(ls, r) # Resample linestring
    '''
    #Removing nodes too close to inner and outer boundary so mesh doesn't go crazy refined (threshold_distance)
    threshold_distance = 1000
    nodes_to_remove = []
    for p1 in ls_resample:
        for p2 in spatial.inner_boundary_poly.exterior.coords:
            p2 = Point(p2)
            if p1.distance(p2) <= threshold_distance:
                nodes_to_remove.append(p1)
        for p3 in spatial.model_boundary_poly.exterior.coords:
            p3 = Point(p3)
            if p1.distance(p3) <= threshold_distance:
                nodes_to_remove.append(p1)

    ls_new = [node for node in ls_resample if node not in nodes_to_remove]'''

    fault_nodes = [(point.x, point.y) for point in ls_resample]  # Convert to list of tuples 
    fault_ls = LineString(fault_nodes)
    faults_gdf = gpd.GeoDataFrame(pd.DataFrame({'geometry': [fault_ls]}), crs=spatial.epsg)
    
    spatial.fault_ls = fault_ls
    spatial.faults_gdf = faults_gdf   
    spatial.fault_nodes = fault_nodes
   
def lakes(spatial):  

    gdf = gpd.read_file('../data/data_shp/Gingin_transect_1_surface water.shp')
    gdf.to_crs(epsg=28350, inplace=True)
    gdf = gpd.clip(gdf, spatial.model_boundary_poly).reset_index(drop=True)
    spatial.lakes_gdf = gdf

def plot_spatial(spatial, extent = None):    # extent[[x0,x1], [y0,y1]]
    
    fig, ax = plt.subplots(figsize = (7,7))
    ax.set_title('Gingin transect 1')
       
    x, y = spatial.model_boundary_poly.exterior.xy
    ax.plot(x, y, '-o', ms = 2, lw = 1, color='black')
    x, y = spatial.inner_boundary_poly.exterior.xy
    ax.plot(x, y, '-o', ms = 2, lw = 0.5, color='black')
    if extent: 
        ax.set_xlim(extent[0][0], extent[0][1])
        ax.set_ylim(extent[1][0], extent[1][1])
        
    spatial.faults_gdf.plot(ax=ax, markersize = 5, color = 'lightblue', zorder=2)
    for node in spatial.fault_nodes: 
        ax.plot(node[0], node[1], 'o', ms = 3, color = 'lightblue', zorder=2)
    
    spatial.lakes_gdf.plot(ax=ax, color = 'darkblue', zorder=2)
    #spatial.chd_east_gdf.plot(ax=ax, markersize = 12, color = 'red', zorder=2)
    #spatial.chd_west_gdf.plot(ax=ax, markersize = 12, color = 'red', zorder=2)
    spatial.obsbore_gdf.plot(ax=ax, markersize = 5, color = 'black', zorder=2)
    spatial.pumpbore_gdf.plot(ax=ax, markersize = 12, color = 'red', zorder=2)
    spatial.geobore_gdf.plot(ax=ax, markersize = 12, color = 'green', zorder=2)

    for x, y, label in zip(spatial.geobore_gdf.geometry.x, spatial.geobore_gdf.geometry.y, spatial.geobore_gdf.ID):
        ax.annotate(label, xy=(x, y), xytext=(2, 2), size = 7, textcoords="offset points")
    for x, y, label in zip(spatial.obsbore_gdf.geometry.x, spatial.obsbore_gdf.geometry.y, spatial.obsbore_gdf.ID):
        ax.annotate(label, xy=(x, y), xytext=(2, 2), size = 7, textcoords="offset points")
    for x, y, label in zip(spatial.pumpbore_gdf.geometry.x, spatial.pumpbore_gdf.geometry.y, spatial.pumpbore_gdf.ID):
        ax.annotate(label, xy=(x, y), xytext=(2, 2), size = 10, textcoords="offset points")

    plt.savefig('../figures/spatial.png')
    