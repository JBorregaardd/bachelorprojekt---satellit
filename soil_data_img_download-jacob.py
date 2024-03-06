#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 23:51:15 2023

@author: mabso
"""

# read all below

import geopandas as gpd
import pandas as pd
import numpy as np
from sentinelhub import SHConfig
from sentinelhub.geometry import BBox, Geometry
from shapely.geometry import Polygon, mapping, MultiPolygon
import json
import shutil
import glob
import datetime
import os

import matplotlib.pyplot as plt
from sentinelhub import CRS, BBox, DataCollection, MimeType, WcsRequest, WmsRequest

from sentinelhub import (
    CRS,
    BBox,
    DataCollection,
    DownloadRequest,
    MimeType,
    MosaickingOrder,
    SentinelHubDownloadClient,
    SentinelHubRequest,
    bbox_to_dimensions,
)

#%% set up using your sentinelhub account

# Set your API credentials
CLIENT_ID = "ca2fd726-2237-4a0e-9ee8-8dbcb0caf3ef"
CLIENT_SECRET = "0vDMaibWNiaWQ9VXaGySAXEszqvOkgP6"

# Configure Sentinel Hub
config = SHConfig()
config.sh_client_id = CLIENT_ID
config.sh_client_secret = CLIENT_SECRET

#%% read in your data

# read in shape files (id number and geometrics shape - lon and lat coordinate pairs)
file_path="/home/mabso/Desktop/lucas_labels/LUCAS_2018_Copernicus/LUCAS_2018_Copernicus_polygons.shp"
df_shape=gpd.read_file(file_path)

# read in the soil file (id number and lots of soil properties per id number)
file_path="/home/mabso/Desktop/soil data analysis/LUCAS-SOIL-2018-data-report-readme-v2/LUCAS-SOIL-2018-v2/LUCAS-SOIL-2018.csv"
df_soil=pd.read_csv(file_path)


#%% set your evalscripts 

# you get these from sentinelhub just copy and paste

# this one is NDWI (a water index)
evalscript_ndwi = """
//VERSION=3
//This script was converted from v1 to v3 using the converter API

//ndwi
var colorRamp1 = [
  	[0, 0xFFFFFF],
  	[1, 0x008000]
  ];
var colorRamp2 = [
  	[0, 0xFFFFFF],
  	[1, 0x0000CC]
  ];

let viz1 = new ColorRampVisualizer(colorRamp1);
let viz2 = new ColorRampVisualizer(colorRamp2);

function evaluatePixel(samples) {
  var val = index(samples.B03, samples.B08);

  if (val < -0) {
    return viz1.process(-val);
  } else {
    return viz2.process(Math.sqrt(Math.sqrt(val)));
  }
}

function setup() {
  return {
    input: [{
      bands: [
        "B03",
        "B08"
      ]
    }],
    output: {
      bands: 3
    }
  }
}
"""

# this is a float ndvi script - you can also compute these manually from raw data
evalscript_all_nvdi = """
   //VERSION=3
    function setup() {
      return{
        input: [{
          bands: ["B04", "B08"]
        }],
        output: {
          id: "default",
          bands: 1,
          sampleType: SampleType.FLOAT32
        }
      }
    }
    function evaluatePixel(sample) {
      let ndvi = (sample.B08 - sample.B04) / (sample.B08 + sample.B04)
      return [ ndvi ]
    }
"""

# this is the screen classifcation map
evalscript_complex="""
//VERSION=3

 function RGBToColor (r, g, b,dataMask){
	return [r/255, g/255, b/255,dataMask];
}

function setup() {
   return {
    input: ["SCL","dataMask"],
    output: { bands: 4 }
  };
}

function evaluatePixel(samples) {
    const SCL=samples.SCL;
    switch (SCL) {
    // No Data (Missing data) - black   
    case 0: return RGBToColor (0, 0, 0,samples.dataMask);
        
    // Saturated or defective pixel - red 
    case 1: return RGBToColor (255, 0, 0,samples.dataMask);

    // Topographic casted shadows ("Dark features/Shadows" for data before 2022-01-25) - very dark grey
    case 2: return RGBToColor (47,  47,  47,samples.dataMask);
        
    // Cloud shadows - dark brown
    case 3: return RGBToColor (100, 50, 0,samples.dataMask);
        
    // Vegetation - green
    case 4: return RGBToColor (0, 160, 0,samples.dataMask);
        
    // Not-vegetated - dark yellow
    case 5: return RGBToColor (255, 230, 90,samples.dataMask);
        
    // Water (dark and bright) - blue
    case 6: return RGBToColor (0, 0, 255,samples.dataMask);
    
    // Unclassified - dark grey
    case 7: return RGBToColor (128, 128, 128,samples.dataMask);
    
    // Cloud medium probability - grey
    case 8: return RGBToColor (192, 192, 192,samples.dataMask);
        
    // Cloud high probability - white
    case 9: return RGBToColor (255, 255, 255,samples.dataMask);
    
    // Thin cirrus - very bright blue
    case 10: return RGBToColor (100, 200, 255,samples.dataMask);
        
    // Snow or ice - very bright pink
    case 11: return RGBToColor (255, 150, 255,samples.dataMask);

    default : return RGBToColor (0, 0, 0,samples.dataMask);  
    }
}
"""

# this is an integer version of the NDVI (less precison)
evalscript_nvdi = """
   //VERSION=3

function setup() {
   return {
       input: ["B04", "B08", "dataMask"],
       output: {bands: 4}
     };
}



function evaluatePixel(sample) {
   let ndvi = index(sample.B08, sample.B04);
   let imgVals = null;
 
   if (ndvi < -0.5) imgVals = [0.05, 0.05, 0.05];
   else if (ndvi < -0.2) imgVals = [0.75, 0.75, 0.75];
   else if (ndvi < -0.1) imgVals = [0.86, 0.86, 0.86];
   else if (ndvi < 0) imgVals = [0.92, 0.92, 0.92];
   else if (ndvi < 0.025) imgVals = [1, 0.98, 0.8];
   else if (ndvi < 0.05) imgVals = [0.93, 0.91, 0.71];
   else if (ndvi < 0.075) imgVals = [0.87, 0.85, 0.61];
   else if (ndvi < 0.1) imgVals = [0.8, 0.78, 0.51];
   else if (ndvi < 0.125) imgVals = [0.74, 0.72, 0.42];
   else if (ndvi < 0.15) imgVals = [0.69, 0.76, 0.38];
   else if (ndvi < 0.175) imgVals = [0.64, 0.8, 0.35];
   else if (ndvi < 0.2) imgVals = [0.57, 0.75, 0.32];
   else if (ndvi < 0.25) imgVals = [0.5, 0.7, 0.28];
   else if (ndvi < 0.3) imgVals = [0.44, 0.64, 0.25];
   else if (ndvi < 0.35) imgVals = [0.38, 0.59, 0.21];
   else if (ndvi < 0.4) imgVals = [0.31, 0.54, 0.18];
   else if (ndvi < 0.45) imgVals = [0.25, 0.49, 0.14];
   else if (ndvi < 0.5) imgVals = [0.19, 0.43, 0.11];
   else if (ndvi < 0.55) imgVals = [0.13, 0.38, 0.07];
   else if (ndvi < 0.6) imgVals = [0.06, 0.33, 0.04];
   else imgVals = [0, 0.27, 0];  
   
   imgVals.push(sample.dataMask)

   return imgVals
}
"""
# sampling month script
def samp_month_fun(mon_num):
    
    if mon_num=='01':
        out="2018-01-01", "2018-01-31"
    
    elif mon_num=='02':
        out="2018-02-01", "2018-02-28"
        
    elif mon_num=='03':
         out="2018-03-01", "2018-03-31"
    
    elif mon_num=='04':
        out="2018-04-01", "2018-04-30"
    
    elif mon_num=='05':
        out="2018-05-01", "2018-05-31"
    
    elif mon_num=='06':
        out="2018-06-01", "2018-06-30"
    
    elif mon_num=='07':
        out="2018-07-01", "2018-07-31"
    
    elif mon_num=='08':
        out="2018-08-01", "2018-08-31"
    
    elif mon_num=='09':
        out="2018-09-01", "2018-09-30"
    
    elif mon_num=='10':
        out="2018-10-01", "2018-10-31"
    
    elif mon_num=='11':
        out="2018-11-01", "2018-11-30"
        
    elif mon_num=='12':
        out="2018-12-01", "2018-12-31"
    
    return out


ndvi_evalscript = """
//VERSION=3

function setup() {
  return {
    input: [
      {
        bands: [
          "B04",
          "B08",
          "dataMask"
        ]
      }
    ],
    output: [
      {
        id: "ndvi",
        bands: 1
      },
      {
        id: "dataMask",
        bands: 1
      }
    ]
  }
}

function evaluatePixel(samples) {
    return {
      ndvi: [index(samples.B08, samples.B04)],
      dataMask: [samples.dataMask]
    };
}
"""
#%% This is the downloading part of the script


# convert shape file data to int
df_shape["POINT_ID"]=df_shape["POINT_ID"].astype(int)

# it might be a good idea to choose only the ids that are associated with crops instead of downloading everything like i have
# so match these two: 
# idx=np.where(df_soil["LC0_Desc"]=="Cropland")[0] # index of crops
# crop_ids=np.array(df_soil.iloc[idx]["POINTID"]) # point ids which matches the df_shape


# set i from 0 to how ever many shapes files there are (you might need to do multiple downloads sessions if your account runs out)
for i in range(18725,19000):
    print(i)
    try:
        # get soild id
        soil_id=df_soil["POINTID"][i]
        samp_date=df_soil["SURVEY_DATE"][i]
        
        # extract polygon and match id name
        col_idx=np.where(df_shape["POINT_ID"]==soil_id)[0]
        id_name=df_shape.iloc[col_idx]['POINT_ID']
        
        # Check if col_idx is empty (length is 0)
        if len(col_idx) == 0:
        # Skip this iteration and move to the next
            continue
        
        polygon_str = df_shape['geometry'][col_idx]
        
        coords_list = []
        
        for polygon in df_shape['geometry'][col_idx]:
            # Extract the exterior coordinates of the polygon
            exterior_coords = list(polygon.exterior.coords)
        
            # Append the exterior coordinates to the list
            coords_list.append(exterior_coords)
        
        polygon = Polygon(coords_list[0])
        
        # Convert the Polygon into a MultiPolygon with a single polygon
        multi_polygon = MultiPolygon([polygon])
        
        # Convert the MultiPolygon to GeoJSON format
        geojson_data = json.dumps(mapping(multi_polygon))
        
        multi_polygon_geometry = json.loads(geojson_data) 
        
        import datetime
        
        dates=samp_month_fun(samp_date[3]+samp_date[4])
        end=dates[1]
        start=dates[0]
        
        end = end.split('-')
        year = int(end[0])
        month = int(end[1])
        day = int(end[2])
        
        end=datetime.datetime(year, month, day)
        
        start = start.split('-')
        year = int(start[0])
        month = int(start[1])
        day = int(start[2])
        
        start=datetime.datetime(year, month, day)
        
        n_chunks=30
        tdelta = (end - start) / n_chunks
        edges = [(start + i * tdelta).date().isoformat() for i in range(n_chunks)]
        #slots = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
        slots = [(edges[i], edges[i]) for i in range(len(edges) - 1)]
        print(slots)
    
        # first evalscript to download from
        def get_NDWI_request(time_interval):
            return SentinelHubRequest(
                evalscript=evalscript_all_nvdi, # set your evalscript to what you want
                data_folder="/media/mabso/Data/downloads/soil_sen2_downloads", # change to where you want your data output
                input_data=[
                    SentinelHubRequest.input_data(
                        data_collection=DataCollection.SENTINEL2_L2A,
                        time_interval=time_interval,
                        mosaicking_order=MosaickingOrder.LEAST_CC,
                        other_args={'processing': {'upsampling': 'BILINEAR', 'downsampling':'BILINEAR'}},
                        maxcc=0.25 # maxium cloud cover
                    )
                ],
                responses=[SentinelHubRequest.output_response("default", MimeType.TIFF)],
                bbox=None,
                geometry=Geometry(multi_polygon_geometry, CRS.WGS84),
                size=(256, 256 ),  # Adjust the size of image as needed
                config=config,
            )
        # second evalscript to download from (optional)
        def get_complex_request(time_interval):
            return SentinelHubRequest(
                evalscript=evalscript_complex, # set your second evalscript to what you want
                data_folder="/media/mabso/Data/downloads/soil_sen2_downloads", # change to where you want your data output
                input_data=[
                    SentinelHubRequest.input_data(
                        data_collection=DataCollection.SENTINEL2_L2A,
                        time_interval=time_interval,
                        mosaicking_order=MosaickingOrder.LEAST_CC,
                        other_args={'processing': {'upsampling': 'BILINEAR', 'downsampling':'BILINEAR'}},
                        maxcc=0.25
                    )
                ],
                responses=[SentinelHubRequest.output_response("default", MimeType.TIFF)],
                bbox=None,
                geometry=Geometry(multi_polygon_geometry, CRS.WGS84),
                size=(256, 256),  # Adjust the size as needed
                config=config,
            )
        
        
        # Create a list of requests
        list_of_requests = [
            get_NDWI_request(slot) for slot in slots
        ] + [
            get_complex_request(slot) for slot in slots  # remove one of these if only wanting to download one thing at at time - usally a good idea to download multiple to speed up
        ]
        
        list_of_requests = [request.download_list[0] for request in list_of_requests]
        
        # Download data with multiple threads
        data = SentinelHubDownloadClient(config=config).download(list_of_requests, max_threads=5)
        
        # note if you only choose to download one at a time you may need to do some major editing down below
        
        # remove all files with blanks
        import numpy as np
        import tifffile as tiff
        
        def are_all_pixels_zero(image):
            return np.all(image == 0)
        
        def process_folder(folder_path):
            subfolders_to_delete = []
            for root, subfolders, files in os.walk(folder_path):
                images_have_only_zeros = True
                for file in files:
                    if file.lower().endswith('.tiff'):
                        image_path = os.path.join(root, file)
                        image = tiff.imread(image_path)
                        if not are_all_pixels_zero(image):
                            images_have_only_zeros = False
                            break
                
                if images_have_only_zeros:
                    subfolders_to_delete.append(root)
        
            for subfolder in subfolders_to_delete:
                print(f"Deleting subfolder: {subfolder}")
                # Uncomment the line below to actually delete the subfolder and its contents
                # shutil.rmtree(subfolder)
                
            return subfolders_to_delete
        if __name__ == "__main__":
            root_folder = "/media/mabso/Data/downloads/soil_sen2_downloads"
            XX=process_folder(root_folder)
            ''
        for i in range(1,len(XX)):
            shutil.rmtree(XX[i])
            
        # move files
        import os
        import shutil
        import tifffile as tiff
        import json
        from datetime import datetime
        
        def process_folder(folder_path):
            image_dict = {}  # To store matching image and mask paths
        
            for root, _, files in os.walk(folder_path):
                for file in files:
                    if file.lower().endswith('.json'):
                        json_path = os.path.join(root, file)
                        with open(json_path) as json_file:
                            json_data = json.load(json_file)
        
                        start_time = json_data["request"]["payload"]["input"]["data"][0]["dataFilter"]["timeRange"]["from"]
                        end_time = json_data["request"]["payload"]["input"]["data"][0]["dataFilter"]["timeRange"]["to"]
                        
                        for tiff_file in files:
                            if tiff_file.lower().endswith('.tiff') or tiff_file.lower().endswith('.tif'):
                                tiff_path = os.path.join(root, tiff_file)
                                image_dict.setdefault(start_time, []).append(tiff_path)
        
            for start_time, files in image_dict.items():
                if len(files) == 2:
                    image_path = max(files, key=lambda x: os.path.getsize(x))
                    mask_path = min(files, key=lambda x: os.path.getsize(x))
        
                    # Convert start_time to a more filesystem-friendly format
                    start_time_dt = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%SZ")
                    formatted_start_time = start_time_dt.strftime("%Y-%m-%d_%H-%M-%S")
                    
                    target_subfolder = os.path.join(folder_path, formatted_start_time)
                    os.makedirs(target_subfolder, exist_ok=True)
                    
                    # Rename and move the image and mask files
                    sr_image_path = os.path.join(target_subfolder, "sr.tiff") # name of first file
                    scl_mask_path = os.path.join(target_subfolder, "scl.tiff") # name of second files
                    shutil.move(image_path, sr_image_path)
                    shutil.move(mask_path, scl_mask_path)
        
                    print(f"Transferred to subfolder: {target_subfolder}")
                    
                    # Delete the original subfolder
                 #   shutil.rmtree(os.path.dirname(image_path))
                    
            print("Processing completed.")
        
        if __name__ == "__main__":
            root_folder = "/media/mabso/Data/downloads/soil_sen2_downloads" # set to your root folder
            process_folder(root_folder)
        
        # remove all non date time files
        from datetime import datetime
        
        def is_datetime_format(string, format):
            try:
                datetime.strptime(string, format)
                return True
            except ValueError:
                return False
        
        def delete_non_datetime_folders(root_folder):
            for folder in os.listdir(root_folder):
                folder_path = os.path.join(root_folder, folder)
                if os.path.isdir(folder_path) and not is_datetime_format(folder, "%Y-%m-%d_%H-%M-%S"):
                    print(f"Deleting folder: {folder_path}")
                    # Uncomment the line below to actually delete the folder and its contents
                    shutil.rmtree(folder_path)
        
        if __name__ == "__main__":
            root_folder = "/media/mabso/Data/downloads/soil_sen2_downloads" # set to your root folder
            delete_non_datetime_folders(root_folder)
        
        import glob
        import rasterio
        
        download_root_directory = "/media/mabso/Data/downloads/soil_sen2_downloads"# set to your root folder
        
        original_file_path_l=[]
        new_file_path_l=[]
        
        inputs = glob.glob(download_root_directory+'/*/*', recursive=True)
        
        # Transpose the dimensions of each downloaded image file and write to new files
        #for root, _, files in os.walk(download_root_directory):
        for file_name in inputs:
            if file_name.endswith("scl.tiff"):  # Make sure to check the correct file extension
                original_file_path = os.path.join(download_root_directory, file_name)
                new_file_path = os.path.join(file_name, file_name.replace("scl.tiff","scl1.tiff"))
                
                # Open the raster file
                with rasterio.open(original_file_path) as src:
                    data = src.read()
                    print("Original data shape:", data.shape)
                    meta = src.meta.copy()  # Copy the metadata
                    
                # Transpose dimensions from (4, h, w) to (h, w, 4)
                data_transposed = data.transpose(1, 2, 0)
                print("Transposed data shape:", data_transposed.shape)
                # Reorder bands to match the original order (B01, B02, B03, B04)
                #data_reordered = data_transposed[:, :, [0, 1, 2, 3]]
                
                # Update the metadata
              #  meta['count'] = 4  # Update the band count
                
                # Save the transposed data to the new file
                with rasterio.open(new_file_path, 'w', **meta) as dst:
                    data_transposed = np.rollaxis(data_transposed, axis=2)
                    dst.write(data_transposed)  # Explicitly specify indexes
                
                    
                print(f"Saved transposed data to {new_file_path}")
                
                original_file_path_l.append(original_file_path)
                
                new_file_path_l.append(new_file_path)
                # Rename the scl1.tiff to scl.tiff
              #  os.rename(new_file_path, original_file_path)
        
              #  print(f"Renamed {new_file_path} to {original_file_path}")
        
        for i in range(0,len(original_file_path_l)):
            os.remove(original_file_path_l[i])
            os.rename(new_file_path_l[i], original_file_path_l[i])
        
        
        import os
        from datetime import datetime
        
        download_root_directory = "/media/mabso/Data/downloads/soil_sen2_downloads" # set to your root folder
        inva_list=[]
        
        for root, dirs, files in os.walk(download_root_directory):
            for folder in dirs:
                try:
                    # Parse the original folder name
                    original_date = datetime.strptime(folder, "%Y-%m-%d_%H-%M-%S")
                    # Create the new folder name with the format %Y-%m-%d
                    new_folder_name = original_date.strftime("%Y-%m-%d")
                    # Construct the paths
                    original_folder_path = os.path.join(root, folder)
                    new_folder_path = os.path.join(root, new_folder_name)
                    # Rename the folder
                    os.rename(original_folder_path, new_folder_path)
                    print(f"Renamed folder: {original_folder_path} to {new_folder_path}")
                except ValueError:
                    print(f"Skipped invalid folder name: {folder}")
                    inva_list.append(i)
        
        
        # read all folder names
        root_dir=download_root_directory
        inputs = glob.glob(root_dir+ '/*', recursive=True)
        dates_array=np.sort(inputs)
        
        # Target date string
        #samp_date_str = '06-07-18'
        
        # Convert the target date string to a datetime object
        samp_date = datetime.strptime(samp_date, '%d-%m-%y')
        
        # Function to extract the date from a path string
        def extract_date_from_path(path):
            return datetime.strptime(path.split('/')[-1], '%Y-%m-%d')
        
        # Convert the array of dates to datetime objects
        date_objects = np.array([extract_date_from_path(date_str) for date_str in dates_array])
        
        # Calculate the time differences between the target date and the array of dates
        time_diffs = np.abs(date_objects - samp_date)
        
        # Find the index of the date with the smallest time difference
        closest_date_index = np.argmin(time_diffs)
        
        # Get the closest date from the array
        closest_date = dates_array[closest_date_index]
        
        print("Closest date:", closest_date)
        
        import re
        
        match = re.search(r'\d{4}-\d{2}-\d{2}', closest_date)
        if match:
            closest_date = match.group()
        else:
            closest_date = None
        
        
        # List all files in the source directory
        source_directory=download_root_directory
        entries = os.listdir(source_directory)
        
        # Iterate through the entries and delete files and directories that don't contain the closest date
        for entry in entries:
            entry_path = os.path.join(source_directory, entry)
            
            if os.path.isdir(entry_path) and closest_date:
                # Check if it's a directory and closest_date is defined
                if closest_date not in entry:
                    # If the closest date is not in the directory name, remove it along with its contents
                    shutil.rmtree(entry_path)
        
        # List all files and directories in the source directory
        entries = os.listdir(source_directory)
        
        
        # Iterate through the entries and process directories
        for entry in entries:
            entry_path = os.path.join(source_directory, entry)
            
            if os.path.isdir(entry_path):
                # Check if it's a directory
                
                if closest_date != samp_date:
            # If the closest date is different from samp_date, include it in the directory name
                    new_dir_name = f"{id_name.iloc[0]}_{closest_date}_close"
                else:
            # Otherwise, use samp_date as is
                    new_dir_name = f"{id_name.iloc[0]}_{samp_date}"
        
                
                # Rename the directory
                new_dir_path = os.path.join(source_directory, new_dir_name)
                os.rename(entry_path, new_dir_path)
                
                # Move the renamed directory to the specified folder
                target_folder = '/media/mabso/Data/downloads/lucas_soil_data/NVDI_floats' # set to your folder name
                target_path = os.path.join(target_folder, new_dir_name)
                shutil.move(new_dir_path, target_path)
    except ValueError:
         print(f"Skipped invalid name: {soil_id}")
         
#%% read in what you have downloaded 
import math
# find the ids that matches crops 
idx=np.where(df_soil["LC0_Desc"]=="Cropland")[0]

# get the point ids
crop_ids=np.array(df_soil.iloc[idx]["POINTID"])

#crop_ids=np.array(df_soil["POINTID"])

# match with NDVI images
imgs=glob.glob("/media/mabso/Data/downloads/lucas_soil_data/NVDI_floats/*/sr.tiff",recursive=True)
#%% Extract ID numbers and dates of images
id_numbers = []
date_l=[]

for img_path in imgs:
    # Extract the folder name from the path
    folder_name = os.path.basename(os.path.dirname(img_path))
    
    # Split the folder name by underscore
    parts = folder_name.split('_')
    
    # The ID number is the first part before the underscore
    id_number = parts[0]
    date=parts[1]
    
    # Append the ID number to the list
    id_numbers.append(int(id_number))
    date_l.append(date)

id_numbers=np.array(id_numbers)

#%% match idx and id numbers of the crops
idx=np.where(np.isin(id_numbers,crop_ids))[0]
date_crops_l=np.array(date_l)[idx] # date list
id_crops_l=np.array(id_numbers)[idx] # id number list
img_crops_l1=np.array(imgs)[idx] # list of images based on id number

#%% plot your results to see what you have downloaded 
import tifffile as tiff

test=tiff.imread(imgs[6])

# plot what you have found
plt.imshow(test)


#%%
import tifffile as tiff
min_l=[]
max_l=[]
median_l=[]
point_id_l=[]
mean_l=[]
date_l=[]



for i in range(0,len(img_crops_l1)):
    test=tiff.imread(img_crops_l1[i])
 
    test = np.where(test >= 0.0, test, np.nan) # this is only relevant for NDVI values
    
    min_value = np.nanmin(test)
    max_value = np.nanmax(test)
    median_value=np.nanmedian(test)
    mean_value=np.nanmean(test)
    std_value=np.nanstd(test)
    pixel_count=np.sum(~np.isnan(test))
    point_id=id_crops_l[i]
    date=date_crops_l[i]
    
    min_l.append(min_value)
    max_l.append(max_value) 
    median_l.append(median_value)
    mean_l.append(mean_value)
    point_id_l.append(point_id)
    date_l.append(date)
    
data = {
    'Min': min_l,
    'Max': max_l,
    'Median': median_l,
    'Mean': mean_l,
    'Point_id': point_id_l,
    'date': date_l
}    
    
df = pd.DataFrame(data)

#%% merge soil dataframe with NDVI frame
idx_l=[]
for i in range(0,len(df["Point_id"])):
    idx=np.where(df["Point_id"][i]==np.array(df_soil["POINTID"]))[0]
    idx_l.append(int(idx))

df_soil_subset=df_soil.iloc[idx_l]

merged_df = pd.merge(df, df_soil_subset, left_on='Point_id', right_on='POINTID')

#%% export CSV
csv_file_path = "/home/mabso/Desktop/lucas_annoations/lucas_soil_ndvi_new8_postive.csv"

# Write the DataFrame to a CSV file
merged_df.to_csv(csv_file_path, index=False)
