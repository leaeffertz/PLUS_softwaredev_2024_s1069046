"""About this program

This file contains all the code neccessary to create and display a cloud mask with s2cloudless, 
while process of what I was doing is in the Jupyter notebook, along with the functions being called.
"""

# importing the earth engine module and folium (to later display the result)
import ee
import folium

# authenticate earth engine
ee.Authenticate()

# initialize the library I use
ee.Initialize(project="ee-psd")

#set parameters
AOI = ee.Geometry.Point(8.642, 49.877)
START_DATE = "2024-06-03"
END_DATE = "2024-06-04"
CLOUD_FILTER = 60
CLD_PRB_THRESH = 50
NIR_DRK_THRESH = 0.15
CLD_PRJ_DIST = 1
BUFFER = 50

def get_s2_sr_cld_col(aoi, start_date, end_date):
    """
    Retrieves filtered collections of sentinel 2 images and joins the collection with 
    a cloud propability dataset to better detect clouds
    
    Parameters:
        aoi - ee.Geometry : Coordinates for area of interest
        start_date - str : start date for collection (inclusive)
        end_date - str : end date for colection (exclusive)
    
    Returns:
        ee.ImageCollection: filtered collection of Sentinel-2 images joined with cloud
        probability data
    """
    # Add sentinel 2 imagery
    s2_sr_col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lte("CLOUDY_PIXEL_PERCENTAGE", CLOUD_FILTER)))
    # Add cloud probability data
    s2_cloudless_col = (ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
        .filterBounds(aoi)
        .filterDate(start_date, end_date))

    # Join both datasets by index
    return ee.ImageCollection(ee.Join.saveFirst("s2cloudless").apply(**{
        "primary": s2_sr_col,
        "secondary": s2_cloudless_col,
        "condition": ee.Filter.equals(**{
            "leftField": "system:index",
            "rightField": "system:index"
        })
    }))
# Store created collection in variable
s2_sr_cld_col_eval = get_s2_sr_cld_col(AOI, START_DATE, END_DATE)

def add_cloud_bands(img):
    """
    Adds data as bands from image collection to the new image,
    data with a cloud probability above the threshold.
    
    Parameters:
        img - ee.Image : input image that will be mask
    
    Returns:
        img - ee.Image : original image, with bands added
    """

    # Get cloud probability from the collection we created earlier
    cld_prb = ee.Image(img.get("s2cloudless")).select("probability")

    # Rename items in dataset to clouds, if they are above the set threshold 
    is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename("clouds")

    # Add both image layers as bands
    return img.addBands(ee.Image([cld_prb, is_cloud]))



def add_shadow_bands(img):
    """
    Identifies potential cloud shadows by pixel value and adds them to mask
    
    Parameters:
        img - ee.Image : input image, that will be mask
    
    Returns:
        img - ee.Image : mask, now with added shadow and dark px bands
    """


    # Select everything that is not water (6) from scene classification values band
    not_water = img.select("SCL").neq(6)

    # Identify dark pixels on the NIR band and name them accordingly
    SR_BAND_SCALE = 1e4
    dark_pixels = img.select("B8").lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename("dark_pixels")

    # Determine the mean of the solar azimuth angle (angle between proj. of sun rays), 
    # to get direction where clouds are projected (UTM projection assumed)
    shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get("MEAN_SOLAR_AZIMUTH_ANGLE")));

    # Project cloud shadows, specified by input (CLD_PRJ_DIST)
    cld_proj = (img.select("clouds").directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject(**{"crs": img.select(0).projection(), "scale": 100})
        .select("distance")
        .mask()
        .rename("cloud_transform"))

    # Classify shadows by multiplying areas where there are cloud projections identified and dark pixels
    shadows = cld_proj.multiply(dark_pixels).rename("shadows")

    # Add all of the above shadow elements to the mask
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

def add_cld_shdw_mask(img):
    """
    Combine both the identified clouds and shadows into one mask
    
    Parameters:
        img - ee.Image : input image, to which the cloud- and shadow mask is added
    
    Returns:
        ee.Image : input image, now with mask applied
    """

    # Run functions to bands
    img_cloud = add_cloud_bands(img)
    img_cloud_shadow = add_shadow_bands(img_cloud)

    # Combine both band layers, if there are neither cloud nor shadow, set value to 0
    is_cld_shdw = img_cloud_shadow.select("clouds").add(img_cloud_shadow.select("shadows")).gt(0)

    # Add buffer to reduce noise (small areas falsely classified as clouds), 20 m scale for performance and return the created mask applied to the image
    is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(BUFFER*2/20)
        .reproject(**{"crs": img.select([0]).projection(), "scale": 20})
        .rename("cloudmask"))

    return img_cloud_shadow.addBands(is_cld_shdw)



# Display Earth Engine image tiles in folium map
def add_ee_layer(self, ee_image_object, vis_params, name, show=True, opacity=1, min_zoom=0):
    """
    Adds ee-layer to folium, in order to display result in interactive map
    
    Parameters:
        self - folium.Map : folium instance, to which layer will be added
        ee_image_object - ee.Image : ee-image object
        vis_params - dict: visualization parameters
        name - str : layer name
        show - bool : set initial visibility of the layer
        opacity - float: layer opacity
        min_zoom - int: minimum zoom level for layer visibility
    """

    map_id_dict = ee.Image(ee_image_object).getMapId(vis_params)
    folium.raster_layers.TileLayer(
        tiles=map_id_dict["tile_fetcher"].url_format,
        attr="Map Data &copy; <a href='https://earthengine.google.com/'>Google Earth Engine</a>",
        name=name,
        show=show,
        opacity=opacity,
        min_zoom=min_zoom,
        overlay=True,
        control=True
        ).add_to(self)

folium.Map.add_ee_layer = add_ee_layer

def display_cloud_layers(col):
    """
    Display the previously created image masks in an interactive folium map
    
    Parameters:
        col - ee.ImageCollection : Sentinel-2 image collection to display
    """
        
    # create mosaic of previously detected clouds
    img = col.mosaic()

    # Select layers for display and updtae values
    clouds = img.select("clouds").selfMask()
    shadows = img.select("shadows").selfMask()
    dark_pixels = img.select("dark_pixels").selfMask()
    probability = img.select("probability")
    cloudmask = img.select("cloudmask").selfMask()
    cloud_transform = img.select("cloud_transform")

    # Create map object with folium
    center = AOI.centroid(10).coordinates().reverse().getInfo()
    m = folium.Map(location=center, zoom_start=12)

    # Add layers to map
    m.add_ee_layer(img,
                   {"bands": ["B4", "B3", "B2"], "min": 0, "max": 2500, "gamma": 1.1},
                   "S2 image", True, 1, 9)
    m.add_ee_layer(probability,
                   {"min": 0, "max": 100},
                   "probability (cloud)", False, 1, 9)
    m.add_ee_layer(clouds,
                   {"palette": "e056fd"},
                   "clouds", False, 1, 9)
    m.add_ee_layer(cloud_transform,
                   {"min": 0, "max": 1, "palette": ["white", "black"]},
                   "cloud_transform", False, 1, 9)
    m.add_ee_layer(dark_pixels,
                   {"palette": "orange"},
                   "dark_pixels", False, 1, 9)
    m.add_ee_layer(shadows, {"palette": "yellow"},
                   "shadows", False, 1, 9)
    m.add_ee_layer(cloudmask, {"palette": "orange"},
                   "cloudmask", True, 0.5, 9)

    # Add layer controls
    m.add_child(folium.LayerControl())

# Display map
s2_sr_cld_col_eval_disp = s2_sr_cld_col_eval.map(add_cld_shdw_mask)

#Test docstrings
#if __name__ == "__main__":
    #help(get_s2_sr_cld_col)
