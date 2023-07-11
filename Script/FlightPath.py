#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:18:35 2019

@author: italo

X/longitude/easting/columns/width/horizontal/left-right/across path, strip, track
Y/latitude/northing/rows/height/vertical/top-bottom/along path, strip, track
Z/altitude/zenith-nadir
"""

import os
from osgeo import gdal, ogr, osr
import numpy as np
from scipy import constants


class CameraUAV():
    '''
    
    Camera Imagery Footprint UAV and Direct Georeferencing
    
    A class used for modeling the internal and external orientation 
    of an camera (optical sensor) with a digital frame on a UAV platform, 
    and compute various relationships between them.
    
    Parameters
    ----------
    xSensor : float
        decimal number for width of the sensor size in millimeters
    ySensor : float
        decimal number for height of the sensor size in millimeters
    xPixels : int
        integer number of pixels recorded for image sensor in width or across path, strip, track
    yPixels : int
        integer number of pixels recorded for image sensor in height or along path, strip, track
    focalLen : float
        decimal number for focal length of the lens in millimeters
        
    Attributes
    ----------
    As a convention, we append an underscore (_) to attributes 
    that are not created upon the initialization of the object.   
    
    
    Methods
    -------
    

    '''
        
        
    def __init__(self, xSensor, ySensor, xPixels, yPixels, focalLen): 

        self.xSensor = xSensor
        self.ySensor = ySensor
        self.xPixels = xPixels
        self.yPixels = yPixels
        self.focalLen = focalLen

    
    def pixelPitch_(self, view=False):
        '''
        Pixel pitch is the distance from the center of one pixel
        to the center of the next pixel measured in micrometers (µm)
        
        The ratio of the sensor size to the number of recorded pixels,
        a ratio of the pixel sides of one is generally assumed.
        
        
        
        Parameters
        ----------
        xSensor : float
            decimal number for width of the sensor size in millimeters
        ySensor : float
            decimal number for height of the sensor size in millimeters
        xPixels : int
            integer number of pixels recorded for image sensor in width
        yPixels : int
            integer number of pixels recorded for image sensor in height
            
        Returns
        -------
        xPitch : float
            decimal number for the pixel x length in microns
        yPitch : float
            decimal number for the pixel y length in microns
            
        '''
        # Convert pixels microns to millimeters
        xPitch = np.round((self.xSensor / self.xPixels) * 1000, 2) 
        yPitch = np.round((self.ySensor / self.yPixels) * 1000, 2) 
        pxRatio = xPitch / yPitch
        
        if view == True:
            print(' Results '.center(79, "*"))
            print(f'Pixel pitch X (um): {xPitch} \nPixel pitch Y (um): {yPitch}')
            print(f'Pixel side ratio is: {pxRatio}')
            print(79 * "=")
            
        return xPitch, yPitch
    
    
    def FoV_(self, view=False):
        '''
        Compute lens Field of View like angles (AOV)
        

        Returns
        -------
        HAOV : float
            a decimal number to Horizontal Angle of View in degree
        VAOV : float
            a decimal number to Vertical Angle of View in degree
        DAOV : float
            a decimal number to Diagonal Angle of View in degree
                
        '''
                
        # Calculamos el arco tangente de la razon entre en ancho del sensor y la distancia focal
        HAOV = np.round(2 * np.degrees(np.arctan(self.xSensor / (2 * self.focalLen))), 1) # xview
        VAOV = np.round(2 * np.degrees(np.arctan(self.ySensor / (2 * self.focalLen))), 1) # yview
        # Elements_of_Photogrametry_and_GIS, equation 3-1
        DAOV = np.round(2 * np.degrees(np.arctan(np.sqrt(self.xSensor**2 + self.ySensor**2) / (2 * self.focalLen))), 1)
        
        if view == True:
            print(' Results '.center(79, "*"))
            print(f'Horizontal Field of View (°): {HAOV} ')
            print(f'Vertical Field of View (°): {VAOV} ')
            print(f'Diagonal Field of View (°): {DAOV} ')
            if DAOV > 0 and DAOV <= 75:
                print('Normal angle (up to 75°)')
            elif DAOV > 75 or DAOV < 100:
                print('Wide angle (75° to 100°)')
            else:
                print('Superwide angle (greater than 100°)')
            print(79 * "=")
        
        return HAOV, VAOV, DAOV

    
    def GSD(self, altitude, view=False):
        '''
        Compute Ground Sampling Distance (GSD) 
        is the distance between two consecutive pixel centers measured on the ground.
         
        Parameters
        ----------
        xPitch : float
            a decimal number for the pixel x length in microns
        yPitch : float
            a decimal number for the pixel y length in microns
        focalLen : float
            a decimal number focal length of lens in millimeters
        altitude : float
            a decimal number of height of flight on the ground in meter
  
        Returns
        -------
        GSD : float
            a decimal number to Ground Sampling Distance (centimeters/pixel)

        '''
        
        self.altitude = altitude
        xPitch, yPitch = self.pixelPitch_()
        
        # Ground Sampling Distance width (centimeters/pixel)
        GSDw = (self.xSensor * self.altitude * 100) / (self.focalLen * self.xPixels)
        # Ground Sampling Distance height (milimeters/pixel)
        GSDh = (yPitch * (self.altitude / self.focalLen)) / 10 # convert to centimeters/pixel

        if view == True:
            print(' Results '.center(79, "*"))
            print(f'Altitude: {self.altitude} \n')
            print(f'GSD X (cm/px): {GSDw} \nGSD Y (cm/pxl): {GSDh}')
            print(79 * "=")
        
        return GSDw, GSDh
    

    def gFoV(self, altitude, view=False):
        '''
        Compute Ground Field of View
        i.e. image sides length
         
        Parameters
        ----------

        altitude : float
            a decimal number of height of flight on the ground in meter

            
        Returns
        -------
        Dh : float
            Field of View Horizontal or width of single image footprint on the ground (meters)
        Dv : float
            Field of View Vertical or height of single image footprint on the ground (meters)
        
        '''

        self.altitude = altitude
        GSDw, GSDh = self.GSD(self.altitude)
        
        # Scale Computation: Convert pixels microns to millimeters * Convert GSD cm to mm
        # scale = round((self.pxSize/1000)/(GSD*10), 1)
        # width of single image footprint on the ground (meters)
        Dh = np.round((GSDw * self.xPixels) / 100, 2)
        # height of single image footprint on the ground (meters)
        Dv = np.round((self.ySensor * self.altitude) / self.focalLen, 2)
        # Area of imagery ground coverage 
        area = Dh * Dv
        
        if view == True:
            print(' Results '.center(79, "*"))
            print(f'Altitude: {self.altitude} \n')
            print(f'Area (m**2): {area} \n')
            print(f'Imagery Orthogonal Width (m): {Dh} \nImagery Orthogonal Height (m): {Dv}')
            print(79 * "=")
        
        return Dh, Dv
    
    def cornersCoordinates(self, altitude, E, N):
        
        '''
        Coordinates of ground field of view 
        with camera oriented in a nadir view and terrain perfectly planar
        
        Parameters
        ----------

        Exterior orientation of the image:

        altitude : float
            a decimal number of height of flight on the ground in meter 
        E: float
            a decimal number to easting coordinate in UTM
        N: float
            a decimal number to northing coordinate in UTM

            
        Returns
        -------

                
        '''
        
        self.altitude = altitude
        self.E = E
        self.N = N
        
        # image sides length 
        Dh, Dv = self.gFoV(self.altitude)
        
        # We will define a rectangle by the coordinates of its four corners
        # with x- and y-coordinates of the four corners, referred to as 
        # "UpperLeft", "LowerLeft", "LowerRight" and "UpperRight".
        # Corner up right
        UR_x = self.E + (Dh / 2)
        UR_y = self.N + (Dv / 2)
        # Corner up left
        UL_x = self.E - (Dh / 2)
        UL_y = self.N + (Dv / 2)
        # Corner down left
        LL_x = self.E - (Dh / 2)
        LL_y = self.N - (Dv / 2)
        # Corner down right
        LR_x = self.E + (Dh / 2)
        LR_y = self.N - (Dv / 2)
        
        return UR_x, UR_y, UL_x, UL_y, LL_x, LL_y, LR_x, LR_y 
    
    
    def footprintProjection(self, altitude, E, N, phi, omega, kappa):
        
        '''
        The area that is visible in one photo
        
        Funcion for direct georeferencing through the image footprint projection 
        (IFP) method, which utilizes collinearity equations on each image individually.
        (Lisein, 2013)
       
        
        Parameters
        ----------
        
        Exterior orientation of the image:
        E: float
            a decimal number to east coordinate
        N: float
            a decimal number to west coordinate
        altitude : float
            a decimal number to height of flight in meter
        phi : float
            a decimal number to x-axis GimbalRoll in degrees
        omega : float
            a decimal number to y-axis GimbalPitch in degrees      
        kappa : float
            a decimal number to z-axis GimbalYaw in degrees
            
        Interior Orientation of the camera
        focalLen : float
            a decimal number focal length of lens in millimeters
        xSensor : float
            a decimal number to width of sensor in millimeters
        ySensor : float
            a decimal number to height of sensor in millimeters
        
        
        Returns
        -------
                       
        '''
        
        self.E = E
        self.N = N
        self.altitude = altitude
        self.phi = phi
        self.omega = omega
        self.kappa = kappa
              
        # Convert Units
        Ls = self.xSensor / 1000 # Sensor length [m] 
        Hs = self.ySensor / 1000 # Sensor heigth [m]
        f = self.focalLen / 1000 # Focal lenght [m]
        # Coordinate of the four corners of the sensor in the image coordinate system
        # Corner up right
        x1 = Ls / 2
        y1 = Hs / 2
        # Corner up left
        x2 = -Ls / 2
        y2 = Hs / 2
        # Corner down left
        x3 = -Ls / 2
        y3 = -Hs / 2
        # Corner down right
        x4 = Ls / 2
        y4 = -Hs / 2
        
        # Rotational matrix to transform the image coordinate system to world coordinate system 
        m11 = np.cos(np.radians(self.phi)) * np.cos(np.radians(self.kappa))
        m12 = -np.cos(np.radians(self.phi)) * np.sin(np.radians(self.kappa))
        m13 = np.sin(np.radians(self.phi))
        m21 = np.cos(np.radians(self.omega)) * np.sin(np.radians(self.kappa)) + np.sin(np.radians(self.omega)) * np.sin(np.radians(self.phi)) * np.cos(np.radians(self.kappa))
        m22 = np.cos(np.radians(self.omega)) * np.cos(np.radians(self.kappa)) - np.sin(np.radians(self.omega)) * np.sin(np.radians(self.phi)) * np.sin(np.radians(self.kappa))
        m23 = -np.sin(np.radians(self.omega)) * np.cos(np.radians(self.phi))
        m31 = np.sin(np.radians(self.omega)) * np.sin(np.radians(self.kappa)) - np.cos(np.radians(self.omega)) * np.sin(np.radians(self.phi)) * np.cos(np.radians(self.kappa))
        m32 = np.sin(np.radians(self.omega)) * np.cos(np.radians(self.kappa)) + np.cos(np.radians(self.omega)) * np.sin(np.radians(self.phi)) * np.sin(np.radians(self.kappa))
        m33 = np.cos(np.radians(self.omega)) * np.cos(np.radians(self.phi))
        
        # The collinearity equations
        # Project the 4 corners of the sensor in the world coordinate system
        X1 = -self.altitude * ((m11 * x1 + m21 * y1 - m31 * f) / (m13 * x1 + m23 * y1 - m33 * f)) + E
        Y1 = -self.altitude * ((m12 * x1 + m22 * y1 - m32 * f) / (m13 * x1 + m23 * y1 - m33 * f)) + N
        X2 = -self.altitude * ((m11 * x2 + m21 * y2 - m31 * f) / (m13 * x2 + m23 * y2 - m33 * f)) + E
        Y2 = -self.altitude * ((m12 * x2 + m22 * y2 - m32 * f) / (m13 * x2 + m23 * y2 - m33 * f)) + N
        X3 = -self.altitude * ((m11 * x3 + m21 * y3 - m31 * f) / (m13 * x3 + m23 * y3 - m33 * f)) + E
        Y3 = -self.altitude * ((m12 * x3 + m22 * y3 - m32 * f) / (m13 * x3 + m23 * y3 - m33 * f)) + N
        X4 = -self.altitude * ((m11 * x4 + m21 * y4 - m31 * f) / (m13 * x4 + m23 * y4 - m33 * f)) + E
        Y4 = -self.altitude * ((m12 * x4 + m22 * y4 - m32 * f) / (m13 * x4 + m23 * y4 - m33 * f)) + N
        
        return X1, Y1, X2, Y2, X3, Y3, X4, Y4 # Corners : up right, up left, down left, down right 
        
    def polygonFootprint(self, vectorPath, crs, UR_x, UR_y, UL_x, UL_y, LL_x, LL_y, LR_x, LR_y):
        
        '''
        Save footprint in vector polygon GeoPackage
        
        Parameters
        ----------
            
        vectorPath : str
            a string with path, name and extention (.gpkg)
        '''
        self.vectorPath = vectorPath
        self.crs = crs
        self.UR_x = UR_x
        self.UR_y = UR_y
        self.UL_x = UL_x
        self.UL_y = UL_y
        self.LL_x = LL_x
        self.LL_y = LL_y
        self.LR_x = LR_x
        self.LR_y = LR_y
        
        # Instead of a polygon being made up of a list of vertices, like a line, 
        # they’re made of rings, which are made of vertices
        # A simple polygon only has one ring, 
        # To create polygon, you need a ring object and then add it to the polygon.
        # Vertices need to be added in order, but the direction around the perimeter 
        # can vary depending on the format you want to use to store the data.
        
        # For example, shapefiles specify that the outer rings are in clockwise order.
        
        # The first and last vertices must have the same coordinates so they close the ring. 
        # or you can call CloseRings on the ring or polygon after adding all vertices
        
        # Create a ring and add vertices
        ring = ogr.Geometry(ogr.wkbLinearRing)
        # starts with the upper left vertex and traverses the perimeter in counter-clockwise direction.
        ring.AddPoint(self.UL_x,self.UL_y) # ULC
        ring.AddPoint(self.UR_x,self.UR_y) # URC
        ring.AddPoint(self.LR_x,self.LR_y) # LRC
        ring.AddPoint(self.LL_x,self.LL_y) # LLC
        ring.CloseRings() # Close all rings in the polygon
        # Create a polygon and add the ring
        polygon = ogr.Geometry(ogr.wkbPolygon)
        polygon.AddGeometry(ring)
        
        # Creating spatial reference objects
        # you need to provide it when you create a new layer because 
        # you have no function to add an SRS to an existing layer.
        # standard EPSG code or a PROJ.4 string
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(self.crs) # The zero returned means that the SRS was imported successfully. 
        srs.GetAttrValue('PROJCS')
        
        # Create a new GeoPackage with OGR
        # select and set up the driver for our gpkg-file creation
        driver = ogr.GetDriverByName('GPKG')
            
        # Remove output vector if it already exists
        if os.path.exists(self.vectorPath):
            driver.DeleteDataSource(self.vectorPath)
                        
        # Create a new layer in data source GeoPackage
        datasource = driver.CreateDataSource(self.vectorPath)
        
        # CreateLayer(name, [srs], [geom_type], [options])
        layer = datasource.CreateLayer("footprintRPAS", srs, geom_type=ogr.wkbPolygon)
        
        # Define the attributes of the layer
        field_defn = ogr.FieldDefn('Area', ogr.OFTReal)
        field_defn.SetWidth(20)
        field_defn.SetPrecision(2) #added line to set precision
        layer.CreateField(field_defn)
        
        # Create a new feature (polygon) in the layer
        feature = ogr.Feature(layer.GetLayerDefn())
        
        # Set the geometry and attribute values for the feature
        area = polygon.GetArea() 
        print('Area (m**2) : ', area)
        feature.SetGeometry(polygon)
        feature.SetField('Area', area)
        
        # Add the feature to the layer
        layer.CreateFeature(feature)

        # # Create centroids of polygon
        centroid = polygon.Centroid()
        centroid_x = centroid.GetX()
        centroid_y = centroid.GetY()
        print('Centroids coordinates x,y : ', centroid_x , centroid_y)

        # Save and close everything
        datasource = layer = feature = None

        return centroid_x , centroid_y

    def GCP2imagery(self, rasterPath, crs, UR_x, UR_y, UL_x, UL_y, LL_x, LL_y, LR_x, LR_y):
        '''
        Add georeference to imagery 
        from known coordinates for specific pixel offsets.
        
        Parameters
        ----------
            
        filename: str
            a string with path, name and extention (.tif)
            
        '''
                
        self.rasterPath = rasterPath
        self.crs = crs
        self.UR_x = UR_x
        self.UR_y = UR_y
        self.UL_x = UL_x
        self.UL_y = UL_y
        self.LL_x = LL_x
        self.LL_y = LL_y
        self.LR_x = LR_x
        self.LR_y = LR_y
                
        ds = gdal.Open(self.rasterPath, gdal.GA_Update)
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(self.crs)
                
        gcps = [gdal.GCP(self.UL_x, self.UL_y, 0, 0, 0), # gdal.GCP([x], [y], [z], [pixel-column], [line-row], 
                gdal.GCP(self.UR_x, self.UR_y, 0, self.xPixels, 0),
                gdal.GCP(self.LR_x, self.LR_y, 0, self.xPixels, self.yPixels),
                gdal.GCP(self.LL_x, self.LL_y, 0, 0, self.yPixels)]
                
        ds.SetGCPs(gcps, sr.ExportToWkt())
        ds = None
        
        # Geotransform
        gt = gdal.GCPsToGeoTransform(gcps)
        
        print('Origin:', (gt[0], gt[3])) # GeoTransform item: 0 = Origin x coordinate, 3 = Origin y coordinate
        print('Pixel size:', (gt[1], gt[5])) # GeoTransform item: 1 = Pixel width, 5 = Pixel height (negative)
        print('X pixel rotation:', (gt[2])) # GeoTransform item: 2 = x pixel rotation (0° if image is north up)
        print('Y pixel rotation', (gt[4])) # GeoTransform item: 4 = y pixel rotation (0° if image is north up)


class LidarUAV:
    ''' 
    A class based on 
    Flight Planning for LiDAR-Based UAS Mapping Applications paper and
    Airborne laser scanning: basic relations and formulas
    
    Based on assumptions like constant flying speed and  flying height, 
    the terrain is flat, vertical scan, roll and pitch angle are zero

    ...

    Attributes
    ----------
    vfov: float
        along-track scanning field of view 
    hfov: float
        cross-track scanning field of view 


    Methods
    -------
    
    '''

    
    def __init__(self, vfov, hfov, prf):
                
        '''
        Parameters
        ----------
        vfov: float
            along-track scanning field of view
        hfov: float
            cross-track scanning field of view, 
            scanning field of view out of the offered 360 FOV
        prf : int
            pulse repetition frequency
                    
        '''
        
        self.vfov = vfov
        self.hfov = hfov
        self.prf = prf
        
        
    def pointCloudDensityEstimate(self, altitude, velocity, overlap=50):
        '''
        Compute Point Cloud Density Estimate

        Parameters
        ----------
        altitude : TYPE
            DESCRIPTION.
        velocity : TYPE
            DESCRIPTION.

        Returns
        -------
        d : float
            Average point density (points / m**2)

        '''
        
        self.altitude = altitude
        self.velocity = velocity
        self.overlap = overlap
        
        L, W = self.gFoV(self.altitude)
        
        if self.overlap > 0:
            d = self.prf / (360/self.hfov) / (L * self.velocity) * (1 / (1 - (self.overlap / 100)))
        else: 
            print('No valid option')
        
        return d
        
    
    def gFoV(self, altitude, view=False):
        '''
        Compute footprint sides length
        
        footprint based on a LiDAR sensor mounted with a nadir orientation
        project area dimensions defined as a rectangle with a length L and width W
         
        Parameters
        ----------

        altitude : float
            a decimal number of height of flight on the ground in meter

            
        Returns
        -------
        W : float
            along track scanning width
            length of single footprint on the ground (meters)
        L : float
            across track swath width of the scanning
            width of single footprint on the ground (meters)
        
        '''

        self.altitude = altitude

        L = 2 * np.tan(np.radians(self.hfov / 2)) * self.altitude
        W = 2 * np.tan(np.radians(self.vfov / 2)) * self.altitude
        
        return L, W
    
    def cornersCoordinates(self, altitude, E, N):
        
        '''
        Coordinates of ground field of view 
        with camera oriented in a nadir view and terrain perfectly planar
        
        Parameters
        ----------

        Exterior orientation of the image:
        
        E: float
            a decimal number to easting coordinate in UTM
        N: float
            a decimal number to northing coordinate in UTM
        altitude : float
            a decimal number of height of flight on the ground in meter 
            
        Returns
        -------

                
        '''
        
        self.altitude = altitude
        self.E = E
        self.N = N
        
        
        # Footprint sides length 
        L, W = self.gFoV(self.altitude)
        
        # We will define a rectangle by the coordinates of its four corners
        # with x- and y-coordinates of the four corners, referred to as 
        # "UpperLeft", "LowerLeft", "LowerRight" and "UpperRight".
        # Corner up right
        UR_x = self.E + (L / 2)
        UR_y = self.N + (W / 2)
        # Corner up left
        UL_x = self.E - (L / 2)
        UL_y = self.N + (W / 2)
        # Corner down left
        LL_x = self.E - (L / 2)
        LL_y = self.N - (W / 2)
        # Corner down right
        LR_x = self.E + (L / 2)
        LR_y = self.N - (W / 2)
        
        return UR_x, UR_y, UL_x, UL_y, LL_x, LL_y, LR_x, LR_y 
    
    
    def rotateCorners(self, kappa, UR_x, UR_y, UL_x, UL_y, LL_x, LL_y, LR_x, LR_y):
        '''
        Rotate vertex coordinate
        
        Parameters
        ----------
        
        
            
        '''
                
        self.UR_x = UR_x
        self.UR_y = UR_y
        self.UL_x = UL_x
        self.UL_y = UL_y
        self.LL_x = LL_x
        self.LL_y = LL_y
        self.LR_x = LR_x
        self.LR_y = LR_y
        self.kappa = kappa
        
        ox = (self.UL_x + self.UR_x) / 2
        oy = (self.UR_y + self.LR_y) / 2
        
        x1r = ox + np.cos(np.radians(self.kappa)) * (self.UR_x - ox) + np.sin(np.radians(self.kappa)) * (self.UR_y - oy)
        y1r = oy + -np.sin(np.radians(self.kappa)) * (self.UR_x - ox) + np.cos(np.radians(self.kappa)) * (self.UR_y - oy)

        x2r = ox + np.cos(np.radians(self.kappa)) * (self.UL_x - ox) + np.sin(np.radians(self.kappa)) * (self.UL_y - oy)
        y2r = oy + -np.sin(np.radians(self.kappa)) * (self.UL_x - ox) + np.cos(np.radians(self.kappa)) * (self.UL_y - oy)

        x3r = ox + np.cos(np.radians(self.kappa)) * (self.LL_x - ox) + np.sin(np.radians(self.kappa)) * (self.LL_y - oy)
        y3r = oy + -np.sin(np.radians(self.kappa)) * (self.LL_x - ox) + np.cos(np.radians(self.kappa)) * (self.LL_y - oy)

        x4r = ox + np.cos(np.radians(self.kappa)) * (self.LR_x - ox) + np.sin(np.radians(self.kappa)) * (self.LR_y - oy)
        y4r = oy + -np.sin(np.radians(self.kappa)) * (self.LR_x - ox) + np.cos(np.radians(self.kappa)) * (self.LR_y - oy)
        
        return x1r,y1r,x2r,y2r,x3r,y3r,x4r,y4r
    


    def polygonFootprint(self, vectorPath, crs, UR_x, UR_y, UL_x, UL_y, LL_x, LL_y, LR_x, LR_y):
        
        '''
        Save footprint in vector polygon GeoPackage
        
        Parameters
        ----------
            
        vectorPath : str
            a string with path, name and extention (.gpkg)
        '''
        self.vectorPath = vectorPath
        self.crs = crs
        self.UR_x = UR_x
        self.UR_y = UR_y
        self.UL_x = UL_x
        self.UL_y = UL_y
        self.LL_x = LL_x
        self.LL_y = LL_y
        self.LR_x = LR_x
        self.LR_y = LR_y
        
        # Instead of a polygon being made up of a list of vertices, like a line, 
        # they’re made of rings, which are made of vertices
        # A simple polygon only has one ring, 
        # To create polygon, you need a ring object and then add it to the polygon.
        # Vertices need to be added in order, but the direction around the perimeter 
        # can vary depending on the format you want to use to store the data.
        
        # For example, shapefiles specify that the outer rings are in clockwise order.
        
        # The first and last vertices must have the same coordinates so they close the ring. 
        # or you can call CloseRings on the ring or polygon after adding all vertices
        
        # Create a ring and add vertices
        ring = ogr.Geometry(ogr.wkbLinearRing)
        # starts with the upper left vertex and traverses the perimeter in counter-clockwise direction.
        ring.AddPoint(self.UL_x,self.UL_y) # ULC
        ring.AddPoint(self.UR_x,self.UR_y) # URC
        ring.AddPoint(self.LR_x,self.LR_y) # LRC
        ring.AddPoint(self.LL_x,self.LL_y) # LLC
        ring.CloseRings() # Close all rings in the polygon
        # Create a polygon and add the ring
        polygon = ogr.Geometry(ogr.wkbPolygon)
        polygon.AddGeometry(ring)
        
        # Creating spatial reference objects
        # you need to provide it when you create a new layer because 
        # you have no function to add an SRS to an existing layer.
        # standard EPSG code or a PROJ.4 string
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(self.crs) # The zero returned means that the SRS was imported successfully. 
        srs.GetAttrValue('PROJCS')
        
        # Create a new GeoPackage with OGR
        # select and set up the driver for our gpkg-file creation
        driver = ogr.GetDriverByName('GPKG')
            
        # Remove output vector if it already exists
        if os.path.exists(self.vectorPath):
            driver.DeleteDataSource(self.vectorPath)
                        
        # Create a new layer in data source GeoPackage
        datasource = driver.CreateDataSource(self.vectorPath)
        
        # CreateLayer(name, [srs], [geom_type], [options])
        layer = datasource.CreateLayer("footprintRPAS", srs, geom_type=ogr.wkbPolygon)
        
        # Define the attributes of the layer
        field_defn = ogr.FieldDefn('Area', ogr.OFTReal)
        field_defn.SetWidth(20)
        field_defn.SetPrecision(2) #added line to set precision
        layer.CreateField(field_defn)
        
        # Create a new feature (polygon) in the layer
        feature = ogr.Feature(layer.GetLayerDefn())
        
        # Set the geometry and attribute values for the feature
        area = polygon.GetArea() 
        print('Area (m**2) : ', area)
        feature.SetGeometry(polygon)
        feature.SetField('Area', area)
        
        # Add the feature to the layer
        layer.CreateFeature(feature)

        # # Create centroids of polygon
        centroid = polygon.Centroid()
        centroid_x = centroid.GetX()
        centroid_y = centroid.GetY()
        print('Centroids coordinates x,y : ', centroid_x , centroid_y)

        # Save and close everything
        datasource = layer = feature = None

        return centroid_x , centroid_y



class UAV():
    '''
    A class for Stereoscopic Coverage of Survey Area                    

    Each photograph overlaps the next photograph in direction of the flightline 
    by approximately 60 % (forward overlap or endlap), 
    while adjacent lines overlap by 20 % – 30 % (sidelap).

    
    
    '''
    def __init__(self, Fw, Fh): 

        self.Fw = Fw # footprint length width
        self.Fh = Fh # footprint length height
        
    def endlap(self, pe):
        
        '''
        Overlap in flight direction, in percent

        The required air base or distance between exposure stations
        is dependent on the dimensions of the image footprint and the desired endlap.
        
        Parameters
        ----------

        pe : float
            endlap percent (%) between adjacent images

        '''
        
        self.pe = pe
        
        el = self.Fh * (1 - self.pe / 100) # air base (m), endlap (along-track)
        
        print(el)
        
        return el
    
    
    def sidelap(self, pe):
        '''
        Overlap between fight lines, in percent

        Separation distance between the flight strips
        
        
        Parameters
        ----------

        pe : float
            sidelap percent (%) between adjacent scanning strips.

        
        '''
        
        self.pe = pe
        
        sl = self.Fw * (1 - self.pe / 100) # air base (m), sidelap (cross-track)
        
        print(sl)
        
        return sl
    
    def intervalExposure(self, Bh, Bw, Vg):
        '''
        Time interval between exposures (s)
        
        For platforms with approximately constant ground velocity Vg

        Parameters
        ----------
        Bh : float
            air base (m)
        Bw : float
            separation distance between the flight strips (m)
        Vg : float
            flying speed, constant ground velocity (m/s)

        Returns
        -------
        None.

        '''
        
        self.Bh = Bh
        self.Bw = Bw        
        self.Vg = Vg
        
        Th = self.Bh / self.Vg # time interval between exposures (s)
        Th
        
        Tw = self.Bw / self.Vg # time interval between exposures (s)
        Tw
        
        return Th, Tw
    
    def cogo(self, E1, N1, distance, azimuth) :
        
        '''
        COGO, an acronym for coordinate geometry, 
        calculates locations using distances and bearings from known reference points.
        
        Create waypoints for flight path with COGO using bearing and distance

        Parameters
        ----------
        E1 : float
            Easting coordinate in UTM.
        N1 : float
            Norting coordinate in UTM.
        distance : float
            Distance in meters.
        azimuth : float
            Azimuth in degrees.

        Returns
        -------
        E2 : TYPE
            DESCRIPTION.
        N2 : TYPE
            DESCRIPTION.

        '''
        
        self.E1 = E1
        self.N1 = N1
        self.distance = distance
        self.azimuth = azimuth
        
        dE = self.distance * np.sin(np.radians(self.azimuth))
        dN = self.distance * np.cos(np.radians(self.azimuth))
        
        print(dE, dN)
        E2 = E1 + dE
        N2 = N1 + dN
        
        return E2, N2
    
    def curvedTurns(self, s, r, Eo, No, azimuth, dr = 'CW') :
        '''
        Curved Turn Waypoint Missions
        
        Based on angle measures of regular polygons
        
        Parameters
        ----------
        s : int
            number of sides (half of sides of the regular polygon)
        r : float
            radii of curved turn
        Eo : float
            easting start coordinate
        No : float
            northing start coordinate   
        azimuth: float
            direction of flight line
        dr : float
            direction of rotation
            Any rightward motion in a circular fashion is known as clockwise, 
            often denoted as CW, 
            the circular motion towards the left is counterclockwise or anticlockwise, 
            often denoted as CCW and ACW, respectively.
            

        Returns
        -------
        None.

        '''
        
        self.Eo = Eo
        self.No = No
        self.s = s
        self.r = r
        self.azimuth = azimuth
        S = self.s * 2
        
        siarp = (S - 2) * 180 # Sum of Interior Angles of a Regular Polygon
        iarp =  siarp / S # Interior Angles of a Regular Polygon
        earp = 360 / S # Exterior Angle of a Regular Polygon
        
        apothem = self.r * np.cos(constants.pi / S)
        #apothem = self.r * np.cos(np.radians(180) / self.n)
        side = 2 * apothem * np.tan(constants.pi / S)
        #side = 2 * self.r * np.sin(np.radians(180) / self.n) 
        
        print(siarp, iarp, earp, apothem, side)
        
        # Generate coordinates for curved turn
        ex, nx = [], []
        
        e, n = self.Eo, self.No
        
        angle = self.azimuth
        
        if dr == 'CW' :
            for i in range(1, int(self.s) + 1):
                print(i)
                if i == 1 :
                    angle = angle + (earp/2)
                    print('angle', angle)
                    e, n = self.cogo(e, n, side, angle)
                    print(e,n)
                    ex.append(e)
                    nx.append(n)
                    
                else :
                    angle = angle + earp
                    print('angle', angle)
                    e, n = self.cogo(e, n, side, angle)
                    print(e,n)
                    ex.append(e)
                    nx.append(n)
                    
        elif dr == 'CCW' :
            for i in range(1, int(self.s) + 1):
                print(i)
                if i == 1 :
                    angle = angle - (earp/2)
                    print('angle', angle)
                    e, n = self.cogo(e, n, side, angle)
                    print(e,n)
                    ex.append(e)
                    nx.append(n)
                    
                else :
                    angle = angle - earp
                    print('angle', angle)
                    e, n = self.cogo(e, n, side, angle)
                    print(e,n)
                    ex.append(e)
                    nx.append(n)
        else:
            print('No valid option')
                
        return ex, nx
    
    def flightPath(self, vectorPath, crs, Ex, Nx):
        '''
        Save flight path in vector line GeoPackage
        

        Parameters
        ----------
        vectorPath : str
            a string with path, name and extention (.gpkg)
        crs : int
            spatial reference systems in EPGS code
        Ex : list
            easting coordinates
        Nx : list
            norting coordinates

        Returns
        -------
        None.

        '''
        
        self.vectorPath = vectorPath
        self.crs = crs
        self.Ex = Ex
        self.Nx = Nx
        
        # Create a line geometry
        line = ogr.Geometry(ogr.wkbLineString)
        
        for e, n in zip(self.Ex, self.Nx): # Add points
            print(e, n)
            line.AddPoint(e, n)  
            
        # Creating spatial reference objects
        # standard EPSG code or a PROJ.4 string
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(self.crs) # The zero returned means that the SRS was imported successfully. 
        
        # Create a new GeoPackage
        driver = ogr.GetDriverByName('GPKG')        
        
        # Remove output vector if it already exists
        if os.path.exists(self.vectorPath):
            driver.DeleteDataSource(self.vectorPath)
                        
        # Create a new layer in data source GeoPackage
        datasource = driver.CreateDataSource(self.vectorPath)
        
        # Create a new layer in the GeoPackage        
        layer = datasource.CreateLayer('flightPath', srs, ogr.wkbLineString)
        
        # Define the attributes of the layer
        field_defn = ogr.FieldDefn('Length', ogr.OFTReal)
        field_defn.SetWidth(20)
        field_defn.SetPrecision(2) #added line to set precision
        layer.CreateField(field_defn)
        
        # Create a new feature (line) in the layer
        feature = ogr.Feature(layer.GetLayerDefn())        
        
        # Set the geometry and attribute values for the feature
        length = line.Length()
        print("Longitud : (m)", length)        
        feature.SetGeometry(line)
        feature.SetField('Length', length)
        
        # Add the feature to the layer
        layer.CreateFeature(feature)

        # Save and close everything
        datasource = layer = feature = None
