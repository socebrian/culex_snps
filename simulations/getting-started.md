# SLiM simulations for training data for dispersal distances analyses: first cx. pipiens models
### 27/01/2025

Some notes about the files:
* make_mortality.R: script created by Peter to obtain the survival curve of adults Cx.pipiens. Used for the first draft of the model (approx two survival peaks on may and september)
* culex.slim: first draft of Cx. pipiens population model 
* river_greyscale.png: spatial map with rivers used in the original paper (find ref)
* survival.csv: numeric value of survival for each day of the year obteined from make_mortality.R
* survival.pdf: plot of survival from make-mortality.R

First step was to obtain a png greyscale image of the water bodies on my study areas. I used a .shp file that contained a grid os 2km*2km of Seville, Huelva, Malaga, Cadiz and Cordoba. In each cell it is represented the percentage that correspond to water bodies. I trnasformed this layer to a .tif raster using QGIS rasterize tool. Then I transformed the .tif to .png in greyscale using <gdal>

```bash
gdal_translate -scale 0 100 0 255 -ot Byte -of PNG raster-waterbody-grid2.tif raster-waterbody-grid2.png
```
This script transforms the scale from my .tif to the -8 bite scale of png greyscale.



