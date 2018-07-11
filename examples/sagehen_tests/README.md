# Tests Models

There are currently four models that have been developed from the Sagehen model to test different aspects of the stream delineation portion of the scripts.

## active_fill

Define the study area to include cells that are well outside the basin.
These cells will be made active and filled so that they flow to the outlet point.

## clipped_elev

Clip the input elevation grid so that all inactive cells have no elevation data.
Check that the model still correctly identifies the outlet point(s).

## closed_basin

The DEM was modified to make Independence Lake a closed basin and a "SWALE" point was defined in the middle of the lake.

## closed_basin_b

This is the same DEM as the other closed basin example but stream "OUTLET" point was removed so that the "SWALE" point in the middle of Independence Lake is the only outlet.  This test was added because the scripts were not running when the model did not include an "OUTLET" point.
