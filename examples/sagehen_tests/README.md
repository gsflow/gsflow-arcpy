# Sagehen Tests Models

There are currently six models that have been developed from the Sagehen model to test different aspects of the scripts.

## active_fill

Define the study area to include cells that are well outside the basin.
These cells will be made active and filled so that they flow to the outlet point.

## clipped_elev

Clip the input elevation grid so that all inactive cells have no elevation data.
Check that the model still correctly identifies the outlet point(s).

## lake

Add a single lake to the model.  The elevations in the lake have been manually adjusted to the same value in the provided DEM.

## precip_zones

## terminal_basin

The DEM was modified to "dam" the eastern portion of the watershed and the outlet point was changed to a "SWALE" type to simulate a closed/terminal basin.  This test was added because the scripts were not running when the model did not include an "OUTLET" point.

## terminal_lake

This is the same DEM as the other terminal basin example but a lake was included around the "SWALE" point.  This test was added to ensure that a SWALE point in a lake works the same as one outside.
