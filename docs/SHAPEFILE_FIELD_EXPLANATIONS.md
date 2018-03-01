Python scripts for GSFLOW -- HRU Parameters									
hru_parameters	explanation/equation	script	documentation source	data source	comment				
ORIG_FID	"HRU numbering, starting at 0, from the top left of the model grid"	hru_parameters.py							
HRU_ID	"HRU numbering, starting at 1, from the top left of the model grid"	hru_parameters.py							
HRUTYPE_IN	Integer values designating HRU cells: 0=inactive 1=active 2=lake	hru_parameters.py		model domain and lake polygons	"Initial hru type, does not change"				
HRU_TYPE	Integer values designating HRU cells: 0=inactive 1=active 2=lake	hru_parameters.py		model domain and lake polygons	"Can change, based on outermost active cells' slope"				
HRU_AREA	Area of each grid cell in acres	hru_parameters.py	Weasel	HRU grid					
HRU_ROW	Row number from top left of model grid	hru_parameters.py							
HRU_COL	Column number from top left of model grid	hru_parameters.py							
HRU_X	UTM X location at centroid of model grid cell	hru_parameters.py			based on specified datum				
HRU_Y	UTM Y location at centroid of model grid cell	hru_parameters.py			based on specified datum				
HRU_LAT	Latitude in Decimal Degrees at centroid of each model grid cell	hru_parameters.py			based on specified datum				
HRU_LON	Longitude in Decimal Degrees at centroid of each model grid cell	hru_parameters.py			based on specified datum				
DEM_MEDIAN	Median of original DEM values within each model grid cell	dem_parameters.py		DEM	based on specified DEM				
DEM_MEAN	Mean of original DEM values within each model grid cell	dem_parameters.py		DEM					
DEM_MIN	Minimum of original DEM values within each model grid cell	dem_parameters.py		DEM					
DEM_MAX	Maximum of original DEM values within each model grid cell	dem_parameters.py		DEM					
DEM_ADJ	Final HRU elevation to be used by CRT and by the PRMS input file	dem_parameters.py		DEM	"Default = DEM_FLOWAC , but typically undergoes manual adjustments to finalize stream network"				
DEM_FLOWAC	Weighs elevations of cells with streams heavier when up-scaling original DEM	dem_parameters.py	ArcHydro	DEM	Allows stream elevations to be better represented when upscaling a DEM to the model grid cell size				
DEM_SUM	Cumulative sum of flow accumulation values	dem_parameters.py		DEM	Zonal statistics				
DEM_COUNT	Flow accumulation count per HRU after a low pass filter?	dem_parameters.py	ArcHydro	DEM	Zonal statistics				
DEM_FEET	DEM_ADJ converted to feet	dem_parameters.py		DEM					
DEM_ASPECT	"Aspect at the original DEM resolution, then averaged within each HRU"	dem_parameters.py	Weasel	DEM					
DEM_SLP_D	"Slope (degrees) at the original DEM resolution, then averaged within each HRU"	dem_parameters.py		DEM					
DEM_SLP_R	"Slope (radians) at the original DEM resolution, then averaged within each HRU"	dem_parameters.py		DEM					
DEM_SLP_P	"Slope (percent) at the original DEM resolution, then averaged within each HRU"	dem_parameters.py		DEM					
DEM_SINK4	Internal algorithm to compute 4-way sinks 	dem_2_streams.py	Algorithm for 'flood_fill' in support_functions.py	DEM	Necessary to be compatible with MODFLOW				
DEM_SINK8	ArcGIS algorithm to compute 8-way sinks	dem_2_streams.py	ArcHydro	DEM					
CRT_ELEV	Value that CRT uses for cascade routing: DEM_ADJ plus CRT_FILL	dem_2_streams.py	CRT	HRU grid					
CRT_FILL	"Amount that DEM_ADJ will need to change to avoid swales, computed "	dem_2_streams.py	CRT	HRU grid					
LAKE_ID	Individual ID number assigned to all HRU cells  within each lake	dem_2_streams.py			Each lake is given a unique value				
LAKE_AREA	Cumulative area of cells that are 50% or more contained within a lake ploygon	dem_2_streams.py							
ISEG	"Numbering of gridded stream segments, unique for each segment"	dem_2_streams.py	SFR2		"Segment refers to a stream section, multiple grid cells long"				
IRUNBOUND	"Numbering of contributing HRUs to each stream segment, based on ISEG numbering"	dem_2_streams.py	SFR2		Negative values assigned to lake cells				
FLOW_DIR	Flow direction values indicate the direction of the steepest decent from an HRU	dem_2_streams.py	ArcHydro						
KRCH	An integer value equal to the layer number of the cell containing the stream reach	stream_parameters.py	SFR2 						
IRCH	Row number of grid cell containing the stream reach 	stream_parameters.py	SFR2 						
JRCH	Column number of grid cell containing the stream reach 	stream_parameters.py	SFR2						
IREACH	Sequential numbering of each stream reach (single grid cell) within each stream segment	stream_parameters.py	SFR2		Starts with 1 at the farthest upstream reach				
RCHLEN	Length of channel within a stream reach (single cell)	stream_parameters.py	SFR2		Can be greater than model cell dimensions due to meandering nature of streams				
MAXREACH	An integer value for each stream segment equal to its number of reaches	stream_parameters.py	SFR2		1 HRU = 1 reach				
OUTSEG	An integer value for a stream segment equal to the ISEG number to which it contributes	stream_parameters.py	SFR2		"If the segment ends in a lake, the value is the negative of the LAKE_ID"				
IUPSEG	An integer value of the upstream segment from which water is diverted (if applicable)	stream_parameters.py	SFR2		"If diverted from a lake, the value is the negative of the LAKE_ID from which it flows"				
SUB_BASIN	Integer values for sub-watershed scale drainages	stream_parameters.py			Sub-basins can be designated with user defined pour points or stream gages				
SEG_BASIN	Corresponds to the ISEG of the lowest segment ?	stream_parameters.py			Negative values assigned for lake segments				
OUTFLOWHRU	Outflow cells that exit the model domain to inactive cells?	stream_parameters.py							
STRM_TOP	Elevation of each stream cell set equal to DEM_ADJ - 1	stream_parameters.py							
STRM_SLOPE	Slope of stream cell	?			I can't find where this gets computed				
PPT_ZONE	Designated precipitation zone(s) based on shapefile provided by user (if applicable)	ppt_ratio_parameters.py			Default ppt zone =1  (no shapefile necessary of only using 1 zone)				
TOPO_INDEX	?				"This field is in the field list, but not in any script"				
DEPLCRV	?	dem_paramerters.py			Depletion curve set to 1 for all active cells?				
JH_TMAX	"Highest mean monthly temp, used  for Jenson-Haise PET  "	PRISM							
JH_TMIN	Mean minimum temp for the same month where highest mean temp was found	PRISM							
JH_COEF	Uses PRISM values... Modified Jenson-Haise to compute PET                                                      JH_COEF = 27.5 - [0.25 x (? sat vapor pressure of jh_tmax and Jh_tmin)] - hru_elev/1000	support_functions.py							
SNAREA_THR	"Maximum SWE threshold, below which the snow-covered-area curve is applied"	dem_parameters.py							
TMAX_ADJ	Adjustment factor based on ratio between PRISM data and climate station data	prms_template_fill.py?							
TMIN_ADJ	Adjustment factor based on ratio between PRISM data and climate station data	prms_template_fill.py?							
IMPERV_PCT	Impervious percentage assigned to each HRU	impervious_parameters.py	Markstrom et al PRMS IV	 Impervious Land Cover Dataset (NLCD)					
CAREA_MAX	Maximum portion of HRU area contributing to runoff (decimal fraction)	impervious_parameters.py	Markstrom et al PRMS IV	 Impervious Land Cover Dataset (NLCD)	Is this just from cells with impervious percent >0?				
COV_TYPE	"Integer values representing coverage type (0=bare soil, 1=grasses, 2=shrubs, 3=desiduous trees, 4=coniferous)"	veg_parameters.py	Markstrom et al PRMS IV	LANDFIRE					
COVDEN_SUM	"Major vegetation type summer coverage density, baeed on cov_type?"	veg_parameters.py	Markstrom et al PRMS IV	LANDFIRE					
COVDEN_WIN	"Major vegetation type winter coverage density, baeed on cov_type?"	veg_parameters.py	Markstrom et al PRMS IV	LANDFIRE					
RAD_TRNCF	Transmission coefficinet for short_wave radiation through winter veg canopy  rad_trncf = 0.9917 * Exp(-2.7557 * (covden_win))	veg_parameters.py	Markstrom et al PRMS IV						
SNOW_INTCP	Storage capacity for snow interception based on cov_type	veg_parameters.py	Markstrom et al PRMS IV						
SRAIN_INTC	Storage capacity for summer rain interception based on cov_type	veg_parameters.py	Markstrom et al PRMS IV						
WRAIN_INTC	Storage capacity for winter rain basedon cov_type and covden_win	veg_parameters.py	Markstrom et al PRMS IV						
ROOT_DEPTH	"Bottom of unsaturated zone, computed from root_depth remap file"	veg_parameters.py	Markstrom et al PRMS IV						
SOIL_TYPE	"Sand, silt, or clay, base on soil database percentages within each cell"	soil_parameters.py	Markstrom et al PRMS IV	SSURGO					
MOIST_INIT	Initial value of available water in capillary reservoir for each HRU	soil_parameters.py	Markstrom et al PRMS IV						
MOIST_MAX	Maximum available water holding capacity (AWC * root depth)	soil_parameters.py	Markstrom et al PRMS IV	SSURGO					
RECHR_INIT	Initial storage for soil recharge zone for each HRU	soil_parameters.py	Markstrom et al PRMS IV	must be < or = MOIST_INIT					
RECHR_MAX	Maximum storage for soil recharge zone for each HRU	soil_parameters.py	Markstrom et al PRMS IV						
SSR2GW_RATE	Linear coefficient in equation to route water from gravity reservoir to the groundwater reservoir (Ksat/sat_threshold)	soil_parameters.py	Markstrom et al PRMS IV	SSURGO 	Ksat must be in/day and sat_threshold = (%sand/100) * moist_max				
SLOWCOEF_L	Linear coefficient in equation to route gravity-reservoir storage down slope  (0.1 *  (Ksat sin??)/L_HRU) 	soil_parameters.py	Markstrom et al PRMS IV		?? = hru slope and L_HRU = length of hru cell				
SLOWCOEF_S	Non-linear coefficient in equation to route gravity reservoir storage down slope (0.9 *  (Ksat sin??)/(L_HRU*sat_threshold)) 	soil_parameters.py	Markstrom et al PRMS IV		*provide link to the pdf with derivation of these eqns*				
FASTCOEF_L	Linear coefficient in equation to route preferential flow storage down slope	soil_parameters.py	Markstrom et al PRMS IV						
FASTCOEF_S	Non-linear coefficient in equation to route preferential flow storage down slope	soil_parameters.py	Markstrom et al PRMS IV						
