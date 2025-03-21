# CATFLOW preprocessing tool

Pre- and postprocessing workflows for CATFLOW using functions from the [Catflow-R-Package](https://github.com/CATFLOW/Catflow-R-Package).

The implementation follows the [Tool Specification](https://vforwater.github.io/tool-specs/) for reusable research software using Docker.

## Tools

### make_representative_hillslope

**Title:** Create a CATFLOW geometry file from .tif files using the method of representative hillslopes (Loritz et al 2017).

**Description:**
Hillslope wizard developed by Ralf Loritz (2015)
Delination of 2D catenas from hillslopes. The idea behind this code is based on a approach from Francke et al. 2006 (Automated 
catena-based discretization of landscapes for the derivation of hydrological modelling 
units) and Cochrane and Flanagan (2003). The idea behind this approach is to select every unique distance and their corrosponding 
elevation from a river inside a hillslope. If more than one raster cell have the same 
distance the elevation is calculated by a weighted mean using the flow accumulation. The 
hillslope width can be estimated by the frequency of the appearance of a single distance. 
In this present version, we assume a fixed horizontal and varying vertical grid spacing for the hillslope element as in Manoj J et al (2024). This would be revised in the
upcoming versions.

**Parameters:**
- `hillslope_id`: Integer ID of the hillsope from hillslope.tif for calculating the geometry. If entire basin is to be used for the representative hillslope, give -1. (default: -1)
- `no_flow_area`: Percentage of no flow area with almost no slope within the area of interest. (default: 0.30)
- `min_cells`: Minimum number of cells within a hillslope. (default: 10)
- `hill_type`: Hillslope type. (1) constant thickness (default), (2) cake-shape, (3) variable thickness with spline approximation of lower boundary. Refer CATFLOW manual for more details. (default: constant)
- `depth`: Thickness of soil profile. (default: 2.1)

**Data:**
- `flow_accumulation`: Flow accumulation .tif file.
- `hillslopes`: .tif file for hillslopes.
- `elev2river`: .tif file for elevation to river.
- `dist2river`: .tif file for the distance to river.
- `filled_dem`: Filled Digital elevation model (DEM).
- `aspect`: .tif file for aspect.
- `river_id`: .tif file for river network.

### define_run_printouts

**Title:** Creates the printout times for the Model run

**Description:**
Writes a file with printout times for a CATFLOW simulation.
This file defines how frequently outputs are saved and dispalyed during the run.

**Parameters:**
- `start.time`: Start date for the simulation ("%d.%m.%Y %H:%M:%S")
- `end.time`: End date for the simulation ("%d.%m.%Y %H:%M:%S")
- `interval`: Time interval between printout times
- `time.unit`: Time units of printout times (Currently only tested for hours and seonds)
  - `hourly`
  - `seconds`
- `flag`: Flag controlling the amount of output at printout times, eventually a vector (1: dump all; 0: dump for surface nodes) (default: 1)

### write_multipliers

**Title:** Creates the multiplier files.

**Description:**
Writes the files ksmulto and thsmult.
Usually the factor is 1. but can be optionally used to have other scaling factors for saturated hydrualic conductivity and soil theta.

**Parameters:**
- `fac_kst`: Scaling factor for saturated hydrualic conductivity (default: 1)
- `fac_ths`: Scaling factor for soil theta (default: 1)

**Data:**
- `geometry`: .geo file for the hillslope geometry



