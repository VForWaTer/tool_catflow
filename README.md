# CATFLOW preprocessing tool

This is an internal development branch for Manoj J et al (2026).

CATFLOW is a physically-based hydrological model for 2‑D hillslopes. The basic modeling unit is a 2‑D hillslope discretized by curvilinear orthogonal coordinates (vertical and downslope); the third dimension is represented via a variable width perpendicular to the slope line. Soil water dynamics are simulated using the Richards equation solved with an implicit mass-conservative Picard iteration (Celia et al., 1990). The model can simulate unsaturated and saturated subsurface flow (no separate groundwater routine), various soil hydraulic parametrizations (e.g., van Genuchten–Mualem), diffusion-wave overland flow (Saint‑Venant), infiltration and saturation excess runoff, re‑infiltration, lateral subsurface flow and return flow. Multiple hillslopes can be connected by a river network for catchment modeling. CATFLOW does not simulate snow or frozen soil. The model has been used in many studies (e.g., Zehe et al., 2005, 2010, 2014; Wienhöfer & Zehe, 2014).

This repository provides containerised preprocessing workflows to create CATFLOW input (representative hillslope geometry and auxiliary files) from raster and optional vector data. It follows the VForWaTer tool-specs for reusable tools.

## Tools

### make_representative_hillslope
Create a CATFLOW geometry file from raster inputs using the representative‑hillslope method (Loritz et al., 2017).

Parameters
- hillslope_id (integer): ID in hillslopes.tif to process; use -1 to process the whole basin. (default: -1)
- no_flow_area (float): Fraction of upper hillslope considered “no flow” when truncating profile. (default: 0.30)
- min_cells (integer): Minimum number of unique distance values required to build a representative profile. (default: 10)
- hill_type (enum): 'constant' | 'cake' | 'variable' (default: constant)
- depth (float): Profile thickness (m). (default: 2.1)
- constant_width (boolean): Force a constant width for the generated geometry. (default: true)
- min_area (integer): Minimum area (m^2) for a hillslope to be processed. (default: 10000)
- freedom (integer): Spline freedom parameter used for smoothing. (default: 10)

Data (inputs)
- flow_accumulation (file, required): Flow accumulation raster (.tif)
- hillslopes (file, required): Hillslope raster (.tif)
- elev2river (file, required): Elevation‑to‑river raster (.tif)
- dist2river (file, required): Distance‑to‑river raster (.tif)
- filled_dem (file, required): Filled DEM (.tif)
- aspect (file, required): Aspect raster (.tif)
- river_id (file, required): River/stream id raster (.tif)
- soil (file, optional): Soil categorical raster (.tif) — used to generate soil.dat
- landuse (file, optional): Land‑use categorical raster (.tif) — used to compute land‑use distributions
- hillslopes_vect (file, optional): Vector GeoPackage (.gpkg) with hillslope attributes (e.g., new precipitation IDs)

Outputs (examples)
- /out/rep_hill.geo — CATFLOW geometry file
- /out/geom.Rds — geometry object for downstream tools
- /out/CATFLOW/in/soil/soil.dat — soil bands (if soil provided)
- /out/CATFLOW/in/landuse/surface.pob — surface nodes / landuse mapping
- plots in /out/plots/ (energy distribution, geometry, landuse distributions, etc.)

### define_run_printouts
Create CATFLOW printout times.

Parameters
- start.time (string): "%d.%m.%Y %H:%M:%S"
- end.time (string): "%d.%m.%Y %H:%M:%S"
- interval (integer): print interval
- time.unit (enum): 'hourly' | 'seconds'
- flag (integer): output flag (1: full, 0: surface nodes, ...)

## Usage

1. Place input files in `/in`.
2. Edit `in/input.json` to specify parameters and data.
3. Build and run the Docker image (or run the script in an R environment with required packages).

Example input.json (make_representative_hillslope)
```json
{
  "make_representative_hillslope": {
    "parameters": {
      "hillslope_id": 137,
      "no_flow_area": 0.30,
      "min_cells": 10,
      "hill_type": "constant",
      "depth": 2.1,
      "min_area": 10000,
      "freedom": 10,
      "constant_width": false
    },
    "data": {
      "flow_accumulation": "/in/flow_accumulation.tif",
      "hillslopes": "/in/hillslopes.tif",
      "hillslopes_vect": "/in/hillslopes.gpkg",
      "elev2river": "/in/elevation.tif",
      "dist2river": "/in/distance.tif",
      "filled_dem": "/in/fill_DEM.tif",
      "aspect": "/in/aspect.tif",
      "river_id": "/in/streams.tif",
      "soil": "/in/soils.tif",
      "landuse": "/in/reproject_lulc_2016.tif"
    }
  }
}
```

## Running with Docker

Build:
```sh
docker build -t catflow .
```

Run (example):
```sh
docker run --rm -v "${PWD}/local/in:/in" -v "${PWD}/local/out:/out" -e TOOL_RUN=make_representative_hillslope catflow
```

## Development
- Code in `/src`.
- Update `tool.yml` and `in/input.json` when changing tool inputs.
- Tests and CI configured in the repository.

## License
MIT — see LICENSE.

## References
- Zehe et al., 2005, 2010, 2014
- Wienhöfer & Zehe, 2014
- Celia et al., 1990
- Loritz et al., 2017
- Francke et al., 2006
- Cochrane & Flanagan, 2003
- CATFLOW R Package: https://github.com/CATFLOW/Catflow-R-Package



