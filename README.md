# CATFLOW Preprocessing Tool

[![Docker Image CI](https://github.com/VForWaTer/tool_catflow/actions/workflows/docker-image.yml/badge.svg)](https://github.com/VForWaTer/tool_catflow/actions/workflows/docker-image.yml)
[![DOI](https://zenodo.org/badge/610682357.svg)](https://zenodo.org/badge/latestdoi/610682357)

**About CATFLOW:**  
CATFLOW is a physically-based hydrological model for simulating water flow and transport in hillslopes and catchments. The basic modeling unit is a 2-D hillslope, discretized by curvilinear orthogonal coordinates in the vertical and downslope directions; the third dimension is represented via a variable width of the slope perpendicular to the slope line at each node. Soil water dynamics are simulated based on the Richards equation in the pressure-based form and numerically solved using an implicit mass conservative Picard iteration (Celia et al., 1990). The model can simulate unsaturated and saturated subsurface flow and hence has no separate groundwater routine. Soil hydraulic functions following van Genuchten–Mualem are commonly used, though several other parametrizations are possible. Overland flow is simulated using the diffusion wave approximation of the Saint-Venant equation and explicit upstreaming. The hillslope module can simulate infiltration excess runoff, saturation excess runoff, re-infiltration of surface runoff, lateral water flow in the subsurface as well as return flow. For catchment modeling several hillslopes can be interconnected by a river network for collecting and routing their runoff contributions, i.e., surface runoff or subsurface flow leaving the hillslope, to the catchment outlet.
The CATFLOW model has been successfully used and specified in numerous studies ([Zehe et al. 2001](https://doi.org/10.1016/S0378-3774(99)00083-9), [Loritz et al. 2017](https://doi.org/10.5194/hess-21-1225-2017) and [Manoj J et al. 2024](https://doi.org/10.1029/2023WR036420))

This tool provides containerized workflows for generating CATFLOW input files, especially representative hillslope geometries, from raster and vector geospatial data.  
It follows the [Tool Specification](https://vforwater.github.io/tool-specs/) for reusable research software using Docker.

**Docker Image Updates:**  
You can always pull the latest released version of this tool from the GitHub Container Registry using:
```sh
docker pull ghcr.io/vforwater/tool_catflow:latest
```
or rebuild locally to get the newest features and fixes.

---

## Table of Contents

- [Tools](#tools)
- [Usage](#usage)
- [Input and Output](#input-and-output)
- [Running with Docker](#running-with-docker)
- [Development](#development)
- [License](#license)
- [References](#references)

---

## Tools

### make_representative_hillslope

**Description:**  
Creates a CATFLOW geometry file from .tif files using the method of representative hillslopes [Loritz et al. 2017](https://doi.org/10.5194/hess-21-1225-2017) and [Manoj J et al. 2024](https://doi.org/10.1029/2023WR036420). Delineates 2D catenas from hillslopes based on raster and vector input data. 

**Parameters:**
- `hillslope_id` (integer): Integer ID of the hillslope from hillslopes.tif for calculating the geometry. If entire basin is to be used for the representative hillslope, give -1.
- `no_flow_area` (float): Percentage of no flow area with almost no slope within the area of interest. (default: 0.30)
- `min_cells` (integer): Minimum number of unique rounded distance values to be considered for the hillslope geometry. (default: 10)
- `hill_type` (enum): Hillslope type. Options: `constant`, `cake`, `variable`. (default: constant)
- `depth` (float): Thickness of soil profile. (default: 2.1)
- `constant_width` (boolean): If true, use a constant width for the hillslope geometry. If false, use varying width. (default: true)
- `min_area` (integer): Minimum area (in square meters) required for a hillslope to be considered valid. (default: 10000)
- `freedom` (integer): Degree of freedom for the spline function used in hillslope geometry calculations. (default: 10)

**Data:**
- `flow_accumulation` (file): Flow accumulation .tif file.
- `hillslopes` (file): .tif file for hillslopes.
- `elev2river` (file): .tif file for elevation to river.
- `dist2river` (file): .tif file for the distance to river.
- `filled_dem` (file): Filled Digital elevation model (DEM).
- `aspect` (file): .tif file for aspect.
- `river_id` (file): .tif file for river network.
- `soil` (file, optional): .tif file for soil properties.

---


## Usage

1. **Prepare Input Files:**  
   Place your input files (e.g., raster TIFFs) in the `/in` directory.

2. **Configure Input JSON:**  
   Edit the `input.json` file to specify parameters and data for the desired tool.

   Example for hillslope geometry:
   ```json
   {
     "make_representative_hillslope": {
       "parameters": {
         "hillslope_id": 223,
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
         "elev2river": "/in/elevation.tif",
         "dist2river": "/in/distance.tif",
         "filled_dem": "/in/fill_DEM.tif",
         "aspect": "/in/aspect.tif",
         "river_id": "/in/streams.tif",
         "soil": "/in/soils.tif"
       }
     }
   }
   ```

3. **Install Docker and Build the Container:**

   - **Install Docker:**  
     Download and install Docker Desktop from [https://www.docker.com/products/docker-desktop/](https://www.docker.com/products/docker-desktop/) and follow the installation instructions for your operating system.

   - **Build the Docker image:**  
     Open a terminal or command prompt in the project directory and run:
     ```sh
     docker build -t catflow . 
     ```

4. **Run with Docker:**  
   Use the following command to run the tool (replace `catflow` with your Docker image name):

   **Note:**  
   Make sure to mount your local input and output directories to the container's `/in` and `/out` directories using the `-v` option.  
   For example, if your local input files are in `local/in` and you want outputs in `local/out`, use:

   For **PowerShell**:
   ```powershell
   docker run --rm -it `
     -v "${PWD}/local/in:/in" `
     -v "${PWD}/local/out:/out" `
     -e TOOL_RUN=make_representative_hillslope `
     catflow
   ```

   For **Command Prompt**:
   ```cmd
   docker run --rm -it ^
     -v "%cd%/local/in:/in" ^
     -v "%cd%/local/out:/out" ^
     -e TOOL_RUN=make_representative_hillslope ^
     catflow
   ```

   *(Adjust the paths as needed for your environment.)*

---

## Input and Output

- **Input Directory (`/in`):**  
  Place all required input files here.  
  Example: `/in/flow_accumulation.tif`, `/in/hillslopes.tif`, etc.

- **Output Directory (`/out`):**  
  All generated files will be saved here.  
  Example: `/out/rep_hill.geo`, `/out/plots/energy_distribution.pdf`, etc.

- **Input JSON Example:**
  ```json
  {
    "make_representative_hillslope": {
      "parameters": {
        "hillslope_id": 223,
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
        "elev2river": "/in/elevation.tif",
        "dist2river": "/in/distance.tif",
        "filled_dem": "/in/fill_DEM.tif",
        "aspect": "/in/aspect.tif",
        "river_id": "/in/streams.tif",
        "soil": "/in/soils.tif"
      }
    }
  }
  ```

---

## Running with Docker

- **Build the Docker image:**
  ```sh
  docker build -t catflow .
  ```

- **Run the container:**
  (See [Usage](#usage) for platform-specific examples.)

---

## Development

- Source code is located in the `/src` directory.
- To add new tools or modify existing ones, edit the R files in `/src` and update `tool.yml` as needed.
- Contributions are welcome! Please open issues or pull requests on GitHub.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## References

- [CATFLOW R Package](https://github.com/CATFLOW/Catflow-R-Package)
- [VForWaTer Tool Specification](https://github.com/VForWaTer/tool-specs)
- Zehe, E., Maurer, T., Ihringer, J., and Plate, E.: Modeling water flow and mass transport in a loess catchment, Phys. Chem. Earth, Part B Hydrol. Ocean. Atmos., 26, 487–507, https://doi.org/10.1016/S1464-1909(01)00041-7, 2001.

