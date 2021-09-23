# CMAQ Volcano Emissions

---
    author: Barron H. Henderson
    date: 2021-09-23
---

This is a processor that generates Volcano SO2 from two types of OMI derived[1]
databases:

1. First, it can use the GEOS-Chem Volcano inputs to create CMAQ emissions. 
2. It can read the OMI Volcano eruption file directly.[1]

The GEOS-Chem processor follows the HEMCO extension. The HEMCO extension reads
rc files that were originally developed for GEOS-5. These rc files use the
latest and greatest Carn database[2] for both eruptive and degassing emissions.
Where year-specific data is not available, degassing is provided using the last
5-years as a climatology and eruptions are omitted.

The OMI Volcano eruption file can also be read directly. This is useful because
the GEOS-Chem inputs include climatalogical degassing, but can missing
eruptions that have occured rescently.

## How To

Open MakeVolcanoEmis.ipynb. Edit the dom and the m3dtmpl and g2dtmpl to point
to MCIP files with HT (g2dtmpl) and ZF (m3dtmpl). Then, run the cells. Two
types of output will be created. First, eruption only files from OMI L4.
Second,

## Temporal allocation

* The emissions are assumed to be constant for each UTC day.

## Vertical allocation

* Degassing and eruptions are emitted in the layer that contains the volcano
  altitude. This follows Itahashi et al. 2021 and allows for degassing at the
  volcano "surface" to be in a nonsurface cell. This is particularly useful
  for volcanos with peaks that are not representative of the grid cell surface
  level.
* The vertical allocation follows GEOS-5 by emitting in the top 1/3 of the
  plume cloud top height. I can imagine wanting to split some mass to the lower
  layers. 
* To identify the vertical layers, the layer heights are read in from MCIP
  meteorology. The meteorology could be day-specific, but I used a monthly
  ensemble. 

## Speciation

Currently, none is done. All emissions are stored as SO2_ERUPT or SO2_DEGAS.

## Uncertainties

Obviously, there are plenty of uncertainties. For example, the effective plume
heights are specified for eruptive volcanoes – but the allocation is a
simplistic. Also, at coarse resolutions its fine to put the whole volcano in
one grid cell – but what if the volcano is near a 4km grid cell edge… So I
expect finer scales to need empirical adjustment.

## Annotated Bibliography

```
.
|-- README.md
|   # This file.
|-- MakeVolcanoEmis.ipynb
|   # Driver for creating CMAQ Volcano emissions
|-- PlotVolcanoEmis.ipynb
|   # Plotting Jupyter Notebook
|-- volcano2ioapi.py
|   # Definition of classes that read rc or text files
|   # and allocate to vertical levels
|-- input
|   `-- MSVOLSO2L4_20210621.txt
|       # Eruptive volcanos
|-- mcip
|   |-- 108NHEMI2
|   `-- 27HI1
|       # Each domain you want to run will need MCIP files that contain
|       # a measure of height in meters. That will inclue HT from GRIDCRO2D
|       # and ZF or ZH from METCRO3D
`-- output
    |-- 108NHEMI2
    `-- 27HI1
        |-- so2_volcanic_emissions_Carns.%Y%m%d.27HI1.nc
        |   # Results include SO2_ERUPT and SO2_DEGAS
        |   # The user should confirm that SO2_ERUPT has non-zero
        |   # values on eruption days. If not, the values from 
        |   # MSVOLSO2L4_%Y-%m-%d.27HI1.nc can be merged.
        `-- MSVOLSO2L4_%Y-%m-%d.27HI1.nc
            # Eruptive emissions read directly from the text file
            # on OMI's website
```

## References

[1] https://so2.gsfc.nasa.gov/measures.html
[2] Carn, S. A., Yang, K., Prata, A. J., & Krotkov, N. A. (2015). Extending the long-term record of volcanic SO 2 emissions with the Ozone Mapping and Profiler Suite nadir mapper: OMPS volcanic SO2 measurements. Geophysical Research Letters, 42(3), 925–932. https://doi.org/10.1002/2014GL062437
[3] Itahashi, S., Mathur, R., Hogrefe, C., Napelenok, S. L., & Zhang, Y. (2021). Incorporation of volcanic SO&lt;sub&gt;2&lt;/sub&gt; emissions in the Hemispheric CMAQ (H-CMAQ) version 5.2 modeling system and assessing their impacts on sulfate aerosol over the Northern Hemisphere. Geoscientific Model Development, 14(9), 5751–5768. https://doi.org/10.5194/gmd-14-5751-2021