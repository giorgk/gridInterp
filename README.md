# gridInterp
Implementing gridded interpolation with data types common in groundwater modeling

# Introduction
This is a header only c++ library that interpolates gridded formatted data.

To use it just add the following include and you are good to go
```
#include "gridInterp.h"
```

All the methods and classes are defined into the `GRID_INTERP` namespace.
The main class is the `GRID_INTERP::interp`

One can defined 1-2-3D interpolation objects as follows:
```
GRID_INTERP::interp<1> oneD;
GRID_INTERP::interp<2> twoD;
GRID_INTERP::interp<3> threeD;
```
These are just empty containers. You can populate the with data using the related methods.

# Interpolation types
The library offers 2 interpolation modes and 2 interpolation methods
## Interpolation methods
1. **Linear** interpolates linearly between the values 
2. **Nearest** returns the value of the nearest coordinate

## Interpolation modes
1. **Point** Assumes that the values are defined on the coordinates
2. **Cell** Assumes that the values are defined at the cell centers and the coordinates correspond to the interface betwee nthe cells. In Cell mode the coordinates have to be the number of values + 1

## Layer interpolation (Not yet available)
This interpolation type is implemented for 3D only. 
In groundwater modelling it is common to have layers with varying elevation. 






## How to configure VSCODE for Simple C++
I hope the following link wont break soon.
https://code.visualstudio.com/docs/cpp/config-linux
