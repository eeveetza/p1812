# MATLAB/Octave Implementation of Recommendation ITU-R P.1812-6

This is a development branch that uses matrices to define simultaneous different path profiles in the function call. 

This code repository contains a MATLAB/Octave software implementation of Recommendation ITU-R P.1812 with a path-specific propagation prediction method for point-to-area terrestrial services in the frequency range 30 MHz to 6000 MHz.  

The reference version of this code (as approved by ITU-R Working Party 3K) is published by Study Group 3 on [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx).

The following table describes the structure of the folder `./matlab/` containing the MATLAB/Octave implementation of Recommendation ITU-R P.1812.

| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`tl_p1812_matr.m`                | MATLAB function implementing Recommendation ITU-R P.1812-6          |
|`validate_p1812_matr.m`          | MATLAB script used to validate the implementation of Recommendation ITU-R P.1812-6 in `tl_p1812_matr.m`             |
|`./validation_profiles/`    | Folder containing a proposed set of terrain profiles and inputs for validation of MATLAB implementation (or any other software implementation) of this Recommendation |
|`./validation_results/`	   | Folder containing all the results written during the transmission loss computations for the set of terrain profiles defined in the folder `./validation_profiles/` |
|`./private/`   |             Folder containing the functions called by `tl_p1812.m` and `validate_p1812.m`|

## Function Call

The function `tl_p1812_matr` can be called

1. by invoking only the required input arguments including latitude of path centre:
~~~
[Lb, Ep] = tl_p1812_matr(f, p, d, h, R, Ct, zone, htg, hrg, pol, 'phi_path', phi_path);
~~~

2. by invoking only the required input arguments including latitude/longitude of Tx/Rx:
~~~
[Lb,Ep] = tl_p1812_matr(f, p, d, h, R, Ct, zone, htg, hrg, pol, ...
                        'phi_t',phi_t,'phi_r',phi_r,'lam_t',lam_t,'lam_r',lam_r);
~~~
3. or by invoking both required and optional parameters. Optional parameters are specified using comma separated pairs of Name, Value arguments in any order:
~~~
[Lb, Ep] = tl_p1812_matr(f, p, d, h, R, Ct, zone, htg, hrg, pol, 'phi_path', phi_path, 'DN', DN, 'N0', N0);
~~~

## Required input arguments of function `tl_p1812_matr`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `f`               | scalar double | GHz   | 0.03 ≤ `f` ≤ 6 | Frequency   |
| `p         `      | scalar double | %     | 1 ≤ `p` ≤ 50 | Time percentage for which the calculated basic transmission loss is not exceeded |
| `d`               | matrix double MxN| km    | ~0.25 ≤ `max(d)` ≤ ~3000 | Terrain profile distances (in the ascending order from the transmitter)|
| `h`          | matrix double MxN| m (asl)   |   | Terrain profile heights |
| `R`           | matrix double MxN   | m      |              |  Representative clutter heights |
| `Ct`           | matrix int  MxN  |       |  1 - Water/sea, 2 - Open/rural, 3 - Suburban, 4 - Urban/trees/forest, 5 - Dense urban             |  Array of representative clutter types. If empty or all zeros, the default clutter type used is Open/rural |
| `zone`           | matrix int   MxN |       | 1 - Sea, 3 - Coastal land, 4 - Inland             |  Radio-climatic zone types |
| `htg`           | scalar double    | m      |   1 ≤ `htg`  ≤ 3000          |  Tx antenna height above ground level |
| `hrg`           | scalar double    | m      |   1 ≤ `hrg`  ≤ 3000          |  Rx antenna height above ground level |
| `pol`           | scalar int    |       |   `pol`  = 1, 2          |  Polarization of the signal: 1 - horizontal, 2 - vertical |

`d`, `h`, `R`, `Ct`, and `zone` can be matrices with as many rows M as there are profiles (each row containing information about one profile). 
M is the number of paths with N profile points. Restriction: Each path has the same number of profile points N.

The required parameters related to path center can be provided in two ways:

* The latitude of the path center is provided either directly

 Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `phi_path`           | array double Mx1   | deg      |   -80 ≤ `phi_path`  ≤ 80          |  Latitude of the path centre |

* or the latitudes and longitudes of the Tx and Rx are provided, in which case the path center will be computed by `tl_p1812_matrix`. 

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `phi_t`           | array double Mx1    | deg      |   -80 ≤ `phi_t`  ≤ 80          |  Latitude of Tx station |
| `phi_r`           | array double Mx1   | deg      |   -80 ≤ `phi_r`  ≤ 80          |  Latitude of Rx station |
| `lam_t`           | array double Mx1   | deg      |   -180 ≤ `lam_t`  ≤ 180          |  Longitude of Tx station |
| `lam_r`           | array double Mx1   | deg      |   -180 ≤ `lam_r`  ≤ 180          |  Longitude of Rx station |


 

## Optional input arguments of function `tl_p1812`
| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `pL`           | scalar double    | %      |   1 ≤ `pL`  ≤ 99          |  Location percentage for which the calculated basic tranmission loss is not exceeded. Default is 50%. |
| `sigmaL`           | scalar double    | dB      |             |  location variability standard deviations computed using stdDev.m according to §4.8 and §4.10; the value of 5.5 dB used for planning Broadcasting DTT; Default: 0 dB |
| `Ptx`           | scalar double    | kW      |   `Ptx` > 0          |  Tx power, default value is 1 kW |
| `DN`            | array double Mx1   | N-units/km      | `DN`> 0           | The average radio-refractive index lapse-rate through the lowest 1 km of the atmosphere at the path-center. It can be derived from an appropriate map. Default: 45. |
| `N0`           | array double   Mx1 | N-units      |             | The sea-level surface refractivity at the path-centre. It can be derived from an appropriate map. Default: 325.|
| `dct`           | array double Mx1   | km      |   `dct` ≥ 0          |  Distance over land from the Tx antenna to the coast along the great-circle interference path. Default: 500 km. Set to zero for a terminal on a ship or sea platform|
| `dcr`           | array double  Mx1  | km      |   `dcr` ≥ 0          |  Distance over land from the Rx antenna to the coast along the great-circle interference path. Default: 500 km. Set to zero for a terminal on a ship or sea platform|
| `flag4`           | scalar int    |       |             |  If `flag4`= 1, the alternative method from Attachment 4 to Annex 1 is used to calculate `Lbulls` without using terrain profile. Default: 0 |
| `debug`           | scalar int    |       |             |  If `debug`= 1, the results are written in log files. Default: 0. |
| `fid_log`           | scalar int    |       |     Only used if `debug`= 1        |  File identifier of the log file opened for writing outside the function. If not provided, a default file with a filename containing a timestamp will be created. |


 
## Outputs ##

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `Lb`    | array double Mx1 | dB    | Basic transmission loss |
| `Ep`    | array double Mx1 | dB(uV/m)    | Electric field strength |

## Software Versions
The code was tested and runs on:
* MATLAB versions 2017a and 2020a
* Octave version 6.1.0

## References

* [Recommendation ITU-R P.1812](https://www.itu.int/rec/R-REC-P.1812/en)

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)
