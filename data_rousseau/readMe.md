# Turbulent flows over steep rough permeable beds - Dataset, plots and scripts

This archive contains experimental data used to obtain the figures shown in the paper entitled *Experimental investigation into turbulent free-surface flows over a steep permeable bed by Rousseau G. & Ancey C.*

## CSV data

Classical CSV files are generated for each case (A1, A2, A3, A4, B1, B2, B3, B4, B5).  The following spatial averaged profiles are stored :

- ***Z_80*** : wall-normal position with the origin located at the porosity equal to 80%
- ***Porosity*** : porosity profile
- ***Ux [m/s]*** : stream-wise velocity
- ***Uz [m/s]*** : wall-normal velocity
- ***utx [m/s]*** : stream-wise turbulence intensity
- ***utz [m/s]*** : wall-normal turbulence intensity
- ***tautxz [N/m2]*** : turbulent stress
- ***upzs [m/s]*** : stream-wise disturbance intensity
- ***upzs [m/s]*** : wall-normal disturbance intensity
- ***tausxz [N/m2]***] : dispersive stress

## JSON database

A global database is stored into a `JSON` file including the above spatially averaged profiles in addition the flow parameters for each case (Reynolds numbers, mean diameter, optical conditions, etc...). `JSON` files can be imported in most computer language. `python` converts it into a *dictionary* while `Matlab` converts `JSON` into structures. Mathematica native `Import` function also supports the `JSON` file format. Python scripts that are used to generate the figures are shared in the `scripts_plot` folder and use the `JSON` file.

 For each case (A1, A2, A3, A4, B1, B2, B3, B4, B5) data structure in the `JSON` file is organized as follows:

- profiles
  - z_80 [m],
  - porosity
  - ux [m/s]
  - uz [m/s]
  - utx [m/s]
  - utz [m/s]
  - tautxz [N/m2]
  - udx [m/s]
    - - udz [m/s]
  - tausxz [N/m2]
- flow characteristics
  - H [m]
  - Ub [m/s]
  - porM
  - q [m2/s]
  - Us [m/s]
  - Uss [m/s]
  - uetp [m/s]
  - uet [m/s]
  - Re
  - Re0
  - Reh
  - Rep
  - ReRL
  - ReK
  - Fr
  - K_{kozeny-carman} [m2/s]
  - d_p [m]
  - i (slope in meter per meter)
  - uetM [m/s]
  - u* Voermans [m/s]
  - Roughness Crest index
  - Free Surface index
  - Origin index
  - Experimental parameters
    - slope in deg
    - small beads diameter [m]
    - large beads diameter [m]
    - flow per unit width [m2/s]
    - Flow  [L/s]
    - slope in %
- Optical system properties:
  - laser scan velocity [m/s]: 0.002
  - exposure time [s]: 0.0005
  - camera Type: Basler ac2040-180kc
  - frame per seconds: 420
  - focal: 35 mm
  - frame dimensions :  700 X 1496
  - number of frames collected: 6930

- Physico-Chemical properties:
  - Ethanol
  - Benzyl Alcohol
  - isoMix

The L cases shared by Joey J. Voermans and described in the article *Voermans, J. J., M. Ghisalberti, and G. N. Ivey. "The variation of flow and turbulence across the sediment–water interface." Journal of Fluid Mechanics 824 (2017): 413-437.* are also contained in the data set.

With the following keys for each cases (L11, L12, L13, L14, L15):

- profiles
  - z_80 [m] ('z_80' key is set for convenience but the origin of the profiles in *Voermans et al. (2017)* is established on the porosity inflection)
  - ux [m/s]
  - tautxz [N/m2]
  - tausxz [N/m2]
  - utx [m/s]
  - utz [m/s]
  - udx [m/s]
  - udz [m/s]
- uetV  [m/s]
- delta [m]
- deltab [m]

Warning from JJ Voermans on data by e-mail the 25/10/2018):
*Just a note on the data: in most cases there is a tiny discontinuity near the top of the boundary layer, which corresponds to the edge of the horizontally oriented camera (landscape). Some caution with using the L15 data where we used the maximum achievable bulk velocity with our experimental setup. We noticed in this experiment a thin layer at the surface that remained stagnant, causing an additional boundary layer from the surface downwards. This remained restricted to the upper part of the boundary layer (see for instance the Reynolds stress increasing again at the top), so lower proportion is still usable.*

## Import JSON database

On python just type the following script to import data

```python
import json
# make sure that you  type the correct path
f = open('full_path_to/data_Rousseau_JFM.json') 
data=json.load(f)
f.close()
# to display the possible keys type :
data['A1'].keys()
# to call the porosity profile for instance type
data['A1']['profiles']['Porosity']
```

## Python script for figures

In the `scripts_plots` folder you may successfully run the allFigures.py file by executing:

```bash
python allFigures.py
```

or for individual figures

```bash
python figureXX.py
```

Make sure that you have all dependencies (`matplotlib`,`json`, numpy, `opyf` [*optional*]).
Latex rendering is disabled by default, but you can uncomment the line in the allFigures.py file


## Raw video and codes to obtain double averaged profiles

An external DropBox folder can be downloaded [here](https://www.dropbox.com/sh/hf9l1h2e7zojjz5/AAAPvVqjaHIccjb8gTIoxkKua?dl=0) containing a raw video and the scripts to obtain the double averaged profiles using the PIVRIMS procedure.

## Obtaining the 3D velocity field data

The PIV-RIMS procedure collected first and second order flow statistics in a three-dimensional domain. The volume of data collected is substantially high (about 1Gb per case). If you are interested to obtain it, please contact us : gauthier.rousseau@gmail.com


