precursorHDF5
=============
[![OpenFOAM 2.3.1](https://img.shields.io/badge/OpenFOAM-2.3.1-blue)](https://openfoam.org/download/2-3-1-source/)
[![OpenFOAM 2.3.x](https://img.shields.io/badge/OpenFOAM-2.3.x-blue)](https://github.com/OpenFOAM/OpenFOAM-2.3.x)

**Please note**  
For now, the code is only valid in `v2.3.1`. If you want to adapt it in newer versions, please start from the native [timeVayingMappedFixedValue](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/fields/fvPatchFields/derived/timeVaryingMappedFixedValue/timeVaryingMappedFixedValueFvPatchField.H) bc, and add the new features of precursorHDF5 to it.

If you are using this bc in your work, please consider citing [this paper](https://doi.org/10.1080/10618562.2024.2370802) :

Zhang, Teng, Jinghua Li, Yingwen Yan, and Yuxin Fan. “Influence of LES Inflow Conditions on Simulations of a Piloted Diffusion Flame.” International Journal of Computational Fluid Dynamics, July 4, 2024, 1–15. https://doi.org/10.1080/10618562.2024.2370802.

This boundary condition is adapted from [timevaryingmappedhdf5fixedvalue](https://gitlab.com/chalmers-marine-technology/timevaryingmappedhdf5fixedvalue).

## Features
- Turbulent inlet velocity library is stored in **one** [**HDF5**](https://www.hdfgroup.org/HDF5/) file.
- No need to do the precursor simulation on the fly.
- Rescaling the turbulent inlet velocity based on experimental data or user-defined profiles.
- Support the **recycling usage** of the precursor's library.
- Other features are similar to the **timeVaryingMappedFixedValue**.

## Compatibility
- For **OpenFOAM v2.3.1** and **OpenFOAM-2.3.x**, use the `master` branch.
  ```
  git checkout master
  ```
- For **OpenFOAM v7**, use the `of7` branch.
  ```
  git checkout of7
  ```

## Schematic
![figure3](https://github.com/TimoLin/precursorHDF5/assets/7792396/4a9dcc5d-681f-452a-9a6a-c317e89a6b48)


## Prerequisite
### HDF5 library

Here is a brief tutorial of how to intall HDF5 lib:  
```sh
# Install path
cd $HOME/software
mkdir hdf5
cd hdf5
# Get hdf5-1.8.3 binary release
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.3/bin/linux-x86_64/5-1.8.3-linux-x86_64-shared.tar.gz
# Unpack 
tar xzf 5-1.8.3-linux-x86_64-shared.tar.gz
rm 5-1.8.3-linux-x86_64-shared.tar.gz
# Make a symbolic link "latest"
ln -s 5-1.8.3-linux-x86_64-shared latest
# Redeploy script
cd latest/bin
./h5redeploy
```
Set `hdf5` path in `.zshrc` or `.bashrc`:
```sh
export HDF5_DIR=$HOME/software/hdf5/latest
export LD_LIBRARY_PATH=$HDF5_DIR/lib:$LD_LIBRARY_PATH
```
### Install the precursorHDF5 bc
Clone this repo to your local drive and compile, for example:
```sh
git clone git@github.com:TimoLin/precursorHDF5.git $WM_PROJECT_USER_DIR/precursorHDF5
cd $WM_PROJECT_USER_DIR/precursorHDF5
wmake
```

## Workflow
### Generate the hdf5 file
Follow this Python script [tVMHDF5FV](https://github.com/TimoLin/pyScriptFoam/tree/master/inletTurb/tVMHDF5FV).  
### Use this bc
The boundary condition expects that the hdf5 file will contain three datasets.

One for the points, of shape N by 3, where N is the number of points you have in the sampling slice.
The first column contains the x coordinates, the second the y coordinates, and the third the z coordinates.

One for the time-values that the data is provided for, of shape nTimeValues by 1, where nTimeValues is simply the number of time-values that you are providing the data for.

One for the data itself, of shape nTimeValues by N by nComponents, where nComponents depends on what type of field one is dealing with.
For a vector field nComponents is 3, for instance.
The ordering of the data should be in agreement with the ordering of the points, i.e. the first value is expected to correspond to the first point etc.

One should provide the following parameters in the definition of the boundary type.

hdf5FileName -- name the hdf5 file.

hdf5PointsDatasetName -- name of the dataset containg the points.

hdf5SampleTimesDatasetName -- name of the dataset containing the sample times.

hdf5FieldValuesDatasetName -- name of the dataset containing the values of the field .

For example, one could have the following in the 0/U file.
```cpp
inlet
{
    type            precursorHDF5;
    setAverage      false;
    offset          (0 0 0);
    perturb         0.0;
    mapMethod       nearest; //or planarInterpolation;
    recycling       true;
    hdf5FileName    "dbTest.hdf5";
    hdf5PointsDatasetName    "points";
    hdf5SampleTimesDatasetName    "times";
    hdf5FieldValuesDatasetName    "velocity";
}
```

Note: For the *hdf5FileName*, the file is under the case's root folder, i.e.:
```sh
├── 0
├── constant   
├── system
└── dbTest.hdf5
```

Add the following line to `controlDict` when using it:
```cpp
libs ("libPrecursorHDF5.so");
```

## Current limitations
Currently the boundary condition only works with vector fields!
