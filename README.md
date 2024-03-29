precursorHDF5
=============
[![OpenFOAM 2.3.1](https://img.shields.io/badge/OpenFOAM-2.3.1-blue)](https://openfoam.org/download/2-3-1-source/)
[![OpenFOAM 2.3.x](https://img.shields.io/badge/OpenFOAM-2.3.x-blue)](https://github.com/OpenFOAM/OpenFOAM-2.3.x)


This boundary condition is adapted from [timevaryingmappedhdf5fixedvalue](https://gitlab.com/chalmers-marine-technology/timevaryingmappedhdf5fixedvalue).

This modification makes the boundary condition read all the data from a single file, in HDF5 format.More about this file format can be found [at the official web-site](https://www.hdfgroup.org/HDF5/).

A recycling usage of the precursor's library is introduced. This can be triggered using `recycling   true;` in the bc file.

## Prerequisite
### HDF5 library

Intallation note:  
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
### Install
Clone this repo to your local drive and compile, say: 
```sh
git clone git@github.com:TimoLin/precursorHDF5.git $WM_PROJECT_USER_DIR/precursorHDF5
cd $WM_PROJECT_USER_DIR/precursorHDF5
wmake
```

## Usage
### Generate the hdf5 file
Follow this Python script [tVMHDF5FV](https://github.com/TimoLin/pyScriptFoam/tree/master/inletTurb/tVMHDF5FV).  
### Use this bc
The boundary condition expects that the hdf5 file will contain three datasets.

One for the points, of shape N by 3, where N is the number of points you have available.
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
    mapMethod       nearest; //planarInterpolation;
    recycling       true;
    hdf5FileName    "dbTest.hdf5";
    hdf5PointsDatasetName    "points";
    hdf5SampleTimesDatasetName    "times";
    hdf5FieldValuesDatasetName    "velocity";
}
```

Add the following line to `controlDict` when using it:
```cpp
libs ("libPrecursorHDF5.so");
```

## Current limitations
Currently the boundary condition only works with vector fields!
