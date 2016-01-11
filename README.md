# timeVaryingMappedHDF5FixedValue #

This code is a modification of the timeVaryingMappedFixedValue boundary condition found in OpenFOAM.
This boundary condition allows to prescribe the value of the field by explicitly providing the values at a number of points.
The boundary condition will then use interpolation in order to obtain the values at the face centers.
Moreover, the boundary condition allows to interpolate between the prescribed values in time.

In the original implementation the data is put into constant/boundaryData/patchName/timeValue.
That is, for each time-value a separate directory is created.
When there is a need to provide data for thousands of different time-values this approach becomes very hard to deal with.

This modification makes the boundary condition read all the data from a single file, in HDF5 format.
More about this file format can be found [at the official web-site](https://www.hdfgroup.org/HDF5/).
The main advantage is that API for reading and writing HDF5 files exist for many languages (python, Matlab, C, Fortran...).

### Installation ###

The main prerequisite is that you have some implementation of the hdf5 library present on your system. 
The compilation goes well at least with version 1.8.13 and 1.8.14 of the library.
The HDF5_DIR environmental variable, pointing to the installation directory of hdf5 should be available at compile time.

Nothing else should be required, just run wmake.

The boundary condition has so far been only used in OpenFOAM version 2.3.1.

### Usage ###
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


```
#!c++

inlet
{
    type            timeVaryingMappedHDF5FixedValue;
    setAverage      false;
    offset          (0 0 0);
    perturb         0.0;
    hdf5FileName    "dbTest.hdf5";
    hdf5PointsDatasetName    "points";
    hdf5SampleTimesDatasetName    "time";
    hdf5FieldValuesDatasetName    "velocity";
}
```

### Current limitations ###
Currently the boundary condition only works with vector fields!