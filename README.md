# timeVaryingMappedHDF5FixedValue #

This code is a modification of the timeVaryingMappedFixedValue boundary condition found in OpenFOAM.
This boundary condition allows to prescribe the value of the field by explicitly providing the values at a number of points.
The boundary condition will then use interpolation in order to obtain the values at the face centers.
Moreover, the boundary condition allows to interpolate between the prescribed values in time.

In the original implementation the data is put into constant/boundaryData/patchName/timeValue.
That is, for each time-value a separate directory is created.
When there is a need to provide data for thousands of different time-values this approach becomes very hard to deal with.

This modification makes the boundary condition read all the data from a single file, in HDF5 format.

### Installation ###

The main prerequisite is that you have some implementation of the hdf5 library present on your system.
The HDF5_DIR environmental variable, pointing to the installation directory of hdf5 should be available at compile time.

Nothing else should be required, just run wmake.

### Usage ###
The boundary condition expects that the hdf5 file will contain three datasets: one for the points, one for the time-values that the data is provided for and one for the data itself.
One should therefore provide the following parameters in the definition of the boundary type.

hdf5FileName -- name the hdf5 file

hdf5PointsDatasetName -- name of the dataset containg the points

hdf5SampleTimesDatasetName -- name of the dataset containing the sample times

hdf5FieldValuesDatasetName -- name of the dataset containing the values of the field 

For example, one could have the following in the 0/U file.

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