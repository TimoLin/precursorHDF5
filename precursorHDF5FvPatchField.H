/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::precursorHDF5FvPatchField

Group
    grpInletBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary conditions interpolates the values from a set of supplied
    points in space and time.  Supplied data should be specified in
    constant/boundaryData/\<patchname\> where:
    - points : pointField with locations
    - ddd    : supplied values at time ddd
    The default mode of operation (mapMethod planarInterpolation) is
    to project the points onto a plane (constructed from the first threee
    points) and construct a 2D triangulation and finds for the face centres
    the triangle it is in and the weights to the 3 vertices.

    The optional mapMethod nearest will avoid all projection and
    triangulation and just use the value at the nearest vertex.

    Values are interpolated linearly between times.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        setAverage   | flag to activate setting of average value | yes |
        perturb      | perturb points for regular geometries | no | 1e-5
        fieldTableName | alternative field name to sample | no| this field name
        mapMethod    | type of mapping | no | planarInterpolation
        recycling    | flag to activate recycling of inlet data library| yes/false|
        hdf5FileName | name the hdf5 file  | --
        hdf5PointsDatasetName      | name of the dataset containg the points | --
        hdf5SampleTimesDatasetName | name of the dataset containing the sample times | --
        hdf5FieldValuesDatasetName | name of the dataset containing the
                                     values of the field | --
    \endtable

    /verbatim
    myPatch
    {
        type            precursorHDF5;
        setAverage      false;
        //perturb       0.0;
        //fieldTableName samples;
        mapMethod       nearest;
        recycling       false;
        hdf5FileName    "db.hdf5";
        hdf5PointsDatasetName    "points";
        hdf5SampleTimesDatasetName    "times";
        hdf5FieldValuesDatasetName    "velocity";
        
    }
    /endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    precursorHDF5FvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef precursorHDF5FvPatchField_H
#define precursorHDF5FvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "FixedList.H"
#include "instantList.H"
#include "pointToPointPlanarInterpolation.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class precursorHDF5FvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class precursorHDF5FvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- If true adjust the mapped field to maintain average value
        bool setAverage_;

        //- Fraction of perturbation (fraction of bounding box) to add
        scalar perturb_;

        //- Interpolation scheme to use
        word mapMethod_;

        //- If true recycling using the inlet data library
        bool recycling_;

        //- 2D interpolation (for 'planarInterpolation' mapMethod)
        autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- List of boundaryData time directories
        instantList sampleTimes_;

        //- Current starting index in sampleTimes
        label startSampleTime_;

        //- Interpolated values from startSampleTime
        Field<Type> startSampledValues_;

        //- If setAverage: starting average value
        Type startAverage_;

        //- Current end index in sampleTimes
        label endSampleTime_;

        //- Interpolated values from endSampleTime
        Field<Type> endSampledValues_;

        //- If setAverage: end average value
        Type endAverage_;

        //- Time varying offset values to interpolated data
        autoPtr<Function1<Type>> offset_;

        //- Name of the hdf5 file
        word hdf5FileName_;

        //- Name of the hdf5 dataset to read the points form
        word hdf5PointsDatasetName_;

        //- Name of the hdf5 dataset to read the sample times from
        word hdf5SampleTimesDatasetName_;
        
        //- Name of the hdf5 dataset to read the field values from
        word hdf5FieldValuesDatasetName_;

public:

    //- Runtime type information
    TypeName("precursorHDF5");


    // Constructors

        //- Construct from patch and internal field
        precursorHDF5FvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        precursorHDF5FvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given precursorHDF5FvPatchField
        //  onto a new patch
        precursorHDF5FvPatchField
        (
            const precursorHDF5FvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        precursorHDF5FvPatchField
        (
            const precursorHDF5FvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new precursorHDF5FvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        precursorHDF5FvPatchField
        (
            const precursorHDF5FvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new precursorHDF5FvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return startSampledValues
            const Field<Type> startSampledValues()
            {
                 return startSampledValues_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Utility functions

            //- Find boundary data inbetween current time and interpolate
            void checkTable();


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "precursorHDF5FvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
