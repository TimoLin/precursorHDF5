/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "precursorHDF5FvPatchField.H"
#include "Time.H"
#include "AverageIOField.H"

#include "hdf5.h"
#include "hdf5_hl.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
precursorHDF5FvPatchField<Type>::
precursorHDF5FvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    fieldTableName_(iF.name()),
    setAverage_(false),
    perturb_(0),
    recycling_(false),
    mapperPtr_(NULL),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(pTraits<Type>::zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(pTraits<Type>::zero),
    offset_(),
    hdf5FileName_("dummy.hdf5"),
    hdf5PointsDatasetName_("points"),
    hdf5SampleTimesDatasetName_("times"),
    hdf5FieldValuesDatasetName_("velocity")
{}


template<class Type>
precursorHDF5FvPatchField<Type>::
precursorHDF5FvPatchField
(
    const precursorHDF5FvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    recycling_(ptf.recycling_),
    mapperPtr_(NULL),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(pTraits<Type>::zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(pTraits<Type>::zero),
    offset_
    (
        ptf.offset_.valid()
      ? ptf.offset_().clone().ptr()
      : NULL
    ),
    hdf5FileName_(ptf.hdf5FileName_),
    hdf5PointsDatasetName_(ptf.hdf5PointsDatasetName_),
    hdf5SampleTimesDatasetName_(ptf.hdf5SampleTimesDatasetName_),
    hdf5FieldValuesDatasetName_(ptf.hdf5FieldValuesDatasetName_)
{}


template<class Type>
precursorHDF5FvPatchField<Type>::
precursorHDF5FvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    fieldTableName_(iF.name()),
    setAverage_(readBool(dict.lookup("setAverage"))),
    perturb_(dict.lookupOrDefault("perturb", 1e-5)),
    mapMethod_
    (
        dict.lookupOrDefault<word>
        (
            "mapMethod",
            "planarInterpolation"
        )
    ),
    recycling_(readBool(dict.lookup("recycling"))),
    mapperPtr_(NULL),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(pTraits<Type>::zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(pTraits<Type>::zero),
    offset_(DataEntry<Type>::New("offset", dict))
{

    if
    (
        mapMethod_ != "planarInterpolation"
     && mapMethod_ != "nearest"
    )
    {
        FatalIOErrorIn
        (
            "precursorHDF5FvPatchField<Type>::\n"
            "precursorHDF5FvPatchField\n"
            "(\n"
            "    const fvPatch&\n"
            "    const DimensionedField<Type, volMesh>&\n"
            "    const dictionary&\n"
            ")\n",
            dict
        )   << "mapMethod should be one of 'planarInterpolation'"
            << ", 'nearest'" << exit(FatalIOError);
    }


    dict.readIfPresent("fieldTableName", fieldTableName_);
    dict.readIfPresent("hdf5FileName", hdf5FileName_);
    dict.readIfPresent("hdf5PointsDatasetName", hdf5PointsDatasetName_);
    dict.readIfPresent("hdf5SampleTimesDatasetName", hdf5SampleTimesDatasetName_);
    dict.readIfPresent("hdf5FieldValuesDatasetName", hdf5FieldValuesDatasetName_);

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
    }
    else
    {
        // Note: we use evaluate() here to trigger updateCoeffs followed
        //       by re-setting of fvatchfield::updated_ flag. This is
        //       so if first use is in the next time step it retriggers
        //       a new update.
        this->evaluate(Pstream::blocking);
    }
}


template<class Type>
precursorHDF5FvPatchField<Type>::
precursorHDF5FvPatchField
(
    const precursorHDF5FvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    recycling_(ptf.recycling_),
    mapperPtr_(NULL),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    startAverage_(ptf.startAverage_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    endAverage_(ptf.endAverage_),
    offset_
    (
        ptf.offset_.valid()
      ? ptf.offset_().clone().ptr()
      : NULL
    ),
    hdf5FileName_(ptf.hdf5FileName_),
    hdf5PointsDatasetName_(ptf.hdf5PointsDatasetName_),
    hdf5SampleTimesDatasetName_(ptf.hdf5SampleTimesDatasetName_),
    hdf5FieldValuesDatasetName_(ptf.hdf5FieldValuesDatasetName_)
{}


template<class Type>
precursorHDF5FvPatchField<Type>::
precursorHDF5FvPatchField
(
    const precursorHDF5FvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(NULL),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    startAverage_(ptf.startAverage_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    endAverage_(ptf.endAverage_),
    offset_
    (
        ptf.offset_.valid()
      ? ptf.offset_().clone().ptr()
      : NULL
    ),
    hdf5FileName_(ptf.hdf5FileName_),
    hdf5PointsDatasetName_(ptf.hdf5PointsDatasetName_),
    hdf5SampleTimesDatasetName_(ptf.hdf5SampleTimesDatasetName_),
    hdf5FieldValuesDatasetName_(ptf.hdf5FieldValuesDatasetName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void precursorHDF5FvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    if (startSampledValues_.size())
    {
        startSampledValues_.autoMap(m);
        endSampledValues_.autoMap(m);
    }
    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void precursorHDF5FvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const precursorHDF5FvPatchField<Type>& tiptf =
        refCast<const precursorHDF5FvPatchField<Type> >(ptf);

    startSampledValues_.rmap(tiptf.startSampledValues_, addr);
    endSampledValues_.rmap(tiptf.endSampledValues_, addr);

    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void precursorHDF5FvPatchField<Type>::checkTable()
{
    hid_t hdf5DataBase = H5Fopen(hdf5FileName_.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT);
    herr_t hdf5Status;

    // Initialise
    if (mapperPtr_.empty())
    {

        hsize_t nPoints[2];

        // Get the size of the points dataset
        H5LTget_dataset_info(hdf5DataBase,
                             hdf5PointsDatasetName_.c_str(),
                             nPoints, 
                             NULL,
                             NULL);

        // Create an array to read in the points
        double points[nPoints[0]][nPoints[1]];

        // Open the points dataset
        hid_t pointsDataset = H5Dopen(hdf5DataBase,
                                      hdf5PointsDatasetName_.c_str(),
                                      H5P_DEFAULT);

        // Read the points dataset
        hdf5Status = H5Dread(pointsDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                H5P_DEFAULT, points);

        if (hdf5Status < 0)
        {
            Info << "HDF5 error when opening points dataset" << endl;
        }

        // Close the points dataset
        hdf5Status = H5Dclose(pointsDataset);

        if (hdf5Status < 0)
        {
            Info << "HDF5 error when closing points dataset" << endl;
        }

        // Create the samplePoints field
        pointField samplePoints(nPoints[0]);

        forAll(samplePoints, pointI)
        {
            samplePoints[pointI][0] = points[pointI][0];
            samplePoints[pointI][1] = points[pointI][1];
            samplePoints[pointI][2] = points[pointI][2];
            //Info<<pointI<<":"<<samplePoints[pointI][0]<<","<<samplePoints[pointI][1]<<","<<samplePoints[pointI][2]<<endl;
        }


        Info << "precursorHDF5FvPatchField :"
             << " Read " << samplePoints.size() << " points" << endl;

        if (debug)
        {
            Info << samplePoints << endl;
        }

        // tbd: run-time selection
        bool nearestOnly =
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                 this->patch().patch().faceCentres(),
                perturb_,
                nearestOnly
            )
        );

        // Read the times for which data is available

        hsize_t     nTimes[2];
        hdf5Status = H5LTget_dataset_info(hdf5DataBase,
                                          hdf5SampleTimesDatasetName_.c_str(),
                                          nTimes,
                                          NULL,
                                          NULL);

        if (hdf5Status < 0)
        {
            Info << "HDF5 error getting info on time dataset" << endl;
        }

        scalar times[nTimes[0]];

        hdf5Status = H5LTread_dataset_double(hdf5DataBase,
                                             hdf5SampleTimesDatasetName_.c_str(),
                                             times);

        if (hdf5Status < 0)
        {
            Info << "HDF5 error reading time dataset" << endl;
        }

        sampleTimes_.resize(nTimes[0]);

        forAll(sampleTimes_, timeI)
        {
            sampleTimes_[timeI] = instant(times[timeI]);
        }

        if(false)
        {
            Info << sampleTimes_ << endl;
        }
    }

    // Find current time in sampleTimes
    label lo = -1;
    label hi = -1;
    
    scalar period = sampleTimes_[sampleTimes_.size()-1].value()-sampleTimes_[0].value();
    
    scalar dbTime = this->db().time().value();
    
    // Take remainder
    scalar seekTime;
    //Info<<"Test:"<<recycling_<<endl;
    if (recycling_)
    {
        seekTime= dbTime - int(dbTime/period)*period;
    }
    else
    {
        seekTime = dbTime;
    }

    if (seekTime < this->db().time().deltaTValue())
    {   
        // Avoid numerical error
        seekTime = 0.0;

        startSampleTime_ = -1;
        endSampleTime_ = -1;
    }

    bool foundTime = mapperPtr_().findTime
    (
        sampleTimes_,
        startSampleTime_,
        seekTime,//this->db().time().value(),
        lo,
        hi
    );
    

    if (true)
    {
        Info<<" Debug Recycling:"
            <<" Current time:"<<dbTime
            <<" Current deltaT:"<<this->db().time().deltaTValue()
            <<",seekTimeeekTime:"<<seekTime
            <<",found:"<<foundTime
            <<",lo:"<<lo<<",hi:"<<hi
            <<",startSampleTime:" <<startSampleTime_
            <<",endSampleTime:"<<endSampleTime_<<"\n"; 
    }

    if (!foundTime)
    {
        FatalErrorIn
        (
            "precursorHDF5FvPatchField<Type>::checkTable()"
        )   << "Cannot find starting sampling values for current time "
            << this->db().time().value() << nl
            << "Have sampling values for times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "In directory "
            <<  this->db().time().constant()/"boundaryData"/this->patch().name()
            << "\n    on patch " << this->patch().name()
            << " of field " << fieldTableName_
            << exit(FatalError);
    }


    // Update sampled data fields.

    hsize_t nVelocity[3];

    // Get the size of the velocity dataset
    hdf5Status = H5LTget_dataset_info(hdf5DataBase,
                                      hdf5FieldValuesDatasetName_.c_str(),
                                      nVelocity, 
                                      NULL,
                                      NULL);

    if (hdf5Status < 0)
    {
        Info << "HDF5 error gettng the dimensions of the velocity dataset" << endl;
    }

    // Create an array to read in the points
    double velocity[nVelocity[1]][nVelocity[2]];

    // Open the velocity dataset
    hid_t velocityDataset = H5Dopen(hdf5DataBase,
                                    hdf5FieldValuesDatasetName_.c_str(),
                                    H5P_DEFAULT);

    // Selection parameters
    hsize_t offset[3]; //Only the first one is used
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;

    hsize_t count[3]; // Time, Points, Dims
    count[0] = 1;
    count[1] = nVelocity[1];
    count[2] = nVelocity[2];

    hsize_t nVelocitySlice[2];
    nVelocitySlice[0] = nVelocity[1];
    nVelocitySlice[1] = nVelocity[2];

    // Create a new simple dataspace by the given rank and dims
    hid_t velocityMemspace = H5Screate_simple(2, nVelocitySlice, NULL);  

    // Return a copy of the dataspace of the given dataset
    hid_t velocityDataspace = H5Dget_space(velocityDataset);

    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;

        if (startSampleTime_ == endSampleTime_)
        {
            // No need to reread since are end values
            if (true)
            {
                Pout<< "checkTable : Setting startValues to (already read) "
                    << sampleTimes_[startSampleTime_].name()
                    << endl;
            }
            startSampledValues_ = endSampledValues_;
            startAverage_ = endAverage_;
        }
        else
        {
            if (true)
            {
                Pout<< "checkTable : Reading startValues from "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[lo].name()
                    << endl;
            }
            
            
            offset[0] = lo;

            hdf5Status = H5Sselect_hyperslab(velocityDataspace, H5S_SELECT_SET, offset, NULL,
                                             count, NULL);

            if (hdf5Status < 0)
            {
                Info << "HDF5 error selecting velocity slice." << endl;
            }
            
            // Read in the slice
            hdf5Status = H5Dread(velocityDataset, H5T_NATIVE_DOUBLE, velocityMemspace, velocityDataspace,
                                 H5P_DEFAULT, velocity);


            if (hdf5Status < 0)
            {
                Info << "HDF5 error reading velocity slice" << endl;
            }

            if(false)
            {
                for(int i=0; i<nVelocity[0]*0.005; i++)
                {
                    Info << scalar(velocity[i][0]) << endl;
                }
            }

            // Reread values and interpolate
            AverageIOField<Type> vals
            (
                IOobject
                (
                   
                    "dummy", //fieldTableName_,
                    this->db().time().constant(),
                    "boundaryData"
                   /this->patch().name()
                   /sampleTimes_[startSampleTime_].name(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                    //false
                ),
                nVelocity[1]
            );

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorIn
                (
                    "precursorHDF5FvPatchField<Type>::"
                    "checkTable()"
                )   << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << vals.objectPath() << exit(FatalError);
            }
            
            if (this->db().template 
                    foundObject<AverageIOField<vector > >("dummy") )
            {
                Info << "Assigning read velocity values" << endl;
                AverageIOField<vector> & field(const_cast<AverageIOField<vector> &>(this->db().template
                            lookupObject<AverageIOField<vector> >("dummy")));
                forAll(vals, i)
                {
                    field[i][0] = velocity[i][0];
                    field[i][1] = velocity[i][1];
                    field[i][2] = velocity[i][2];
                }
                //Info<<"Check table @ lo:"<<vals[100]<<"-("<<velocity[100][0]<<" "<<velocity[100][1]<<" "<<velocity[100][2]<<")\n";
            }

            startAverage_ = vals.average();
            startSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }

    if (hi != endSampleTime_)
    {
        endSampleTime_ = hi;

        if (endSampleTime_ == -1)
        {
            // endTime no longer valid. Might as well clear endValues.
            if (true)
            {
                Pout<< "checkTable : Clearing endValues" << endl;
            }
            endSampledValues_.clear();
        }
        else
        {
            if (true)
            {
                Pout<< "checkTable : Reading endValues from "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[endSampleTime_].name()
                    << endl;
            }

            offset[0] = hi;

            hdf5Status = H5Sselect_hyperslab(velocityDataspace, H5S_SELECT_SET, offset, NULL, 
                                              count, NULL);
            
            if (hdf5Status < 0)
            {
                Info << "HDF5 error selecting velocity slice." << endl;
            }

            // Read in the slice
            hdf5Status = H5Dread(velocityDataset, H5T_NATIVE_DOUBLE, velocityMemspace, velocityDataspace,
                    H5P_DEFAULT, velocity);
            

            if (hdf5Status < 0)
            {
                Info << "HDF5 error reading velocity slice." << endl;
            }

            // Reread values and interpolate
            AverageIOField<Type> vals
            (
                IOobject
                (
                    "dummy", //fieldTableName_,
                    this->db().time().constant(),
                    "boundaryData"
                   /this->patch().name()
                   /sampleTimes_[endSampleTime_].name(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                    //false
                ),
                nVelocity[1]
            );

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorIn
                (
                    "precursorHDF5FvPatchField<Type>::"
                    "checkTable()"
                )   << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << vals.objectPath() << exit(FatalError);
            }

            if (this->db().template 
                    foundObject<AverageIOField<vector> >("dummy"))
            {
                Info << "Assigning read velocity values" << endl;
                AverageIOField<vector> & field(const_cast<AverageIOField<vector> &>(this->db().template 
                            lookupObject<AverageIOField<vector> >("dummy")));
                forAll(vals, i)
                {
                    field[i][0]= velocity[i][0];
                    field[i][1]= velocity[i][1];
                    field[i][2]= velocity[i][2];
                }
                //Info<<"Check table @ hi:"<<vals[100]<<"-("<<velocity[100][0]<<" "<<velocity[100][1]<<" "<<velocity[100][2]<<")\n";
            }

            endAverage_ = vals.average();
            endSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }

    // Close the dataset and dataspace and file
    H5Dclose(velocityDataset);
    H5Sclose(velocityDataspace);
    H5Sclose(velocityMemspace);
    H5Fclose(hdf5DataBase);
}


template<class Type>
void precursorHDF5FvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


    checkTable();

    // Interpolate between the sampled data

    Type wantedAverage;

    if (endSampleTime_ == -1)
    {
        // only start value
        if (true)
        {
            Pout<< "updateCoeffs : Sampled, non-interpolated values"
                << " from start time:"
                << sampleTimes_[startSampleTime_].name() << nl;
        }

        this->operator==(startSampledValues_);
        wantedAverage = startAverage_;
    }
    else
    {
        scalar period = sampleTimes_[sampleTimes_.size()-1].value()-sampleTimes_[0].value();
        scalar dbTime = this->db().time().value();

        // Take remainder
        scalar seekTime;        
        if (recycling_)
        {
            seekTime = dbTime - int(dbTime/period)*period;
        }
        else
        {
            seekTime = dbTime;
        }

        if (seekTime < this->db().time().deltaTValue())
        {
            seekTime = 0.0;
        }

        scalar start = sampleTimes_[startSampleTime_].value();
        scalar end = sampleTimes_[endSampleTime_].value();

        scalar s = (seekTime - start)/(end - start);

        if (true)
        {
            Pout<< "updateCoeffs : Sampled, interpolated values"
                << " between start time:"
                << sampleTimes_[startSampleTime_].name()
                << " and end time:" << sampleTimes_[endSampleTime_].name()
                << " for seek time:"<< seekTime
                << " with weight:" << s << endl;
        }

        this->operator==((1 - s)*startSampledValues_ + s*endSampledValues_);
        wantedAverage = (1 - s)*startAverage_ + s*endAverage_;
    }

    // Enforce average. Either by scaling (if scaling factor > 0.5) or by
    // offsetting.
    if (setAverage_)
    {
        const Field<Type>& fld = *this;

        Type averagePsi =
            gSum(this->patch().magSf()*fld)
           /gSum(this->patch().magSf());

        if (debug)
        {
            Pout<< "updateCoeffs :"
                << " actual average:" << averagePsi
                << " wanted average:" << wantedAverage
                << endl;
        }

        if (mag(averagePsi) < VSMALL)
        {
            // Field too small to scale. Offset instead.
            const Type offset = wantedAverage - averagePsi;
            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " offsetting with:" << offset << endl;
            }
            this->operator==(fld + offset);
        }
        else
        {
            const scalar scale = mag(wantedAverage)/mag(averagePsi);

            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " scaling with:" << scale << endl;
            }
            this->operator==(scale*fld);
        }
    }

    // apply offset to mapped values
    const scalar t = this->db().time().timeOutputValue();
    this->operator==(*this + offset_->value(t));

    if (debug)
    {
        Pout<< "updateCoeffs : set fixedValue to min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this) << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void precursorHDF5FvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("setAverage") << setAverage_ << token::END_STATEMENT << nl;
    if (perturb_ != 1e-5)
    {
        os.writeKeyword("perturb") << perturb_ << token::END_STATEMENT << nl;
    }

    if (fieldTableName_ != this->dimensionedInternalField().name())
    {
        os.writeKeyword("fieldTableName") << fieldTableName_
            << token::END_STATEMENT << nl;
    }

    if
    (
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        )
    )
    {
        os.writeKeyword("mapMethod") << mapMethod_
            << token::END_STATEMENT << nl;
    }

    os.writeKeyword("recycling") << recycling_ << token::END_STATEMENT << nl;

    offset_->writeData(os);

    os.writeKeyword("hdf5FileName") << hdf5FileName_ << token::END_STATEMENT <<nl;

    os.writeKeyword("hdf5PointsDatasetName") << hdf5PointsDatasetName_
        << token::END_STATEMENT <<nl;

    os.writeKeyword("hdf5SampleTimesDatasetName") << hdf5SampleTimesDatasetName_
        << token::END_STATEMENT <<nl;

    os.writeKeyword("hdf5FieldValuesDatasetName") << hdf5FieldValuesDatasetName_
        << token::END_STATEMENT <<nl;

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

//// ************************************************************************* //
