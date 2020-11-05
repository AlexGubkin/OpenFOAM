/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::adaptiveTotalPressureFvPatchScalarField

Description
    This boundary condition provides a adaptiveTotalPressure condition,
    calculated as:

        \f[
            Q = Q_{0} + Q_{p} + s*Q_{t}
        \f]

    where
    \vartable
        s      | single scalar value [units]
        Q_{0}  | single scalar value [units]
        Q_{p}  | scalar field across patch [units]
        Q_{t}  | scalar function of time [units]
    \endtable

Usage
    \table
        Property   | Description                  | Req'd? | Default
        scalarData | single scalar value          | yes    |
        data       | single scalar value          | yes    |
        fieldData  | scalar field across patch    | yes    |
        timeVsData | scalar function of time      | yes    |
        wordData   | word, eg name of data object | no     | wordDefault
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       adaptiveTotalPressure;
        scalarData -1;
        data       1;
        fieldData  uniform 3;
        timeVsData table (
                             (0 0)
                             (1 2)
                          );
        wordName   anotherName;
        value      uniform 4; // optional initial value
    }
    \endverbatim

SourceFiles
    adaptiveTotalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveTotalPressure_H
#define adaptiveTotalPressure_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class adaptiveTotalPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class adaptiveTotalPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

//         //- Inlet integral flow rate
//         autoPtr<Function1<scalar>> flowRateMin_;
// 
//         //- Inlet integral flow rate
//         autoPtr<Function1<scalar>> flowRateMax_;

        //- Inlet integral flow rate
        autoPtr<Function1<scalar>> criticalFlowRate_;

        //- Is volumetric?
        bool volumetric_;

        //- Switcher for increase/decrease pressure
        bool switcher_;

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if necessary
        word rhoName_;

        //- Name of the compressibility field used to calculate the wave speed
        word psiName_;

        //- Characteristic time of pressure changing
        scalar tc_;
        
        //- Switching time
        scalar ts_;

        //- Heat capacity ratio
        scalar gamma_;

        //- Rho initialisation value (for start; if value not supplied)
        scalar rhoInlet_;

        //- Patch averaged pressure
        scalar pAveragedp_;

//         //- High level of total pressure
//         scalarField p0High_;
// 
//         //- Low level of total pressure
//         scalarField p0Low_;

        //- Table of time vs total pressure, including the bounding treatment
        autoPtr<Function1<scalar>> p0High_;

        //- Table of time vs total pressure, including the bounding treatment
        autoPtr<Function1<scalar>> p0Low_;

    // Private Member Functions

        //- Return current time
        scalar t() const;


public:

    //- Runtime type information
    TypeName("adaptiveTotalPressure");


    // Constructors

        //- Construct from patch and internal field
        adaptiveTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        adaptiveTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        adaptiveTotalPressureFvPatchScalarField
        (
            const adaptiveTotalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        adaptiveTotalPressureFvPatchScalarField
        (
            const adaptiveTotalPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new adaptiveTotalPressureFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        adaptiveTotalPressureFvPatchScalarField
        (
            const adaptiveTotalPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new adaptiveTotalPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& UName()
            {
                return UName_;
            }

            //- Return the name of the flux field
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of the flux field
            //  to allow adjustment
            word& phiName()
            {
                return phiName_;
            }

            //- Return the name of the density field
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of the density field
            //  to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return the name of the compressibility field
            const word& psiName() const
            {
                return psiName_;
            }

            //- Return reference to the name of the compressibility field
            //  to allow adjustment
            word& psiName()
            {
                return psiName_;
            }

//             //- Return the heat capacity ratio
//             scalar gamma() const
//             {
//                 return gamma_;
//             }
// 
//             //- Return reference to the heat capacity ratio to allow adjustment
//             scalar& gamma()
//             {
//                 return gamma_;
//             }
// 
//             //- Return the rho initialisation value
//             scalar rhoInlet() const
//             {
//                 return rhoInlet_;
//             }
// 
//             //- Return reference to the rho initialisation value to allow adjustment
//             scalar& rhoInlet()
//             {
//                 return rhoInlet_;
//             }

//             //- Return the high level of total pressure
//             const scalarField& p0High() const
//             {
//                 return p0High_;
//             }
// 
//             //- Return reference to the high level of total pressure to allow adjustment
//             scalarField& p0High()
//             {
//                 return p0High_;
//             }
// 
//             //- Return the low level of total pressure
//             const scalarField& p0Low() const
//             {
//                 return p0Low_;
//             }
// 
//             //- Return reference to the total pressure to allow adjustment
//             scalarField& p0Low()
//             {
//                 return p0Low_;
//             }


//         // Mapping functions
// 
//             //- Map (and resize as needed) from self given a mapping object
//             //  Used to update fields following mesh topology change
//             virtual void autoMap(const fvPatchFieldMapper&);
// 
//             //- Reverse map the given fvPatchField onto this fvPatchField
//             //  Used to reconstruct fields
//             virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Inherit updateCoeffs from fixedValueFvPatchScalarField
            using fixedValueFvPatchScalarField::updateCoeffs;

            //- Update the coefficients associated with the patch field
            //  using the given patch total pressure and velocity fields
            virtual void updateCoeffs
            (
//                 const scalarField& p0Highp,
//                 const scalarField& p0Lowp,
                const vectorField& Up
            );

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
