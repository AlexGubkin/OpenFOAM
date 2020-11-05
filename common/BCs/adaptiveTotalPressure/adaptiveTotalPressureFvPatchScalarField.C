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

\*---------------------------------------------------------------------------*/

#include "adaptiveTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::adaptiveTotalPressureFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveTotalPressureFvPatchScalarField::
adaptiveTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
//     flowRateMin_(),
//     flowRateMax_(),
    criticalFlowRate_(),
    volumetric_(false),
    switcher_(true),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    psiName_("none"),
    tc_(0.0),
    ts_(0.0),
    gamma_(0.0),
    rhoInlet_(0.0),
    pAveragedp_(0.0),
    p0High_(),
    p0Low_()
{}


Foam::adaptiveTotalPressureFvPatchScalarField::
adaptiveTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    psiName_(dict.lookupOrDefault<word>("psi", "none")),
    tc_(dict.lookupOrDefault<scalar>("tc", 1e-6)),
    gamma_(psiName_ != "none" ? readScalar(dict.lookup("gamma")) : 1),
    rhoInlet_(dict.lookupOrDefault<scalar>("rhoInlet", -vGreat)),
    p0High_(Function1<scalar>::New("p0High", dict)),
    p0Low_(Function1<scalar>::New("p0Low", dict))
{
    switcher_ = true;

    ts_ = t();

    if (dict.found("criticalVolumetricFlowRate")/*dict.found("volumetricFlowRateMin") && dict.found("volumetricFlowRateMax")*/)
    {
        volumetric_ = true;
//         flowRateMin_ = Function1<scalar>::New("volumetricFlowRateMin", dict);
//         flowRateMax_ = Function1<scalar>::New("volumetricFlowRateMax", dict);
        criticalFlowRate_ = Function1<scalar>::New("criticalVolumetricFlowRate", dict);
//         rhoName_ = "rho";
    }
    else if (dict.found("criticalMassFlowRate")/*dict.found("massFlowRateMin") && dict.found("massFlowRateMax")*/)
    {
        volumetric_ = false;
//         flowRateMin_ = Function1<scalar>::New("massFlowRateMin", dict);
//         flowRateMax_ = Function1<scalar>::New("massFlowRateMax", dict);
        criticalFlowRate_ = Function1<scalar>::New("criticalMassFlowRate", dict);
//         rhoName_ = word(dict.lookupOrDefault<word>("rho", "rho"));
    }
    else
    {
//         FatalIOErrorInFunction
//         (
//             dict
//         )   << "Please supply either 'volumetricFlowRateMin' and 'volumetricFlowRateMax' or"
//             << " 'massFlowRateMin', 'massFlowRateMax' and 'rho'" << exit(FatalIOError);
        FatalIOErrorInFunction
        (
            dict
        )   << "Please supply either 'criticalVolumetricFlowRate' or"
            << " 'criticalMassFlowRate' and 'rho'" << exit(FatalIOError);
    }

    pAveragedp_ = p0Low_->value(t());

    fvPatchField<scalar>::operator==(p0Low_->value(t()));

//     evaluate(Pstream::commsTypes::blocking);

//     fixedValueFvPatchScalarField::evaluate();

    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
    */
}


Foam::adaptiveTotalPressureFvPatchScalarField::
adaptiveTotalPressureFvPatchScalarField
(
    const adaptiveTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper, false), // Don't map
    criticalFlowRate_(ptf.criticalFlowRate_, false),
    volumetric_(ptf.volumetric_),
    switcher_(ptf.switcher_),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    tc_(ptf.tc_),
    ts_(ptf.ts_),
    gamma_(ptf.gamma_),
    rhoInlet_(ptf.rhoInlet_),
    pAveragedp_(ptf.pAveragedp_),
    p0High_(ptf.p0High_, false),
    p0Low_(ptf.p0Low_, false)
{
    patchType() = ptf.patchType();

    // Set the patch pressure to the current total pressure
    // This is not ideal but avoids problems with the creation of patch faces
    fvPatchScalarField::operator==(p0Low_->value(t()));
}


Foam::adaptiveTotalPressureFvPatchScalarField::
adaptiveTotalPressureFvPatchScalarField
(
    const adaptiveTotalPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    criticalFlowRate_(ptf.criticalFlowRate_, false),
    volumetric_(ptf.volumetric_),
    switcher_(ptf.switcher_),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    tc_(ptf.tc_),
    ts_(ptf.ts_),
    gamma_(ptf.gamma_),
    rhoInlet_(ptf.rhoInlet_),
    pAveragedp_(ptf.pAveragedp_),
    p0High_(ptf.p0High_, false),
    p0Low_(ptf.p0Low_, false)
{}


Foam::adaptiveTotalPressureFvPatchScalarField::
adaptiveTotalPressureFvPatchScalarField
(
    const adaptiveTotalPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    criticalFlowRate_(ptf.criticalFlowRate_, false),
    volumetric_(ptf.volumetric_),
    switcher_(ptf.switcher_),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    tc_(ptf.tc_),
    ts_(ptf.ts_),
    gamma_(ptf.gamma_),
    rhoInlet_(ptf.rhoInlet_),
    pAveragedp_(ptf.pAveragedp_),
    p0High_(ptf.p0High_, false),
    p0Low_(ptf.p0Low_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void Foam::adaptiveTotalPressureFvPatchScalarField::autoMap
// (
//     const fvPatchFieldMapper& m
// )
// {
//     fixedValueFvPatchScalarField::autoMap(m);
//     m(p0High_, p0High_);
//     m(p0Low_, p0Low_);
// }
// 
// 
// void Foam::adaptiveTotalPressureFvPatchScalarField::rmap
// (
//     const fvPatchScalarField& ptf,
//     const labelList& addr
// )
// {
//     fixedValueFvPatchScalarField::rmap(ptf, addr);
// 
//     const adaptiveTotalPressureFvPatchScalarField& tiptf =
//         refCast<const adaptiveTotalPressureFvPatchScalarField>(ptf);
// 
//     p0High_.rmap(tiptf.p0High_, addr);
//     p0Low_.rmap(tiptf.p0Low_, addr);
// }


void Foam::adaptiveTotalPressureFvPatchScalarField::updateCoeffs
(
//     const scalarField& p0Highp,
//     const scalarField& p0Lowp,
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }

//     const scalar t = db().time().timeOutputValue();

    scalar p0Highp = p0High_->value(t());
    scalar p0Lowp = p0Low_->value(t());

    const scalar criticalFlowRate = criticalFlowRate_->value(t());

//     const surfaceScalarField&
//         phi = db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

//     const volScalarField& p =
//         db().lookupObject<volScalarField>(this->internalField().name());

    fvPatchField<scalar> p =
            patch().lookupPatchField<volScalarField, scalar>(this->internalField().name());

//     label patchI = this->patch().index();
// 
//     const surfaceScalarField& pp = p.boundaryField()[patchi];

    if (internalField().dimensions() == dimPressure)
    {
        const fvPatchField<scalar>& rho =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

//         scalarField p(this->patchInternalField());

        if (psiName_ == "none")
        {
            // Variable density and low-speed compressible flow

            p += 0.5*rho*(1.0 - pos0(phip))*magSqr(Up);

            const scalar calculatedFlowRate =
                (volumetric_) ? -gSum(phip/rho) : -gSum(phip);

            if (calculatedFlowRate <= criticalFlowRate)
            {
                if (switcher_)
                {
                    switcher_ = false;
    
                    ts_ = t();

                    pAveragedp_ =
                        gSum(this->patch().magSf() * p)/gSum(this->patch().magSf());
                }

                operator==
                (
                    min
                    (
                        pAveragedp_ + (p0Highp - p0Lowp)/tc_*(t() - ts_),
                        p0Highp
                    )
                    -0.5*rho*(1.0 - pos0(phip))*magSqr(Up)
                );
            }
            else
            {
                if (!switcher_)
                {
                    switcher_ = true;

                    ts_ = t();

                    pAveragedp_ =
                        gSum(this->patch().magSf() * p)/gSum(this->patch().magSf());
                }

                operator==
                (
                    max
                    (
                        pAveragedp_ - (p0Highp - p0Lowp)/tc_*(t() - ts_),
                        p0Lowp
                    )
                    -0.5*rho*(1.0 - pos0(phip))*magSqr(Up)
                );
            }
        }
        else
        {
            // High-speed compressible flow

            const fvPatchField<scalar>& psip =
                patch().lookupPatchField<volScalarField, scalar>(psiName_);

            if (gamma_ > 1)
            {
                scalar gM1ByG = (gamma_ - 1)/gamma_;

                p *=
                    pow
                    (
                        (1.0 + 0.5*psip*gM1ByG*(1.0 - pos0(phip))*magSqr(Up)),
                        1.0/gM1ByG
                    );

//                 const scalar pAveragedp_ =
//                     gSum(this->patch().magSf() * p)
//                     / gSum(this->patch().magSf());

                const scalar calculatedFlowRate =
                    (volumetric_) ? -gSum(phip/rho) : -gSum(phip);

//                 Info<< "calculatedFlowRate = " << calculatedFlowRate
//                     << ((volumetric_) ? " [m^3 / s]" : " [kg / s]") << nl
//                     << "criticalFlowRate = " << criticalFlowRate 
//                     << ((volumetric_) ? " [m^3 / s]" : " [kg / s]") << nl;

                operator==
                (
                    (calculatedFlowRate <= criticalFlowRate)
                    ?
                        min(pAveragedp_ + 1e-5*pAveragedp_,p0Highp)
                        /pow
                        (
                            (1.0 + 0.5*psip*gM1ByG*(1.0 - pos0(phip))*magSqr(Up)),
                            1.0/gM1ByG
                        )
                    :
                        max(pAveragedp_,p0Lowp)
                        /pow
                        (
                            (1.0 + 0.5*psip*gM1ByG*(1.0 - pos0(phip))*magSqr(Up)),
                            1.0/gM1ByG
                        )
                );
            }
            else
            {
                p *= (1.0 + 0.5*psip*(1.0 - pos0(phip))*magSqr(Up));

//                 const scalar pAveragedp_ =
//                     gSum(this->patch().magSf() * p)
//                     / gSum(this->patch().magSf());

                const scalar calculatedFlowRate =
                    (volumetric_) ? -gSum(phip/rho) : -gSum(phip);

//                 Info<< "calculatedFlowRate = " << calculatedFlowRate
//                     << ((volumetric_) ? " [m^3 / s]" : " [kg / s]") << nl
//                     << "criticalFlowRate = " << criticalFlowRate 
//                     << ((volumetric_) ? " [m^3 / s]" : " [kg / s]") << nl;

                operator==
                (
                    (calculatedFlowRate <= criticalFlowRate)
                    ?
                        min(pAveragedp_ + 1e-5*pAveragedp_,p0Highp)
                        /(1.0 + 0.5*psip*(1.0 - pos0(phip))*magSqr(Up))
                    :
                        max(pAveragedp_,p0Lowp)
                        /(1.0 + 0.5*psip*(1.0 - pos0(phip))*magSqr(Up))
                );
            }
        }

    }
    else if (internalField().dimensions() == dimPressure/dimDensity)
    {
        // Incompressible flow

        // Use constant density
        if (!volumetric_ && rhoInlet_ < 0)
        {
            FatalErrorInFunction
                << "Did not find constant density 'rhoInlet' specified"
                << exit(FatalError);
        }

        p += 0.5*(1.0 - pos0(phip))*magSqr(Up);

//         const scalar pAveragedp_ =
//             gSum(this->patch().magSf() * p)
//             / gSum(this->patch().magSf());

        const scalar calculatedFlowRate =
            (volumetric_) ? -gSum(phip) : -rhoInlet_*gSum(phip);

//             Info<< "calculatedFlowRate = " << calculatedFlowRate
//                 << ((volumetric_) ? " [m^3 / s]" : " [kg / s]") << nl
//                 << "criticalFlowRate = " << criticalFlowRate 
//                 << ((volumetric_) ? " [m^3 / s]" : " [kg / s]") << nl;

        operator==
        (
            (calculatedFlowRate <= criticalFlowRate)
            ?
                min(pAveragedp_ + 1e-5*pAveragedp_,p0Highp)
                - 0.5*(1.0 - pos0(phip))*magSqr(Up)
            :
                max(pAveragedp_,p0Lowp)
                - 0.5*(1.0 - pos0(phip))*magSqr(Up)
        );

    }
    else
    {
        FatalErrorInFunction
            << " Incorrect pressure dimensions " << internalField().dimensions()
            << nl
            << "    Should be " << dimPressure
            << " for compressible/variable density flow" << nl
            << "    or " << dimPressure/dimDensity
            << " for incompressible flow," << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::adaptiveTotalPressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs
    (
//         p0High(),
//         p0Low(),
        patch().lookupPatchField<volVectorField, vector>(UName())
    );
}

void Foam::adaptiveTotalPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, criticalFlowRate_());
    writeEntry(os, "tc", tc_);
    if (!volumetric_)
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntryIfDifferent<scalar>(os, "rhoInlet", -vGreat, rhoInlet_);
    }
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry(os, "rho", rhoName_);
    writeEntry(os, "psi", psiName_);
    writeEntry(os, "gamma", gamma_);
    writeEntry(os, p0High_());
    writeEntry(os, p0Low_());
//     writeEntry(os, "p0High", p0High_);
//     writeEntry(os, "p0Low", p0Low_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adaptiveTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
