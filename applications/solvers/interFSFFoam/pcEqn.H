// fcf.storePrevIter();
fcfFilt.storePrevIter();

// const scalar cPhicFilterCoeff(alphaControls.lookupOrDefault<scalar>("cPhicFilterCoeff", 0.1));
const scalar cPhicThresholdCoeff(alphaControls.lookupOrDefault<scalar>("cPhicThresholdCoeff", 0.01));

// label pcRefCell = 0;
// scalar pcRefValue = 0.0;

// gradPcf =
//     fvc::interpolate(w1)*fvc::interpolate(fvc::grad(pc));

// gradPcf =
//     fvc::interpolate(fvc::grad(pc));

gradPc =
    fvc::grad(pc);

fcf =
    fvc::interpolate(sigmaPtr->sigma())
   *fvc::interpolate(w2*kappa)
   /fvc::interpolate(w2)*deltaf;

fc = fvc::reconstruct(fcf*mesh.magSf());

// fcfFilt =
//     deltaf/(mag(deltaf) + mixture.deltaN())
//    *(
//         fcfFilt.oldTime()
//       + cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf()/mesh.magSf())
//     );

fcfFilt =
    mag(deltaf)/(mag(deltaf) + mixture.deltaN())
   *fvc::interpolate(gradPc - (gradPc & nHatv)*nHatv) & (mesh.Sf()/mesh.magSf());

fcfFilt.relax();

fcFilt = fvc::reconstruct(fcfFilt*mesh.magSf());

fcf -= fcfFilt;

// fcf.relax();

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
     == fvc::div(fcf*mesh.magSf())
    );

    pcEqn.setReference
    (
        pressureReference.refCell(),
        getRefCellValue(pc, pressureReference.refCell())
    );

//     pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();
    pc.correctBoundaryConditions();

    if (pimple.finalNonOrthogonalIter())
    {
        const dimensionedScalar dimFactor
        (
            dimensionSet(dimForce/dimLength),
            1.0
        );
//         phic = pcEqn.flux();
        phic = (fcf - fvc::snGrad(pc))*mesh.magSf();
//         phicFilt =
//             pos(phic)*mag(phic)
//           - (1.0 - pos(phic))*mag(phic);
//         phic = mag(deltaf)/(mag(deltaf) + mixture.deltaN())*fvc::snGrad(pc)*mesh.magSf();
//         phic = fvc::snGrad(pc)*mesh.magSf();

        const dimensionedScalar phicAvg
        (
            dimensionSet(dimForce/dimLength),
//             average(mag(fcf))/(average(pos(mag(fcf))) + SMALL)
//             average(mag(fcf))
            cPhicThresholdCoeff*gSum(mag(fcf.field()*mesh.magSf()))
           /gSum(SMALL + pos(mag(fcf.field()) - SMALL))
        );

        Info<< "phicAvg = " << phicAvg.value() << nl;

        phic -=
            max
            (
                min
                (
                    phic,
                    phicAvg
                ),
              - phicAvg
            );

//         phic -= cPhicThresholdCoeff*phicFilt;
    }
}
