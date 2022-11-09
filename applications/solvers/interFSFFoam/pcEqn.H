fcf.storePrevIter();
fcfFilt.storePrevIter();

const scalar cPhicFilterCoeff(alphaControls.lookupOrDefault<scalar>("cPhicFilterCoeff", 0.1));
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
fcFilt = gradPc - (gradPc & nHatv)*nHatv;

fcf.relax();

// fcFilt.relax();

// fcfFilt =
//     mag(deltaf)/(mag(deltaf) + mixture.deltaN())
//    *(
//         fcfFilt.oldTime()
//       + cPhicFilterCoeff*(fvc::interpolate(fcFilt) & (mesh.Sf()/mesh.magSf()))
//     );

fcfFilt =
   mag(deltaf)/(mag(deltaf) + mixture.deltaN())
  *cPhicFilterCoeff*fvc::interpolate(fcFilt) & (mesh.Sf()/mesh.magSf());

fcfFilt.relax();

// fcFilt = fvc::reconstruct(fcfFilt*mesh.magSf());

// fcf -= fcfFilt;

// fcf.relax();

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
     == fvc::div((fcf - fcfFilt)*mesh.magSf())
    );

    pcEqn.setReference
    (
        pcPressureReference.refCell(),
        getRefCellValue(pc, pcPressureReference.refCell())
    );

//     pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();
//     pc.correctBoundaryConditions();

    if (pimple.finalNonOrthogonalIter())
    {
        pc += dimensionedScalar
        (
            "pc",
            pc.dimensions(),
            pcPressureReference.refValue()
          - getRefCellValue(pc, pcPressureReference.refCell())
        );

//         phic = pcEqn.flux();
        phic = (fcf - fcfFilt - fvc::snGrad(pc))*mesh.magSf();

        const dimensionedScalar fcfAvg
        (
            dimensionSet(dimForce/dimVolume),
            average(mag(fcf.field()))
        );

        Info<< "fcfAvg = " << fcfAvg.value() << nl;

        phic =
            phic
          - max
            (
                min
                (
                    phic,
                    cPhicThresholdCoeff*fcfAvg*mesh.magSf()
                ),
              - cPhicThresholdCoeff*fcfAvg*mesh.magSf()
            );
    }
}
