const scalar cPhicFilterCoeff(alphaControls.lookupOrDefault<scalar>("cPhicFilterCoeff", 0.1));

label pcRefCell = 0;
scalar pcRefValue = 0.0;

gradPcf =
    fvc::interpolate(fvc::grad(pc));

fcf =
    fvc::interpolate(sigmaPtr->sigma())
   *fvc::interpolate(w2*kappa)
   /fvc::interpolate(w2)*deltaf;

fcf -=
    deltaf/(mag(deltaf) + mixture.deltaN())
   *(
        fcf.oldTime()
      + cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf()/mesh.magSf())
    )

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
     == fvc::div(fcf*mesh.magSf())
    );

    pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();

//     if (pimple.finalNonOrthogonalIter())
//     {
//         phic -= pcEqn.flux();
//         phic -=
//             max
//             (
//                 min
//                 (
//                     phic,
//                     cPhicThresholdCoeff*fcAvg*mesh.magSf()
//                 ),
//               - cPhicThresholdCoeff*fcAvg*mesh.magSf()
//             );
//     }
}
