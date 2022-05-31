label pcRefCell = 0;
scalar pcRefValue = 0.0;
const scalar cPhicFilterCoeff(alphaControls.lookupOrDefault<scalar>("cPhicFilterCoeff", 0.1));
const scalar cPhicRelaxCoeff(alphaControls.lookupOrDefault<scalar>("cPhicRelaxCoeff", 0.7));
const scalar cPhicThresholdCoeff(alphaControls.lookupOrDefault<scalar>("cPhicThresholdCoeff", 0.01));

gradPcf =
    fvc::interpolate(w1)*fvc::interpolate(fvc::grad(pc));

phic =
    fvc::interpolate(sigmaPtr->sigma())
   *fvc::interpolate(w2*kappa)
   /fvc::interpolate(w2)*deltaf*mesh.magSf();

// const dimensionedScalar fcAvg
// (
//     average(mag(phic/mesh.magSf()))
//    /(average(pos(mag(phic))) + SMALL)
// );
// 
// Info<< "fcAvg = " << fcAvg.value() << nl;
   
phicFilt =
    phicFilt.oldTime()
  + cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf());

phic -=
    phicFilt;

phic =
    cPhicRelaxCoeff*phicPrevIter
  + (1.0 - cPhicRelaxCoeff)*phic;

phicPrevIter = phic;

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
     == fvc::div(phic)
    );

    pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        gradPcf =
            fvc::interpolate(w1)*fvc::interpolate(fvc::grad(pc));

        phic =
            fvc::interpolate(sigmaPtr->sigma())
           *fvc::interpolate(w2*kappa)
           /fvc::interpolate(w2)*deltaf*mesh.magSf();

        phicFilt =
            phicFilt.oldTime()
          + cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf());

        phic -= phicFilt;

//         phic -= fvc::snGrad(pc)*mesh.magSf();

//         phic =
//             phic
//           - max
//             (
//                 min
//                 (
//                     phic,
//                     cPhicThresholdCoeff*fcAvg*mesh.magSf()
//                 ),
//               - cPhicThresholdCoeff*fcAvg*mesh.magSf()
//             );

//         gradSigmaf =
//             fvc::interpolate(fvc::grad(sigmaPtr->sigma()));
// 
//         phict =
//             mag(gradAlphaShf)*(gradSigmaf - (gradSigmaf & nHatfv)*nHatfv) & mesh.Sf();
    }
}
