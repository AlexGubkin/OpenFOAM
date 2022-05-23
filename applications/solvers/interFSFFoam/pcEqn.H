label pcRefCell = 0;
scalar pcRefValue = 0.0;
const scalar cPhicThresholdCoeff(alphaControls.lookupOrDefault<scalar>("cPhicThresholdCoeff", 0.01));
const scalar cPhicRelaxCoeff(alphaControls.lookupOrDefault<scalar>("cPhicRelaxCoeff", 0.9));

// pc.storeOldTime();
// pc.oldTime();

gradPcf =
    deltaf/(mag(deltaf) + mixture.deltaN())*fvc::interpolate(fvc::grad(pc));

phic =
    fvc::interpolate(sigmaPtr->sigma())
   *fvc::interpolate(w2*kappa)
   /fvc::interpolate(w2)*deltaf*mesh.magSf();

// phic -=
//     (gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf();

phicFilt =
    cPhicRelaxCoeff*phicFilt
  + (1.0 - cPhicRelaxCoeff)
   *((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf());

phic -=
    phicFilt;

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
        ==
        fvc::div(phic)
    );

    pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
//         gradPcf =
//             mag(deltaf)/(mag(deltaf) + mixture.deltaN())*fvc::interpolate(fvc::grad(pc));
// 
//         phic -=
//             (gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf();
// 
//         phicFilt =
//             cPhicRelaxCoeff*phicFilt + (1.0 - cPhicRelaxCoeff)*(gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf();

        phic =
            mag(deltaf)/(mag(deltaf) + mixture.deltaN())*fvc::snGrad(pc)*mesh.magSf();

        const dimensionedScalar phicThreshold
        (
            dimForce/dimLength,
            cPhicThresholdCoeff
           *gSum(static_cast<surfaceScalarField>(mag(phic)*mesh.magSf()))
           /(gSum(static_cast<surfaceScalarField>(pos(mag(phic))*mesh.magSf())) + SMALL)
        );

        Info<< "phicThreshold = " << phicThreshold.value() << nl;

        phic =
            phic
          - max
            (
                min
                (
                    phic,
                    phicThreshold
                ),
              - phicThreshold
            );
            
        gradSigmaf =
            fvc::interpolate(fvc::grad(sigmaPtr->sigma()));

        phict =
            mag(gradAlphaShf)*(gradSigmaf - (gradSigmaf & nHatfv)*nHatfv) & mesh.Sf();

//         tangentGradSigmaf =
//             mag(gradAlphaShf)*(gradSigmaf - (gradSigmaf & nHatfv)*nHatfv) & mesh.Sf()/mesh.magSf();

//         surfaceTensionForce =
//             (
//                 phic
//               - max
//                 (
//                     min
//                     (
//                         phic,
//                         phicThreshold
//                     ),
//                   - phicThreshold
//                 )
//             )/mesh.magSf()
//           + tangentGradSigmaf;
    }
}
