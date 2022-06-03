label pcRefCell = 0;
scalar pcRefValue = 0.0;
const scalar cPhicFilterCoeff(alphaControls.lookupOrDefault<scalar>("cPhicFilterCoeff", 0.1));
// const scalar cPhicRelaxCoeff(alphaControls.lookupOrDefault<scalar>("cPhicRelaxCoeff", 0.7));
const scalar cPhicThresholdCoeff(alphaControls.lookupOrDefault<scalar>("cPhicThresholdCoeff", 0.01));

// gradPcf =
//     fvc::interpolate(w1)*fvc::interpolate(fvc::grad(pc));

gradPcf =
    mag(deltaf)/(mag(deltaf) + mixture.deltaN())*fvc::interpolate(fvc::grad(pc));

phic =
    fvc::interpolate(sigmaPtr->sigma())
   *fvc::interpolate(w2*kappa)
   /fvc::interpolate(w2)*deltaf*mesh.magSf();
   
// phicFilt =
//     phicFilt.oldTime()
//   + cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf());

phicFilt =
    cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf());

phic -= phicFilt;

// phic =
//     cPhicRelaxCoeff*phicPrevIter
//   + (1.0 - cPhicRelaxCoeff)*phic;
// 
// phicPrevIter = phic;

phic.relax();

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc)
     == fvc::div(phic)
    );

    pcEqn.setReference(pcRefCell, pcRefValue);
    pcEqn.solve();

//     gradPcf =
//         deltaf/(mag(deltaf) + mixture.deltaN())*fvc::interpolate(fvc::grad(pc));

//     gradPcf =
//         fvc::interpolate(w1)*fvc::interpolate(fvc::grad(pc));
// 
//     phicFilt =
//         cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf());

//     phic -= phicFilt;

    if (pimple.finalNonOrthogonalIter())
    {
//         gradPcf =
//             fvc::interpolate(w1)*fvc::interpolate(fvc::grad(pc));

        gradPcf =
            mag(deltaf)/(mag(deltaf) + mixture.deltaN())*fvc::interpolate(fvc::grad(pc));

        phic =
            fvc::interpolate(sigmaPtr->sigma())
           *fvc::interpolate(w2*kappa)
           /fvc::interpolate(w2)*deltaf*mesh.magSf();

//         phicFilt =
//             phicFilt.oldTime()
//           + cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf());

        phicFilt =
            cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf());

        phic -= phicFilt;

//         phic -= fvc::snGrad(pc)*mesh.magSf();

//         phic -= phicFilt;

        const dimensionedScalar fcAvg
        (
            average(mag(phic/mesh.magSf()))
           /(average(pos(mag(phic))) + SMALL)
        );

        Info<< "fcAvg = " << fcAvg.value() << nl;

        phic =
            phic
          - max
            (
                min
                (
                    phic,
                    cPhicThresholdCoeff*fcAvg*mesh.magSf()
                ),
              - cPhicThresholdCoeff*fcAvg*mesh.magSf()
            );

        pc.relax();

//         gradSigmaf =
//             fvc::interpolate(fvc::grad(sigmaPtr->sigma()));
// 
//         phict =
//             mag(gradAlphaShf)*(gradSigmaf - (gradSigmaf & nHatfv)*nHatfv) & mesh.Sf();
    }
}
