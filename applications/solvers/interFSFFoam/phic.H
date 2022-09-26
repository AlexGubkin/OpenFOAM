const scalar cPhicFilterCoeff(alphaControls.lookupOrDefault<scalar>("cPhicFilterCoeff", 0.1));
const scalar cPhicThresholdCoeff(alphaControls.lookupOrDefault<scalar>("cPhicThresholdCoeff", 0.01));

phic =
    fvc::interpolate(sigmaPtr->sigma())
   *fvc::interpolate(w2*kappa)
   /fvc::interpolate(w2)*deltaf*mesh.magSf();

const dimensionedScalar fcAvg
(
    average(mag(phic/mesh.magSf()))
   /(average(pos(mag(phic))) + SMALL)
);

Info<< "fcAvg = " << fcAvg.value() << nl;

// gradPcf =
//     fvc::interpolate(w1)*fvc::interpolate(fvc::grad(pc));

gradPcf =
    fvc::interpolate(fvc::grad(pc));

phicFilt =
    deltaf/(mag(deltaf) + mixture.deltaN())
   *(
        phicFilt.oldTime()
      + cPhicFilterCoeff*((gradPcf - (gradPcf & nHatfv)*nHatfv) & mesh.Sf())
    );

phic -= phicFilt;

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

phic.relax();
