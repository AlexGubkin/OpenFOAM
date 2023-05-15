    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + advectionFactor*fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevTau(rho, U)
      + fvm::Sp(resistanceFactor*mixture.mu(), U)
     ==
       phaseChange.SU(rho, rhoPhi, U)
     + fvModels.source(rho, U)
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                  + fvc::snGrad(pc)
                ) * mesh.magSf()
              + phic
            )
        );

        fvConstraints.constrain(U);
    }
