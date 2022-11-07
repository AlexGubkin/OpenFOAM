fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U) - fvm::Sp(contErr, U)
  + MRF.DDt(rho, U)
  + turbulence.divDevTau(U)
 ==
    fvModels.source(rho, U)
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

    K = 0.5*magSqr(U);
}
