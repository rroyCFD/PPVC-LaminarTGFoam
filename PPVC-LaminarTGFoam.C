/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    PPVC-LaminarFoam

Description
    Transient solver for incompressible flow using a pressure projection,
    rotational velocity-correction strategy.
    Following the paper "A robust and accurate outflow boundary condition for
    incompressible flow simulations on severely-truncated unbounded domains" by
    Dong et al.

    Work for constant viscosity, laminar flow
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
//#include "turbulenceModel.H"
#include "IFstream.H"
#include "OFstream.H"
//#include "LESfilter.H"
#include "IOmanip.H" // for input/ouput format control
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H" 
    #include "readTransportProperties.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"    
    
    #include "TaylorGreenFiles/readAndDeclareVariables.H"    
    #include "TaylorGreenFiles/createErrorFields.H"

    #include "TaylorGreenFiles/initialize.H"    
    #include "TaylorGreenFiles/errorNorm.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    scalar gamma0; // temoral discretization coeffs
    dimensionedScalar dt = runTime.deltaT();
    
    //U0.correctBoundaryConditions();     
       
    while (runTime.loop())
    {   
        Info<< "Time = " << runTime.timeName()<< endl;
        #include "CourantNo.H"
        #include "setDeltaT.H" 
            
        // Start with one first order backward step
        if(runTime.value() == dt.value())
        {
            Info<< "\n first time step: 1st order backward" 
                    "temporal discretization applied\n" << endl;
            gamma0 = 1.0;
            Uhat   = U;
            Ustar  = U;           
        } else {
            //Info<< "\nSecond-order backward temporal discretization\n" << endl;
            gamma0 = 1.5;
            Uhat   = 2.*U - 0.5*U0; 
            Ustar  = 2.*U - U0;      
        }
        
        // Update in every time step
        dt = runTime.deltaT();  // update deltaT  
        U0 = U;
        
        Info<< "Updating W^(n+1) velocity.... " << endl;
        W.correctBoundaryConditions(); // update W =U^(n+1) at dirichlet boundaries
        
        // Pressure Poisson for new pressure----------------------------------------//
        phi = (fvc::interpolate(U) & mesh.Sf());       

        // Non-orthogonal pressure corrector loop
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (      
            fvm::laplacian(p) 
            == 
            fvc::div(
                Uhat/dt 
                - (fvc::div(phi, U) - fvc::div(phi)*U) 
                )
            );
            pEqn.setReference(pRefCell, pRefValue);

            if (nonOrth == nNonOrthCorr)
            {
                pEqn.solve(mesh.solver("pFinal"));
            } 
            else
            {
                pEqn.solve(mesh.solver("p"));
            }
        }
              
        // update p at all boundaries
        Info << "Final Pressure Update....." << endl;
        p.correctBoundaryConditions();
        
        //Info<< "Calcuting phiStar for pEqu:..." << endl;
        surfaceScalarField phiStar = (fvc::interpolate(Ustar) & mesh.Sf()); 

        // Helmholz equation for velocity ------------------------------------------// 
        fvVectorMatrix UEqn
            (
              fvm::Sp(gamma0/dt, U)
            - fvm::laplacian(nu, U)  
            == 
            (
              Uhat/dt 
            - (fvc::div(phiStar, Ustar) - fvc::div(phiStar)*Ustar)
            - fvc::grad(p)
            )

        );
        UEqn.solve();
        
        #include "continuityErrs.H"
        
        #include "TaylorGreenFiles/errorNorm.H"
        #include "TaylorGreenFiles/globalProperties.H"

        runTime.write();      
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //