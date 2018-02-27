// Reynolds number
dimensionedScalar Re = Uinit*L/nu;

// Intialize Velocity field
volScalarField x = mesh.C().component(vector::X);
volScalarField y = mesh.C().component(vector::Y);
volScalarField z = mesh.C().component(vector::Z); 

// Reference analytical pressure at pRefCell
const scalar xPRefCell = U.mesh().C()[pRefCell].x();
const scalar yPRefCell = U.mesh().C()[pRefCell].y();
const scalar zPRefCell = U.mesh().C()[pRefCell].z();

dimensionedScalar pRefValueAna("pRefValueAna", dimensionSet(0,2,-2,0,0), 0);


if(TaylorGreenDimension =="2D")
{
    Ustart =  Uinit* ( vector(1,0,0) * sin(x/L) * cos(y/L) 
                                    - vector(0,1,0) * cos(x/L) * sin(y/L) 
                                    + vector(0,0,1) * scalar(0.)) ;

    pstart = sqr(Uinit)/4 * (cos(2.*x/L) + cos(2.*y/L));

    if(pRefOn) {
        pRefValueAna  = sqr(Uinit)* 0.25*(
                                      Foam::cos(2.*xPRefCell/L.value())
                                    + Foam::cos(2.*yPRefCell/L.value()) );
    }
} else if (TaylorGreenDimension == "3D")
{
    Ustart = Uinit* ( vector(1,0,0) * sin(x/L) * cos(y/L) * cos(z/L)
                                   - vector(0,1,0) * cos(x/L) * sin(y/L) * cos(z/L)
                                   + vector(0,0,1) * scalar(0.)) ;

    pstart = sqr(Uinit)/16 * (cos(2.*x/L) + cos(2.*y/L)) 
                                          * (cos(2.*y/L) + 2);

    if(pRefOn) {
        pRefValueAna  = sqr(Uinit)/16 * ( Foam::cos(2.*xPRefCell/L.value())
                                        + Foam::cos(2.*yPRefCell/L.value()) )
                                      * ( Foam::cos(2.*zPRefCell/L.value()) + 2.0); 
    }                                          
}

// Initialize velocity and pressur field at time 0
if (runTime.timeName() == "0")
{
    Info<< "\nInitializing Taylor-Greem Vortex flow with Re = "<< Re.value() << "\n" << endl;
    U = Ustart;
    U.correctBoundaryConditions();
    Info<< "Writing field U" << endl;
    U.write();

    p = pstart - pRefValueAna;
    p.correctBoundaryConditions();
    Info<< "Writing field p" << endl;
    p.write();
}
