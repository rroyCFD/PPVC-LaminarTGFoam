// create error fields----------------------------//
// if(TaylorGreenDimension =="2D")
// {   
    #include "vector.H"

    volVectorField Uerror
    (
        IOobject
        (
            "Uerror",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE // AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Uerror",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
    );

    volScalarField perror
    (
        IOobject
        (
            "perror",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE // AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("perror",dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0)
    );


    dimensionedScalar ULinfErr = max(mag(Uerror));          // infinity norm
    dimensionedScalar UL2err   = sqrt(sum(magSqr(Uerror))); // 2nd norm

    dimensionedScalar pLinfErr = max(mag(perror));          // infinity norm
    dimensionedScalar pL2err   = sqrt(sum(magSqr(perror))); // 2nd norm

    volVectorField Ustart
    (
        IOobject
        (
            "Ustart",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE // AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Ustart",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
    );

    volScalarField pstart
    (
        IOobject
        (
            "pstart",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE // AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pstart",dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0)
    );
// }
