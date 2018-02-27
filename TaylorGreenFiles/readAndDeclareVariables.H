// Read the dimension of the case
    const word TaylorGreenDimension
    (
        transportProperties.lookupOrDefault<word>("TaylorGreenDimension", "2D")   
    );
    Info << "TaylorGreenDimension: " << TaylorGreenDimension << endl;
    //  Read initial Velocity V_0
    dimensionedScalar Uinit
    (
        // "Uinit",
        // dimVelocity,    
        // transportProperties.lookupOrDefault<scalar>("Uinit", 1.0)
        transportProperties.lookup("Uinit")
    );
    Info << "Velocity Scale: " << Uinit.value() << endl;

    //  Read length scale
    dimensionedScalar L
    (
        // "L",
        // dimLength,
        // transportProperties.lookupOrDefault<scalar>("L", 1.0)
        transportProperties.lookup("L")
    );
    Info << "Length Scale: " << L.value() << endl;

    // Declare variables -------------------------------------------- //    
    dimensionedScalar epsilon
    (
        "epsilon",
        dimensionSet(0, 0, 0, 0, 0),
        0.0
    );
    
    dimensionedScalar Ek
    (
        "Ek",
        dimensionSet(0, 0, 0, 0, 0),
        0.0
    );
    
    // characteristic time
    dimensionedScalar tc
    (
        "tc",
        dimensionSet(0, 0, 1, 0, 0),
        0.0
    );
    tc = L/Uinit;

    // create output file
    OFstream globalPropertiesFile
    (
      runTime.path()/("TaylorGreenVortexProperties.dat_"+runTime.timeName())
    );

    OFstream globalPropertiesFilePar
    (
      runTime.path()/"../"/("TaylorGreenVortexPropertiesParallel.dat_"+runTime.timeName())
    );
    
    // decide on, which output file to use
    if (Pstream::parRun() && Pstream::master())
    { 
        // write to file
        if(TaylorGreenDimension =="2D")
        {
            globalPropertiesFilePar << "time" << "\t" << "Ek" << "\t" << "epsilon\t\t" 
            << "L2(U)\t" << "L2(p)\t" << "Linf(U)\t" << "Linf(p)" << endl;
        }
        else if (TaylorGreenDimension =="3D")
        {
            globalPropertiesFilePar << "time" << "\t" << "Ek" << "\t" << "epsilon" << endl;
        }
        else {
            globalPropertiesFilePar << "TG Dimension mismatch: check for dimension in the transportproperties file" << endl;
        }
    } 
    else {
        // write to file
        if(TaylorGreenDimension =="2D")
        {
            globalPropertiesFile << "time" << "\t" << "Ek" << "\t" << "epsilon\t\t" 
            << "L2(U)\t" << "L2(p)\t" << "Linf(U)\t" << "Linf(p)" << endl;
        }
        else if (TaylorGreenDimension =="3D")
        {
            globalPropertiesFile << "time" << "\t" << "Ek" << "\t" << "epsilon" << endl;
        }
        else {
            globalPropertiesFile << "TG Dimension mismatch: check for dimension in the transportproperties file" << endl;
        }
    }

    


    