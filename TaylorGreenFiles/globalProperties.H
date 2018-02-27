// Dissipation rate normalized
//epsilon =  (0.5*nu*average( magSqr(fvc::grad(U) + fvc::grad(U)().T())))/(Uinit*Uinit*Uinit/L);
epsilon = 0.5*nu*magSqr(fvc::grad(U) + fvc::grad(U)().T())().weightedAverage(mesh.V())/(Uinit*Uinit*Uinit/L);

// Kinetic enery normalized
//Ek = 0.5 * average(magSqr(U)) / (Uinit*Uinit);
Ek = 0.5*magSqr(U)().weightedAverage(mesh.V())/ (Uinit*Uinit);

if (Pstream::parRun() && Pstream::master())
{
    if(TaylorGreenDimension =="2D")
    {   
        globalPropertiesFilePar 
        << runTime.value()/tc.value() << "\t" << Ek.value() << "\t" << epsilon.value() 
        << "\t" << UL2err.value() << "\t" << pL2err.value() 
        << "\t" << ULinfErr.value() << "\t" << pLinfErr.value()<< endl;
    }
    else if (TaylorGreenDimension =="3D")
    {
        globalPropertiesFilePar 
        << runTime.value()/tc.value() << "\t" << Ek.value() << "\t" << epsilon.value()<< endl;
    }
} 
else {
    if(TaylorGreenDimension =="2D")
    {   
        globalPropertiesFile 
        << runTime.value()/tc.value() << "\t" << Ek.value() << "\t" << epsilon.value() 
        << "\t" << UL2err.value() << "\t" << pL2err.value() 
        << "\t" << ULinfErr.value() << "\t" << pLinfErr.value()<< endl;
    }
    else if (TaylorGreenDimension =="3D")
    {
        globalPropertiesFile 
        << runTime.value()/tc.value() << "\t" << Ek.value() << "\t" << epsilon.value()<< endl;
    }    
}
  


