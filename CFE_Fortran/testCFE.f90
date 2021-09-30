! TEST PROGRAM - READ THE HEADER OF F_CFE.F90 FOR DETAILS

PROGRAM testCFE
    use CFEmod
    REAL pc,msun,myr,surfMW,qMW,omegaMW,cfe(1:4)
    
    pc=3.1e16 !parsec in meters
    Msun=2.e30 !solar mass in kg
    Myr=3.16e13 !million years in seconds
    surfMW=9.3*msun/pc**2. !surface density
    qMW=2.0 !Toomre Q
    omegaMW=0.026/Myr !angular velocity
    cfe=f_cfe(surfMW,qMW,omegaMW) !CFE
    PRINT*,cfe(1)*100. 
END PROGRAM


