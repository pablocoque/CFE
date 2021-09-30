; IDL FUNCTIONS FOR CALCULATING THE CLUSTER FORMATION EFFICIENCY (LOCAL VERSION)
; Copyright (C) 2012  Diederik Kruijssen
;
; NAME:
;       F_CFElocal (function)
;
; PURPOSE:
;       Calculate the fraction of star formation occurring in bound clusters,
;       i.e. the cluster formation efficiency or CFE
;
; TERMS OF USE:
;       If you use this routine while preparing a paper, please cite the paper
;       in which this model was presented:
;       Kruijssen, J. M. D., 2012, MNRAS 426, 3008
;
; CALLING SEQUENCE:
;       f_cfelocal,rholoc,sigmaloc,csloc,sflaw,qvir,tsn,tview,surfGMC,ecore,beta0,radfb
;
; INPUT PARAMETERS:
;       NOTE: ALL INPUT SHOULD BE IN SI UNITS!
;       rholoc - local gas volume density (required)
;       sigmaloc - local 1D gas velocity dispersion (required)
;       csloc - local gas sound speed (required)
;       NOTE: ALL INPUT BELOW IS OPTIONAL
;       WHEN SETTING ONE PARAMETER, ALL PRECEDING ONES SHOULD BE SET MANUALLY TOO
;       sflaw - star formation law: 0=Elmegreen (2002, default)
;                                   1=Krumholz & McKee (2005)
;       qvir - giant molecular cloud virial parameter (default=1.3)
;       tsn - time of the first supernova (default=3.*Myr)
;       tview - time at which CFE is determined (default=10.*Myr)
;       surfGMC - giant molecular cloud surface density (default=100.*Msun/pc^2.)
;       ecore - maximum (protostellar core) star formation efficiency (default=.5)
;       beta0 - turbulent-to-magnetic pressure ratio (default=1.d10)
;       radfb - feedback: 0=supernova, 1=radiative, 2=both (default=0)
;
; OUTPUT:
;       cfearray - array with four elements, the 1st being the CFE, the 2nd the
;       naturally bound fraction of star formation, the third the fraction of
;       bound star formation that survives the cruel cradle effect, and the 4th
;       the fraction of all star formation that survives the cruel cradle effect
;
; EXAMPLE:
;       Calculate the CFE in the solar neighbourhood
;
;       IDL> msun=2.d30
;       IDL> pc=3.1d16
;       IDL> rhoMW=0.03*msun/pc^3.
;       IDL> sigmaMW=7.d3
;       IDL> csMW=0.2d3
;       IDL> cfelocal=f_cfelocal(rhoMW,sigmaMW,csMW) 
;       IDL> print,cfelocal(0)*100.
;              8.8275440
;
; FILE STRUCTURE: 
;       Auxiliary functions come first, the actual CFE function is located at
;       the bottom of this file
;
; REVISION HISTORY:
;       Written by Diederik Kruijssen, August 2012
;

function f_tff,rho ;free-fall time as a function of density
    G=6.67d-11 ;gravitational constant
    tff=sqrt(3.*!pi/32./G/rho) ;free-fall time
    return,tff
end

function f_mach,sigmaloc,csloc ;Mach number as a function of velocity dispersion and sound speed
    mach=sigmaloc/csloc ;Mach number
    return,mach
end

function f_sigrho,mach,beta0 ;dispersion of overdensity PDF as a function of Mach number and magnetic pressure ratio, based on e.g. Padoan & Nordlund (2011)
    b=0.5 ;constant
    sig=sqrt(alog(1.+3.*b^2.*mach^2.*beta0/(beta0+1.))) ;dispersion
    return,sig
end

function f_xcrit,qvir,mach ;critical overdensity for star formation in the KM05 sSFR_ff as a function of Mach number and cloud virial ratio, from Krumholz & McKee (2005)
    phix=1.12 ;constant
    x=!pi^2.*phix^2./15.*qvir*mach^2. ;critical overdensity
    return,x
end

function f_sfrff,qvir,mach,beta0,ecore,sflaw ;specific star formation rate per free-fall time (sSFR_ff) for Elmegreen (2002, sflaw=0) or Krumholz & McKee (2005, sflaw=1)
    phit=1.91 ;constant
    xcrit=f_xcrit(qvir,mach) ;critical overdensity for star formation in the KM05 model, see above
    sigrho=f_sigrho(mach,beta0) ;overdensity PDF dispersion, see above
    if sflaw then f=.5*ecore/phit*(1.+erf((-2.*alog(xcrit)+sigrho^2.)/(2.^1.5*sigrho))) else f=0.012 ;specific star formation rate per free-fall time
    return,f
end

function f_dpdx,x,mulnx,sig ;overdensity PDF of the ISM as a function of overdensity x and its logarithmic mean and dispersion
    dpdx=1./(sqrt(2.*!pi*sig^2.)*x)*exp(-.5*(alog(x)-mulnx)^2./sig^2.) ;overdensity PDF
    return,dpdx
end

function f_surfg,rholoc,sigmaloc ;rough relation between gas surface density and volume density and velocity dispersion, assuming an equilibrium disc
    G=6.67d-11 ;gravitational constant
    phiP=3. ;constant
    surfg=sqrt(2.*rholoc*sigmaloc^2./(!pi*G*phiP)) ;gas surface density
    return,surfg
end

function f_fstar,rholoc,sigmaloc,csloc,x,ecore,beta0,qvir,tsn,tview,surfGMC,sflaw,radfb ;naturally bound fraction of star formation
    G=6.67d-11 ;gravitational constant
    sigSB=5.67d-8 ;Stefan-Boltzmann constant
    c=299792458. ;speed of light
    phifb=1.6d-5 ;feedback efficiency
    kappa0=2.4d-5 ;opacity constant
    psi=.3 ;light-to-mass ratio
    phitrap=.2 ;trapping ratio
    surfg=f_surfg(rholoc,sigmaloc) ;estimate of gas surface density
    surffb=max([surfGMC,surfg]) ;surface density on which radiative feedback acts
    mach=f_mach(sigmaloc,csloc) ;Mach number, see above
    sfrff=f_sfrff(qvir,mach,beta0,ecore,sflaw) ;specific star formation rate per free-fall time, see above
    rhog=x*rholoc ;gas volume density
    tff=f_tff(rhog) ;free-fall time, see above
    ;Star Formation Efficiencies
    if radfb eq 0 or radfb eq 2 then efb=0.5*sfrff*tsn/tff*(1.+sqrt(1.+4.*tff*sigmaloc^2./(phifb*sfrff*tsn^2.*x))) else efb=1. ;SN feedback
    einc=sfrff*tview/tff ;star formation is incomplete/still ongoing
    if radfb gt 0 then efbrad=2.*sigSB/(phitrap*kappa0^2.*psi*surffb^3.)*(sqrt(1.+2.*!pi*c*G*phitrap*kappa0^2.*surffb^4./(1.*sigSB))-1.) else efbrad=1. ;radiative feedback
    epsilons=[ecore,efb,efbrad,einc] ;SFEs for [maximum,SNfeedback,radiativefeedback,incomplete]
    fstar=min(epsilons) ;local SFE is the minimum of those
    return,fstar
end

function f_integrate,xsurv,mulnx,sig,rholoc,sigmaloc,csloc,ecore,beta0,qvir,tsn,tview,surfGMC,cce,sflaw,radfb ;obtain fractions from integrating the overdensity PDFs
    nx=1000 ;number of integration steps (checked to be sufficient for convergence)
    xmin1=exp(mulnx-5.*sig) ;minimum overdensity
    xmax=exp(mulnx+10.*sig) ;maximum overdensity
    if cce ge 1 and xsurv lt xmin1 then xmin1=xsurv ;if calculating the cruel cradle effect and critical overdensity below minimum, then adjust
    if cce ge 1 and xsurv gt xmax then xmax=xsurv ;if calculating the cruel cradle effect and critical overdensity below maximum, then adjust
    xarr=xmin1*(xmax/xmin1)^((dindgen(nx)+0.5)/nx) ;integration array
    f1=0. ;initialize integral
    for i=0,nx-1 do begin ;denominator integral
        dx=xarr(i)*((xmax/xmin1)^(1./(2.*nx))-(xmax/xmin1)^(-1./(2.*nx))) ;step size
        xg=xarr(i) ;overdensity
        fstar=f_fstar(rholoc,sigmaloc,csloc,xg,ecore,beta0,qvir,tsn,tview,surfGMC,sflaw,radfb) ;local SFE
        bound=fstar/ecore ;local bound fraction
        if cce eq 0 then bound=1. ;if not calculating the cruel cradle effect but the naturally bound fraction of SF, the denominator should contain all SF
        if cce eq 2 then bound=1. ;if calculating the cruel cradle effect with respect to all SF, the denominator should contain all SF
        int=bound*fstar*xarr(i) ;integral part 1
        dpdx=f_dpdx(xg,mulnx,sig) ;overdensity PDF, i.e. integral part 2
        f1+=int*dpdx*dx ;integral
    endfor
    if cce ge 1 then xmin2=xsurv else xmin2=xmin1 ;if calculating the cruel cradle effect set minimum overdensity to critical overdensity
    xarr=xmin2*(xmax/xmin2)^((dindgen(nx)+0.5)/nx) ;integration array
    f2=0. ;initialize integral
    for i=0,nx-1 do begin ;numerator integral
        dx=xarr(i)*((xmax/xmin2)^(1./(2.*nx))-(xmax/xmin2)^(-1./(2.*nx))) ;step size
        xg=xarr(i) ;overdensity
        fstar=f_fstar(rholoc,sigmaloc,csloc,xg,ecore,beta0,qvir,tsn,tview,surfGMC,sflaw,radfb) ;local SFE
        bound=fstar/ecore ;local bound fraction
        if cce eq 2 then bound=1.  ;if calculating the cruel cradle effect with respect to all SF, the numerator should contain all SF
        int=bound*fstar*xarr(i) ;integral part 1
        dpdx=f_dpdx(xg,mulnx,sig) ;overdensity PDF, i.e. integral part 2
        f2+=int*dpdx*dx ;integral
    endfor
    if f1 ne 0 then frac=f2/f1 else frac=0. ;numerator divided by denominator
    return,frac
end

function f_phit,qvir,x ;ratio of encounter timescale to energy dissipation timescale as a function of cloud virial ratio and overdensity
    phit=3.1*sqrt((qvir/1.3)*(x/1.d4)) ;ratio of encounter timescale to energy dissipation timescale
    return,phit
end

function f_phiad,qvir,x ;adiabatic correction as a function of cloud virial ratio and overdensity
    phiad=exp(-2.*f_phit(qvir,x)) ;adiabatic correction
    return,phiad
end

function f_xcce,sigmaloc,surfGMC,qvir,tview ;critical overdensity to remain bound despite the cruel cradle effect (also see Kruijssen et al. 2011)
    eta=2.*1.305*3.*!pi/64. ;for Plummer
    G=6.67d-11 ;gravitational constant
    g_close=1.5 ;close encounter correction
    phish=2.8 ;higher-order energy loss correction
    rh2r2av=.25 ;for Plummer
    f=0.7 ;fraction of injected energy that is used for unbinding the region
    ;SOLVE implicit relation for x_cce
    xmin=1.d-4 ;minimum x
    xmax=1.d8 ;maximum x
    nx=101 ;length of x array
    niter=0 ;number of elapsed iterations
    itermax=10 ;maximum number of iterations
    xfit=xmin^2./xmax ;initialisation of fitted x
    accuracy=1.d-6 ;desired logarithmic accuracy
    xfit0=xfit*accuracy ;initialisation of previously fitted x
    while abs(alog10(xfit/xfit0)) gt accuracy and niter lt itermax do begin ;while iteration does not give desired convergence, do
        xfit0=xfit ;previously fitted x
        xarr=10.^(dindgen(nx)/(nx-1.)*alog10(xmax/xmin)+alog10(xmin)) ;x array
        xarr2=dblarr(nx) ;right-hand side of equation
        for j=0,nx-1 do begin ;for all x
            phiad=f_phiad(qvir,xarr(j)) ;adiabatic correction, see above
            xarr2(j)=87.5*sqrt(!pi)*f*g_close*G*phish*surfGMC*phiad*tview/sigmaloc ;right-hand side of equation
        endfor
        jfit=min(where(abs(xarr-xarr2) eq min(abs(xarr-xarr2)))) ;index where x equals right-hand side of equation
        xfit=xarr(jfit) ;solution for x_cce
        xmin=xarr(jfit-1) ;new minimum x
        xmax=xarr(jfit+1) ;new maximum x
        niter+=1 ;increase number of elapsed iterations by 1
    endwhile
    if niter eq itermax then begin ;if we stopped due to reaching maximum number of iterations, then stop
        print,' no convergence, increase nx or itermax in the f_xcce subroutine'
        stop
    endif
    return,xfit
end

function f_cfelocal,rholoc,sigmaloc,csloc,sflaw,qvir,tsn,tview,surfGMC,ecore,beta0,radfb ;FUNCTION TO CALCULATE THE CLUSTER FORMATION EFFICIENCY
    ;SET CONSTANTS
    pc=3.086d16 ;parsec in meters
    Msun=1.989d30 ;solar mass in kg
    Myr=1.d6*86400.*365.25 ;million years in seconds
    ;SET UNSPECIFIED PARAMETERS
    if n_elements(sflaw) eq 0 then sflaw=0 ;if star formation law is not specified, set to Elmegreen (2002) - NOTE: set to 1 for Krumholz & McKee (2005)
    if n_elements(qvir) eq 0 then qvir=1.3 ;giant molecular cloud virial parameter
    if n_elements(tsn) eq 0 then tsn=3.*Myr ;time of the first supernova
    if n_elements(tview) eq 0 then tview=10.*Myr ;time at which CFE is determined
    if n_elements(surfGMC) eq 0 then surfGMC=100.*Msun/pc^2. ;giant molecular cloud surface density
    surfg=f_surfg(rholoc,sigmaloc) ;estimate of gas surface density
    surfGMC=max([surfGMC,surfg])
    if n_elements(ecore) eq 0 then ecore=0.5 ;maximum (protostellar core) star formation efficiency
    if n_elements(beta0) eq 0 then beta0=1.d10 ;if turbulent-to-magnetic pressure ratio is not specified, set to turbulent-only
    if n_elements(radfb) eq 0 then radfb=0 ;if feedback is not specified, set to supernovae only - NOTE: set to 1 for radiative only, and to 2 for both
    ;CALCULATE DERIVED PARAMETERS
    mach=f_mach(sigmaloc,csloc) ;Mach number
    sigrho=f_sigrho(mach,beta0) ;dispersion of overdensity PDF
    mulnx=-.5*sigrho^2. ;logarithmic mean of overdensity PDF
    ;CALCULATE F_BOUND
    fbound=f_integrate(0,mulnx,sigrho,rholoc,sigmaloc,csloc,ecore,beta0,qvir,tsn,tview,surfGMC,0,sflaw,radfb) ;naturally bound part of star formation
    ;CALCULATE F_CCE
    xsurv=f_xcce(sigmaloc,surfGMC,qvir,tview) ;critical overdensity to remain bound despite the cruel cradle effect
    fcce=f_integrate(xsurv,mulnx,sigrho,rholoc,sigmaloc,csloc,ecore,beta0,qvir,tsn,tview,surfGMC,1,sflaw,radfb) ;part of bound SF surviving the cruel cradle effect
    fcce2=f_integrate(xsurv,mulnx,sigrho,rholoc,sigmaloc,csloc,ecore,beta0,qvir,tsn,tview,surfGMC,2,sflaw,radfb) ;part of all SF surviving the cruel cradle effect
    ;CALCULATE CFE
    cfearray=[fbound*fcce,fbound,fcce,fcce2] ;array containing the cluster formation efficiency, fbound, fcce, and fcce2 (i.e. fcce with respect to all SF)
    return,cfearray
end


