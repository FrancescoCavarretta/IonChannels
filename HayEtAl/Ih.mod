:Comment :
:Reference : :		Kole,Hallermann,and Stuart, J. Neurosci. 2006

NEURON	{
	SUFFIX Ih
	NONSPECIFIC_CURRENT ihcn
	RANGE gbar
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2) 
	ehcn =  -45.0 (mV)
}

ASSIGNED	{
	v	(mV)
	ihcn	(mA/cm2)
	hInf
	hTau    (ms)
	hAlpha
	hBeta
}

STATE	{ 
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	ihcn = gbar*h*(v-ehcn)
}

DERIVATIVE states	{
	rates()
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	h = hInf
}

FUNCTION efun(z) {
	 if (fabs(z) < 1e-4) {
	    efun = 1 - z/2
	 }else{
	    efun = z/(exp(z) - 1)
         }
}

PROCEDURE rates(){
	UNITSOFF
		hAlpha =  0.001*6.43*11.9*efun((v+154.9)/11.9)
		hBeta  =  0.001*193*exp(v/33.1)
		hInf = hAlpha/(hAlpha + hBeta)
		hTau = 1/(hAlpha + hBeta)
	UNITSON
}
