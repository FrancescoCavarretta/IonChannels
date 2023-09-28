:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca_HVA
	USEION ca READ eca WRITE ica
	RANGE gbar, vhm1, km1, vhh1, kh1, vhm2, km2, vhh2, kh2, a0m, a0h, b0m, b0h
        GLOBAL vshm, pkm1, pkm2, vshh, pkh1, pkh2
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2)

        a0m = 0.209
        a0h = 0.000457
        b0m = 0.94
        b0h = 0.0065
        vhm1 = -27 (mV)
        vhm2 = -75 (mV)
        
        km1 = 3.8 (/mV)
        km2 = 17 (/mV)
        
        vshm = 0 (mV)
        pkm1 = 0
        pkm2 = 0

        vhh1 = -13 (mV)
        vhh2 = -15 (mV)
        
        kh1 = 50 (/mV)
        kh2 = 28 (/mV)
        
        vshh = 0 (mV)
        pkh1 = 0
        pkh2 = 0
}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	mInf    (1)
	mTau    (ms)
	mAlpha  (1)
	mBeta   (1)
	hInf    (1)
	hTau    (ms)
	hAlpha  (1)
	hBeta   (1)

        minf_vh1 (mV)
        minf_k1  (/mV)
        minf_vh2 (mV)
        minf_k2  (/mV)

        hinf_vh1 (mV)
        hinf_k1  (/mV)
        hinf_vh2 (mV)
        hinf_k2  (/mV)
}

STATE	{ 
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	ica = gbar*m*m*h*(v-eca)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
        minf_vh1 = vhm1 + vshm
        minf_k1 = km1 * (1 + pkm1)

        minf_vh2 = vhm2 + vshm
        minf_k2 = km2 * (1 + pkm2)

        hinf_vh1 = vhh1 + vshh
        hinf_k1 = kh1 * (1 + pkh1)

        hinf_vh2 = vhh2 + vshh
        hinf_k2 = kh2 * (1 + pkh2)
        
	rates()
	m = mInf
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
		mAlpha =  a0m*efun(-(v - minf_vh1)/minf_k1)
    		mBeta  =  b0m*exp(-(v - minf_vh2)/minf_k2)
		mInf = mAlpha / (mAlpha + mBeta)
        
		hAlpha =  a0h * exp(-(v - hinf_vh1)/hinf_k1)
		hBeta  =  b0h / (exp(-(v - hinf_vh2)/hinf_k2)+1)
		hInf = hAlpha / (hAlpha + hBeta)
        
		mAlpha =  0.055*3.8*efun(-(v+27)/3.8)
    		mBeta  =  (0.94*exp((-75-v)/17))        
		mTau = 1/(mAlpha + mBeta)
        
		hAlpha =  (0.000457*exp((-13-v)/50))
		hBeta  =  (0.0065/(exp((-v-15)/28)+1))
		hTau = 1/(hAlpha + hBeta)
	UNITSON
}
