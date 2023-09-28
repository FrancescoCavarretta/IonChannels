:Comment :
:Reference : :		Kole,Hallermann,and Stuart, J. Neurosci. 2006

NEURON	{
	SUFFIX Ih
	NONSPECIFIC_CURRENT ihcn
	RANGE gbar, vhh1, kh1, vhh2, kh2, a0h, b0h
        GLOBAL vshh, pkh1, pkh2
        GLOBAL htau_min, htau_max_var, htau_sh, htau_k_var1, htau_k_var2
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2) 
	ehcn =  -45.0 (mV)

        vhh1 = -154.9 (mV)
        vhh2 = 0 (mV)
        
        kh1 = 11.9 (/mV)
        kh2 = 33.1 (/mV)
        
        vshh = 0 (mV)
        pkh1 = 0
        pkh2 = 0
        a0h = 0.076517
        b0h = 0.193

        htau_min = 0 (ms)
        htau_max_var = 0
        htau_sh = 0 (mV)
        htau_k_var1 = 0
        htau_k_var2 = 0
}

ASSIGNED	{
	v	(mV)
	ihcn	(mA/cm2)
	hInf
	hTau    (ms)
	hAlpha
	hBeta

        hinf_vh1 (mV)
        hinf_k1  (/mV)
        hinf_vh2 (mV)
        hinf_k2  (/mV)
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
        hinf_vh1 = vhh1 + vshh
        hinf_k1 = kh1 * (1 + pkh1)

        hinf_vh2 = vhh2 + vshh
        hinf_k2 = kh2 * (1 + pkh2)
        
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
	hAlpha =  a0h * efun((v - hinf_vh1)/hinf_k1)
	hBeta  =  b0h * exp((v - hinf_vh2)/hinf_k2)
	hInf = hAlpha/(hAlpha + hBeta)

        hAlpha =  a0h * efun((v - vhh1 - htau_sh) / (kh1 * (1 + htau_k_var1)))
	hBeta  =  b0h *  exp((v  - vhh2 - htau_sh) / (kh2 * (1 + htau_k_var2)))
	hTau = htau_min + (1 + htau_max_var) / (hAlpha + hBeta)
	UNITSON
}
