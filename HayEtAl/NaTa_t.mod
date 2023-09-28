:Reference :Colbert and Pan 2002

NEURON	{
	SUFFIX NaTa_t
	USEION na READ ena WRITE ina
	RANGE gbar, vhm, km, vhh, kh
        GLOBAL vshm, pkm, vshh, pkh
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2)

        vhm = -40.302 (mV)
        km = 6 (/mV)
        vshm = 0 (mV)
        pkm = 0

        vhh = -66 (mV)
        kh = 6 (/mV)
        vshh = 0 (mV)
        pkh = 0
}

ASSIGNED	{
        celsius (degC)
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	mInf
	mTau    (ms)
	mAlpha
	mBeta
	hInf
	hTau    (ms)
	hAlpha
	hBeta
        
        minf_vh (mV)
        minf_k  (/mV)
        hinf_vh (mV)
        hinf_k  (/mV)
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	ina = gbar*m*m*m*h*(v-ena)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
        minf_vh = vhm + vshm
        minf_k = km * (1 + pkm)
        
        hinf_vh = vhh + vshh
        hinf_k = kh * (1 + pkh)
        
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
  LOCAL qt
  qt = 2.3^((celsius-21)/10)
	
  UNITSOFF
  mInf = 1 / (1 + exp(-(v - minf_vh) / minf_k))
  hInf = 1 / (1 + exp((v - hinf_vh) / hinf_k))
  
  mAlpha = 0.182 * 6 * efun(-(v + 38)/6)
  mBeta = 0.124 * 6 * efun((v + 38)/6)
  mTau = (1/(mAlpha + mBeta))/qt
  
  hAlpha = 0.015 * 6 * efun((v + 66) / 6) 
  hBeta  = 0.015 * 6 * efun(-(v + 66) / 6)
  hTau = (1/(hAlpha + hBeta))/qt
  
  UNITSON
}
