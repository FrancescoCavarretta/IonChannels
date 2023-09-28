:Reference : :		Adams et al. 1982 - M-currents and other potassium currents in bullfrog sympathetic neurones
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX Im
	USEION k READ ek WRITE ik
	RANGE gbar, vhm, km
        GLOBAL vshm, pkm
        GLOBAL mtau_min, mtau_max_var, mtau_sh, mtau_k_var1, mtau_k_var2
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2)

        vhm = -35 (mV)
        km = 5 (/mV)
        vshm = 0 (mV)
        pkm = 0

        mtau_min = 0 (ms)
        mtau_max_var = 0
        mtau_sh = 0 (mV)
        mtau_k_var1 = 0
        mtau_k_var2 = 0
}

ASSIGNED	{
        celsius (degC)
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	mInf
	mTau    (ms)
	mAlpha
	mBeta

        minf_vh (mV)
        minf_k  (/mV)
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	ik = gbar*m*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
        minf_vh = vhm + vshm
        minf_k = km * (1 + pkm)
                
	rates()
	m = mInf
}

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((celsius-21)/10)

  UNITSOFF
  mInf = 1 / (1 + exp(-(v - minf_vh) / minf_k))
  
  mAlpha = 3.3e-3 * exp( (v - (-35 + mtau_sh)) / (10 * (1 + mtau_k_var1)))
  mBeta  = 3.3e-3 * exp(-(v - (-35 + mtau_sh)) / (10 * (1 + mtau_k_var2)))
  mTau = ( mtau_min + (1 + mtau_max_var) / (mAlpha + mBeta) ) / qt
  UNITSON
}
