:Comment : LVA ca channel. Note: mtau is an approximation from the plots
:Reference : :		Avery and Johnston 1996, tau from Randall 1997
:Comment: shifted by -10 mv to correct for junction potential
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX Ca_LVAst
	USEION ca READ eca WRITE ica
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

        vhm = -40 (mV)
        km = 6 (/mV)
        vshm = 0 (mV)
        pkm = 0

        vhh = -90 (mV)
        kh = 6.4 (/mV)
        vshh = 0 (mV)
        pkh = 0
}

ASSIGNED	{
        celsius (degC)
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	mInf
	mTau    (ms)
	hInf
	hTau    (ms)

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
	ica = gbar*m*m*h*(v-eca)
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

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((celsius-21)/10)

  UNITSOFF
  mInf = 1 / (1 + exp(-(v - minf_vh) / minf_k))
  hInf = 1 / (1 + exp((v - hinf_vh) / hinf_k))
  
		mTau = (5.0000 + 20.0000/(1+exp((v + 10 - -25.000)/5)))/qt
		hTau = (20.0000 + 50.0000/(1+exp((v + 10 - -40.000)/7)))/qt
  UNITSON
}
