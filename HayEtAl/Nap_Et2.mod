:Comment : mtau deduced from text (said to be 6 times faster than for NaTa)
:Comment : so I used the equations from NaT and multiplied by 6
:Reference : Modeled according to kinetics derived from Magistretti & Alonso 1999
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX Nap_Et2
	USEION na READ ena WRITE ina
	RANGE gbar
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2)
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
  mInf = 1.0/(1+exp(-(v + 52.6)/4.6))
  mAlpha = 0.182 * 6 * efun(-(v + 38)/6)
  mBeta = 0.124 * 6 * efun((v + 38)/6)
  mTau = 6*(1/(mAlpha + mBeta))/qt
  hInf = 1.0/(1+exp((v + 48.8)/10))
  hAlpha = 2.88e-6 * 4.63 * efun((v + 17)/4.63)
  hBeta = 6.94e-6 * 2.63 * efun(-(v + 64.4)/2.63)
  hTau = (1/(hAlpha + hBeta))/qt
  UNITSON
}
