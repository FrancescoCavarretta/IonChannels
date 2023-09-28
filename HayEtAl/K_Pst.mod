:Comment : The persistent component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21


NEURON	{
	SUFFIX K_Pst
	USEION k READ ek WRITE ik
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

        vhm = -11 (mV)
        km = 12 (/mV)
        vshm = 0 (mV)
        pkm = 0

        vhh = -64 (mV)
        kh = 11 (/mV)
        vshh = 0 (mV)
        pkh = 0
}

ASSIGNED	{
        celsius (degC)
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
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
	ik = gbar*m*m*h*(v-ek)
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
  
        if(v<-50){
		    mTau =  (1.25+175.03*exp(-(v + 10) * -0.026))/qt
        }else{
            mTau = ((1.25+13*exp(-(v + 10)*0.026)))/qt
        }
		hTau =  (360+(1010+24*(v + 10+55))*exp(-((v + 10+75)/48)^2))/qt
  UNITSON
}
