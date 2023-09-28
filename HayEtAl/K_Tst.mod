:Comment : The transient component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX K_Tst
	USEION k READ ek WRITE ik
	RANGE gbar, vhm, km, vhh, kh
        GLOBAL vshm, pkm, vshh, pkh

        GLOBAL mtau_min, mtau_max_var, mtau_sh, mtau_k_var
        GLOBAL htau_min, htau_max_var, htau_sh, htau_k_var
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2)

        vhm = -10 (mV)
        km = 19 (/mV)
        vshm = 0 (mV)
        pkm = 0

        vhh = -76 (mV)
        kh = 10 (/mV)
        vshh = 0 (mV)
        pkh = 0

        mtau_min = 0 (ms)
        mtau_max_var = 0
        mtau_sh = 0 (mV)
        mtau_k_var = 0

        htau_min = 0 (ms)
        htau_max_var = 0
        htau_sh = 0 (mV)
        htau_k_var = 0
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
	ik = gbar*(m^4)*h*(v-ek)
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

  mTau =  (0.34 + mtau_min + 0.92 * (1 + mtau_max_var) * exp(-(( v - (-81 + mtau_sh)) / (59 * (1 + mtau_k_var))) ^ 2)) / qt
  hTau =  (8    + htau_min + 49   * (1 + htau_max_var) * exp(-(( v - (-83 + htau_sh)) / (23 * (1 + htau_k_var))) ^ 2)) / qt
  UNITSON
}
