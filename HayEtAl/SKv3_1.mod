:Comment :
:Reference : :		Characterization of a Shaw-related potassium channel family in rat brain, The EMBO Journal, vol.11, no.7,2473-2486 (1992)

NEURON	{
	SUFFIX SKv3_1
	USEION k READ ek WRITE ik
	RANGE gbar, vhm, km
        GLOBAL vshm, pkm
        
        RANGE  mtau_max, mtau_k, mtau_k_var
        GLOBAL mtau_max_var, mtau_k_var, mtau_min, mtau_sh
        
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2)
        
        vhm = 18.7 (mV)
        km = 9.7 (/mV)
        vshm = 0 (mV)
        pkm = 0
        
        mtau_min = 0 (ms)
        mtau_max = 4.0 (ms)
        mtau_max_var = 0
        mtau_sh = 0 (mV)
        mtau_k = 44.14 (/mV)
        mtau_k_var = 0
        
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	mInf	(1)
	mTau	(ms)

        minf_vh (mV)
        minf_k  (/mV)
}

STATE	{ 
	m      (1)
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
	UNITSOFF
	mInf =  1 / (1 + exp( -(v - minf_vh) / minf_k ) )
        
	mTau =  mtau_min + mtau_max * (1 + mtau_max_var) / (1 + exp( - (v - (-46.560 + mtau_sh)) / (mtau_k * (1 + mtau_k_var))))
	UNITSON
}
