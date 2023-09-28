: SK-type calcium-activated potassium current
: Reference : Kohler et al. 1996

NEURON {
       SUFFIX SK_E2
       USEION k READ ek WRITE ik
       USEION ca READ cai
       
       RANGE gbar, chm, km, mtau_max
       GLOBAL cshm, pkm

       GLOBAL mtau_max_var
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
}

PARAMETER {
        gbar = .000001 (mho/cm2)
        
        chm = -7.752 (1)
        km = 0.208 (1)
        cshm = 0 (mV)
        pkm = 0

        mtau_max = 1 (ms)
        mtau_max_var = 0
}

ASSIGNED {
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
        cai     (mM)
        
	mInf	(1)
	mTau	(ms)

        minf_ch (mM)
        minf_k  (/mM)
}

STATE {
      m   FROM 0 TO 1
}

INITIAL {
        minf_ch = chm + cshm
        minf_k = km * (1 + pkm)
        
        rates()
        m = mInf
}

BREAKPOINT {
           SOLVE states METHOD cnexp
           ik   =  gbar * m * (v - ek)
}

DERIVATIVE states {
        rates()        
        m' = (mInf - m) / mTau
}

PROCEDURE rates() {
  UNITSOFF
  mTau = mtau_max * (1 + mtau_max_var)
  mInf = 1 / (1 + exp(-(log(cai) - minf_ch) / minf_k))
  UNITSON
}


