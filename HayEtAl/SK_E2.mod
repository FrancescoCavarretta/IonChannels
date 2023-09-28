: SK-type calcium-activated potassium current
: Reference : Kohler et al. 1996

NEURON {
       SUFFIX SK_E2
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gbar
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
}

PARAMETER {
          v            (mV)
          gbar = .000001 (mho/cm2)
          ek           (mV)
          cai          (mM)
}

ASSIGNED {
         mInf
         mTau (ms)
         ik            (mA/cm2)
}

STATE {
      m   FROM 0 TO 1
}

BREAKPOINT {
           SOLVE states METHOD cnexp
           ik   =  gbar * m * (v - ek)
}

DERIVATIVE states {
        rates(cai)
        mTau = 1              
        m' = (mInf - m) / mTau
}

PROCEDURE rates(ca(mM)) {
  mInf = 1 / (1 + exp(-(log(ca) + 7.752) / 0.208))
}

INITIAL {
        rates(cai)
        m = mInf
}
