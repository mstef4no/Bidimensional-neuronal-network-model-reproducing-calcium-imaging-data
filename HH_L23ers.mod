TITLE L23 regular spiking (excitatory) neuron HH channels

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}


NEURON {
	SUFFIX hhers
        NONSPECIFIC_CURRENT ina
        NONSPECIFIC_CURRENT ikdr
        NONSPECIFIC_CURRENT im
	RANGE gnabar, gkdrbar, gmbar, vt,taumax
       RANGE ena, ek
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}


PARAMETER {
	ena     = 50    (mV)
	ek      = -100   (mV)
	vt	= -55	(mV)
	gnabar  = 0.05  (mho/cm2)
        gkdrbar = 0.005 (mho/cm2)
        gmbar = 7e-05 (mho/cm2)
        taumax = 1000 (ms)
}

STATE {
n m h p
}

ASSIGNED {
	v	(mV)
	ina     (mA/cm2)
	ikdr    (mA/cm2)
	im    (mA/cm2)
	ninf	(/)
	minf	(/)
        hinf	(/)
        pinf	(/)
        taun	(ms)
        taum	(ms)
        tauh	(ms)
        taup	(ms)
}


BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gnabar*h*m*m*m*(v-ena)
    ikdr = gkdrbar*n*n*n*n*(v-ek)
    im = gmbar*p*(v-ek)
}



INITIAL {
evaluate_fct(v)
}

DERIVATIVE states {
    evaluate_fct(v)
    n' = (ninf - n) / taun
    m' = (minf - m) / taum
    h' = (hinf - h) / tauh
    p' = (pinf - p) / taup
}


PROCEDURE evaluate_fct(v(mV)) {

LOCAL a,b,v2,x,y

UNITSOFF

       v2 = v-vt

       x=(v2-13)
       y=4
       if(fabs(x/y)>1e-6){
       a = 0.32 * x / ( 1 - exp(-x/y))
       } else {
       a = 0.32 * y/(1-x/(2*y))
       }


       x=(v2-40)
       y=5
      if(fabs(x/y)>1e-6){
       b = -0.28 * x / ( 1 - exp(x/y))
      } else {
      b = -0.28 * y/(-1-x/(2*y))
      }

	minf = a / (a + b)
        taum = 1/(a + b)



	a = 0.128 * exp((17-v2)/18)
	b = 4 / ( 1 + exp((40-v2)/5) )

	hinf = a/(a + b)
	tauh = 1/(a + b)



       x=(v2-15)
       y=5
       if(fabs(x/y)>1e-6){
       a = 0.032 * x / ( 1 - exp(-x/y))
       } else {
       a = 0.032 * y/(1-x/(2*y))
       }


       b = 0.5 * exp((10-v2)/40)

	ninf = a/(a + b)
	taun = 1/(a + b)


	pinf=1/(1+exp(-(v+35)/10))
	taup=taumax/(3.3*exp((v+35)/20)+exp(-(v+35)/20))

UNITSON

}


