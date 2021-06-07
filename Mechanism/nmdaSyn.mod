COMMENT

Author: Mark Cembrowski, 2012

This is an extension of the Exp2Syn class to incorporate NMDA-like properties,
and incorporates some NMDA features from Elena Saftenku, 2001.

First, Exp2Syn is described:

Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

Next, two extensions have been included:
1.  A switch whether to control whether the condutance is on (isOn)
2.  Voltage gating, mimicking Mg block

ENDCOMMENT

NEURON {
	POINT_PROCESS Exp2SynNmda
	RANGE tau1, tau2, e, i, mgBlock, extMgConc
	NONSPECIFIC_CURRENT i
	RANGE isOn

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	e=0	(mV)
	alpha_vspom = -0.062 (/mV) :voltage-dependence of Mg2+ block from Maex and De Schutter 1998
	v0_block = 10 (mV)
	extMgConc = 1 (mM) : external Mg concentration
	isOn = 0
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	mgBlock
	:extMgConc (mM)
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	mgBlock = vspom(v)
	i = isOn*g*mgBlock*(v - e)

}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}

FUNCTION vspom (v(mV))( ){
	vspom=1./(1.+0.2801*extMgConc*exp(alpha_vspom*(v-v0_block)))
}
