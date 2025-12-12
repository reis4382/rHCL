#include <R.h>
#include <math.h>

static double parms[16];  // define parameter inputs
static double forc[1];    // define forcing input

/* Define additional variables to help track parameters and forcings. Note: input order must match */
/* Sleep and wake parameters */
#define mu parms[0] 		  // upper asymptote for sleep pressure: note zero-indexing
#define chi parms[1] 		  // time constant for sleep pressure decay/rise
#define Hzero parms[2] 		// mean level of wake propensity rhythm
#define delta parms[3] 		// separation between thresholds
#define ca_par parms[4] 	// circadian wake propensity amplitude

/* Circadian parameters */
#define tau_c parms[5] 		// circadian period
#define f_par parms[6] 		// correction factor for oscillator
#define G_par parms[7] 		// gain factor determining magnitude of light effect
#define p_par parms[8] 		// light sensitivity
#define k_par parms[9] 		// determines relative effect of light on oscillator variables
#define little_b_par parms[10] 	// sensitivity modulation factor
#define gamma parms[11] 	// stiffness of van der Pol oscillator
#define alpha_zero parms[12]	  // magnitude of effect of light on fraction of activated photoreceptors
#define beta parms[13]		// decay rate of fraction of activated photoreceptors
#define Izero parms[14]		// scaling factor for light
#define time_scale parms[15]		// time scaling factor

# define Itilde forc[0] 	// Interpolated forcing value

/* initializers */
void parmsc_p(void (* odeparms)(int *, double *))
{
	int N=16;
	odeparms(&N, parms);
}

// Initialize forcings if not using enforced wake periods
void forcc_p(void (* odeforcs)(int *, double *))
{
	int N=1;
	odeforcs(&N, forc);
}

/* derivative function */
// Derivatives if not using enforced wake periods //
// ref: y[0] = h; y[1] = n; y[2] = x; y[3] = y; y[4] = S; Order of input for ode (y*) variables
void derivsc_p(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
	// Gate light (i.e., set to 0 lux) if sleeping
	Itilde = (1 - y[4]) * Itilde; // eq.5: y[4] = 1 if sleeping

	// Auxiliary variables (makes code more readable). Make sure to define variables.
	double beta_hat = G_par * alpha_zero * pow(Itilde / Izero, p_par) * (1 - y[1]); // eq. 7
	double B_par = (1 - little_b_par * y[2]) * (1 - little_b_par * y[3]) * beta_hat; // eq. 10

	/* Derivative formulas */
	// ref: ydot[0] = dhdt; ydot[1] = dndt; ydot[2] = dxdt; ydot[3] = dydt; ydot[4] = dsleepdt

	// dhdt - derivative of sleep pressure
	ydot[0] = (-y[0] + (1 - y[4]) * mu) / chi; // eq. 2

	// dndt - derivative of photoreceptor activation
  // Note addition of leading "60*", which is present in Forger 1999
  // Skeldon 2023 appear to drop this and incorporate it into default parameters
  // (alpha_zero and beta), but that seems like it would introduce an error.
  // Correcting formula here.
	ydot[1] = (60*(alpha_zero * pow(Itilde / Izero, p_par) * (1 - y[1]) - beta * y[1])) / time_scale; // eq. 6


	// dxdt - derivative of x (I believe this is xc in forger 1999)
	ydot[2] = (gamma * (y[2] - (4 * pow(y[2], 3) / 3)) - y[3] * (pow(24 / (f_par * tau_c), 2) + k_par * B_par)) / (12/M_PI * time_scale); // eq. 8

	// dydt - deriative of y (I believe this is x in forger 1999)
	ydot[3] = (y[2] + B_par) / (12/M_PI * time_scale); // eq. 9

	// dsleepdt - "derivative" of sleep state variable (always 0 b/c it doesn't change dynamically, only during root function)
	ydot[4] = 0;

}


/* Root finding function (for when sleep homeostasis crosses threshold)*/
// This function identifies a root
void rootc_p(int *neq, double *t, double *y, int *ng, double *gout, double *out, int *ip)
{

	/* First, calculate wake propensity score */
	// Scalars are taken from appendix of skeldon 2023 paper
	// TODO: Could modify this to allow different values based on user input
	double circ_val = 0.7896 +
	                  -0.3912 * y[2] +
	                  0.7583 * y[3] +
	                  -0.4442 * pow(y[2], 2) +
	                  0.0250 * y[2] * y[3] +
	                  -0.9647 * pow(y[3], 2); // eq. 11

	/* Determine current threshold based on sleep/wake state */
	double h_thresh; // define variable
	// Use conditional logic
	if(y[4] == 0){
		h_thresh = Hzero + 0.5 * delta + ca_par * circ_val; // eq. 3; threshold for sleep if awake
	} else if (y[4] == 1){
		h_thresh = Hzero - 0.5 * delta + ca_par * circ_val; // eq. 4; threshold for wake if asleep
	}

	// return value of current sleep pressure minus threshold
	gout[0] = y[0] - h_thresh;
}

/* Event function - triggers when root is found (switches sleep and wake) */
void eventc_p(int *n, double *t, double *y)
{
  // If currently asleep, swtich to wake, and vice versa
  if (y[4] == 1){
    y[4] = 0; // If asleep, switch to wake
  } else if (y[4] == 0){
    y[4] = 1; // If awake, switch to sleep
  }
}

/* TODO Determine if there is a way to add enforced wake to C code. The root
 * finding function causes problems with the forcing functions - updates do
 * not seem to work properly around events. I cannot figure out why this is.
 * Essentially, the derivatives work, but during the root it seems to search
 * around different time values and causes all enforced wakes to be on after
 * the first root trigger, even when forced wake should not be on.
*/
