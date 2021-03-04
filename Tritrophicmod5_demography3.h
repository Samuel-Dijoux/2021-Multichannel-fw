/*
 Tritrophicmod5_demography3.h
 v3: Replacement of predator functional response by constant additional induced-mortality rate "Mup"						
*/
  
/*
*===========================================================================
*   SECTION 1: PROBLEM DIMENSIONS, NUMERICAL SETTINGS AND MODEL PARAMETERS
*===========================================================================
*/
// Dimension settings: Required
#define POPULATION_NR       2
#define STAGES              3
#define	I_STATE_DIM         2
#define	PARAMETER_NR       22

#define ENVIRON_DIM         5
#define INTERACT_DIM        5 

// Descriptive names of parameters in parameter array (at least two parameters are required)
char  *parameternames[PARAMETER_NR] =
{ "Lb", "Lj", "Lm", "Omega", "Imax", "Rh1", "Rh2", "Nu", "Rm", "Mub", "Rmax1", "Rmax2", "Rho1", "Rho2", "Lv", "A", "Th", "Epsilon", "Delta", "Alpha", "Bsratio", "Mup" };

// Default values of all parameters 
double	parameter[PARAMETER_NR] =
{ 7.0, 110.0, 300.0, 9.0E-6, 1.0E-4, 1.5E-5, 1.5E-5, 0.006, 0.003, 0.01, 5E-6, 3E-4, 0.1, 0.2, 27.0, 5000.0, 0.1, 0.5, 0.01, 0.0, 1.2, 0.9};

// Aliases definitions for all istate variables
#define AGE(n)                 istate[n][0]
#define LENGTH(n)              istate[n][1]

// Aliases definitions for all environment variables
#define R1                  E[0]
#define R2                  E[1]
#define P                   E[2]
#define B1                  E[3]
#define B2                  E[4]

// Aliases definitions for all parameters
#define LB                  parameter[ 0]  // Default: 7 mm				        Length at birth.
#define LJ                  parameter[ 1]  // Default: 110 mm			        Length at maturation
#define LM                  parameter[ 2]  // Default: 300 mm			        Maximum Length.

#define OMEGA               parameter[ 3]  // Default: 9.0E-6 g/mm^3	  	Proportiononality constant.

#define IMAX                parameter[ 4]  // Default: 1.0E-4 g/day/mm?		Proportiononality constant.
#define RH1                 parameter[ 5]  // Default: 1.5E-5 g/L		      Half-saturation constant of the resource R1.
#define RH2                 parameter[ 6]  // Default: 1.5E-5 g/L		      Half-saturation constant of the resource R2.

#define NU                  parameter[ 7]  // Default: 0.006 /day		      Growth rate.
#define RM                  parameter[ 8]  // Default: 0.003 /day/mm?	    Proportiononality constant. 

#define MUB                 parameter[ 9]  // Default: 0.01 /day		      Mortality rate.

#define RMAX1               parameter[10]  // Default: 3.0E-4 g/L		      Carrying capacity K of R1.
#define RMAX2               parameter[11]  // Default: 3.0E-4 g/L		      Carrying capacity K of R2.
#define RHO1                parameter[12]  // Default: 0.1 /day
#define RHO2				        parameter[13]  // Default: 0.2 /day	

#define LV                  parameter[14]  // Default: 27 mm			        Length treshold of predation vulnerability
#define A                   parameter[15]  // Default: 5000.0 L/day		  	Attack rate.
#define TH                  parameter[16]  // Default: 0.1 day/g		      Handling time.
#define EPSILON             parameter[17]  // Default: 0.5				        Conversion efficiency.
#define DELTA               parameter[18]  // Default: 0.01 /day		      Predator Mortality rate.

#define ALPHA               parameter[19]  // Default: 0				          Resource consumption partitioning (Alpha = 0 or 1: specific consumers; 0<Alpha<1: generalists consumers).
#define BSRATIO             parameter[20]  // Default: 1.2				        Body size ratio btw B1 and B2.
#define MUP                 parameter[21]  // Default: 0				          Predation induced Mortality.
/*
  *===========================================================================
  * 	SECTION 2: DEFINITION OF THE INDIVIDUAL LIFE HISTORY
*===========================================================================
  */
  
  /*
  * Specify the number of states at birth for the individuals in all structured
* populations in the problem in the vector BirthStates[].
*/
  
  void SetBirthStates(int BirthStates[POPULATION_NR], double E[])
{
  BirthStates[0] = 1;
  BirthStates[1] = 1;
  return;
  }


/*
  * Specify all the possible states at birth for all individuals in all
* structured populations in the problem. BirthStateNr represents the index of
* the state of birth to be specified. Each state at birth should be a single,
* constant value for each i-state variable.
*
  * Notice that the first index of the variable 'istate[][]' refers to the
* number of the structured population, the second index refers to the
* number of the individual state variable. The interpretation of the latter
* is up to the user.
*/
  
  void StateAtBirth(double *istate[POPULATION_NR], int BirthStateNr, double E[])
{
  AGE(0) = 0.0;
  AGE(1) = 0.0;
  LENGTH(0) = LB;
  LENGTH(1) = LB * BSRATIO;
  return;
  }


/*
  * Specify the threshold determining the end point of each discrete life
* stage in individual life history as function of the i-state variables and
* the individual's state at birth for all populations in every life stage.
*
* Notice that the first index of the variable 'istate[][]' refers to the
* number of the structured population, the second index refers to the
* number of the individual state variable. The interpretation of the latter
* is up to the user.
*/

void IntervalLimit(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
double limit[POPULATION_NR])
{
  switch (lifestage[0])
  {
  case 0:
  limit[0] = LENGTH(0) - LV;
  break;
  case 1:
  limit[0] = LENGTH(0) - LJ;
  break;
  }
  
  switch (lifestage[1])
  {
  case 0:
  limit[1] = LENGTH(1) - LV;
  break;
  case 1:
  limit[1] = LENGTH(1) - (BSRATIO*LJ);
  break;
  }
  return;
}


/*
* Specify the development of individuals as a function of the i-state
* variables and the individual's state at birth for all populations in every
* life stage.
*
  * Notice that the first index of the variables 'istate[][]' and 'development[][]'
* refers to the number of the structured population, the second index refers
* to the number of the individual state variable. The interpretation of the
* latter is up to the user.
*/
  
  void Development(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
                   double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
                   double development[POPULATION_NR][I_STATE_DIM])
{
  development[0][0] = 1.0;
  development[1][0] = 1.0;
  development[0][1] = NU*(LM*((1-ALPHA)*R1/(R1 + RH1) + ALPHA*R2/(R2+RH2)) - LENGTH(0));
  development[1][1] = NU*((BSRATIO*LM)*((1-ALPHA)*R2/(R2 + RH2) + ALPHA*R1/(R1+RH1)) - LENGTH(1));
  
  return;
  }


/*
  * Specify the possible discrete changes (jumps) in the individual state
* variables when ENTERING the stage specified by 'lifestage[]'.
*
  * Notice that the first index of the variable 'istate[][]' refers to the
* number of the structured population, the second index refers to the
* number of the individual state variable. The interpretation of the latter
* is up to the user.
*/
  
  void DiscreteChanges(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
                       double *birthstate[POPULATION_NR], int BirthStateNr, double E[])
{
  return;
  }


/*
  * Specify the fecundity of individuals as a function of the i-state
* variables and the individual's state at birth for all populations in every
* life stage.
*
* The number of offspring produced has to be specified for every possible
* state at birth in the variable 'fecundity[][]'. The first index of this
* variable refers to the number of the structured population, the second
* index refers to the number of the birth state.
* Notice that the first index of the variable 'istate[][]' refers to the
* number of the structured population, the second index refers to the
*
* number of the individual state variable. The interpretation of the latter
* is up to the user.
*/

void Fecundity(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
double *fecundity[POPULATION_NR])
{
  for(int i=0; i<2; i++)	{
  fecundity[i][0] = 0.0;
  if (lifestage[i] == 2) {
  if (i==0)
  fecundity[i][0] = RM*((1-ALPHA)*R1/(R1 + RH1) + ALPHA*R2/(R2 + RH2)) * LENGTH(i)*LENGTH(i);
  else if (i==1)
  fecundity[i][0] = RM*((1-ALPHA)*R2/(R2 + RH2) + ALPHA*R1/(R1 + RH1)) * LENGTH(i)*LENGTH(i);
  else 
  printf("Population index i outside bounds\n");
  }
  }
  return;
}


/*
* Specify the mortality of individuals as a function of the i-state
* variables and the individual's state at birth for all populations in every
* life stage.
*
  * Notice that the first index of the variable 'istate[][]' refers to the
* number of the structured population, the second index refers to the
* number of the individual state variable. The interpretation of the latter
* is up to the user.
*/
  
  void Mortality(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
                 double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
                 double mortality[POPULATION_NR])
{
  for(int i=0; i<2; i++)  {
    if (lifestage[i] == 0)
      mortality[i] = MUB + MUP;
    else
      mortality[i] = MUB;
  }
  return;
  }


/*
  *===========================================================================
  * 	SECTION 3: FEEDBACK ON THE ENVIRONMENT
*===========================================================================
  */
  
  /*
  * For all the integrals (measures) that occur in interactions of the
* structured populations with their environments and for all the integrals
* that should be computed for output purposes (e.g. total juvenile or adult
                                               * biomass), specify appropriate weighing function dependent on the i-state
* variables, the individual's state at birth, the environment variables and
* the current life stage of the individuals. These weighing functions should
* be specified for all structured populations in the problem. The number of
* weighing functions is the same for all of them.
*
* Notice that the first index of the variables 'istate[][]' and 'impact[][]'
* refers to the number of the structured population, the second index of the
* variable 'istate[][]' refers to the number of the individual state variable,
* while the second index of the variable 'impact[][]' refers to the number of
* the interaction variable. The interpretation of these second indices is up
* to the user.
*/

void Impact(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
double impact[POPULATION_NR][INTERACT_DIM])
{
  
  impact[0][0] = IMAX*(1-ALPHA)*(R1/(R1 + RH1))*LENGTH(0)*LENGTH(0); 
  impact[1][0] = IMAX*(ALPHA)*(R1/(R1 + RH1))*LENGTH(1)*LENGTH(1);
  impact[0][1] = IMAX*(ALPHA)*(R2/(R2 + RH2))*LENGTH(0)*LENGTH(0); 
  impact[1][1] = IMAX*(1-ALPHA)*(R2/(R2 + RH2))*LENGTH(1)*LENGTH(1);
  
  for(int i=0; i<2; i++) {
  switch (lifestage[i])	{
  case 0: /* Small juveniles */
  impact[i][2] = OMEGA*LENGTH(i)*LENGTH(i)*LENGTH(i);
  impact[i][3] = 0;
  impact[i][4] = 0;
  break;
  case 1: /* Non-vulnerable juveniles */
  impact[i][2] = 0;
  impact[i][3] = OMEGA*LENGTH(i)*LENGTH(i)*LENGTH(i);
  impact[i][4] = 0;
  break;
  case 2: /* Adults */
  impact[i][2] = 0;
  impact[i][3] = 0;
  impact[i][4] = OMEGA*LENGTH(i)*LENGTH(i)*LENGTH(i);
  break;
  }
  }
  return;
}


/*
* Specify the type of each of the environment variables by setting
* the entries in EnvironmentType[ENVIRON_DIM] to PERCAPITARATE, GENERALODE
* or POPULATIONINTEGRAL based on the classification below:
*
* Set an entry to PERCAPITARATE if the dynamics of E[j] follow an ODE and 0
* is a possible equilibrium state of E[j]. The ODE is then of the form
* dE[j]/dt = P(E,I)*E[j], with P(E,I) the per capita growth rate of E[j].
* Specify the equilibrium condition as condition[j] = P(E,I), do not include
* the multiplication with E[j] to allow for detecting and continuing the
* transcritical bifurcation between the trivial and non-trivial equilibrium.
*
* Set an entry to GENERALODE if the dynamics of E[j] follow an ODE and 0 is
* NOT an equilibrium state of E. The ODE then has a form dE[j]/dt = G(E,I).
* Specify the equilibrium condition as condition[j] = G(E,I).
*
* Set an entry to POPULATIONINTEGRAL if E[j] is a (weighted) integral of the
* population distribution, representing for example the total population
* biomass. E[j] then can be expressed as E[j] = I[p][i]. Specify the
* equilibrium condition in this case as condition[j] = I[p][i].
*
* Notice that the first index of the variable 'I[][]' refers to the
* number of the structured population, the second index refers to the
* number of the interaction variable. The interpretation of the latter
* is up to the user. Also notice that the variable 'condition[j]' should
* specify the equilibrium condition of environment variable 'E[j]'.
*/

/* DB correction: number of GENERALODE/PERCAPITARATE/POPULATIONINTEGRAL specs must equal the length of the condition[] vector */
const int EnvironmentType[ENVIRON_DIM] = {GENERALODE, GENERALODE, PERCAPITARATE, POPULATIONINTEGRAL, POPULATIONINTEGRAL}; 

void EnvEqui(double E[], double I[POPULATION_NR][INTERACT_DIM],
double condition[ENVIRON_DIM])
{
  
  condition[0] = RHO1*(RMAX1 - R1) - (I[0][0] + I[1][0]);
  condition[1] = RHO2*(RMAX2 - R2) - (I[0][1] + I[1][1]); 
  condition[2] = EPSILON*A*(I[0][2] + I[1][2])/(1+A*TH*(I[0][2] + I[1][2])) - DELTA;
  condition[3] = I[0][2];
  condition[4] = I[1][2];
  
  return;
}


/*==============================================================================*/

