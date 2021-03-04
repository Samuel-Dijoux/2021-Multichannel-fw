/*
  Tritrophicmod5_demography2.h
 Extension of Tritrophicmod5_demography.h with the exclusion of Environmental and Interaction variables
 
 Extension of Tritrophicmod5_fin.h - Header file specifying the elementary life-history functions of 
 the Multichannel food web tri-trophic model analyzed in:
 S. Dijoux & D.S. Boukal, 2021. Community structure and collapses in multi-channel food webs:
 role of consumer body sizes and mesohabitat productivities.
 
 
 Model extension of the initial Tritrophic.h 
 A.M. de Roos & L. Persson, 2002. Size-dependent life-history
 traits promote catastrophic collapses of top predators.
 Proc. Natl. Acad. Sciences 99(20): 12907-12912.		
*/
  
/*
*===========================================================================
*   SECTION 1: PROBLEM DIMENSIONS, NUMERICAL SETTINGS AND MODEL PARAMETERS
*===========================================================================
*/
// Dimension settings: Required
#define POPULATION_NR       1
#define STAGES              2
#define	I_STATE_DIM         2
#define	PARAMETER_NR       17

// Numerical sttings: Optionnal
#define MIN_SURVIVAL       1E-20  // Survival at which individual is considered dead
#define ALLOWNEGATIVE      0      // Permissive negative values for the solution

// Descriptive names of parameters in parameter array (at least two parameters are required)
char  *parameternames[PARAMETER_NR] =
{ "Lb", "Lj", "Lm", "Omega", "Imax", "Rh", "Nu", "Rm", "Mub", "R", "Rho", "Lv", "A", "Th", "Epsilon", "Delta", "Bsratio" };

// Default values of all parameters 
double	parameter[PARAMETER_NR] =
{ 7.0, 110.0, 300.0, 9.0E-6, 1.0E-4, 1.5E-5, 0.006, 0.003, 0.01, 3E-4, 0.2, 27.0, 5000.0, 0.1, 0.5, 0.01, 0.5 };

// Aliases definitions for all istate variables
#define AGE                 istate[0][0]
#define LENGTH              istate[0][1]

// Aliases definitions for all parameters
#define LB                  parameter[ 0]  // Default: 7 mm				        Length at birth.
#define LJ                  parameter[ 1]  // Default: 110 mm			        Length at maturation
#define LM                  parameter[ 2]  // Default: 300 mm			        Maximum Length.

#define OMEGA               parameter[ 3]  // Default: 9E-6 g/mm^3	  	  Proportiononality constant.

#define IMAX                parameter[ 4]  // Default: 1E-4 g/day/mm?		  Proportiononality constant.
#define RH                  parameter[ 5]  // Default: 1.5E-5 g/L		      Half-saturation constant of the resource R.

#define NU                  parameter[ 6]  // Default: 0.006 /day		      Growth rate.
#define RM                  parameter[ 7]  // Default: 0.003 /day/mm?	    Proportiononality constant. 

#define MUB                 parameter[ 8]  // Default: 0.01 /day		      Mortality rate.

#define R                   parameter[ 9]  // Default: 3E-4 g/L		        Fixed resource R.
#define RHO                 parameter[10]  // Default: 0.2 /day

#define LV                  parameter[11]  // Default: 27 mm			        Length treshold of predation vulnerability
#define A                   parameter[12]  // Default: 5000 L/day		  	  Attack rate.
#define TH                  parameter[13]  // Default: 0.1 day/g		      Handling time.
#define EPSILON             parameter[14]  // Default: 0.5				        Conversion efficiency.
#define DELTA               parameter[15]  // Default: 0.01 /day		      Predator Mortality rate.

#define BSRATIO             parameter[16]  // Default: 1.2				        Body size ratio btw B1 and B2.

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
  AGE = 0.0;
  LENGTH = LB * BSRATIO;
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
    limit[0] = LENGTH - (BSRATIO*LJ);
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
  development[0][1] = NU*(BSRATIO*LM*R/(R+RH) - LENGTH);
  
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
  
  if (lifestage[0] == 1) 
    fecundity[0][0] = RM*R/(R+RH)* LENGTH*LENGTH;
    else
    fecundity[0][0] = 0.0;
  
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
  mortality[0] = MUB;
  return;
  }


/*==============================================================================*/

