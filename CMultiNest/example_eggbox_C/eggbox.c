#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <sys/stat.h>

#include "multinest.h"

/*From the source folder we include the following files*/
#include "../../DirectDetectionCalc/XENON/binning_Xe.h"
#include "../../DirectDetectionCalc/source/coeffs_eft.h"
#include "../../DirectDetectionCalc/source/difRateGen.h"
#include "../../DirectDetectionCalc/XENON/Loglike_Xe.h"
#include "../../DirectDetectionCalc/GERMANIUM/Loglike_Ge.h"
#include "../../DirectDetectionCalc/Argon/Loglike_Ar.h"
#include "../../DirectDetectionCalc/source/halo.h"


void RemoveSpaces(char* source)
{
    char* i = source;
    char* j = source;
    while(*j != 0)
    {
        *i = *j++;
        if(*i != ' ')
            i++;
    }
    *i = 0;
}

struct things_to_pass { char * root; char * data_points_Xe ; char * data_points_Ge; char * data_points_Ar; double exposureGe; double Eth_Ge; double Emax_Ge; int Nbins_Ge; double exposureXe; double Eth_Xe; double Emax_Xe; int Nbins_Xe; double exposureAr; double Eth_Ar; double Emax_Ar; int Nbins_Ar;};

struct halo_params struct_halo = { "SHM", 544, 220., 0., 90., 245., 220., 0., 2, 15., 151., 365., "halo_table.dat"};

//================================================================
//			Parameter Ranges
//================================================================

#define mchi_MIN 6.0
#define mchi_MAX 1000.0
//#define mchi_MIN 10.0
//#define mchi_MAX 10.0
double mchi_Cube;

#define sigsi_MIN 1.0E-12
#define sigsi_MAX 1.0E-8
double sigsi_Cube;

#define sigan_MIN 1.0E-4
#define sigan_MAX 1.0E0
double sigan_Cube;

#define C1_MIN 1.0E-6
#define C1_MAX 1.0E-3
//#define C1_MIN 5.0E-5
//#define C1_MAX 5.0E-5
double C1_Cube;
double C1_pCube;
double C1_nCube;

#define C3_MIN 1.0E-3
#define C3_MAX 1.0E0
double C3_Cube;

#define C4_MIN 1.E-3
#define C4_MAX 1.E0
//#define C4_MIN 1.E-1
//#define C4_MAX 1.E-1
double C4_Cube;

#define C7_MIN 1.0E0
#define C7_MAX 1.0E3
double C7_Cube;

#define C8_MIN 1.0E-2
#define C8_MAX 1.0E0
//#define C8_MIN 1.0E-1
//#define C8_MAX 1.0E-1
double C8_Cube;
double C8_pCube;
double C8_nCube;

#define C9_MIN 1.0E-1
#define C9_MAX 1.0E2
//#define C9_MIN 1.5E0
//#define C9_MAX 1.5E0
double C9_Cube;


#define C10_MIN 1.E-2
#define C10_MAX 1.E1
//#define C10_MIN 2.0E0
//#define C10_MAX 2.0E0
double C10_Cube;

#define C11_MIN 1.E0
#define C11_MAX 1.E3
double C11_Cube;

#define C_M_SI_MIN 1.E-6
#define C_M_SI_MAX 1.E-3
//#define C_M_SI_MIN 2.E-4
//#define C_M_SI_MAX 2.E-4
double C_M_SI_Cube;

#define C_MV_SI_MIN 1.E-4
#define C_MV_SI_MAX 1.E-2
//#define C_MV_SI_MIN 2.E-4
//#define C_MV_SI_MAX 2.E-4
double C_MV_SI_Cube;

#define C_pi_SI_MIN 1.E-4
#define C_pi_SI_MAX 1.E-1
//#define C_pi_SI_MIN 6.E-3
//#define C_pi_SI_MAX 6.E-3
double C_pi_SI_Cube;

double C1_isoscalar_Cube;
double C1_isovector_Cube;

#define C11_MIN 1.E0
#define C11_MAX 1.E3
double C11_Cube;

#define rho_MIN 0.2
#define rho_MAX 0.6
double rho_Cube;

#define vesc_MIN 478
#define vesc_MAX 610
double vesc_Cube;

#define v0_MIN 170
#define v0_MAX 290
double v0_Cube;

#define k_MIN 0.5
#define k_MAX 3.5
double k_Cube;



//================================================================
//			Coupling on/off
//================================================================
int mass_analysis = 1;
int sigsi_analysis = 0;
int sigan_analysis = 0;
int C1_analysis = 1; // (1-Switched on) (0-Switched off)
int C1_Isospin_analysis = 0;
int C3_analysis = 0;
int C4_analysis = 0; // (1-Switched on) (0-Switched off)
int C7_analysis = 0;
int C8_analysis = 0; // (1-Switched on) (0-Switched off)
int C8_Isospin_analysis = 0;
int C9_analysis = 0;
int C10_analysis = 0;
int C10_Isospin_analysis = 0;
int C11_analysis = 0;
int C_M_SI_analysis  = 0;
int C_MV_SI_analysis  = 0;
int C_pi_SI_analysis  =0;
int C1_isoscalar_analysis = 0;
int C1_isovector_analysis = 0;

//================================================================
//           HALO UNCERTAINTIES (rho, vesc, v0, k)
//================================================================

 int HALO_UNCERTAINTIES = 1; // (1-Switched on) (0-Switched off)

//================================================================
//================================================================

double gaussian_loglike(double x, double median, double sigma){
return -pow(x-median,2.)/(2*sigma*sigma);
}





/******************************************** loglikelihood routine ****************************************************/

// Now an example, sample an egg box likelihood

// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
// context						void pointer, any additional information
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//
// Output arguments
// lnew 						= loglikelihood

void LogLike(double *Cube, int *ndim, int *npars, double *lnew, void *context)
{
	double LL = 0.;
	/*double chi = 1.0;
	int i;
	for(i = 0; i < *ndim; i++)
	{
		double x = Cube[i]*10.0*M_PI;
		chi *= cos(x/2.0);
		Cube[i] = x;
	}
	*lnew = powf(chi + 2.0, 5.0);*/

    //unpack context
    struct things_to_pass * pass_vars
    = (struct things_to_pass *)context;
    char * root_passed = (pass_vars->root);
    //printf("\n %s \n", root_passed);

    char * data_points_passed_Xe = (pass_vars->data_points_Xe);
    char * data_points_passed_Ge = (pass_vars->data_points_Ge);
    char * data_points_passed_Ar = (pass_vars->data_points_Ar);
    double exposure_Ge = (pass_vars-> exposureGe);
    double E_th_Ge = (pass_vars-> Eth_Ge);
    double E_max_Ge = (pass_vars-> Emax_Ge);
    int N_bins_Ge = (pass_vars->Nbins_Ge);
    double exposure_Xe = (pass_vars-> exposureXe);
    double E_th_Xe = (pass_vars-> Eth_Xe);
    double E_max_Xe = (pass_vars-> Emax_Xe);
    int N_bins_Xe = (pass_vars->Nbins_Xe);
    double exposure_Ar = (pass_vars-> exposureAr);
    double E_th_Ar = (pass_vars-> Eth_Ar);
    double E_max_Ar = (pass_vars-> Emax_Ar);
    int N_bins_Ar = (pass_vars->Nbins_Ar);
    double rho_MOCK = 0.4;
 	double vesc_MOCK = 544.;
 	double v0_MOCK = 220.;
 	double k_MOCK = 2.;
 	double ve_MOCK=220.;
 	char* halo_table_path[32];
	struct halo_params struct_halo = { "SHM", 544, 220., 0., 90., 245., 220., 0., 2, 15., 151., 365., "halo_table.dat"};

    set_coeffs_gen_SI ( );

	//-------------Clustered Params-----------------------------------------
	if( mass_analysis==1){
	Cube[0] = exp(log(mchi_MIN)+Cube[0]*log(mchi_MAX/mchi_MIN));
	mchi_Cube = Cube[0];
	set_any_coeffs(mchi_Cube, 0);
	printf("mchi %.5E | ", mchi_Cube);

	}
	if (mass_analysis==0){
	mchi_Cube = 10.0;
	}
	if( sigsi_analysis==1){
	Cube[mass_analysis] = exp(log(sigsi_MIN)+Cube[mass_analysis]*log(sigsi_MAX/sigsi_MIN));
	sigsi_Cube = Cube[mass_analysis];
	set_coeffs_matching_Standard( mchi_Cube, sigsi_Cube, 0.0);
    printf("sigsi %.5E |", sigsi_Cube);
	}

    if( sigan_analysis==1){
        Cube[mass_analysis] = exp(log(sigan_MIN)+Cube[mass_analysis]*log(sigan_MAX/sigan_MIN));
        sigan_Cube = Cube[mass_analysis];
        set_coeffs_anapole( mchi_Cube, sigan_Cube);
        printf("sigan %.5E |", sigan_Cube);
    }

	if(C1_analysis==1){
        Cube[mass_analysis] = exp(log(C1_MIN)+Cube[mass_analysis]*log(C1_MAX/C1_MIN));
        C1_Cube = Cube[mass_analysis];
        set_any_coeffs(Cube[mass_analysis], 1);
        printf("C1 %.5E |", Cp(1));
    }

	if(C1_analysis==0){
		//C1_Cube = Inputs_mock[1] = 1.E-99;
		C1_Cube = 0.0;
	}

	if(C1_Isospin_analysis==1){
        Cube[1] = exp(log(C1_MIN)+Cube[1]*log(C1_MAX/C1_MIN));
	Cube[2] = exp(log(C1_MIN)+Cube[2]*log(C1_MAX/C1_MIN));
	//printf("Cp_O1 %.5E |", Cube[1]);
	set_given_coeff(Cube[1], 1, "p");
	//printf("Cn_O1 %.5E |", Cube[2]);
	set_given_coeff(Cube[2], 1, "n");
        }


	if(C3_analysis==1){
        Cube[C1_analysis + mass_analysis] = exp(log(C3_MIN)+Cube[C1_analysis + mass_analysis]*log(C3_MAX/C3_MIN));
        C3_Cube = Cube[C1_analysis + 1];
	printf("C3 %.5E |", C3_Cube);
	set_given_coeff(C3_Cube, 3, "p");
	set_given_coeff(C3_Cube, 3, "n");
        }

	if(C3_analysis==0){
		//C1_Cube = Inputs_mock[1] = 1.E-99;
		C3_Cube = 0.0;
	}
	if(C4_analysis==1){
        Cube[C1_analysis + C3_analysis + 1] = exp(log(C4_MIN)+Cube[C1_analysis + C3_analysis + 1]*log(C4_MAX/C4_MIN));
        C4_Cube = Cube[C1_analysis + C3_analysis + 1];
	printf("C4 %.5E |", C4_Cube);
	set_given_coeff(C4_Cube, 4, "p");
	set_given_coeff(C4_Cube, 4, "n");
        }

	if(C4_analysis==0){
		//C4_Cube = Inputs_mock[4] = 1.E-99;
		C4_Cube = 1.E-99;
	}

	if(C7_analysis==1){
        Cube[C1_analysis + C3_analysis + C4_analysis + 1] = exp(log(C7_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + 1]*log(C7_MAX/C7_MIN));
        C7_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + 1];
	printf("C7 %.5E |", C7_Cube);
	set_given_coeff(C7_Cube, 7, "p");
	set_given_coeff(C7_Cube, 7, "n");
        }

	if(C7_analysis==0){
		//C8_Cube = Inputs_mock[8] = 1.E-99;
		C7_Cube = 1.E-99;
	}

	if(C8_analysis==1){
        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + mass_analysis] = exp(log(C8_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + mass_analysis]*log(C8_MAX/C8_MIN));
        C8_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + mass_analysis];
	printf("C8 %.5E |", C8_Cube);
	set_given_coeff(C8_Cube, 8, "p");
	set_given_coeff(C8_Cube, 8, "n");
        }

	if(C8_analysis==0){
		//C8_Cube = Inputs_mock[8] = 1.E-99;
		C8_Cube = 1.E-99;
	}

	if(C8_Isospin_analysis==1){
        Cube[3] = exp(log(C8_MIN)+Cube[3]*log(C8_MAX/C8_MIN));
	Cube[4] = exp(log(C8_MIN)+Cube[4]*log(C8_MAX/C8_MIN));
	//printf("Cp_O8 %.5E |", Cube[3]);
	set_given_coeff(Cube[3], 8, "p");
	//printf("Cn_O8 %.5E |", Cube[4]);
	set_given_coeff(Cube[4], 8, "n");
        }

	if(C9_analysis==1){
        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + 1] = exp(log(C9_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + 1]*log(C9_MAX/C9_MIN));
        C9_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + 1];
	printf("C9 %.5E |", C9_Cube);
	set_given_coeff(C9_Cube, 9, "p");
	set_given_coeff(C9_Cube, 9, "n");
        }

	if(C9_analysis==0){
		//C8_Cube = Inputs_mock[8] = 1.E-99;
		C9_Cube = 1.E-99;
	}

	if(C10_analysis==1){
        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + mass_analysis] = exp(log(C10_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + mass_analysis]*log(C10_MAX/C10_MIN));
        C10_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + mass_analysis];
	//printf("C10 %.5E |", C10_Cube);
        set_any_coeffs(Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + mass_analysis], 10);
        printf("C10 %.5E |", Cp(10));
    }

	if(C10_analysis==0){
		//C8_Cube = Inputs_mock[8] = 1.E-99;
		C10_Cube = 1.E-99;
	}

	if(C10_Isospin_analysis==1){
        Cube[3] = exp(log(C10_MIN)+Cube[3]*log(C10_MAX/C10_MIN));
	Cube[4] = exp(log(C10_MIN)+Cube[4]*log(C10_MAX/C10_MIN));
	//printf("Cp_O10 %.5E |", Cube[3]);
	set_given_coeff(Cube[3], 10, "p");
	//printf("Cn_O10 %.5E |", Cube[4]);
	set_given_coeff(Cube[4], 10, "n");
        }

	if(C11_analysis==1){
        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + 1] = exp(log(C11_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + 1]*log(C11_MAX/C11_MIN));
        C11_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + 1];
	printf("C11 %.5E |", C11_Cube);
	set_given_coeff(C11_Cube, 11, "p");
	set_given_coeff(C11_Cube, 11, "n");
        }

	if(C11_analysis==0){
		//C8_Cube = Inputs_mock[8] = 1.E-99;
		C11_Cube = 1.E-99;
	}

	if(C_M_SI_analysis==1){
        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + mass_analysis] = exp(log(C_M_SI_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis+mass_analysis]*log(C_M_SI_MAX/C_M_SI_MIN));

        C_M_SI_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis];
	set_coeffs_gen(C_M_SI_Cube, "S");
	printf("C_M_SI_Cube %.5E|", C_scalar_SI(1));
        }

	if(C_MV_SI_analysis==1){
        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + mass_analysis + C_M_SI_analysis] = exp(log(C_MV_SI_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis+mass_analysis+ C_M_SI_analysis]*log(C_MV_SI_MAX/C_MV_SI_MIN));

        C_MV_SI_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + + C_M_SI_analysis];
	set_coeffs_gen(C_MV_SI_Cube, "V");
	printf("C_MV_SI_Cube %.5E|", C_vector_SI(1));
        }

	if(C_pi_SI_analysis==1){
        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C_M_SI_analysis + C_MV_SI_analysis + mass_analysis] = exp(log(C_pi_SI_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis+ C_M_SI_analysis + C_MV_SI_analysis +mass_analysis]*log(C_pi_SI_MAX/C_pi_SI_MIN));

        C_pi_SI_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + C_M_SI_analysis + mass_analysis];
	set_coeffs_gen(C_pi_SI_Cube, "Pi");
	printf("C_pi_SI_Cube %.5E|", C_scalar_SI(3));
        }

	if(C1_isoscalar_analysis==1 && C1_isovector_analysis == 0 ){
        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis] = exp(log(C_M_SI_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis+mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis]*log(C_M_SI_MAX/C_M_SI_MIN));

        C1_isoscalar_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis];

	set_given_coeff(C1_isoscalar_Cube, 1, "p");
	set_given_coeff(C1_isoscalar_Cube, 1, "n");
	printf("C1_isoscalar %.5E|", (Cp(1)+Cn(1))/2);
        }

	if(C1_isoscalar_analysis==0 && C1_isovector_analysis == 1){
        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis] = exp(log(C_MV_SI_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis+mass_analysis+ C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis]*log(C_MV_SI_MAX/C_MV_SI_MIN));

        C1_isovector_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis];
	set_given_coeff(C1_isovector_Cube, 1, "p");
	set_given_coeff(-C1_isovector_Cube, 1, "n");
	printf("C_isovector %.5E|", (Cp(1)-Cn(1))/2);
        }

	if(C1_isoscalar_analysis==1 && C1_isovector_analysis == 1){

        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis] = exp(log(C_M_SI_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis+mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis]*log(C_M_SI_MAX/C_M_SI_MIN));

        Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis + C1_isoscalar_analysis] = exp(log(C_MV_SI_MIN)+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis+mass_analysis+ C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis + C1_isoscalar_analysis]*log(C_MV_SI_MAX/C_MV_SI_MIN));


	C1_isoscalar_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis];


        C1_isovector_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis + C1_isoscalar_analysis];


	set_given_coeff(C1_isoscalar_Cube + C1_isovector_Cube, 1, "p");
	set_given_coeff(C1_isoscalar_Cube - C1_isovector_Cube, 1, "n");
	printf("C_isoscalar %.5E, %.5E|", (Cp(1)+Cn(1))/2, C1_isoscalar_Cube);
	printf("C_isovector %.5E, %.5E|", (Cp(1)-Cn(1))/2, C1_isovector_Cube);
        }


	if(HALO_UNCERTAINTIES==1){

  		Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis] = rho_MIN+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis]*(rho_MAX-rho_MIN);

  		rho_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis];

  		Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis+ 1] = vesc_MIN+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis+ 1]*(vesc_MAX-vesc_MIN);

  		vesc_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis+ 1];

  		Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis+ 2] = v0_MIN+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis+ 2]*(v0_MAX-v0_MIN);

  		v0_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis+ 2];

  		Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis+ 3] = k_MIN+Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis+ 3]*(k_MAX-k_MIN);

  		k_Cube = Cube[C1_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + C10_analysis + C11_analysis + mass_analysis + C_MV_SI_analysis + C_M_SI_analysis + C_pi_SI_analysis+ 3];

		struct halo_params struct_halo = { "Lisanti", vesc_Cube, v0_Cube, 0., 90., 245., 220., k_Cube, 2, 15., 151., 365., "halo_table.dat"};
    		//{ char * profile; double vesc; double v0; double beta; double vt; double vc; double ve; double k; int i ; double ve0; double t0; double T; char* halo_path;};
//struct difcros_params struct_difcros_Ge = {"Ge", 1., 100., 0.5, 220., "All", 300.*365.};
		define_and_write_halo_path("halo_table_temp.dat", struct_halo.profile, struct_halo.vesc, struct_halo.v0, struct_halo.beta, struct_halo.vt, struct_halo.vc, struct_halo.ve, struct_halo.k, struct_halo.i);
		snprintf(halo_table_path, sizeof halo_table_path, "halo_table_temp.dat", halo_table_path);
		printf("rho_Cube %.5E | vesc_Cube %.5E | v0_Cube %.5E | k_Cube %.5E|", rho_Cube, vesc_Cube, v0_Cube, k_Cube);

  }

  if(HALO_UNCERTAINTIES==0){
  rho_Cube = rho_MOCK;
  vesc_Cube = vesc_MOCK;
  v0_Cube = v0_MOCK;
  k_Cube = k_MOCK;
  snprintf(halo_table_path, sizeof halo_table_path, "../DirectDetectionCalc/halo_table.dat", halo_table_path);
  //printf("\n halo path %s \n", halo_table_path);

  }
	struct difcros_params struct_difcros_Xe = {"Xe", 1., mchi_Cube, 0.5, 220., "All", 5.6e6};
	struct difcros_params struct_difcros_Ge = {"Ge", 1., mchi_Cube, 0.5, 220., "All", 5.6e6};
	struct difcros_params struct_difcros_Ar = {"Ar", 1., mchi_Cube, 0.5, 220., "All", 5.6e6};
	LL += Loglike_Xe_difrate(halo_table_path, rho_Cube, & struct_difcros_Xe, exposure_Xe, E_th_Xe, E_max_Xe, struct_halo.t0, N_bins_Xe, data_points_passed_Xe) + gaussian_loglike(rho_Cube, rho_MOCK, (rho_MOCK-rho_MIN)/2.) + gaussian_loglike(vesc_Cube, vesc_MOCK, (vesc_MOCK-vesc_MIN)/2.) + gaussian_loglike(v0_Cube, v0_MOCK, (v0_MOCK-v0_MIN)/2.);
  //LL += Loglike_Xe_difrate(halo_table_path, rho_Cube, & struct_difcros_Xe, exposure_Xe, E_th_Xe, E_max_Xe, struct_halo.t0, N_bins_Xe, data_points_passed_Xe);
	//LL += Loglike_Xe_difrate_generalised_SI(0.4, & struct_difcros_Xe, exposure_Xe, E_th_Xe, E_max_Xe, struct_halo.t0, N_bins_Xe, data_points_passed_Xe);
	//LL += Loglike_Ar_difrate_generalised_SI(0.4, & struct_difcros_Ar, exposure_Ar, E_th_Ar, E_max_Ar, struct_halo.t0, N_bins_Ar, data_points_passed_Ar);
	printf( "LogLi %.5e \n", LL);
	*lnew = LL;
}

/***********************************************************************************************************************/




/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value from the default (non-INS) mode
// INSlogZ						= log evidence value from the INS mode
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, void *context)
{
	// convert the 2D Fortran arrays to C arrays


	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns

	int i, j;

	double postdist[*nSamples][*nPar + 2];
	for( i = 0; i < *nPar + 2; i++ )
		for( j = 0; j < *nSamples; j++ )
			postdist[j][i] = posterior[0][i * (*nSamples) + j];



	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column

	double pLivePts[*nlive][*nPar + 1];
	for( i = 0; i < *nPar + 1; i++ )
		for( j = 0; j < *nlive; j++ )
			pLivePts[j][i] = physLive[0][i * (*nlive) + j];
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{

    char experimental_card_Ge[256];



    double exposure_Ge, E_th_Ge, E_max_Ge;
    int N_bins_Ge;

    snprintf(experimental_card_Ge, sizeof experimental_card_Ge, "%s/GERMANIUM.dat", argv[2]);
    FILE* exp_card_Ge;
    if (!(exp_card_Ge = fopen(experimental_card_Ge, "r"))){
        printf("\nopening experiment card failed... make one!\n");
    }
    fscanf(exp_card_Ge, "%lf %lf %lf %i", &exposure_Ge, &E_th_Ge, &E_max_Ge, &N_bins_Ge);

    char experimental_card_Xe[256];


    double exposure_Xe, E_th_Xe, E_max_Xe;
    int N_bins_Xe;

    snprintf(experimental_card_Xe, sizeof experimental_card_Xe, "%s/XENON.dat", argv[2]);
    FILE* exp_card_Xe;
    if (!(exp_card_Xe = fopen(experimental_card_Xe, "r"))){
        printf("\nopening experiment card failed... make one!\n");
    }
    fscanf(exp_card_Xe, "%lf %lf %lf %i", &exposure_Xe, &E_th_Xe, &E_max_Xe, &N_bins_Xe);

    char experimental_card_Ar[256];


    double exposure_Ar, E_th_Ar, E_max_Ar;
    int N_bins_Ar;

    snprintf(experimental_card_Ar, sizeof experimental_card_Ar, "%s/ARGON.dat", argv[2]);
    FILE* exp_card_Ar;
    if (!(exp_card_Ar = fopen(experimental_card_Ar, "r"))){
        printf("\nopening experiment card failed... make one!\n");
    }
    fscanf(exp_card_Ar, "%lf %lf %lf %i", &exposure_Ar, &E_th_Ar, &E_max_Ar, &N_bins_Ar);

    /*char experimental_card_F[256];


    double exposure_F, E_th_F, E_max_F;
    int N_bins_F;

    snprintf(experimental_card_F, sizeof experimental_card_F, "%s/FLOURINE.dat", argv[2]);
    FILE* exp_card_F;
    if (!(exp_card_F = fopen(experimental_card_F, "r"))){
        printf("\n F opening experiment card failed... make one!\n");
    }
    fscanf(exp_card_F, "%lf %lf %lf %i", &exposure_F, &E_th_F, &E_max_F, &N_bins_F); */


// set the MultiNest sampling parameters
	int D = mass_analysis + sigsi_analysis + sigan_analysis + C1_analysis + 2*C1_Isospin_analysis + C3_analysis + C4_analysis + C7_analysis + C8_analysis + 2*C8_Isospin_analysis + C9_analysis + C10_analysis + 2*C10_Isospin_analysis + C11_analysis + C_M_SI_analysis + C_pi_SI_analysis + C_MV_SI_analysis + C1_isoscalar_analysis + C1_isovector_analysis + 4*HALO_UNCERTAINTIES;

	int i;

	// set the MultiNest sampling parameters


	int IS = 1;					// do Nested Importance Sampling?

	int mmodal = 0;					// do mode separation?

	int ceff = 0;					// run in constant efficiency mode?

	int nlive = 1000;				// number of live points

	double efr = 0.8;				// set the required efficiency

	double tol = 0.1;				// tol, defines the stopping criteria

	int ndims = D;					// dimensionality (no. of free parameters)

	int nPar = D;					// total no. of parameters including free & derived parameters

	int nClsPar = D;				// no. of parameters to do mode separation on

	int updInt = 1000;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations

	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored

	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)

	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(i = 0; i < ndims; i++) pWrap[i] = 0;

    char root[100] ;		// root for output files

    snprintf(root, sizeof root, "%s", argv[1]);
    RemoveSpaces(root);

	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock

	int fb = 1;					// need feedback on standard output?

	int resume = 0;					// resume from a previous job?

	int outfile = 1;				// write output files?

	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization

	double logZero = -DBL_MAX;			// points with loglike < logZero will be ignored by MultiNest

	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied

    struct things_to_pass pass_vars = {root, argv[3], argv[4], argv[5], exposure_Ge, E_th_Ge, E_max_Ge, N_bins_Ge, exposure_Xe, E_th_Xe, E_max_Xe, N_bins_Xe, exposure_Ar, E_th_Ar, E_max_Ar, N_bins_Ar};
	void *context = &pass_vars;				// not requ typically a 10 am meeting is difficult to attend on timeired by MultiNest, any additional information user wants to pass



	// calling MultiNest
	printf("------ Starting time of the scan ------ \n");
    	time_t timer;
    	char buffer[25];
	struct tm* tm_info;
    	time(&timer);
    	tm_info = localtime(&timer);
	strftime(buffer, 25, "%Y:%m:%d %H:%M:%S", tm_info);
	puts(buffer);
	printf("-------------------------------------\n");

	double start_multinest = clock();

	run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
	logZero, maxiter, LogLike, dumper, context);


	double end_multinest = clock();
	double cpu_time_used_multinest = ((double)(end_multinest - start_multinest)) / (CLOCKS_PER_SEC*60.0*60.0);
	printf("Time Elapsed for MultiNest %f hours \n", cpu_time_used_multinest);
	printf("------ Ending time of the scan ------ \n");
	time_t timerend;
	char bufferend[25];
	struct tm* tm_info_end;
	time(&timerend);
    	tm_info_end = localtime(&timerend);
	strftime(bufferend, 25, "%Y:%m:%d %H:%M:%S", tm_info_end);
    	puts(bufferend);
}

/***********************************************************************************************************************/
