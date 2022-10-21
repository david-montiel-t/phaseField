#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
    // Constructor
    customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {
	 
  		 std::cout << "time scale tau0 = " << tau0 << " s";
  		 std::cout << ", Dimensionless solute diffusion coefficient: D = " << D << " W0^2/tau0" << std::endl;
//		 std::cout << "dt should be less than dx^2/(4*D) = " << 
  		 std::cout << ", cos_theta0 = " << cos_theta0;
		 std::cout << ", sin_theta0 = " << sin_theta0;
		 std::cout << ", a2_minus = "   << a2_minus;
		 std::cout << ", thermal_gradient_dimless = "   << thermal_gradient_dimless;
		 std::cout << ", Vtilde = "   << Vtilde;
		 std::cout << ", ltilde = "   << ltilde;
		 std::cout << std::endl; 
	 };

    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

private:
	#include "../../include/typeDefs.h"

	const userInputParameters<dim> userInputs;

	// Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the RHS of the governing equations for all other equations (in equations.h)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set the LHS of the governing equations (in equations.h)
	void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set postprocessing expressions (in postprocess.h)
	#ifdef POSTPROCESS_FILE_EXISTS
	void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
					const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
	#endif

	// Function to set the nucleation probability (in nucleation.h)
	#ifdef NUCLEATION_FILE_EXISTS
	double getNucleationProbability(variableValueContainer variable_value, double dV) const;
	#endif

	// ================================================================
	// Methods specific to this subclass
	// ================================================================


	// ================================================================
	// Model constants specific to this subclass
	// ================================================================


  double eps4 = userInputs.get_model_constant_double("eps4");
  double Omega = userInputs.get_model_constant_double("Omega");
  double k = userInputs.get_model_constant_double("k");
  double d0inW = userInputs.get_model_constant_double("d0inW");
  double W = userInputs.get_model_constant_double("W");
  double tau = userInputs.get_model_constant_double("tau");
  double cl0 = userInputs.get_model_constant_double("cl0");
  double regval = userInputs.get_model_constant_double("regval");

  double d0 = userInputs.get_model_constant_double("d0"); // capillary length [m]
  double D_liquid = userInputs.get_model_constant_double("D_liquid");
  double V_D_PF = userInputs.get_model_constant_double("V_D_PF");
  double beta0 = userInputs.get_model_constant_double("beta0");

  double m_liquidus                = userInputs.get_model_constant_double("m_liquidus");
  double alloy_nominal_composition = userInputs.get_model_constant_double("alloy_nominal_composition");
  double thermal_gradient          = userInputs.get_model_constant_double("thermal_gradient");
  double pulling_speed             = userInputs.get_model_constant_double("pulling_speed");
  double x0                        = userInputs.get_model_constant_double("x0");

  double theta0                    = userInputs.get_model_constant_double("theta0"); // crystal frame orientation. Parallel to coordinates if set to zero

  double sine_amplitude            = userInputs.get_model_constant_double("sine_amplitude"); // amplitude of the initial planar perturbation 
  double sine_wavelength           = userInputs.get_model_constant_double("sine_wavelength"); // wavelength of the initial planar perturbation 
//  double sine_wave_length = 20.0; // wave length of the initial perturbation
//  double sine_amplitude = 2.0; // amplitude of the initial sinusoidal perturbation 

// Fixed and derived constants
  double a1 = 0.8839;
//  double a2 = 0.6267;
  double at = 1.0/(2.0*std::sqrt(2.0)); // In the anti-trapping current, at is now multiplied with (1-A_trapping*(1-phi)*(1-phi))
  double lambda = a1/d0inW;
//  double D = a2*lambda*W*W/tau; // a more general form of dimensionless diffusion coefficient is used, which allows for kinetic coefficient beta0 > 0

  double W_physical = d0/d0inW; 
  double A_trapping = D_liquid/(W_physical*V_D_PF);
  //double A_trapping = 0;

// Fixed and derived constants
  double K_bar       = 0.0638-0.0505*A_trapping; 
  double F_bar_minus = std::sqrt(2.0)*std::log(2)/2.0 + 3.*std::sqrt(2.0)/4.*A_trapping; 
  double sigma_phi   = 2.0*std::sqrt(2.0)/3.0;
  double J           = 16./15.;
  double a2_minus    = J/sigma_phi*(K_bar + F_bar_minus); // a2 is replaced with a2_minus
  //double a2_minus = 0.6267; 

  // time scale "magnitude" tau0, which is the isotropic version of the time scale tau(n)
  double tau0 = W_physical*lambda/a1 * (beta0 + a1*a2_minus*W_physical/D_liquid);// time scale magnitude tau0, used to compute the dimensionless diffusion coefficient D

  double D = D_liquid*tau0/(W_physical*W_physical); // Dimensionless diffusion coefficient 
  // std::cout << "Dimensionless solute diffusion coefficient: D = " << D << " W0^2/tau0" << std::endl;


  double thermal_gradient_dimless = thermal_gradient*W_physical; // conversion from K/m to K/W0
  double Vtilde = pulling_speed * tau0/W_physical; // conversion m/s to W0/tau0
  double ltilde = m_liquidus*(1-k)*alloy_nominal_composition/thermal_gradient_dimless; // dimensionless thermal length
  //double D = 2.0;
  //
  // let's pre-calculate these trigonometric values, used to rotate the normal to crystal frame
  double pi =3.14159265;  
  // double theta0 = pi/4.0; 
  double cos_theta0 = std::cos(theta0*pi/180.0); 
  double sin_theta0 = std::sin(theta0*pi/180.0); 

	// ================================================================
// python script to check the values 
//
/*
import numpy as np 
D_liquid = 4.4e-9
d0 = 12.17e-9
V_D_PF = 0.2
d0inW = 0.277
beta0=0.0
a1 = 0.8839;

at = 1.0/(2.0*np.sqrt(2.0)) #; // In the anti-trapping current, at is now multiplied with (1-A_trapping*(1-phi)*(1-phi))
Lambda = a1/d0inW;

W_physical = d0/d0inW;
A_trapping = D_liquid/(W_physical*V_D_PF);

K_bar       = 0.0638-0.0505*A_trapping;
F_bar_minus = np.sqrt(2.0)*np.log(2)/2.0 + 3.*np.sqrt(2.0)/4.*A_trapping;
sigma_phi   = 2.0*np.sqrt(2.0)/3.0;
J           = 16./15.;
a2_minus    = J/sigma_phi*(K_bar + F_bar_minus); # // a2 is replaced with a2_minus

tau0 = W_physical*Lambda/a1 * (beta0 + a1*a2_minus*W_physical/D_liquid) # ;// time scale magnitude tau0, used to compute the dimensionless diffusion coefficient 
D = D_liquid*tau0/(W_physical*W_physical) # ; // Dimensionless diffusion coefficient 
D
 */


};
