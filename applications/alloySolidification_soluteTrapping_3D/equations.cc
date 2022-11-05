// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"phi");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "phi,c,xi");
    set_dependencies_gradient_term_RHS(0, "");

    // Variable 1
	set_variable_name				(1,"c");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(1, "phi,c,grad(phi),grad(c)");
    set_dependencies_gradient_term_RHS(1, "phi,c,grad(phi),grad(c)");

	// Variable 2
	set_variable_name				(2,"xi");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,AUXILIARY);

    set_dependencies_value_term_RHS(2, "phi,c,grad(phi)");
    set_dependencies_gradient_term_RHS(2, "grad(phi)");

	 // std::cout << "D = " << D << std::endl; 
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

// The dimensionless solute supersaturation and its derivatives
scalarvalueType phi = variable_list.get_scalar_value(0);
scalargradType phix = variable_list.get_scalar_gradient(0);

// The order parameter and its derivatives
scalarvalueType c = variable_list.get_scalar_value(1);
scalargradType cx = variable_list.get_scalar_gradient(1);

// The auxiliary parameter and its derivatives
scalarvalueType xi = variable_list.get_scalar_value(2);

// --- Setting the expressions for the terms in the governing equations ---

// The azimuthal angle
//scalarvalueType theta;
//for (unsigned i=0; i< phi.n_array_elements;i++){
//	theta[i] = std::atan2(phix[1][i],phix[0][i]);
//}

// Calculation of (outward) interface normal vector
//scalarvalueType normgradn = std::sqrt(phix.norm_square());
scalarvalueType normgradn = phix[0]*phix[0] + phix[1]*phix[1];
if (dim==3)
{
	normgradn += phix[2]*phix[2];
}
normgradn = std::sqrt(normgradn); 

scalargradType normal = phix/(normgradn+constV(regval));

// let's rotate the interface normal to crystal frame
// 2D rotation matrix https://en.wikipedia.org/wiki/Rotation_matrix
scalarvalueType n_x = normal[0]*cos_theta0 - normal[1]*sin_theta0;
scalarvalueType n_y = normal[0]*sin_theta0 + normal[1]*cos_theta0;
scalarvalueType n_z = constV(0.0);
if (dim == 3)
{
	n_z = normal[2]; 
}

  
//The cosine of theta
//scalarvalueType cth = normal_crystal_x; // [0];
//The sine of theta
//scalarvalueType sth = normal_crystal_y; // [1];
scalarvalueType n_x_4 = n_x*n_x*n_x*n_x;
scalarvalueType n_y_4 = n_y*n_y*n_y*n_y;
scalarvalueType n_z_4 = n_z*n_z*n_z*n_z; 


// cth = normal[0];
// sth = normal[1];
//The cosine of 4 theta
//scalarvalueType c4th =sth*sth*sth*sth + cth*cth*cth*cth - constV(6.0)*sth*sth*cth*cth;


// Anisotropic term
scalarvalueType a_c = constV(1.0) - constV(3*eps4     )  + constV(4.0*eps4     )*(n_x_4 + n_y_4 + n_z_4);
scalarvalueType a_k = constV(1.0) + constV(3*epsilon_k)  - constV(4.0*epsilon_k)*(n_x_4 + n_y_4 + n_z_4);
//a_n = (constV(1.0)+constV(eps4)*c4th);
 
  
scalarvalueType c_eq = 0.5*( constV(1.0+k) - constV(1.0-k)*phi );
scalarvalueType inv_c_eq = 1./c_eq; // calculate inverse only once, divison is costly
  
//e^u
//scalarvalueType eu = constV(2.0/cl0)*c/(constV(1.0+k)-constV(1.0-k)*phi);
scalarvalueType eu = constV(1.0/cl0)*c*inv_c_eq;
  
//tau(theta)
//scalarvalueType tau_th=tau*a_n*a_n;
//scalarvalueType tau_th=tau*a_n*a_n*eu;
scalarvalueType tau_th= constV(lambda*W_physical/a1)*a_c * ( constV(beta0)*a_k + constV(a1*a2_minus*W_physical/D_liquid)*a_c*eu ) * tau;
// use e^u prefactor correction as proposed originally in Echebarria 2004 
// (and implemented e.g. in https://doi.org/10.1016/j.jcrysgro.2019.125418 )
  
//dphi/dt
scalarvalueType dphidt=xi/tau_th;

// q(phi) term
//scalarvalueType q_phi = (constV(1.0)-phi)/(constV(1.0+k)-constV(1.0-k)*phi);
scalarvalueType q_phi = (constV(1.0)-phi)*constV(0.5)*inv_c_eq; // /(constV(1.0+k)-constV(1.0-k)*phi);

// Antitrapping term
//scalargradType j_at = -constV(at*W*cl0*(1.0-k))*eu*dphidt*normal;
//scalargradType j_at = -constV(at)*(constV(1.0)-constV(A_trapping)*(constV(1.0)-phi)*(constV(1.0)-phi))*constV(W*cl0*(1.0-k))*eu*dphidt*normal;
scalargradType j_at   = -constV(at*W*cl0*(1.0-k))*(constV(1.0)-constV(A_trapping)*(constV(1.0)-phi)*(constV(1.0)-phi))*eu*dphidt*normal;
//scalargradType j_at = -constV(at*W*cl0*(1.0-k))                                                                     *eu*dphidt*normal;
// Note to self: Deal.II uses a data structure called VectoriedArray. Constants (some times) need to converted into this data structure via constV( ... )
    
// Diffusion term 1
scalargradType diff_term_1 = constV(D)*q_phi*cx;
 
// Diffusion term 2
//scalargradType diff_term_2 = constV(D)*q_phi*c*constV(1.0-k)*phix/(constV(1.0+k)-constV(1.0-k)*phi);
scalargradType diff_term_2 = constV(D)*q_phi*c*constV(1.0-k)*phix*0.5*inv_c_eq; // /(constV(1.0+k)-constV(1.0-k)*phi);

//Diffusion term
scalargradType diff_term = diff_term_1 + diff_term_2;

// Define required equations
scalarvalueType eq_phi = (phi+constV(userInputs.dtValue)*xi/tau_th);

scalarvalueType eq_c = c;

scalargradType eqx_c =  - constV(userInputs.dtValue)*(diff_term - j_at);

// --- Submitting the terms for the governing equations ---

// Terms for the equation to evolve the order parameter
variable_list.set_scalar_value_term_RHS(0,eq_phi);

// Terms for the equation to evolve the concentration
variable_list.set_scalar_value_term_RHS(1,eq_c);
variable_list.set_scalar_gradient_term_RHS(1,eqx_c);
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

 // --- Getting the values and derivatives of the model variables ---


// The order parameter and its derivatives
 scalarvalueType phi = variable_list.get_scalar_value(0);
 scalargradType phix = variable_list.get_scalar_gradient(0);
  
// The concentraton
scalarvalueType c = variable_list.get_scalar_value(1);

// --- Setting the expressions for the terms in the governing equations ---

  // Calculation of interface normal vector
  //scalarvalueType normgradn_sq = phix.norm_square();
  //scalarvalueType = inv_normgradn_sq = 1./(normgradn_sq + constV(regval));

  // scalarvalueType normgradn = std::sqrt(phix.norm_square());
  scalarvalueType normgradn = phix[0]*phix[0] + phix[1]*phix[1];
  if (dim==3)
  {
  	normgradn += phix[2]*phix[2];
  }
  normgradn = std::sqrt(normgradn); 

  scalargradType normal = phix/(normgradn+constV(regval));

  // 2D rotation matrix https://en.wikipedia.org/wiki/Rotation_matrix
  scalarvalueType n_x = normal[0]*cos_theta0 - normal[1]*sin_theta0;
  scalarvalueType n_y = normal[0]*sin_theta0 + normal[1]*cos_theta0;
  scalarvalueType n_z = constV(0.0);
  if (dim == 3)
  {
   	n_z = normal[2]; 
  }
  
  //The cosine of theta
  //scalarvalueType cth = normal_crystal_x; //[0];
  //The sine of theta
  //scalarvalueType sth = normal_crystal_y; // [1];

  scalarvalueType n_x_2 = n_x*n_x;
  scalarvalueType n_y_2 = n_y*n_y;
  scalarvalueType n_z_2 = n_z*n_z;
  scalarvalueType n_x_4 = n_x_2*n_x_2;
  scalarvalueType n_y_4 = n_y_2*n_y_2;
  scalarvalueType n_z_4 = n_z_2*n_z_2;

  // cth = normal[0];
  // sth = normal[1];
  // jotain menee pieleen tässä, tulee erroria. Debggaus ...

  //The cosine of 4 theta
  //scalarvalueType c4th = sth*sth*sth*sth + cth*cth*cth*cth - constV(6.0)*sth*sth*cth*cth;
  //The sine of 4 theta
  //scalarvalueType s4th = constV(4.0)*sth*cth*cth*cth - constV(4.0)*sth*sth*sth*cth;


  // Anisotropic term
 // scalarvalueType a_n;
 scalarvalueType a_n = constV(1.0)-constV(3*eps4) + constV(4.0)*(n_x_4 + n_y_4 + n_z_4);
 //a_n = (constV(1.0)+constV(epsilon)*std::cos(constV(4.0)*(theta)));
  //a_n = (constV(1.0)+constV(eps4)*c4th);
  
//gradient energy coefficient, its derivative and square
 scalarvalueType a_d;
 a_d = constV(16.0*eps4)*normgradn*a_n;
//  a_d = -constV(4.0)*constV(eps4)*s4th;

// The anisotropy term that enters in to the  equation for xi
 scalargradType aniso;
 aniso[0] = a_n*a_n*phix[0] + n_x*(n_x_2 - (n_x_4 + n_y_4)); 
 aniso[1] = a_n*a_n*phix[1] + n_y*(n_y_2 - (n_x_4 + n_y_4)); 
 if (dim==3)
 {
   aniso[2] = a_n*a_n*phix[2] + n_z*(n_z_2 - (n_x_4 + n_y_4)); 
 }
 //aniso[0] = constV(W*W)*(a_n*a_n*phix[0]-a_n*a_d*phix[1]);
 //aniso[1] = constV(W*W)*(a_n*a_n*phix[1]+a_n*a_d*phix[0]);

//Calculation of value term
  scalarvalueType fprime = -phi + phi*phi*phi;

  scalarvalueType c_eq = 0.5*( constV(1.0+k) - constV(1.0-k)*phi );
  scalarvalueType inv_c_eq = 1./c_eq; // calculate inverse only once, divison is costly
  
  //e^u
  // scalarvalueType eu = constV(2.0/cl0)*c/(constV(1.0+k)-constV(1.0-k)*phi);
  scalarvalueType eu = constV(1.0/cl0)*c*inv_c_eq; // /(constV(1.0+k)-constV(1.0-k)*phi);

 // dimensionless temperature changes
  scalarvalueType x =q_point_loc[0]; // The x-component
  scalarvalueType t_n =constV(this->currentTime); // The time
  scalarvalueType tep =((x-constV(x0)-Vtilde*t_n)/ltilde);
  
// scalarvalueType dforceterm = -constV(lambda/(1.0-k))*(constV(1.0)-phi*phi)*(constV(1.0)-phi*phi)*(eu-constV(1.0));
  scalarvalueType dforceterm  = -constV(lambda/(1.0-k))*(constV(1.0)-phi*phi)*(constV(1.0)-phi*phi)*(eu-constV(1.0) + constV(1.0-k)*tep);

 // Define the terms in the equations
 scalarvalueType eq_xi = (-fprime + dforceterm);

 scalargradType eqx_xi = (-aniso);

  // --- Submitting the terms for the governing equations ---

 variable_list.set_scalar_value_term_RHS(2,eq_xi);
 variable_list.set_scalar_gradient_term_RHS(2,eqx_xi);

}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}
