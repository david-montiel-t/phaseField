// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for
// each function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void
customAttributeLoader::loadVariableAttributes()
{
  // Variable 0 - Order Parameter
  set_variable_name(0, "n");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "n, dndt");
  set_dependencies_gradient_term_RHS(0, "");

  // Variable 1 - Time Derivative of Order Parameter
  set_variable_name                (1, "dndt");
  set_variable_type                (1, SCALAR);
  set_variable_equation_type        (1, AUXILIARY);

  set_dependencies_value_term_RHS(1, "n, grad(n)");
  set_dependencies_gradient_term_RHS(1, "grad(n)");
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time
// dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a
// list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one
// proportional to the test function and one proportional to the gradient of the
// test function. The index for each variable in this list corresponds to the
// index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---
  // Get the value of the order parameter
  scalarvalueType n = variable_list.get_scalar_value(0);

  // Get the time derivative of the order parameter (calculated as AUXILIARY field)
  scalarvalueType dndt = variable_list.get_scalar_value(1);

  // Prevent the order parameter from decreasing (no detwinning)
  //dndt = std::max(dndt, constV(0.0));

  // Calculate the new order parameter value
  scalarvalueType new_n = n + constV(userInputs.dtValue) * dndt;

  // Restrict the new order parameter to be between zero and one
  new_n = std::min(new_n, constV(1.0));
  new_n = std::max(new_n, constV(0.0));

  // --- Submitting the terms for the governing equations ---
  variable_list.set_scalar_value_term_RHS(0,new_n);
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time
// independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are
// not explicit time-dependent equations. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is
// given by "q_point_loc". The function outputs two terms to variable_list --
// one proportional to the test function and one proportional to the gradient of
// the test function. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
    // The order parameter and its derivatives
    scalarvalueType n = variable_list.get_scalar_value(0);
    scalargradType nx = variable_list.get_scalar_gradient(0);

    scalarvalueType n_b; // Bounding n
    n_b = std::min(n, constV(1.0));
    n_b = std::max(n_b, constV(0.0));

    //Scalar value term (double obstacle term)
    scalarvalueType dobs_term = -4.0 * MnV * U * (1.0-2.0 * n_b);

      // Calculating gradient term
    scalargradType grad_term = -MnV * KnV * nx;

    // Calculating term p(|grad(n)|)
    scalarvalueType grad_n_c = constV(1.0/lt);
    scalarvalueType p =  tanh((nx.norm()-grad_n_c)/grad_n_c);
    scalargradType fbrterm = -MnV * mu * p * nx;


    variable_list.set_scalar_value_term_RHS(1, dobs_term);
    variable_list.set_scalar_gradient_term_RHS(1, grad_term + fbrterm);
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. The
// (x,y,z) location of that quadrature point is given by "q_point_loc". The
// function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function -- for the
// left-hand-side of the equation. The index for each variable in this list
// corresponds to the index given at the top of this file. If there are multiple
// elliptic equations, conditional statements should be sed to ensure that the
// correct residual is being submitted. The index of the field being solved can
// be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{}
