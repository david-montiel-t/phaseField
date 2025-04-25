// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                                            [[maybe_unused]] const unsigned int index,
                                            [[maybe_unused]] double            &scalar_IC,
                                            [[maybe_unused]] Vector<double>    &vector_IC)
{
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index

	  // The initial condition is two circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "rad".
  
    //Vertical interface position
    double zint_1;
    double zint_2;
    double Lx = userInputs.domain_size[0];
    double Ly = userInputs.domain_size[1];
    double pi = 3.14159265358979323846;
  
    //Default return value
    scalar_IC = 0;
    // Initial condition for the concentration field
    if (index == 0){
        scalar_IC = U0;
    }
    // Initial condition for the order parameter field
    else if (index == 1) {
      //Short wave perturbation 
      zint_1 = zint0+0.25*A*(2.0+std::cos(kq*pi*(0.5*Lx-p[0])/(0.5*Lx))+std::cos(kq*pi*(0.5*Ly-p[1])/(0.5*Ly)));

      //Head start for corner dendrite
      zint_2 = zint_1+0.125*A*(std::cos(2.0*pi*p[0]/(0.5*Lx))+std::cos(2.0*pi*p[1]/(0.5*Lx)));
      
      scalar_IC = (-std::tanh((p[2]-zint_2)/(sqrt(2.0))));
    }
	  // --------------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


    // -------------------------------------------------------------------------

}
