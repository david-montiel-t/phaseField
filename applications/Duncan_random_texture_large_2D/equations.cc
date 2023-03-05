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

  //Number of order parameters
  unsigned int number_ops_ini = 91;
  
  //Number of empty ops (for nucleation)
  unsigned int number_ops_nuc = 50;
  
  //Number of ops (total)
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;
  
  //Field for dislocation_density
  //Variable zero
  set_variable_name                (0,"rho");
  set_variable_type                (0,SCALAR);
  set_variable_equation_type       (0,EXPLICIT_TIME_DEPENDENT);
  
  set_dependencies_value_term_RHS(0,"rho");
  
  //Fields for grains
  // For the input file 'parameters.in'
  for (unsigned int var_index=1; var_index<=number_ops_ini; var_index++){
    
    // For the input file 'parameters_large_2D.in'
    //for (unsigned int var_index=0; var_index<12; var_index++){
    std::string var_name = "n";
    var_name.append(std::to_string(var_index));
    std::string grad_var_name = "grad("+var_name+")";
    
    set_variable_name        (var_index,var_name);
    set_variable_type        (var_index,SCALAR);
    set_variable_equation_type    (var_index,EXPLICIT_TIME_DEPENDENT);
    
    // For the input file 'parameters.in'
    set_dependencies_value_term_RHS(var_index, var_name);
    set_dependencies_gradient_term_RHS(var_index, grad_var_name);
    
  }
  
  for (unsigned int var_index=number_ops_ini+1; var_index<=number_ops_total; var_index++){
    
    // For the input file 'parameters_large_2D.in'
    //for (unsigned int var_index=0; var_index<12; var_index++){
    std::string var_name = "n";
    var_name.append(std::to_string(var_index));
    std::string grad_var_name = "grad("+var_name+")";
    
    set_variable_name        (var_index,var_name);
    set_variable_type        (var_index,SCALAR);
    set_variable_equation_type    (var_index,EXPLICIT_TIME_DEPENDENT);
    
    // For the input file 'parameters.in'
    set_dependencies_value_term_RHS(var_index, var_name);
    set_dependencies_gradient_term_RHS(var_index, grad_var_name);
    
    set_allowed_to_nucleate      (var_index, true);
    set_need_value_nucleation    (var_index, false);
    
  }

  // Sum of multiplication of 2 order parameters (double product). Variable index: N+1
  set_variable_name                (number_ops_total+1,"mult2Op");
  set_variable_type                (number_ops_total+1,SCALAR);
  set_variable_equation_type        (number_ops_total+1,AUXILIARY);
  
  set_dependencies_value_term_RHS(number_ops_total+1,"n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, n41, n42, n43, n44, n45, n46, n47, n48, n49, n50, n51, n52, n53, n54, n55, n56, n57, n58, n59, n60, n61, n62, n63, n64, n65, n66, n67, n68, n69, n70, n71, n72, n73, n74, n75, n76, n77, n78, n79, n80, n81, n82, n83, n84, n85, n86, n87, n88, n89, n90, n91, n92, n93, n94, n95, n96, n97, n98, n99, n100, n101, n102, n103, n104, n105, n106, n107, n108, n109, n110, n111, n112, n113, n114, n115, n116, n117, n118, n119, n120, n121, n122, n123, n124, n125, n126, n127, n128, n129, n130, n131, n132, n133, n134, n135, n136, n137, n138, n139, n140, n141, mult2Op, mult3Op");
  
  set_allowed_to_nucleate      (number_ops_total+1, false);
  set_need_value_nucleation    (number_ops_total+1, true);
  
  //  Sum of multiplication of 3 order parameters (triple product). Variable index: N+2
  set_variable_name                (number_ops_total+2,"mult3Op");
  set_variable_type                (number_ops_total+2,SCALAR);
  set_variable_equation_type        (number_ops_total+2,AUXILIARY);
  
  set_dependencies_value_term_RHS(number_ops_total+2,"n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, n41, n42, n43, n44, n45, n46, n47, n48, n49, n50, n51, n52, n53, n54, n55, n56, n57, n58, n59, n60, n61, n62, n63, n64, n65, n66, n67, n68, n69, n70, n71, n72, n73, n74, n75, n76, n77, n78, n79, n80, n81, n82, n83, n84, n85, n86, n87, n88, n89, n90, n91, n92, n93, n94, n95, n96, n97, n98, n99, n100, n101, n102, n103, n104, n105, n106, n107, n108, n109, n110, n111, n112, n113, n114, n115, n116, n117, n118, n119, n120, n121, n122, n123, n124, n125, n126, n127, n128, n129, n130, n131, n132, n133, n134, n135, n136, n137, n138, n139, n140, n141, mult2Op, mult3Op");
  set_allowed_to_nucleate      (number_ops_total+2, false);
  set_need_value_nucleation    (number_ops_total+2, true);
  
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

  //Number of order parameters
  unsigned int number_ops_ini = 91;
  
  //Number of empty ops (for nucleation)
  unsigned int number_ops_nuc = 50;
  
  //Number of ops (total)
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;
  
  //N= Number of order parameters
  unsigned int N = number_ops_total;

  dealii::VectorizedArray<double> fnV = constV(0.0);
  dealii::VectorizedArray<double> fsV = constV(0.0);
  dealii::VectorizedArray<double> sum_nsq = constV(0.0);
  scalarvalueType rhocalc = constV(0.0);
  scalarvalueType ni, nj;
  scalargradType nix;
  scalarvalueType rho = variable_list.get_scalar_value(0);

  std::vector<scalarvalueType> value_terms;
  value_terms.resize(N);
  std::vector<scalargradType> gradient_terms;
  gradient_terms.resize(N);
  
  // -------------------------------------------------
  // Nucleation expressions
  // -------------------------------------------------
  //scalarvalueType source_term = constV(0.0);
  //scalarvalueType gamma = constV(1.0);
  std::vector<scalarvalueType> source_term;
  source_term.resize(N);
  std::vector<scalarvalueType> gamma;
  gamma.resize(N);
  
  //Initializing source and gamma terms
  for (unsigned int i=0; i<N; i++){
    source_term[i] = constV(0.0);
    gamma[i] = constV(1.0);
  }
  seedNucleus(q_point_loc, source_term, gamma);
  
  //Calculating rho
  //Sum of squares of order parameters
  for (unsigned int i=0; i<N; i++){
    ni = variable_list.get_scalar_value(i+1);
    sum_nsq += ni*ni;
  }

  for (unsigned int i=0; i<N; i++){
    ni = variable_list.get_scalar_value(i+1);
    rhocalc += ni*ni*rho_i[i]/sum_nsq;
  }
  

  for (unsigned int i=0; i<N; i++){
    ni = variable_list.get_scalar_value(i+1);
    nix = variable_list.get_scalar_gradient(i+1);
    fnV = - ni + ni*ni*ni;
    for (unsigned int j=0; j<N; j++){
      if (i != j){
        nj = variable_list.get_scalar_value(j+1);
        fnV += constV(2.0*alpha) * ni * nj*nj;
      }
    }
    fnV = m0*fnV;
    fsV = constV(2.0*m0*fsbs)*ni*(rho_i[i]-rhocalc)/sum_nsq;
    
    //Uncomment this to remove the dependency of stored energy
    //fsV = constV(0.0);
    if (this->currentTime < userInputs.nucleation_start_time - 0.5*userInputs.dtValue)
    {
      if ((i+1 >= number_ops_ini+1) && (i+1 <= number_ops_total))
      {
        //std::cout << "evolution of n6(t < t_start)" << std::endl;
        value_terms[i] = ni - gamma[i]*constV(userInputs.dtValue*MnV)*fnV;
        gradient_terms[i] = gamma[i]*constV(-userInputs.dtValue*KnV*MnV)*nix;
      }
      else
      {
        value_terms[i] = ni - constV(userInputs.dtValue*MnV)*fnV;
        gradient_terms[i] = constV(-userInputs.dtValue*KnV*MnV)*nix;
      }
    } else if (this->currentTime >= userInputs.nucleation_start_time - 0.5*userInputs.dtValue && this->currentTime <= userInputs.nucleation_end_time)
    {
      if ((i+1 >= number_ops_ini+1) && (i+1 <= number_ops_total))
      {
        //std::cout << "evolution of n6 (t_start < t < t_end)" << std::endl;
        value_terms[i] = ni + source_term[i];
        gradient_terms[i] = nix*constV(0.0);
      }
      else
      {
        value_terms[i] = ni;
        gradient_terms[i] = nix*constV(0.0);
      }
    } else if (this->currentTime <= userInputs.nucleation_end_time + globalSeedingTime)
    {
      if ((i+1 >= number_ops_ini+1) && (i+1 <= number_ops_total))
      {
        //std::cout << "evolution of n6 (t_end < t)" << std::endl;
        value_terms[i] = ni - gamma[i]*constV(userInputs.dtValue*MnV)*(fnV);
        gradient_terms[i] = gamma[i]*constV(-userInputs.dtValue*KnV*MnV)*nix;
      }
      else
      {
        value_terms[i] = ni - constV(userInputs.dtValue*MnV)*(fnV);
        gradient_terms[i] = constV(-userInputs.dtValue*KnV*MnV)*nix;
      }
    } else
    {
      if ((i+1 >= number_ops_ini+1) && (i+1 <= number_ops_total))
      {
        //std::cout << "evolution of n6 (t_end < t)" << std::endl;
        value_terms[i] = ni - gamma[i]*constV(userInputs.dtValue*MnV)*(fnV + fsV);
        gradient_terms[i] = gamma[i]*constV(-userInputs.dtValue*KnV*MnV)*nix;
      }
      else
      {
        value_terms[i] = ni - constV(userInputs.dtValue*MnV)*(fnV + fsV);
        gradient_terms[i] = constV(-userInputs.dtValue*KnV*MnV)*nix;
      }
    }
  }
  for (unsigned int i=0; i<N; i++){
    variable_list.set_scalar_value_term_RHS(i+1,value_terms[i]);
    variable_list.set_scalar_gradient_term_RHS(i+1,gradient_terms[i]);
  }
  // --- Submitting the terms for the governing equations ---
  variable_list.set_scalar_value_term_RHS(0,rhocalc);
}
  
  // =================================================================================
  // seedNucleus: a function particular to this app
  // =================================================================================
  template <int dim,int degree>
  void customPDE<dim,degree>::seedNucleus(const dealii::Point<dim,
                                          dealii::VectorizedArray<double> > & q_point_loc,
                                          std::vector<dealii::VectorizedArray<double> > & source_term,
                                          std::vector<dealii::VectorizedArray<double> > & gamma) const
  {
    
    for (typename std::vector<nucleus<dim> >::const_iterator thisNucleus=this->nuclei.begin(); thisNucleus!=this->nuclei.end(); ++thisNucleus){
      
      //if (thisNucleus->seededTime + thisNucleus->seedingTime > this->currentTime){
      //Site saturation condition
      if (this->currentTime < userInputs.nucleation_end_time + thisNucleus->seedingTime){
        // Calculate the weighted distance function to the order parameter freeze boundary (weighted_dist = 1.0 on that boundary)
        dealii::VectorizedArray<double> weighted_dist = this->weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_freeze_semiaxes(thisNucleus->orderParameterIndex), q_point_loc, thisNucleus->orderParameterIndex);
        
        for (unsigned i=0; i<gamma[0].n_array_elements;i++)
        {
          if (weighted_dist[i] <= 1.0)
          {
            gamma[thisNucleus->orderParameterIndex - 1][i] = 0.0;
            
            // Seed a nucleus if it was added to the list of nuclei this time step
            if (thisNucleus->seedingTimestep == this->currentIncrement)
            {
              // Find the weighted distance to the outer edge of the nucleus and use it to calculate the order parameter source term
              dealii::Point<dim,double> q_point_loc_element;
              for (unsigned int j=0; j<dim; j++)
              {
                q_point_loc_element(j) = q_point_loc(j)[i];
              }
              double r = this->weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_semiaxes(thisNucleus->orderParameterIndex), q_point_loc_element, thisNucleus->orderParameterIndex);
              
              double avg_semiaxis = 0.0;
              
              for (unsigned int j=0; j<dim; j++)
              {
                avg_semiaxis += thisNucleus->semiaxes[j];
              }
              avg_semiaxis /= dim;
              
              source_term[thisNucleus->orderParameterIndex - 1][i] = 0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));
            }
          }
        }
      }
    }
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
  
  //Number of order parameters
  unsigned int number_ops_ini = 91;
  
  //Number of empty ops (for nucleation)
  unsigned int number_ops_nuc = 50;
  
  //Number of ops (total)
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;
  
  //N= Number of order parameters
  unsigned int N = number_ops_total;
  
  scalarvalueType sum_nsq = constV(0.0);
  scalarvalueType ni = constV(0.0);
  scalarvalueType nj = constV(0.0);
  scalarvalueType nk = constV(0.0);
  
  //Calculation of double product and triple product
  scalarvalueType mult2Op_calc = constV(0.0);
  scalarvalueType mult3Op_calc = constV(0.0);
  
  for (unsigned int i = 0; i < number_ops_ini; i++)
  {
    ni = variable_list.get_scalar_value(i+1);
    for (unsigned int j = 0; j < number_ops_ini; j++)
    {
      nj = variable_list.get_scalar_value(j+1);
      if (j > i)
      {
        mult2Op_calc += ni*nj;
      }
      for (unsigned int k = 0; k < number_ops_ini; k++)
      {
        nk = variable_list.get_scalar_value(k+1);
        if(j > i && k > j)
        {
          mult3Op_calc += ni*nj*nk;
        }
      }
    }
  }

  variable_list.set_scalar_value_term_RHS(number_ops_total+1,mult2Op_calc);
  variable_list.set_scalar_value_term_RHS(number_ops_total+2,mult3Op_calc);

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

