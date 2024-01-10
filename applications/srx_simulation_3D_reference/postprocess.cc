// =================================================================================
// Set the attributes of the postprocessing variables
// =================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but for
// the postprocessing expressions. It sets the attributes for each postprocessing
// expression, including its name, whether it is a vector or scalar (only scalars are
// supported at present), its dependencies on other variables and their derivatives,
// and whether to calculate an integral of the postprocessed quantity over the entire
// domain.

void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"feature_ids");
	set_variable_type				(0,SCALAR);

    // For the input file 'parameters.in'
    set_dependencies_value_term_RHS(0,"n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, n41, n42, n43, n44, n45, n46, n47, n48, n49, n50, n51, n52, n53, n54, n55, n56, n57, n58, n59, n60, n61, n62, n63, n64, n65, n66, n67, n68, n69, n70, n71, n72, n73, n74, n75, n76, n77, n78, n79, n80, n81, n82, n83, n84, n85, n86, n87, n88, n89, n90, n91, n92, n93, n94, n95, n96, n97, n98, n99, n100, n101, n102, n103, n104, n105, n106, n107, n108, n109, n110, n111, n112, n113, n114, n115, n116, n117, n118, n119, n120, n121, n122, n123, n124, n125, n126, n127, n128, n129, n130, n131, n132, n133, n134, n135, n136, n137, n138, n139, n140, n141, n142, n143, n144, n145, n146, n147, n148, n149, n150, n151, n152, n153, n154, n155, n156, n157, n158, n159, n160");

    set_dependencies_gradient_term_RHS(0, "");

    set_output_integral         	(0,false);

    // Variable 1
    set_variable_name				(1,"op_ids");
	set_variable_type				(1,SCALAR);

    set_dependencies_value_term_RHS(1, "n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, n41, n42, n43, n44, n45, n46, n47, n48, n49, n50, n51, n52, n53, n54, n55, n56, n57, n58, n59, n60, n61, n62, n63, n64, n65, n66, n67, n68, n69, n70, n71, n72, n73, n74, n75, n76, n77, n78, n79, n80, n81, n82, n83, n84, n85, n86, n87, n88, n89, n90, n91, n92, n93, n94, n95, n96, n97, n98, n99, n100, n101, n102, n103, n104, n105, n106, n107, n108, n109, n110, n111, n112, n113, n114, n115, n116, n117, n118, n119, n120, n121, n122, n123, n124, n125, n126, n127, n128, n129, n130, n131, n132, n133, n134, n135, n136, n137, n138, n139, n140, n141, n142, n143, n144, n145, n146, n147, n148, n149, n150, n151, n152, n153, n154, n155, n156, n157, n158, n159, n160");

    set_dependencies_gradient_term_RHS(1, "");

    set_output_integral         	(1,false);
  
  // Variable 2
  set_variable_name        (2,"sum_sq_op_ids");
  set_variable_type        (2,SCALAR);
  
  set_dependencies_value_term_RHS(2, "n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, n41, n42, n43, n44, n45, n46, n47, n48, n49, n50, n51, n52, n53, n54, n55, n56, n57, n58, n59, n60, n61, n62, n63, n64, n65, n66, n67, n68, n69, n70, n71, n72, n73, n74, n75, n76, n77, n78, n79, n80, n81, n82, n83, n84, n85, n86, n87, n88, n89, n90, n91, n92, n93, n94, n95, n96, n97, n98, n99, n100, n101, n102, n103, n104, n105, n106, n107, n108, n109, n110, n111, n112, n113, n114, n115, n116, n117, n118, n119, n120, n121, n122, n123, n124, n125, n126, n127, n128, n129, n130, n131, n132, n133, n134, n135, n136, n137, n138, n139, n140, n141, n142, n143, n144, n145, n146, n147, n148, n149, n150, n151, n152, n153, n154, n155, n156, n157, n158, n159, n160");
  
  set_dependencies_gradient_term_RHS(2, "");
  
  set_output_integral           (2,false);
  
  // Variable 3
  set_variable_name        (3,"sum_rec_ids");
  set_variable_type        (3,SCALAR);
  
  set_dependencies_value_term_RHS(3, "n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, n41, n42, n43, n44, n45, n46, n47, n48, n49, n50, n51, n52, n53, n54, n55, n56, n57, n58, n59, n60, n61, n62, n63, n64, n65, n66, n67, n68, n69, n70, n71, n72, n73, n74, n75, n76, n77, n78, n79, n80, n81, n82, n83, n84, n85, n86, n87, n88, n89, n90, n91, n92, n93, n94, n95, n96, n97, n98, n99, n100, n101, n102, n103, n104, n105, n106, n107, n108, n109, n110, n111, n112, n113, n114, n115, n116, n117, n118, n119, n120, n121, n122, n123, n124, n125, n126, n127, n128, n129, n130, n131, n132, n133, n134, n135, n136, n137, n138, n139, n140, n141, n142, n143, n144, n145, n146, n147, n148, n149, n150, n151, n152, n153, n154, n155, n156, n157, n158, n159, n160");
  
  set_dependencies_gradient_term_RHS(3, "");
  
  set_output_integral           (3,false);

}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================
// This function is analogous to 'explicitEquationRHS' and 'nonExplicitEquationRHS' in
// equations.h. It takes in "variable_list" and "q_point_loc" as inputs and outputs two terms in
// the expression for the postprocessing variable -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file (for
// submitting the terms) and the index in 'equations.h' for assigning the values/derivatives of
// the primary variables.

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

  //Number of order parameters
  unsigned int number_ops_ini = 118;
  
  //Number of empty ops (for nucleation)
  unsigned int number_ops_nuc = 42;
  
  //Number of ops (total)
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;
  
scalarvalueType ni;

scalarvalueType max_val = constV(-100.0);
scalarvalueType max_op = constV(100.0);
for (unsigned int i=1; i<=number_ops_total; i++){
    ni = variable_list.get_scalar_value(i);

    for (unsigned int v=0; v<ni.size();v++){
        if (ni[v] > max_val[v]){
            max_val[v] = ni[v];
            max_op[v] = i;
        }
    }
}

scalarvalueType feature_ids = constV(-1.0);
for (unsigned int v=0; v<ni.size();v++){
    for (unsigned int g=0; g<this->simplified_grain_representations.size(); g++){

        unsigned int max_op_nonvec = (unsigned int)std::abs(max_op[v]);

        if (this->simplified_grain_representations[g].getOrderParameterId() == max_op_nonvec){
            dealii::Point<dim> q_point_loc_nonvec;
            for (unsigned int d=0;d<dim;d++){
                q_point_loc_nonvec(d) = q_point_loc(d)[v];
            }

            double dist = 0.0;
            for (unsigned int d=0;d<dim;d++){
                dist += (q_point_loc_nonvec(d)-this->simplified_grain_representations[g].getCenter()(d))*(q_point_loc_nonvec(d)-this->simplified_grain_representations[g].getCenter()(d));
            }
            dist = std::sqrt(dist);

            if ( dist < (this->simplified_grain_representations[g].getRadius() + userInputs.buffer_between_grains/2.0) ){
                feature_ids[v] = (double)(this->simplified_grain_representations[g].getGrainId());
            }

        }

    }
}

scalarvalueType sum_n = constV(0.0);
for (unsigned int i=1; i<=number_ops_total; i++){
    ni = variable_list.get_scalar_value(i);
    sum_n += ni;
}
  
for (unsigned int v=0; v<ni.size();v++){
    if (sum_n[v] < 0.01){
        max_op[v] = -1.0;
        feature_ids[v] = -1.0;
    }
}

scalarvalueType sum_n_sq = constV(0.0);
for (unsigned int i=1; i<=number_ops_total; i++){
  ni = variable_list.get_scalar_value(i);
  sum_n_sq += ni*ni;
}

scalarvalueType sum_n_nuc = constV(0.0);
for (unsigned int i=number_ops_ini+1; i<=number_ops_total; i++){
  ni = variable_list.get_scalar_value(i);
  sum_n_nuc += ni;
}
// --- Submitting the terms for the postprocessing expressions ---

pp_variable_list.set_scalar_value_term_RHS(0, feature_ids);
pp_variable_list.set_scalar_value_term_RHS(1, max_op);
pp_variable_list.set_scalar_value_term_RHS(2, sum_n_sq);
pp_variable_list.set_scalar_value_term_RHS(3, sum_n_nuc);

}
