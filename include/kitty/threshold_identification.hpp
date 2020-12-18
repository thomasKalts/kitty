/* kitty: C++ truth table library
 * Copyright (C) 2017-2020  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file threshold_identification.hpp
  \brief Threshold logic function identification

  \author CS-472 2020 Fall students
*/

#pragma once

#include <vector>
#include <iostream>
#include <math.h>
#include <lpsolve/lp_lib.h> /* uncomment this line to include lp_solve */
//#include <lp_lib.h>
#include "traits.hpp"
#include "print.hpp"

#include "static_truth_table.hpp"
#include "dynamic_truth_table.hpp"
#include "properties.hpp"
#include "bit_operations.hpp"
#include "operations.hpp"
#include "cube.hpp"
#include "isop.hpp"

//using namespace std;



namespace kitty
{

/*! \brief Threshold logic function identification

  Given a truth table, this function determines whether it is a threshold logic function (TF)
  and finds a linear form if it is. A Boolean function is a TF if it can be expressed as

  f(x_1, ..., x_n) = \sum_{i=1}^n w_i x_i >= T

  where w_i are the weight values and T is the threshold value.
  The linear form of a TF is the vector [w_1, ..., w_n; T].

  \param tt The truth table
  \param plf Pointer to a vector that will hold a linear form of `tt` if it is a TF.
             The linear form has `tt.num_vars()` weight values and the threshold value
             in the end.
  \return `true` if `tt` is a TF; `false` if `tt` is a non-TF.
*/
template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool is_threshold( const TT& tt, std::vector<int64_t>* plf = nullptr )
{
  std::vector<int64_t> linear_form;
 // std::cout << to_binary(tt);

  /* TODO */
  /* if tt is non-TF: */
  //return false; 

  int variables = tt.num_vars();  //number of variables
  TT new_truth_table = tt;
  int variables_unateness[variables];

  //std::cout << "num vars is" << variables << std::endl;
  

  for (int i = 0; i<variables ; i++) {
  
	  auto F = cofactor1( new_truth_table , i );       //we will use the cofactors to check the unateness
	  auto F_star = cofactor0( new_truth_table , i );
	  
	  //std::cout << "F cofactor is " << to_binary(F) <<std::endl;
	  //std::cout << "F star cofactor is " << to_binary(F_star) <<std::endl;

	  if (implies( F, F_star )) {	// check negative unateness in i
			
		//std::cout << "F is negative unate in " << i <<std::endl;	
		new_truth_table = flip( new_truth_table , i ); //change bit in variable i
		variables_unateness[i] = 1;

			  }

	 else if (implies ( F_star , F )) {		//check positive unate in i
	 
		//std::cout << "F is positive unate in  "  << i   << std::endl; 
		variables_unateness[i] = 0;
	 }

	 else {		// if F is binate in 1 
		//std::cout << "F is biante !" << std::endl; 
	 	return false;
		}
  	}
	
	/*for (int j = 0 ;  j < variables ; j++){
		std::cout << variables_unateness[j] << std::endl; 
	} 	*/																					 

	//std::cout << "truth table new is" << to_binary(new_truth_table) << std::endl;
	
	 auto F_onset = isop(new_truth_table); // onset of new(fliped) truth table
	 auto F_offset = isop( unary_not(new_truth_table)); //offset of new (fliped)  truth table
	 
	 //std::cout << "size of onset is  " << F_onset.size() <<std::endl;
	 //std::cout << "size of offset is  " << F_offset.size() << std::endl; 

	//build ilp

  	lprec *lp;
	int Ncol, *colno = NULL; 
	REAL *row = NULL;
	int ret  = 0;

	Ncol = variables + 1;
	lp = make_lp( 0, Ncol);

	if (lp == NULL)     //couldn't construct a model 
		return false ; 
	
	
	//colno = (int *) malloc(Ncol * sizeof(*colno));
	//row = (REAL *) malloc(Ncol * sizeof(*row));

//	if (condition == 0) {  // naming the variables 
		
		/*	set_col_name(lp, 1, "w1" );
			set_col_name(lp, 2, "w2");
			set_col_name(lp, 3 , "w3");
		    set_col_name(lp, 4, "w4");
			set_col_name(lp, 5, "w5");
			set_col_name(lp, 6, "w6");
			set_col_name(lp, 7, "w7");
			set_col_name(lp, 8 ,"w8");
			set_col_name(lp, 9, "w9");	
		
			set_col_name( lp, variables + 1, "T");*/
			
			if (ret  ==0) {

				colno = (int *) malloc(Ncol * sizeof(*colno));
				row = (REAL *) malloc(Ncol * sizeof(*row));

			if ( (colno == NULL) || (row == NULL))
				return false;  
			}

	if (ret == 0) { 	//every variable is positive	   INITIAL CONSTRAINTS
	
		set_add_rowmode(lp, TRUE);

		for (int i = 0; i<Ncol; i++) {
			
			for (int j = 0; j<Ncol; j++){
				
			colno[j] = j + 1 ;
			
			if (i == j) {
				
				row[j] = 1;
			}
			
			else{
				
				row[j] = 0;
				
			}
			}
			
		add_constraintex(lp, Ncol, row, colno, GE, 0);
		
		}
	
		
		for (long unsigned int i = 0 ; i < F_onset.size() ; i++) {   //constraints for the onset set of the function
			
			auto ONSET = F_onset.at(i); 
			
			for (int k = 0 ; k< variables ; k++) {
				
				auto positive = ONSET.get_mask(k);
				colno[k] = k+1; 
				
				if (positive){
					row[k] = 1;
					
				}
				else  {
					row[k] = 0 ;
					
				}
			}
				
				colno[variables] = Ncol;
				row[variables] = -1;
				
				add_constraintex(lp, Ncol, row, colno, GE , 0);	
		}
			
	}
	
	
	if (ret == 0) {    // constraints with OFFSET set of function 
	
	//set_add_rowmode(lp, TRUE);
	
		for (long unsigned int i = 0 ; i < F_offset.size() ; i++) {
			
			auto OFFSET = F_offset.at(i); 
		
			for (int k = 0 ; k< variables ; k++) {
				
				colno[k] = k + 1; 
				auto negative = OFFSET.get_mask(k); 
				
				if (negative){
					row[k] = 0;
					
				}
				else {
					row[k] = 1 ;
					
				}
			}
				
				colno[variables] = Ncol;
				row[variables] = -1;
				
				add_constraintex(lp, Ncol, row, colno, LE , -1);	
		}
			
	} 
	
//	std::cout << to_binary(F_onset) << std::endl;

	if (ret == 0 ) {  		//objective function    			//objective function is fine 
		set_add_rowmode(lp, FALSE);

		for (int i = 0 ; i< Ncol; i++) {
			
			colno[i] = i+1;
			
			row[i]=1;

		}
		
		if ( !set_obj_fnex( lp, Ncol, row, colno) )
			return false ; 

	}


	if (ret  == 0 ) {
		
		set_minim(lp); //object direction to minimize
		
		//write_LP(lp, stdout);
		set_verbose(lp , IMPORTANT);
		ret = solve(lp);
		 if (ret == OPTIMAL) {
			ret = 0;
		}
		else {
			return false ;	
		} 
	} 


	if (ret == 0) {
		
		get_variables(lp, row);

		for (int i =0 ; i<Ncol ; i++) {
			linear_form.push_back(row[i]); 
		}
	}
	
	/*for (int j = 0 ;  j < linear_form.size() ; j++){
		std::cout << linear_form.at(j) << std::endl; 
	} 	*/

		for (int k = 0; k<variables; k++) {
			if (variables_unateness[k] == 1){
				linear_form.at(k) = -linear_form.at(k);	
				linear_form.at(variables) = linear_form.at(variables) + linear_form.at(k);
			}
		}

	
	if (row != NULL)
		free(row);
	
	if ( colno != NULL) 
		free(colno);	
	

	if ( lp != NULL) {
		delete_lp(lp);
	}
	
  

  /* if tt is TF: */
  /* push the weight and threshold values into `linear_form` */
  if ( plf )
  {
    *plf = linear_form;
  }
  return true;
}

} /* namespace kitty */
