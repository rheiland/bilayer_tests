/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

    cell_defaults.functions.cell_division_function = custom_division_function; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	// cell_defaults.functions.update_phenotype = phenotype_function; 
	// cell_defaults.functions.custom_cell_rule = custom_function; 
	// cell_defaults.functions.contact_function = contact_function; 
    // cell_defaults.functions.cell_division_function = custom_division_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

// void dist_to_BM_function( Cell* pCell, Phenotype& phenotype, double dt )
// { 
//     extern std::vector<double> xctr,yctr;
//     double xdel,ydel,zdel, dist2;

//     // if (pCell->custom_data["dist_BM"] > 0)
//     //     return;

//     // std::cout << __FUNCTION__ << ": t= " << PhysiCell_globals.current_time << std::endl;
//     // std::cout << __FUNCTION__ << ": t= " << PhysiCell_globals.current_time << std::endl;
//     double d2_min = 1.0e6;
//     int npts = xctr.size();
//     for (int idx=0; idx<npts; idx++ )
//     {
//         xdel = pCell->position[0] - xctr[idx];
//         ydel = pCell->position[1] - yctr[idx];
//         zdel = pCell->position[2];
//         dist2 = xdel*xdel + ydel*ydel + zdel*zdel;
//         if (dist2 < d2_min) d2_min = dist2;
//     }
//     double dist = sqrt(d2_min);
//     // std::cout << "         dist= " << dist << std::endl;

//     pCell->custom_data["dist_BM"] = dist;
//     if (dist > parameters.doubles("tube_radius")) pCell->is_movable = false;
//     // return dist; 
//     // return 0.0; 
//     return; 
// }

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	

	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();
	
	return; 
}

void compute_level_set( void )
{
    static int bm_dist_idx = BioFVM::microenvironment.find_density_index( "bm_dist" ); 

    int idx_xmax = microenvironment.mesh.x_coordinates.size() - 1; 
    int idx_ymax = microenvironment.mesh.y_coordinates.size() - 1;
    int idx_zmax = microenvironment.mesh.z_coordinates.size() - 1;

	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 

    int ivox = 0;
    double xval,yval;
    double xctr = 100;
    double yctr = 0;
    for (int idz=0; idz<microenvironment.mesh.z_coordinates.size(); idz++)
        for (int idy=0; idy<microenvironment.mesh.y_coordinates.size(); idy++)
        {
            yval = Ymin + idy * microenvironment.mesh.dy;
            for (int idx=0; idx<microenvironment.mesh.x_coordinates.size(); idx++)
            {
                xval = Xmin + idx * microenvironment.mesh.dx;
                double d2_min = 1.e6;
                for( int i=0; i < (*all_cells).size(); i++ )
                {
                    if ((*all_cells)[i]->type_name == "epi")
                    {

                    double xctr = (*all_cells)[i]->position[0];
                    double yctr = (*all_cells)[i]->position[1];
                    double xdel = xval - xctr;
                    double ydel = yval - yctr;
                    // double d = sqrt(xdel*xdel + ydel*ydel);
                    double d2 = xdel*xdel + ydel*ydel;
                    if (d2 < d2_min)
                        d2_min = d2;
                    }
                }

                BioFVM::microenvironment.density_vector(ivox)[bm_dist_idx] = sqrt(d2_min);
                ivox++;
            }
        }
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; }

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void custom_division_function( Cell* pCell1, Cell* pCell2 )
{ 
    std::cout << __FUNCTION__ << ": parent,child IDs= " << pCell1->ID << ", " << pCell2->ID << "; type= " << pCell1->type_name << std::endl;
    double x1 = pCell1->position[0];
    double y1 = pCell1->position[1];
    std::cout << "   parent x,y= " << x1 << ", " << y1 << std::endl;
	for( int idx=0; idx<pCell1->state.neighbors.size(); idx++ )  // for all j nbrs
	{
		Cell* pC = pCell1->state.neighbors[idx]; 
        if (pC->type != pCell1->type)  // skip if not same type
            continue;

        // compute chord of intersection (if any)
        // radii of cells
        double r2 = pC->phenotype.geometry.radius;
        // centers of cells
        double x2 = (*pC).position[0];
        double y2 = (*pC).position[1];
        std::cout << "nbr ID> " << pC->ID << ">  at x,y= " << x2 << ", " << y2 << std::endl;
        double xdiff = x1-x2;
        double ydiff = y1-y2;
        // double d = sqrt(xdiff*xdiff + ydiff*ydiff);
        pCell2->position[0] = x1 + xdiff/2.0;
        pCell2->position[1] = y1 + ydiff/2.0;
    }
}
// void custom_division_function_v0( Cell* pCell1, Cell* pCell2 )
// { 
//     static int idx_default = find_cell_definition_index("default");
//     static int idx_ctype1 = find_cell_definition_index("ctype1");
//     std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell IDs= " << pCell1->ID << ", " << pCell2->ID << std::endl;

//     // Asymmetric division
//     if (UniformRandom() < 0.5)
//     {
//         pCell2->convert_to_cell_definition( *cell_definitions_by_index[idx_default] ); 
//     }
//     else
//     {
//         pCell2->convert_to_cell_definition( *cell_definitions_by_index[idx_ctype1] ); 
//     }

//     return; 
// }

//--------------------------------
// void custom_function_monolayer( Cell* pCell, Phenotype& phenotype, double dt )
// { 
//     // std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << std::endl;

//     static double pi = 3.1415926535897932384626433832795; 
// 	// static double two_pi = 6.283185307179586476925286766559; 

//     if ((*all_cells).size() < 2)
//         return;

//     double r1 = phenotype.geometry.radius;
//     double r1_2 = r1*r1;
//     double x1 = (*pCell).position[0];
//     double y1 = (*pCell).position[1];
//     double gamma = 0.0;
//     double beta = 0.0;

// 	for( int idx=0; idx<pCell->state.neighbors.size(); idx++ )  // for all j nbrs
// 	{
// 		Cell* pC = pCell->state.neighbors[idx]; 

//         // compute chord of intersection (if any)
//         // radii of cells
//         double r2 = pC->phenotype.geometry.radius;
//         // centers of cells
//         double x2 = (*pC).position[0];
//         double y2 = (*pC).position[1];
//         double xdiff = x1-x2;
//         double ydiff = y1-y2;
//         double d = sqrt(xdiff*xdiff + ydiff*ydiff);
//         if (d < r1+r2)
//         {
//             // phi_ij = ( d_ij*d_ij - r_j*r_j + r_i*r_i ) / (2 * d_ij * r_i)
//             double phi = (d*d - r2*r2 + r1_2 ) / (2 * d * r1);
//             // f_i = 1.0 - 1.0/np.pi * np.sqrt( 1.0 - phi_ij*phi_ij )   # for all nbrs j:  1 - 1/pi * SUM_j (sqrt(1 - phi_ij)
//             gamma += sqrt( 1.0 - phi*phi);
//             // std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << ", gamma= " << gamma << std::endl;
//             beta += acos(phi) - phi * sqrt(1 - phi*phi);

//         }
//     }
//     double gamma_inv = 1.0/pi * gamma;   
//     gamma = 1.0 - gamma_inv;   // free surface fraction
//     // pCell->custom_data["gamma"] = 1.0 - 1.0/pi * gamma;
//     pCell->custom_data["f_i"] = gamma;
//     // pCell->custom_data["gamma_inv"] = gamma_inv;

//     beta = 1.0 - 1.0/pi * beta;
//     // pCell->custom_data["beta"] = 1.0 - 1.0/pi * beta;
//     pCell->custom_data["a_i"] = beta;

//     {
//         set_single_behavior( pCell , "cycle entry" , 0.01124);  // 1/89
//         // if ((beta < 0.9) || (gamma < 0.9))
//         // if ((beta < 0.5) || (gamma < 0.9))
//         if ((beta < beta_threshold) || (gamma < gamma_threshold))
//         {
//             set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
//         }
//         return;
//     }