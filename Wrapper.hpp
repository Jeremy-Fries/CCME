//  Wrapper.hpp
//  Mapping_Elites v1
//
//  Created by Jeremy Fries on 11/16/15.
//  Copyright © 2015 Jeremy Fries. All rights reserved.
//

class Wrapper;

#ifndef Wrapper_hpp
#define Wrapper_hpp

#include <stdio.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <ostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "Individual.hpp"
#include "Map_space.hpp"
#include "Map_Elites.hpp"
#include "ME_Input.hpp"
#include "ME_Output.hpp"
#include "ME_TestBed.hpp"
#include "ME_TestBed_Sim.hpp"
#include "ME_TestBed_NN.hpp"

/*
 
 - This is a TEST Wrapper.
 - Functions in this Wrapper will be for testing the Map_Elite framework and functions.

*/

using namespace std;

class Wrapper{
    friend class Map_Elites;
// --------------------------------------------------
protected:
    Map_Elites ME;
    Map_Elites* pME = &ME;
    ME_Input ME_I;
    ME_Input* pME_I = &ME_I;
    ME_Output ME_O;
    ME_Output* pME_O = &ME_O;
    ME_TestBed_NN NN;
    ME_TestBed_NN* pN = &NN;
    ME_TestBed_Sim Sim;
    ME_TestBed_Sim* pS = &Sim;
    
    int size_of_genome1, size_of_genome2; // Possibly get rid of.
    
// --------------------------------------------------
public:
    int hidden_layer_size;      // Possibly get rid of.
//    vector<double> best_fit;    // Possibly not used.
// --------------------------------------------------
    void initialize_wrapper(int,int);
    void wrapper_sets_I_params(int size1, int size2, double mut_mag1, double mut_mag2, int mut_amo1, int mut_amo2); // Change for ME_Input
//    vector<double>& get_individual_1_IH(vector<double>);    // Possibly get rid of.
//    vector<double>& get_individual_1_HO(vector<double>);    // Possibly get rid of.
//    vector<double>& get_individual_2_IH(vector<double>);    // Possibly get rid of.
//    vector<double>& get_individual_2_HO(vector<double>);    // Possibly get rid of.
// --------------------------------------------------
            // Main Wrapper functions
    void wrapper_runs_sim(vector<double>,vector<double>);   // Not used.
    void phenotype_calculation();
    void fitness_calculation();
    void rand_bin();
    void fill_MAP();
    void mutate_MAP();
    void run_single_individual();   // Not used.
    void print_stuff();
    void print_entire_map_solution();
    void always_last();
    void load_genome1();
    void load_genome2();
    void load_bins();
    void write_from_old_genomes();
    void clear_map();
// --------------------------------------------------
private:
    int isize_1, isize_2, imutate_amount_1, imutate_amount_2;
    double imutate_mag_1,imutate_mag_2;
    int bin1, bin2;
    double fit_rating;
    double phenotype_1, phenotype_2;
    int map_solutions;
    int in_layer_size, out_layer_size;
    vector < vector <double> > Map_of_genome1;
    vector < vector <double> > Map_of_genome2;
    vector < vector<int> > Map_of_bins;     // Not used.
    
};
// ------------------------------------------------------------------------------------------------ ^^ Declarations
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ------------------------------------------------------------------------------------------------ vv Definitions

// -------------------------------------------------------------------------------------------------------------------------------------------------    SETTINGS
// -------------------------------------------------------------------------------------------------------------------------------------------------
            // Initialize Wrapper
void Wrapper::initialize_wrapper(int FILL, int MUTATE){
    
    int states = 7;
    int outs = 2;
    
    hidden_layer_size = 5;
    
    /// PHENOTYPE LIMITS
    //- TODO - CCME carrying amount
    pME->set_map_params(0, 100, 0, 100, 3, 2, FILL, MUTATE);             //-------- To Change Map Settings
    pME->set_carrying_capacity(10);

    // (dim1_min, dim1_max, dim2_min, dim2_max, resolution 1,2, fill generation, mutate generation)
    //pME->display_Map_params();        // TODO - delete and add print()
//    
//    in_layer_size=(states+1)*hidden_layer_size;
//    out_layer_size=(hidden_layer_size+1)*outs;
//    int num_weights=in_layer_size+out_layer_size;
//    
//    wrapper_sets_I_params(num_weights, num_weights, 0.2, 0.2, 10, 10);        //-------- To Change Individual Settings
    
    wrapper_sets_I_params(4, 4, 0.2, 0.2, 10, 10);        //-------- To Change Individual Settings

    // individual_size 1,2, mutate_magnitude 1,2, mutation_amount 1,2)
    // int size1, int size2, double mut_mag1, double mut_mag2, int mut_amo1, int mut_amo2
    
}
// -------------------------------------------------------------------------------------------------------------------------------------------------    SETTINGS
// -------------------------------------------------------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------
            // Fitness Function
void Wrapper::fitness_calculation(){    // Call from input file. %$--ME_Input--$%
    fit_rating=0;
    double r1 = ((double)rand() / RAND_MAX);
    double r2 = ((double)rand() / RAND_MAX);
    fit_rating=r1-r2;
    
    //- TODO - functions dependent on the test.
}
// -----------------------------------------------------------------------------------
            // Phenotypes
void Wrapper::phenotype_calculation(){
    
    phenotype_1 = 0;    // Call from input file. %$--ME_Input--$%
    phenotype_2 = 0;    // Call from input file. %$--ME_Input--$%
    
    //- TODO - functions dependent on the test.
//    cout << endl << "P1,2: " << phenotype_1 << " , " << phenotype_2;// << endl;
  
        // ----------------------------------------------------
    /// for a Binary Selection Parallel Simulated Annealing.....
    double r1 = ((double)rand() / RAND_MAX);
    double r2 = ((double)rand() / RAND_MAX);
    phenotype_1=r1*100;
    phenotype_2=r2*100;
        // ----------------------------------------------------
    
    //    cout << endl << "P1,2: " << phenotype_1 << " , " << phenotype_2;// << endl;
}
// --------------------------------------------------
// Clear Map
void Wrapper::clear_map(){
    ME.Map.clear();
}
// -----------------------------------------------------------------------------------
            // Fill Map
void Wrapper::fill_MAP(){
    int fill_gen = ME.get_fill_generation();
    int fill_round=0;
    
    for (int g=0; g<fill_gen; g++){
        Individual I;
        
        I.set_individual_params(isize_1, isize_2, imutate_mag_1, imutate_mag_2, imutate_amount_1, imutate_amount_2);
        I.build_individual();
        //I.display_individual1();    // comment out
        //I.display_individual2();    // comment out

        //---------------------------------------------------
        //- TODO - Create functions in ME_TestBed_Sim and NN for these functions that replicate them. VVVV
        
//        Sim.initialize_sim();
//        NN.take_input_limits(Sim.currentstate.state_variables_LowLimit, Sim.currentstate.state_variables_UpLimit);
//        NN.take_output_limits(Sim.currentstate.control_LowLimits, Sim.currentstate.control_UpLimits);
//        NN.take_num_hidden_units(this->hidden_layer_size); /// repeated from initialize, but no harm.
//        NN.take_num_controls(Sim.currentstate.num_of_controls); /// repeated from initialize, but to no harm.
//        NN.take_weights(I.get_individual1(), I.get_individual2());
        
        //---------------------------------------------------
        
            /// %%% /// %%% BEGIN SIMULATION LOOP %%% /// %%% ///

            /// %%% /// %%% END SIMULATION LOOP %%% /// %%% ///

      
        fitness_calculation();
        I.set_fit_rating(fit_rating);
        //I.display_fit_rating();

        phenotype_calculation();
        I.set_phenotypes(phenotype_1,phenotype_2);
        //I.display_phenotype1();
        //I.display_phenotype2();
        
        ME.challenger = I;
        ME.place_individual_in_map();
        
        // Carrying Capcity function.
        
        fill_round++;
        
        if(rand()%1000 < 100){
            cout << endl << "FILL Round is ---- " << g << " of " << fill_gen << endl;
        }
        
    }
    cout << "completed " << fill_round << " FILL rounds" << endl;
}
// --------------------------------------------------
            // Mutate Map
void Wrapper::mutate_MAP(){
    int mut_gen = ME.get_mutate_generation();
    int mutation_round=0;
    
    for (int g=0; g<mut_gen; g++){
        
        rand_bin();                             // random bin location
        ME.individual_from_map(bin1, bin2);     // copy Individual's vectors at bin location    // CHECK to see if Quiver will affect anything.
        ME.challenger.mutate1();
        ME.challenger.mutate2();
        
        //---------------------------------------------------
        //- TODO - Create functions in ME_TestBed_Sim and NN for these functions that replicate them. VVVV
//        Sim.initialize_sim();
        //NN.take_weights(ME.challenger.get_individual1(), ME.challenger.get_individual2());
//        NN.take_weights(get_individual_1_IH(ME.challenger.get_individual1()), get_individual_1_HO(ME.challenger.get_individual1()));
//        NN.take_input_limits(Sim.currentstate.state_variables_LowLimit, Sim.currentstate.state_variables_UpLimit);
//        NN.take_output_limits(Sim.currentstate.control_LowLimits,Sim.currentstate.control_UpLimits);
        
        //---------------------------------------------------
        
            /// %%% /// %%% BEGIN SIMULATION LOOP %%% /// %%% ///

            /// %%% /// %%% END SIMULATION LOOP %%% /// %%% ///
        
        
        fitness_calculation();                           // fitness for new Individual
        ME.challenger.set_fit_rating(fit_rating);   // Why is this different than above??????????
        //ME.challenger.display_fit_rating();
        
        phenotype_calculation();
        ME.challenger.set_phenotypes(phenotype_1,phenotype_2);
        //ME.challenger.display_phenotype1();
        //ME.challenger.display_phenotype2();
        ME.place_individual_in_map();
        mutation_round++;
        
        if(rand()%1000 < 100){
            cout << endl << "MUTATE Round is ---- " << g << " of " << mut_gen << endl;
        }
    }
    cout << "completed " << mutation_round << " MUTATION rounds" << endl;
}
// --------------------------------------------------
            // Print stuff
void Wrapper::print_stuff(){
    ME.how_many_full_bins();
    ME.best_fit_bin();
    ME.print_fit_ratings_of_map();
    ME.print_best_occupants_fitness();
    ME.print_all_occupants();
    ME.print_best_parents_fitness();
    ME.print_best_parents_id();
    ME.print_best_full_trace();
    ME.print_heat_map();
    ME.print_corresponding_genome1();
    ME.print_corresponding_genome2();
    ME.print_corresponding_bins();
    //ME.print_CC_vec();
}
// --------------------------------------------------
void Wrapper::print_entire_map_solution(){  // why does this recalculate all individuals in the bins? Can't it just grab the solutions from memory????????
    map_solutions=0;
    
    for(int element=0; element<ME.full_bins.size();element++){
        cout << endl << "Bin ID: " << ME.full_bins.at(element).id << endl;
    
//        Sim.initialize_sim();
        
        // --------------------------------------------- Display
//        cout << endl << "Genome 1 is:  " << endl;
//        ME.full_bins.at(element).current_individual.at(0).display_individual1();
//        cout << endl << endl << "Genome 2 is:  " << endl;
//        ME.full_bins.at(element).current_individual.at(0).display_individual2();
        // ---------------------------------------------
        
        /// model from previous.
        //NN.take_weights(ME.full_bins.at(element).current_individual.at(0).get_individual1() , ME.full_bins.at(element).current_individual.at(0).get_individual2());
        
//        vector<double> ih = get_individual_1_IH(ME.full_bins.at(element).current_individual.at(0).get_individual1());
//        vector<double> ho = get_individual_1_HO(ME.full_bins.at(element).current_individual.at(0).get_individual1());
//        NN.take_weights(ih,ho);
//        NN.take_input_limits(Sim.currentstate.state_variables_LowLimit, Sim.currentstate.state_variables_UpLimit);
//        NN.take_output_limits(Sim.currentstate.control_LowLimits,Sim.currentstate.control_UpLimits);
    
    
        /// %%% /// %%% BEGIN SIMULATION LOOP %%% /// %%% ///
        
        /// %%% /// %%% END SIMULATION LOOP %%% /// %%% ///
        map_solutions++;
        fitness_calculation();
        phenotype_calculation();
    }
    cout << endl << endl << "Number of Solutions OUTPUT: " << map_solutions << endl;
    cout << endl << "Number of Solutions is: " << ME.full_bins.size()  << endl;
}
// --------------------------------------------------
void Wrapper::load_genome1(){
    ifstream co("print_corresponding_genome1.txt"); // Get rid of txt file location. <<this marker is old???????
    
    double read;
    vector<double> apush;
    
    cout << "WE GET INTO THIS FUNCTION " << endl;
    while(co >> read){
        apush.push_back(read);
        cout <<"SIZE OF APUSH IS: " << apush.size() << endl;
        cout <<"SIZE OF GENOME 1 IS: " << Map_of_genome1.size() << endl;
        if(apush.size()>=size_of_genome1){
            Map_of_genome1.push_back(apush);
            cout << "SIZE OF APUSH IS: " << apush.size() << endl;
            apush.clear();
        }
    }
    cout << endl << "Starting size Map of old genome 1 is: " << Map_of_genome1.size() << endl;
    cout << endl << "Starting size of old genome 1 is: " << Map_of_genome1.at(0).size() << endl;
}
// --------------------------------------------------
void Wrapper::load_genome2(){
    ifstream co("print_corresponding_genome2.txt");
    
    double read;
    vector<double> apush;
    
    while(co >> read){
        apush.push_back(read);
        
        if(apush.size()>=size_of_genome2){
            Map_of_genome2.push_back(apush);
            apush.clear();
        }
    }
    cout << endl << "Starting size Map of old genome 2 is: " << Map_of_genome2.size() << endl;
    cout << endl << "Starting size of old genome 2 is: " << Map_of_genome2.at(0).size() << endl;
}
// --------------------------------------------------
//write
void Wrapper::write_from_old_genomes(){
    /// Read from txt file
    
    //-----*********-----//
    int RunCountEach = 1;       //Change this value to run each set of weights multiple times.
    //-----*********-----//
    
    Map_of_genome1.clear();
    Map_of_genome2.clear();
    
    load_genome1();
    cout << endl << "Starting size of old genome 1 is: " << Map_of_genome1.at(0).size() << endl;
    load_genome2();
    cout << endl << "Starting size of old genome 2 is: " << Map_of_genome2.at(0).size() << endl;
    
    int old_placed_counter=0;
    if (Map_of_genome1.size()== Map_of_genome2.size()){
        for(int j=0;j<RunCountEach;j++){
            for (int i=0;i<Map_of_genome1.size();i++){
                Individual I;
                I.set_individual_params(Map_of_genome1.at(i).size(), Map_of_genome2.at(i).size(), imutate_mag_1, imutate_mag_2, imutate_amount_1, imutate_amount_2);
                I.build_individual_1_from_another(Map_of_genome1.at(i));
                I.build_individual_2_from_another(Map_of_genome2.at(i));
        
//                Sim.initialize_sim();
                
//                NN.take_input_limits(Sim.currentstate.state_variables_LowLimit, Sim.currentstate.state_variables_UpLimit);
//                NN.take_output_limits(Sim.currentstate.control_LowLimits, Sim.currentstate.control_UpLimits);
//                NN.take_num_hidden_units(this->hidden_layer_size); /// repeated from initialize, but no harm.
//                NN.take_num_controls(Sim.currentstate.num_of_controls); /// repeated from initialize, but to no harm.
//                NN.take_weights(I.get_individual1(), I.get_individual2());
                
                /// %%% /// %%% BEGIN SIMULATION LOOP %%% /// %%% ///

                /// %%% /// %%% END SIMULATION LOOP %%% /// %%% ///
            

                
                fitness_calculation();
                I.set_fit_rating(fit_rating);
                //I.display_fit_rating();
                
                phenotype_calculation();
                I.set_phenotypes(phenotype_1,phenotype_2);
                //I.display_phenotype1();
                //I.display_phenotype2();
                
                ME.challenger = I;
                ME.place_individual_in_map();
                old_placed_counter++;
            }
        }
    }
    else {
        cout << endl << "Error, Genomes from txt file are NOT the same size." << endl;
    }
    cout << endl << "Quantity of Old Genomes Placed is: " << old_placed_counter << endl;
}
// --------------------------------------------------
            // Wrapper sets Individual
void Wrapper::wrapper_sets_I_params(int size1, int size2, double mut_mag1, double mut_mag2, int mut_amo1, int mut_amo2){
    cout << endl << " size of first layer is: " << size1 << endl;
    cout << endl << " size of last layer is: " << size2 << endl;
    
    isize_1=size1;
    isize_2=size2;
    imutate_mag_1=mut_mag1;
    imutate_mag_2=mut_mag2;
    imutate_amount_1=mut_amo1;
    imutate_amount_2=mut_amo2;
    /// assigns values to actual individual.
    
    size_of_genome1=size1;  // For ME input_genome()
    size_of_genome2=size2;
}
// --------------------------------------------------
            // Random Bin
void Wrapper::rand_bin(){
    bin1=0;
    bin2=0;
    int b1_min=0;
    int b2_min=0;
    int b1_max=ME.get_resolution1();
    int b2_max=ME.get_resolution2();
    bin1=rand()%(b1_max-b1_min)+b1_min;
    bin2=rand()%(b2_max-b2_min)+b2_min;
    //cout << "bin to mutate is: (" << bin1 << "," << bin2 << ")" << endl;
}
// --------------------------------------------------
void Wrapper::always_last(){
//    Sim.myfile.close();
//    Sim.windfile.close();
}
// --------------------------------------------------


#endif /* Wrapper_hpp */
