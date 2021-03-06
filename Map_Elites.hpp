//
//  Map_Elites.hpp
//  Mapping_Elites v1
//
//  Created by Jeremy Fries on 10/20/15.
//  Copyright © 2015 Jeremy Fries. All rights reserved.
//

class Map_Elites;

#ifndef Map_Elites_hpp
#define Map_Elites_hpp

#include <stdio.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <ostream>

#include "Individual.hpp"
#include "Map_space.hpp"
#include "Carrying_Capactiy.hpp"

/*

 - 2D Map_Elites Algorithm
 
*/

using namespace std;

class Map_Elites{
    friend class Wrapper;
    friend class Map_space;
// --------------------------------------------------
protected:
    int resolution1, resolution2;
    double dim1_min, dim2_min, dim3_min;
    double dim1_max, dim2_max, dim3_max;
    int num_spacing1, num_spacing2, num_spacing3;
    double spacing1,  spacing2, spacing3;
    int fill_generation, mutate_generation;
    double center_dist;
    int best_fit_index;
    vector<Map_space> full_bins;
    vector<double> distance_between_centerbin_phenotype;
    vector<double> fit_ratings_in_map;
    vector<double> best_individual1, best_individual2;
// --------------------------------------------------
public:
            // Dimensions of Map
// Min of dim1
    void set_min_dim1(double);
    double get_min_dim1();
    void display_min_dim1();
// Min of dim2
    void set_min_dim2(double);
    double get_min_dim2();
    void display_min_dim2();
// Min of dim3
    void set_min_dim3(double);
    double get_min_dim3();
    void display_min_dim3();
// Max of dim1
    void set_max_dim1(double);
    double get_max_dim1();
    void display_max_dim1();
// Max of dim2
    void set_max_dim2(double);
    double get_max_dim2();
    void display_max_dim2();
// Max of dim3
    void set_max_dim3(double);
    double get_max_dim3();
    void display_max_dim3();
// --------------------------------------------------
            // Resolution
    void set_resolution(int,int);
    int get_resolution1();
    void display_resolution1();
    int get_resolution2();
    void display_resolution2();
// --------------------------------------------------
            // Fill Generation
    void set_fill_generation(int);
    int get_fill_generation();
    void display_fill_generation();
// --------------------------------------------------
            // Mutate Generation
    void set_mutate_generation(int);
    int get_mutate_generation();
    void display_mutate_generation();
// --------------------------------------------------
            // Map Parameters
    void set_map_params(double dim1_min,double dim1_max,double dim2_min,double dim2_max,int resolution1,int resolution2, int fill_generation, int mutate_generation);
    void display_Map_params();
    
    // TODO - write map_prams to txt file
    
// --------------------------------------------------
        // Makes Map
    void initialize_map();
    void place_individual_in_map();
    void individual_from_map(int p1, int p2);
    void build_map_from_old(vector<vector<double>>, vector<vector<double>>, vector<vector<int>>);
// --------------------------------------------------
        // CCME
    void set_carrying_capacity(int);
    int get_carrying_capacity();
    void display_carrying_capacity();
    void find_quiver_size(int);
    void add_CC(int);
    void arrange_CC();
    void calc_current_carrying_capacity();
    void calc_carrying_bins();
    void find_id_0();
    void find_id(int);
    void find_id_end();
    void find_placement(int);
    bool check_CC(int);
    void dynamic_CC(int,int);
    void arrange_test_CC_vec();
    void test_qsize_correspond();
// --------------------------------------------------
        // Create full bin vector
    void create_full_bin();
    void find_center_bin(int p1,int p2);
    void find_pheno_dist_to_center_bin(int p1, int p2, int cb1, int cb2);
// --------------------------------------------------
        // To text file
    void print_fit_ratings_of_map();
    void print_best_occupants_fitness();
    void best_fit_bin();
    void print_all_occupants();
    void print_best_parents_id();
    void print_best_parents_fitness();
    void print_best_full_trace();
    void print_heat_map();
    void print_corresponding_genome1();
    void print_corresponding_genome2();
    void print_corresponding_bins();
    void print_CC_vec();      // TODO - J 9/8
// --------------------------------------------------
    // how many bins are full?
    void how_many_full_bins();
// --------------------------------------------------
    // Get Best Individuals
    vector<double>& get_best_individual1();
    vector<double>& get_best_individual2();
// --------------------------------------------------
    /// LYLY Tracking
    Individual challenger;
    
private:
    vector< vector<Map_space> > Map;
    int current_bin1, current_bin2;
//    vector < vector <double> > Map_of_genome1;
//    vector < vector <double> > Map_of_genome2;
    
    vector <Carrying_Capacity> CC_vec;
    vector <Carrying_Capacity> test_CC_vec;
    int carrying_capacity_size;
    int current_carrying_capacity_size;
    int carrying_bins;
    int id_0, id_0_row, id_0_ele;
    int id_row, id_ele;
    int id_end, id_end_row, id_end_ele;
    int id_element;
};
// ------------------------------------------------------------------------------------------------ ^^ Declarations
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ------------------------------------------------------------------------------------------------ vv Definitions
            // Dimensions of Map
    // Min of dim1
void Map_Elites::set_min_dim1(double min1){
    dim1_min = min1;
}
double Map_Elites::get_min_dim1(){
    return dim1_min;
}
void Map_Elites::display_min_dim1(){
    cout << endl << "Min of dim1 is: " << dim1_min << endl;
}
    // Min of dim2
void Map_Elites::set_min_dim2(double min2){
    dim2_min = min2;
}
double Map_Elites::get_min_dim2(){
    return dim2_min;
}
void Map_Elites::display_min_dim2(){
    cout << endl << "Min of dim2 is: " << dim2_min << endl;
}
    // Min of dim3
void Map_Elites::set_min_dim3(double min3){
    dim3_min = min3;
}
double Map_Elites::get_min_dim3(){
    return dim3_min;
}
void Map_Elites::display_min_dim3(){
    cout << endl << "Min of dim3 is: " << dim3_min << endl;
}
    // Max of dim1
void Map_Elites::set_max_dim1(double max1){
    dim1_max = max1;
}
double Map_Elites::get_max_dim1(){
    return dim1_max;
}
void Map_Elites::display_max_dim1(){
    cout << endl << "Max of dim1 is: " << dim1_max << endl;
}
    // Max of dim2
void Map_Elites::set_max_dim2(double max2){
    dim2_max = max2;
}
double Map_Elites::get_max_dim2(){
    return dim2_max;
}
void Map_Elites::display_max_dim2(){
    cout << endl << "Max of dim2 is: " << dim2_max << endl;
}
    // Max of dim3
void Map_Elites::set_max_dim3(double max3){
    dim3_max = max3;
}
double Map_Elites::get_max_dim3(){
    return dim3_max;
}
void Map_Elites::display_max_dim3(){
    cout << endl << "Max of dim3 is: " << dim3_max << endl;
}
// --------------------------------------------------
            // Set Resolution
    /// determines size of Map_space and how many bins.
void Map_Elites::set_resolution(int r1, int r2) {
    resolution1=r1;
    resolution2=r2;
    num_spacing1=r1;
    num_spacing2=r2;
    //num_spacing3=r3;
    spacing1=(dim1_max-dim1_min)/resolution1;
    spacing2=(dim2_max-dim2_min)/resolution2;
    // spacing3=(dim3_max-dim3_min)/resolution3;
}
int Map_Elites::get_resolution1(){
    return resolution1;
}
void Map_Elites::display_resolution1(){
    cout << endl << "Resolution 1 is: " << resolution1 << endl;
}
int Map_Elites::get_resolution2(){
    return resolution2;
}
void Map_Elites::display_resolution2(){
    cout << endl << "Resolution 2 is: " << resolution2 << endl;
}
// --------------------------------------------------
            // Fill Generation
    /// fills Map with individuals for desired amount.
void Map_Elites::set_fill_generation(int f){
    fill_generation=f;
}
int Map_Elites::get_fill_generation(){
    return fill_generation;
}
void Map_Elites::display_fill_generation(){
    cout << endl << "Fill generation is: " << fill_generation << endl;
}
// --------------------------------------------------
            // Mutate Generation
    /// muatates individuals in Map for desired amount.
void Map_Elites::set_mutate_generation(int m){
    mutate_generation=m;
}
int Map_Elites::get_mutate_generation(){
    return mutate_generation;
}
void Map_Elites::display_mutate_generation(){
    cout << endl << "Mutate generation is: " << mutate_generation << endl;
}
// --------------------------------------------------
            // Set Map Parameters
    /// Itialization function used outside of class.
void Map_Elites::set_map_params(double d1_min, double d1_max, double d2_min, double d2_max, int res1, int res2, int fill_gen, int mutate_gen){
    set_min_dim1(d1_min);
    set_max_dim1(d1_max);
    set_min_dim2(d2_min);
    set_max_dim2(d2_max);
    set_resolution(res1, res2);
    set_fill_generation(fill_gen);
    set_mutate_generation(mutate_gen);
    initialize_map();
}
// --------------------------------------------------
            // Display Map Parameters
void Map_Elites::display_Map_params(){
    cout << endl << "------------------------------- Map Parameters" << endl;
    display_min_dim1();
    display_max_dim1();
    display_min_dim2();
    display_max_dim2();
    display_resolution1();
    display_resolution2();
    display_fill_generation();
    display_mutate_generation();
    cout << endl << "-------------------------------" << endl;
}
// --------------------------------------------------
            // Initialze Map
    /// Builds Map based on parameters, ALSO creates and sets Map_space LB and UB
void Map_Elites::initialize_map(){
    double pre_LB1=dim1_min;
    double pre_LB2=dim2_min;
    int id_number=0;
    Map.reserve(num_spacing1);
    for(int d1=0; d1<num_spacing1; d1++){
        vector<Map_space> Row;
        Row.reserve(num_spacing2);
        for(int d2=0; d2<num_spacing2; d2++){
            Map_space M;
            M.previous_genome1.clear();  // clears occupant vector
            M.previous_genome2.clear();  // clears occupant vector
            M.set_id(id_number);
            M.bin1=d1;
            M.bin2=d2;
            M.set_LB1(pre_LB1);
            M.set_LB2(pre_LB2);
            M.set_center_bin1(pre_LB1+(spacing1/2));            // Sets center of Bin
            M.set_center_bin2(pre_LB2+(spacing2/2));
            double calc_UB1= pre_LB1+spacing1;
            M.set_UB1(calc_UB1);
            double calc_UB2= pre_LB2+spacing2;
            M.set_UB2(calc_UB2);
            M.build_map_space();
            Row.push_back(M);
            pre_LB2+=spacing2;
            id_number++;
        }
        Map.push_back(Row);
        pre_LB1+=spacing1;
        pre_LB2=dim2_min;
    }
    cout << endl << "Map is made" << endl;
}
// --------------------------------------------------
        // Carrying Capacity
void Map_Elites::set_carrying_capacity(int cc){
    carrying_capacity_size=cc;
}
int Map_Elites::get_carrying_capacity(){
    return carrying_capacity_size;
}
void Map_Elites::display_carrying_capacity(){
    cout << endl << "carrying_capacity_size is: " << carrying_capacity_size << endl;
}
// --------------------------------------------------
            // Add CC
void Map_Elites::add_CC(int bid){
    Carrying_Capacity CC;
    CC.bin=bid;
    CC.bin_quiver_size=0;
    CC_vec.push_back(CC);
}
// --------------------------------------------------
            // Arrange CC
    /// arrange CC_vec from high to low quiver_size
bool sort_by_qsize(const Carrying_Capacity& a, const Carrying_Capacity& b){
    return a.bin_quiver_size > b.bin_quiver_size;
}
void Map_Elites::arrange_CC(){
    sort(CC_vec.begin(), CC_vec.end(), sort_by_qsize);
}
// --------------------------------------------------
            // Calc Current Carrying Capacity
    /// calculate amount currently being carried
void Map_Elites::calc_current_carrying_capacity(){
    int cc_size1=0;
    int it=0;
//    for( auto i = CC_vec.begin(); i != CC_vec.end(); i++){
//        cc_size1=+CC_vec.at(it).bin_quiver_size;
//        it++;
//    }
//    current_carrying_capacity_size=cc_size1;
    cout << endl << "CC_vec.size() is " << CC_vec.size() << endl;
    for( int i=0; i < CC_vec.size(); i++){
        cout << endl << "CC.vec at " << it << ", has quiver_size = " << CC_vec.at(it).bin_quiver_size << " | and ID of " << CC_vec.at(it).bin;
        cc_size1+=CC_vec.at(it).bin_quiver_size;
        it++;
    }
    current_carrying_capacity_size=cc_size1;
    
    cout << endl << "Calc'd current CC is: " << current_carrying_capacity_size << endl;
}
// --------------------------------------------------
                // Calc Carrying bins
    /// calc # of bins that are in CC_vec
void Map_Elites::calc_carrying_bins(){
    int it=0;
    for( auto i = CC_vec.begin(); i != CC_vec.end(); i++){
        it++;
    }
    carrying_bins=it;
}
// --------------------------------------------------
            // Find ID 0
    /// finds map_space id of first element in CC_vec
void Map_Elites::find_id_0(){   // id of first element in CC_Vec will be biggest map_space.
//    cout << endl << "IN - find_id_0: "  << endl;
    id_0=CC_vec.at(0).bin;      // ID of most full bin.
    id_0_row=0;
    id_0_ele=0;
    int id_it=0;
    int id_it_row=0;
    int id_it_ele=0;
    for(int d1=0; d1<num_spacing1; d1++){
        int id_it_ele=0;
        for(int d2=0; d2<num_spacing2; d2++){
            if(id_0==id_it){
                break;
            }
            else{
                id_it++;
                id_it_ele++;
            }
        }
        if(id_0==id_it){
            break;
        }
        else{
            id_it_row++;
        }
    }
    id_0_row=id_it_row;
    id_0_ele=id_it_ele;
    //cout << endl << "id_0_row is: " << id_0_row << endl;
    //cout << endl << "id_0_ele is: " << id_0_ele << endl;
//    cout << endl << "OUT - find_id_0: "  << endl;
}
// --------------------------------------------------
            // Find ID
    /// finds row and element value based on the ID
void Map_Elites::find_id(int passed_id){
    id_row=0;
    id_ele=0;
    int id_it=0;
    int id_it_row=0;
    int id_it_ele=0;
    for(int d1=0; d1<num_spacing1; d1++){
        int id_it_ele=0;
        for(int d2=0; d2<num_spacing2; d2++){
            if(passed_id==id_it){
                break;
            }
            else{
                id_it++;
                id_it_ele++;
            }
        }
        if(passed_id==id_it){
            break;
        }
        id_it_row++;
    }
    id_row=id_it_row;
    id_ele=id_it_ele;
}
// --------------------------------------------------
            // Find ID end
    /// Finds the id or bin number for the last element of CC_vec, smallest.
void Map_Elites::find_id_end(){
    id_end=CC_vec.end()->bin;   //  TODO- Check to see if works
    id_end_row=0;
    id_end_ele=0;
    int id_it=0;
    int id_it_row=0;
    int id_it_ele=0;
    for(int d1=0; d1<num_spacing1; d1++){
        int id_it_ele=0;
        for(int d2=0; d2<num_spacing2; d2++){
            if(id_end==id_it){
                break;
            }
            else{
                id_it++;
                id_it_ele++;
            }
        }
        if(id_end==id_it){
            break;
        }
        id_it_row++;
    }
    id_end_row=id_it_row;
    id_end_ele=id_it_ele;
}
// --------------------------------------------------
            // Find Placement
    /// find placement of id in CC_vec
void Map_Elites::find_placement(int passed_id){
    int id_it=0;
//    for(auto i = CC_vec.begin(); i != CC_vec.end(); i++){
//        if(passed_id==CC_vec.at(id_it).bin){
//            id_element=id_it;
//        }
//        else{
//            id_it++;
//        }
//    }
    for(int i = 0; i < CC_vec.size(); i++){
        if(passed_id==CC_vec.at(id_it).bin){
            id_element=id_it;
            break;
        }
        else{
            id_it++;
        }
    }
    cout << endl << "id_element is: " << id_element << endl;
}
// --------------------------------------------------
            // Check CC
    /// check to see if id is already in quiver
bool Map_Elites::check_CC(int passed_id){
    int id_it=0;
    int in_CC=0;
    for(auto i = CC_vec.begin(); i != CC_vec.end(); i++){
        if(passed_id==CC_vec.at(id_it).bin){
            in_CC=1;
            return 0;
        }
        else{
            id_it++;
        }
    }
        return 1;
}
// --------------------------------------------------
            // Dynamic CC
    /// Contorls the carrying capacity for the map_spaces.
void Map_Elites::dynamic_CC(int rval, int eval){      // pass bin id
    cout << endl << "--------------- DYNAMIC IN ---------------";
    
    if(check_CC(Map.at(rval).at(eval).id)==1){  // TODO - construct cc vec
        add_CC(Map.at(rval).at(eval).id);       // TODO- change eval to e_val // correlate
    }
    else{
        cout << endl << "Bin is in CC_vec" << endl;
    }
    
    calc_current_carrying_capacity();  
    
    cout << endl << "CHOOSEN Bin ID is " << Map.at(rval).at(eval).id << endl;
    cout << "ZERO    Bin ID is " << CC_vec.at(0).bin << endl;
    
    if (current_carrying_capacity_size < carrying_capacity_size){  /// if carrying capacity IS NOT met.
        Map.at(rval).at(eval).increase_quiver_size();     // increase quiver size - Map_space
        //cout << endl << "Increase quiver size of Map at " << Map.at(rval).at(eval).id << " to " << Map.at(rval).at(eval).get_quiver_size() << endl;
        find_placement(Map.at(rval).at(eval).id);         // find corresponding CC
        CC_vec.at(id_element).increase_quiver_size();         // increase quiver size - CC
        //cout << endl << "Correspinding CC_vec is at " << CC_vec.at(id_element).bin << " increased to " << CC_vec.at(id_element).bin_quiver_size << endl;
        
        
        cout << endl << "CC  - + : size " << CC_vec.at(id_element).bin_quiver_size;
        cout << endl << "Map - + : size " << Map.at(rval).at(eval).quiver_size << endl;
        arrange_CC();
    }
    else if (Map.at(rval).at(eval).id != CC_vec.at(0).bin){     /// if Carrying capactiy IS met, and NOT in most full Map_Space.
        if(Map.at(rval).at(eval).quiver_size==CC_vec.at(0).bin-1 || Map.at(rval).at(eval).quiver_size>CC_vec.at(0).bin-1){
            cout << endl << "NOT CHANGING. TOO CLOSE. I DON'T WANT TOO, NOOOOOOO!" << endl;
        }
        else{
            find_id_0();
            cout << endl << "BEFORE - CC  - at 0 : element is " << CC_vec.at(0).bin << ", with size " << CC_vec.at(0).bin_quiver_size;
            cout << endl << "BEFORE - Map - at 0 : element is " << Map.at(id_0_row).at(id_0_ele).id << ", with size " << Map.at(id_0_row).at(id_0_ele).quiver_size;

            Map.at(rval).at(eval).increase_quiver_size();                   // increase quiver size of new - Map_space
            CC_vec.at(0).decrease_quiver_size();                                // decrease quiver size of 0, most full - CC
            find_id_0();
//            cout << endl << "ID of Most full is: " << Map.at(id_0_row).at(id_0_ele).id << endl;
            Map.at(id_0_row).at(id_0_ele).decrease_quiver_size();           // decrease quiver size of 0, and delete last element in bin_quiver- Map_space
            find_placement(Map.at(rval).at(eval).id);                       // find corresponding CC
            CC_vec.at(id_element).increase_quiver_size();                       // increase quiver size of new - CC
            arrange_CC();
        
            cout << endl << "AFTER - CC  - 0 : element is " << CC_vec.at(0).bin << ", with size " << CC_vec.at(0).bin_quiver_size;;
            cout << endl << "AFTER - Map - 0 : element is " << Map.at(id_0_row).at(id_0_ele).id << ", with size " << Map.at(id_0_row).at(id_0_ele).quiver_size;
            cout << endl << "AFTER - CC  - + : element is " << CC_vec.at(id_element).bin << ", with size " << CC_vec.at(id_element).bin_quiver_size;
            cout << endl << "AFTER - Map - + : element is " << Map.at(rval).at(eval).id << ", with size " << Map.at(rval).at(eval).quiver_size << endl;
        }
    }
    else{
        cout << endl << "Carrying Capacity of Bin is full" << endl; // Bin is 0, most full.
    }
    // TODO - Run Tests -
//    - 1) Total carrying capacities correspond. (quick operation)
//    - 2) CC and Map_space correspond, bin and quiver size for each. (long operation)
    
    cout << "--------------- DYNAMIC OUT ---------------" << endl;
}
// --------------------------------------------------
            // Place individual
    /// places into corresponding map_space in Map
void Map_Elites::place_individual_in_map(){
    double p1=challenger.get_phenotype1();
    double p2=challenger.get_phenotype2();
    
    /// out of bounds check
    double p1_lower=dim1_min;
    double p2_lower=dim2_min;
    double p1_upper=dim1_max;
    double p2_upper=dim2_max;
    
    double pLB1=dim1_min;
    double pLB2=dim2_min;
    int row_value=0;
    int element_value=0;
    
    /// to find which row phenotype belongs in
    for(int d1=0; d1<num_spacing1; d1++){
        double pUB1= pLB1+spacing1;
        /// under bounds
        if (p1<p1_lower){
            row_value=0;
            break;
        }
        else if (p1>=p1_upper){
            row_value=resolution1-1;
            break;
        }
        /// between range of LB1 and UB1
        else if (pLB1<p1 && p1<=pUB1){
            row_value=d1;
            break;
        }
        else {
            pLB1+=spacing1;
        }
    }
    /// to find which element in row phenotype beolongs in
    for(int d2=0; d2<num_spacing2; d2++){
        double pUB2= pLB2+spacing2;
        // under bounds
        if (p2<p2_lower){
            element_value=0;
            break;
        }
        else if (p2>=p2_upper){
            element_value=resolution2-1;
            break;
        }
        /// between range of LB2 and UB2
        else if (pLB2<p2 && p2<=pUB2){
            element_value=d2;
            break;
        }
        else {
            pLB2+=spacing2;
        }
    }
    
    Map.at(row_value).at(element_value).current_individual.push_back(this->challenger);
    
    //cout << endl << endl << endl << "Placed in row "<< row_value << " column " <<element_value << endl;
    
    // ---------- CCME Main Function ---------- Begin
    dynamic_CC(row_value, element_value);
        // TODO - Run Tests -
            // - 1) Total carrying capacities correspond. (quick operation)
            // - 2) CC and Map_space correspond, bin and quiver size for each. (long operation)
    
    //test_qsize_correspond(); // TODO - Fix, Runs forever.
    
    // ---------- CCME Main Function ---------- End
    
    Map.at(row_value).at(element_value).compare_new_individual(this->challenger);
        /// compare new individual in map space and erase worse
        /// Also, controls Carrying Capacity dynamics in Map_Space.
}
// --------------------------------------------------
            // Print CC Vec
void Map_Elites::print_CC_vec(){
    ofstream myfile;
    myfile.open ("CC_vec.txt");
    for(int i=0; i<CC_vec.size(); i++){
        myfile << CC_vec.at(i).bin << '\t';
        myfile << CC_vec.at(i).bin_quiver_size << '\n';
    }
    myfile.close();
    cout << "CC_vec.txt file created." << endl;
}
// --------------------------------------------------
/*
 CCME Test   // TODO
 
 - 1) Total carrying capacities correspond. (quick operation)
 - 2) CC and Map_space correspond, bin and quiver size for each. (long operation)
 
 */
// --------------------------------------------------
            // Arrange Test CC
bool sort_by_id(const Carrying_Capacity& a, const Carrying_Capacity& b){
    return a.bin < b.bin;
}
void Map_Elites::arrange_test_CC_vec(){
    // arrange CC_vec from first to last id
    sort(test_CC_vec.begin(), test_CC_vec.end(), sort_by_id);
}
// --------------------------------------------------
            // Test Quiver_size Correspond
    /// long test, create a test vector from CC_vec, test each element to determine if quiver_size of CC_vec and Map are the same.
void Map_Elites::test_qsize_correspond(){
    test_CC_vec.clear();
    int it=0;
    for(auto i = CC_vec.begin(); i != CC_vec.end(); i++){
        test_CC_vec.push_back(CC_vec.at(it));
        it++;
    }
    arrange_test_CC_vec();
    
    int it2=0;
    for(auto ii = test_CC_vec.begin(); ii != test_CC_vec.end(); ii++){
        find_id(test_CC_vec.at(it2).bin);   // TODO - flag for debug
        if(test_CC_vec.at(it2).bin_quiver_size != Map.at(id_row).at(id_ele).quiver_size){
            cout << endl << "ERROR 12784902: QUIVER_SIZE NOT ALLIGNED" << endl;
        }
        else{
            cout << endl << "QUIVER_SIZE ALLIGNED" << endl;
        }
        it2++;
    }
}
// --------------------------------------------------
            // Individual from Map
    /// gets an individual from a map_space in Map.
void Map_Elites::individual_from_map(int p1, int p2){
    int row_value=0;
    int element_value=0;
    row_value = p1;
    element_value = p2;
    current_bin1=p1;
    current_bin2=p2;
    if (Map.at(row_value).at(element_value).current_individual.size()>0){
        challenger = Map.at(row_value).at(element_value).current_individual.at(0);
        challenger.home_id = Map.at(row_value).at(element_value).id;
        
        Map.at(row_value).at(element_value).mutate_counter++;           // counter on bin pulled for mutation.
        Map.at(row_value).at(element_value).current_individual.at(0).mutate_counter_individual++;      // counter on Individual pulled for mutation.
    }
    else {
        //----------------------------------------------------------------     YL fix to check ----- VVV
        //while(Map.at(row_value).at(element_value).current_individual.size()<=0){
        
        //            int b1_min=0;
        //            int b2_min=0;
        //            int b1_max=get_resolution1();
        //            int b2_max=get_resolution2();
        //            row_value=rand()%(b1_max-b1_min)+b1_min;
        //            element_value=rand()%(b2_max-b2_min)+b2_min;
        //            //cout << "bin to mutate is: (" << bin1 << "," << bin2 << ")" << endl;
        //
        //        }
        //        challenger = Map.at(row_value).at(element_value).current_individual.at(0);
        //        challenger.home_id = Map.at(row_value).at(element_value).id;
        //----------------------------------------------------------------------
        
        //cout << endl << "Full_bin vector before has  " << full_bins.size() << endl;
        //cout << "Trying to Place in  (" << p1 << " , " << p2 << ")  Bin ID is: " << Map.at(row_value).at(element_value).id << endl;
        
        create_full_bin();
        //cout << endl << "Full_bin vector before has  " << full_bins.size() << endl;
        
        int min_index = min_element(distance_between_centerbin_phenotype.begin(), distance_between_centerbin_phenotype.end()) - distance_between_centerbin_phenotype.begin();
        
        //cout << "min distance is in bin: " << min_index << endl;
        
        challenger = full_bins.at(min_index).current_individual.at(0);
        challenger.home_id = full_bins.at(min_index).id;
        
        full_bins.at(min_index).mutate_counter++;               /// counter for amount bin is pulled for mutation.
        full_bins.at(min_index).current_individual.at(0).mutate_counter_individual++;       // counter for amount Individual pulled for mutation.
    }
}
// --------------------------------------------------
// grab individual for testing

// --------------------------------------------------
            // Create Full Bin
    /// creates a vector of all bins that contain an Individual.
void Map_Elites::create_full_bin(){
    full_bins.clear();
    distance_between_centerbin_phenotype.clear();
    for(int row_value=0; row_value<num_spacing1; row_value++){
        for(int element_value=0; element_value<num_spacing2; element_value++){
            if (Map.at(row_value).at(element_value).full_bin_check()==0){ // look for enum
                full_bins.push_back(Map.at(row_value).at(element_value));       // push_back Map_space that is full
                find_pheno_dist_to_center_bin(row_value, element_value, current_bin1, current_bin2);        // calculates distance from phenotype to center of desired bin
                distance_between_centerbin_phenotype.push_back(center_dist);
                //cout << "-----------------------Bin is: " << full_bins.size() << "  Center distance is: " << center_dist << endl;
                // "  Fit is: " << Map.at(row_value).at(element_value).current_individual.at(0).fit_rating << endl;
            }
        }
    }
    //cout << endl << "Full_bin vector has  " << full_bins.size() << endl;
    //cout << endl << "Distance vector has  " << distance_between_centerbin_phenotype.size() << endl;
}
// --------------------------------------------------
            // Find Phenotype distance to Center of Bin
    /// find distance from Individual in a bin to center of desired bin.
void Map_Elites::find_pheno_dist_to_center_bin(int p1, int p2, int cb1, int cb2){
    center_dist = 0;
    double d1 = Map.at(p1).at(p2).current_individual.at(0).get_phenotype1();
    double d2 = Map.at(p1).at(p2).current_individual.at(0).get_phenotype2();
    //cout << endl << "Pheno 1 is:  " << Map.at(p1).at(p2).current_individual.at(0).get_phenotype1() << " , Pheno 2 is:  " << Map.at(p1).at(p2).current_individual.at(0).get_phenotype2() << endl;
    
    double d1c = Map.at(cb1).at(cb2).center_bin1;
    double d2c = Map.at(cb1).at(cb2).center_bin2;
    //cout << "Center bin 1 is:  " << Map.at(cb1).at(cb2).center_bin1 << " , Center bin 2 is:  " << Map.at(cb1).at(cb2).center_bin2 << endl;
    double d1_sq = (d1c-d1)*(d1c-d1);
    double d2_sq= (d2c-d2)*(d2c-d2);
    center_dist = sqrt(d1_sq + d2_sq);
}
// --------------------------------------------------
            // How many bins are full?
void Map_Elites::how_many_full_bins(){
    int num_of_full_bins=0;
    for(int row_value=0; row_value<num_spacing1; row_value++){
        for(int element_value=0; element_value<num_spacing2; element_value++){
            if (Map.at(row_value).at(element_value).full_bin_check()==0){
                num_of_full_bins++;
            }
        }
    }
    cout << endl << "Final number of full bins: " << num_of_full_bins << endl;
}
// --------------------------------------------------
            // Print fit_ratings of Map
void Map_Elites::print_fit_ratings_of_map(){
    ofstream myfile;
    myfile.open ("fit_ratings_of_map.txt");
    for(int row_value=0; row_value<num_spacing1; row_value++){
        for(int element_value=0; element_value<num_spacing2; element_value++){
            if (Map.at(row_value).at(element_value).full_bin_check()==0){
                myfile << Map.at(row_value).at(element_value).current_individual.at(0).get_fit_rating() << '\n';
            }
        }
    }
    myfile.close();
    cout << endl << "fit_ratings_of_map.txt file created." << endl;
}// --------------------------------------------------
            // Get Best Individuals
vector<double>& Map_Elites::get_best_individual1(){
    return best_individual1;
}
vector<double>& Map_Elites::get_best_individual2(){
    return best_individual2;
}
// --------------------------------------------------
            // Best Fit Genome
void Map_Elites::best_fit_bin(){
    full_bins.clear();
    fit_ratings_in_map.clear();
    best_fit_index=0;
    for(int row_value=0; row_value<num_spacing1; row_value++){
        for(int element_value=0; element_value<num_spacing2; element_value++){
            if (Map.at(row_value).at(element_value).full_bin_check()==0){
                full_bins.push_back(Map.at(row_value).at(element_value));       // push_back Map_space that is full
                fit_ratings_in_map.push_back(Map.at(row_value).at(element_value).current_individual.at(0).get_fit_rating());
            }
        }
    }
    best_fit_index = max_element(fit_ratings_in_map.begin(), fit_ratings_in_map.end()) - fit_ratings_in_map.begin();
    
    best_individual1=  full_bins.at(best_fit_index).current_individual.at(0).get_individual1();
    //best_individual2=  full_bins.at(best_fit_index).current_individual.at(0).get_individual2();           // TODO - Remove Comments for main program--------------------
    
    //full_bins.at(best_fit_index).current_individual.at(0).display_individual1();
    
    cout << endl << "best fit bin is: " << best_fit_index << " with a fitness of " << full_bins.at(best_fit_index).current_individual.at(0).fit_rating << endl;
    //cout << endl << "best fit bin has been accessed: " << full_bins.at(best_fit_index).get_counter() << " times." << endl;
    //cout << endl << "best fit has " << full_bins.at(best_fit_index).previous_fit_rating.size() << " different occupants." << endl;
    //cout << endl << "best fit has deleted " << full_bins.at(best_fit_index).old_counter << " past occupants." << endl;
    //cout << endl << "best fit has deleted " << full_bins.at(best_fit_index).new_counter << " potential occupants." << endl;
}
// --------------------------------------------------
            // Print Best Bin Occupants fit_ratings
void Map_Elites::print_best_occupants_fitness(){
    ofstream myfile;
    myfile.open ("best_occupant_fit_ratings.txt");
    for(int occupant=0; occupant < full_bins.at(best_fit_index).previous_fit_rating.size(); occupant++){
        myfile << full_bins.at(best_fit_index).previous_fit_rating.at(occupant) << '\n';
    }
    myfile.close();
    cout << "best_occupant_fit_ratings.txt file created." << endl;
}
// --------------------------------------------------
            // Print All Occupants
void Map_Elites::print_all_occupants(){
    ofstream myfile;
    myfile.open ("all_occupants_with_fit_ratings.txt");
    for(int element=0; element<full_bins.size();element++){
        myfile << "occupant for bin: " << element << "\t accessed count: " << full_bins.at(element).counter << '\n';
        for(int occupant=0; occupant < full_bins.at(element).previous_fit_rating.size(); occupant++){
            myfile << full_bins.at(element).previous_fit_rating.at(occupant) << '\n';
        }
    }
    myfile.close();
    cout << "all_occupants_with_fit_ratings.txt file created." << endl;
}
// --------------------------------------------------
            // Print Best Parents fitness
void Map_Elites::print_best_parents_fitness(){
    //cout << endl << "Parent fitness vector has " << full_bins.at(best_fit_index).current_individual.at(0).parents_fitness.size() << endl;
    ofstream myfile;
    myfile.open ("best_parents_fitness.txt");
    for(int parent=0; parent < full_bins.at(best_fit_index).current_individual.at(0).parents_fitness.size(); parent++){
        myfile << full_bins.at(best_fit_index).current_individual.at(0).parents_fitness.at(parent) << '\n';
    }
    myfile.close();
    cout << "best_parents_fitness.txt file created." << endl;
}
//--------------------------------------------------
            // Print Best Parents ID
void Map_Elites::print_best_parents_id(){
    //cout << endl << "Parent id vector has " << full_bins.at(best_fit_index).current_individual.at(0).parents_id.size() << endl;
    ofstream myfile;
    myfile.open ("best_parents_id.txt");
    for(int parent=0; parent < full_bins.at(best_fit_index).current_individual.at(0).parents_id.size(); parent++){
        myfile << full_bins.at(best_fit_index).current_individual.at(0).parents_id.at(parent) << '\n';
    }
    myfile.close();
    cout << "best_parents_id.txt file created." << endl;
}
// --------------------------------------------------
void Map_Elites::print_best_full_trace(){
    ofstream myfile;
    myfile.open ("best_tracer.txt");
    for(int parent=0; parent < full_bins.at(best_fit_index).current_individual.at(0).parents_id.size(); parent++){
        myfile << full_bins.at(best_fit_index).current_individual.at(0).parents_id.at(parent) << '\t';
        myfile << full_bins.at(best_fit_index).current_individual.at(0).parents_fitness.at(parent) << '\n';
    }
    myfile.close();
    cout << "best_tracer.txt file created" << endl;
}
// --------------------------------------------------
            // Print Heat Map
void Map_Elites::print_heat_map(){
    ofstream myfile;
    myfile.open ("heat_map.txt");
    for(int element=0; element<full_bins.size();element++){
        myfile << full_bins.at(element).id << '\t';                             // ID of Bin
        myfile << full_bins.at(element).current_individual.at(0).get_phenotype1() << '\t';       // phenotype 1
        myfile << full_bins.at(element).current_individual.at(0).get_phenotype2() << '\t';       // phenotype 2
        myfile << full_bins.at(element).best_fit_rating  << '\t';               // Current fit rating in Bin
        myfile << full_bins.at(element).get_counter()  << '\t';                 // Times accessed
        myfile << full_bins.at(element).get_old_deleted_counter()  << '\t';     // Times new Individual is better
        myfile << full_bins.at(element).get_mutation_counter()  << '\t';        // Times bin has been selected for mutation
        myfile << full_bins.at(element).current_individual.at(0).mutate_counter_individual  << '\n';        // Times Individual has been selected for mutation
//        
//        for (int i = 0; i<full_bins.at(element).current_individual.at(0).genome1.size(); i++){
//        myfile << full_bins.at(element).current_individual.at(0).genome1.at(i) << '
//            
//            /// you can also get the current date and time to go in that file as well.
//        
//        
//        
//        myfile << full_bins.at(element).current_individual.at(0).display_individual2() << '\n';
    }
    myfile.close();
    cout << "heat_map.txt file created." << endl;
}
// --------------------------------------------------
            // Print Read Me
    /// For Heat Map and Corresponding genomes



// TODO



// --------------------------------------------------
            // Print corresponding genome1
void Map_Elites::print_corresponding_genome1(){
    ofstream myfile;
    myfile.open ("print_corresponding_genome1.txt");
    for(int element=0; element<full_bins.size();element++){
        for (int i = 0; i<full_bins.at(element).current_individual.at(0).genome1.size(); i++){
             myfile << full_bins.at(element).current_individual.at(0).genome1.at(i) << '\t';
        }
    myfile << '\n';
    }
    myfile.close();
    cout << "print_corresponding_genome1.txt created." << endl;
}
// --------------------------------------------------
            // Print corresponding genome2
void Map_Elites::print_corresponding_genome2(){
    ofstream myfile;
    myfile.open ("print_corresponding_genome2.txt");
    for(int element=0; element<full_bins.size();element++){
        //myfile << "Bin is " << full_bins.at(element).id  << '\n';
        for (int i = 0; i<full_bins.at(element).current_individual.at(0).genome2.size(); i++){
            myfile << full_bins.at(element).current_individual.at(0).genome2.at(i) << '\t';
        }
    myfile << '\n';
    }
    myfile.close();
    cout << "print_corresponding_genome2.txt created." << endl;
}
// --------------------------------------------------
            // Print corresponding Bins
void Map_Elites::print_corresponding_bins(){
    ofstream myfile;
    myfile.open ("print_corresponding_bins.txt");
    for(int element=0; element<full_bins.size();element++){
        myfile << full_bins.at(element).id  << '\n';
    }
    myfile.close();
    cout << "print_corresponding_bins.txt created." << endl;
}
// --------------------------------------------------

// Initialze Map
/// Builds Map based on old genomes, ALSO creates and sets Map_space LB and UB

// Bring in parameters from old run?? read in from txt file?

// ----------------------Map Parameters MUST BE CONSISTENT------------------
void Map_Elites::build_map_from_old(vector<vector<double>> old_genome1_map, vector<vector<double>> old_genome2_map, vector<vector<int>> bin_id){
    double pre_LB1=dim1_min;
    double pre_LB2=dim2_min;
    int id_number=0;
    int id_equals=0;
    Map.reserve(num_spacing1);
    for(int d1=0; d1<num_spacing1; d1++){
        vector<Map_space> Row;
        Row.reserve(num_spacing2);
        for(int d2=0; d2<num_spacing2; d2++){
            Map_space M;
            M.previous_genome1.clear();  // clears occupant vector
            M.previous_genome2.clear();  // clears occupant vector
            M.set_id(id_number);
            M.bin1=d1;
            M.bin2=d2;
            M.set_LB1(pre_LB1);
            M.set_LB2(pre_LB2);
            double cbin1=pre_LB1+(spacing1/2);
            double cbin2=pre_LB2+(spacing2/2);
            M.set_center_bin1(cbin1);            // Sets center of Bin
            M.set_center_bin2(cbin2);
            double calc_UB1= pre_LB1+spacing1;
            M.set_UB1(calc_UB1);
            double calc_UB2= pre_LB2+spacing2;
            M.set_UB2(calc_UB2);
            M.build_map_space();
            
            if(id_number == bin_id.at(id_equals).at(0)){
                Individual I;
                I.build_individual_1_from_another(old_genome1_map.at(id_number));
                I.build_individual_2_from_another(old_genome2_map.at(id_number));
                /// Place Individual in center of Bin, Phenotypes are not transfered.
                I.set_phenotypes(cbin1,cbin2);
                M.current_individual.push_back(I);
                id_equals++;
            }
            
            Row.push_back(M);
            pre_LB2+=spacing2;
            id_number++;
        }
        Map.push_back(Row);
        pre_LB1+=spacing1;
        pre_LB2=dim2_min;
    }
    cout << endl << "Bins transfered: " << id_equals << endl;
    cout << endl << "Map is made" << endl;
}
// --------------------------------------------------






//void Map_Elites::load_genome1(){
//    ifstream co("print_corresponding_genome1.txt");
//    
//    double read;
//    vector<double> apush;
//    
//    while(co >> read){
//        apush.push_back(read);
//        
//        if(apush.size()>=W.size_of_genome1){
//            Map_of_genome1.push_back(apush);
//            apush.clear();
//        }
//    }
//}
//// --------------------------------------------------
//void Map_Elites::load_genome2(){
//    ifstream co("print_corresponding_genome2.txt");
//    
//    double read;
//    vector<double> apush;
//    
//    while(co >> read){
//        apush.push_back(read);
//        
//        if(apush.size()>=W.size_of_genome2){
//            Map_of_genome2.push_back(apush);
//            apush.clear();
//        }
//    }
//}
// --------------------------------------------------
// read
// write
// run













#endif /* Map_Elites_hpp */
