//
//  Map_space.hpp
//  Mapping_Elites v1
//
//  Created by Jeremy Fries on 10/20/15.
//  Copyright © 2015 Jeremy Fries. All rights reserved.
//

class Map_space;

#ifndef Map_space_hpp
#define Map_space_hpp

#include <stdio.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <ostream>

#include "Individual.hpp"

using namespace std;

class Map_space{
    friend class Wrapper;
    friend class Map_Elites;
    friend class Individual;
protected:
    int full_bin_check();
    int bin1, bin2;
    int mutate_counter=0;
    
public:
// --------------------------------------------------
    void display_space();
    // Lower Bound 1
    void set_LB1(double);
    double get_LB1();
    void display_LB1();
    // Lower Bound 2
    void set_LB2(double);
    double get_LB2();
    void display_LB2();
    // Lower Bound 3
    void set_LB3(double);
    double get_LB3();
    void display_LB3();
    // Upper Bound 1
    void set_UB1(double);
    double get_UB1();
    void display_UB1();
    // Upper Bound 2
    void set_UB2(double);
    double get_UB2();
    void display_UB2();
    // Upper Bound 3
    void set_UB3(double);
    double get_UB3();
    void display_UB3();
// --------------------------------------------------
    // Center of Bin
    void set_center_bin1(double c1);
    double get_center_bin1();
    void display_center_bin1();
    void set_center_bin2(double c2);
    double get_center_bin2();
    void display_center_bin2();
// --------------------------------------------------
    void build_map_space();
// --------------------------------------------------
    void arrange_quiver();
    void set_quiver_size(int);
    void increase_quiver_size();
    void decrease_quiver_size();
// --------------------------------------------------
    void set_best_fit();
    double get_best_fit();
    void compare_new_individual(Individual challenger);
    int get_counter();
    void display_counter();
    int get_new_deleted_counter();
    void display_new_deleted_counter();
    int get_old_deleted_counter();
    void display_old_deleted_counter();
    int get_mutation_counter();
    void display_mutation_counter();
// --------------------------------------------------
    void set_stat_run(int);
    int get_stat_run();
// --------------------------------------------------
    void set_id(int);
    int get_id();
// --------------------------------------------------
private:
    int counter=0;
    int old_counter=0;
    int new_counter=0;
    
    int stat_run;
    int id;
    
    double LB1, UB1, LB2, UB2, LB3, UB3;
    double center_bin1, center_bin2;
    double best_fit_rating;
    vector<double> LB_vec;
    vector<double> UB_vec;
    vector<Individual> current_individual;
    vector< vector <double> > previous_genome1,previous_genome2;
    vector <double> previous_fit_rating;
    
    vector<Individual> bin_quiver;
    int quiver_size=0;
    double low_quiver_fitness;
    double high_quiver_fitness;
};
// ------------------------------------------------------------------------------------------------ ^^ Declarations
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ------------------------------------------------------------------------------------------------ vv Definitions
// Display bounds of Map_space
void Map_space::display_space(){
    display_LB1();
    display_UB1();
    display_LB2();
    display_UB2();
    display_LB3();
    display_UB3();
}
// Lower Bound 1
void Map_space::set_LB1(double L1){
    LB1 = L1;
}
double Map_space::get_LB1(){
    return LB1;
}
void Map_space::display_LB1(){
    cout << endl << "LB1 is: " << LB1 << endl;
}
// Lower Bound 2
void Map_space::set_LB2(double L2){
    LB2 = L2;
}
double Map_space::get_LB2(){
    return LB2;
}
void Map_space::display_LB2(){
    cout << endl << "LB2 is: " << LB2 << endl;
}
// Lower Bound 3
void Map_space::set_LB3(double L3){
    LB3 = L3;
}
double Map_space::get_LB3(){
    return LB3;
}
void Map_space::display_LB3(){
    cout << endl << "LB3 is: " << LB3 << endl;
}
// Upper Bound 1
void Map_space::set_UB1(double U1){
    UB1 = U1;
}
double Map_space::get_UB1(){
    return UB1;
}
void Map_space::display_UB1(){
    cout << endl << "UB1 is: " << UB1 << endl;
}
// Upper Bound 2
void Map_space::set_UB2(double U2){
    UB2 = U2;
}
double Map_space::get_UB2(){
    return UB2;
}
void Map_space::display_UB2(){
    cout << endl << "UB2 is: " << UB2 << endl;
}
// Upper Bound 3
void Map_space::set_UB3(double U3){
    UB3 = U3;
}
double Map_space::get_UB3(){
    return UB3;
}
void Map_space::display_UB3(){
    cout << endl << "UB3 is: " << UB3 << endl;
}
// --------------------------------------------------
// Center of Bin
void Map_space::set_center_bin1(double c1){
    center_bin1 = c1;
}
double Map_space::get_center_bin1(){
    return center_bin1;
}
void Map_space::display_center_bin1(){
    cout << endl << "center bin 1 is: " << center_bin1 << endl;
}
void Map_space::set_center_bin2(double c2){
    center_bin2 = c2;
}
double Map_space::get_center_bin2(){
    return center_bin2;
}
void Map_space::display_center_bin2(){
    cout << endl << "center bin 2 is: " << center_bin2 << endl;
}
// --------------------------------------------------
// Build Map Space
// Builds lower bound and uppper bound vector of Map_space, used to deterine where Map_space lies in Map.
void Map_space::build_map_space(){
    LB_vec.push_back(LB1);
    LB_vec.push_back(LB2);
    //LB_vec.push_back(LB3);
    
    UB_vec.push_back(UB1);
    UB_vec.push_back(UB2);
    //UB_vec.push_back(UB3);
}
// --------------------------------------------------
// If later a Map_space can hold many Individuals, best fitness rating map_space will be first element. Currently Map_space can hold only one Individual, so that will also lie in the first element. This function may not be used in the case of only one Individual.
void Map_space::set_best_fit(){
    best_fit_rating=current_individual.at(0).get_fit_rating();
}
double Map_space::get_best_fit(){
    return best_fit_rating;
}
// --------------------------------------------------
void Map_space::set_quiver_size(int q){
    quiver_size=q;
}
// --------------------------------------------------
void Map_space::increase_quiver_size(){
    ++quiver_size;
}
// --------------------------------------------------
void Map_space::decrease_quiver_size(){
    --quiver_size;
    // delete last element in bin_quiver here?
    bin_quiver.erase(bin_quiver.end());
}
// --------------------------------------------------
bool sort_by_fit(const Individual& a, const Individual& b){
    return a.fit_rating > b.fit_rating;
}
// --------------------------------------------------
void Map_space::arrange_quiver(){
// arrange quiver from high to low fitness
    sort(bin_quiver.begin(), bin_quiver.end(),sort_by_fit);
    
    cout << endl << "quiver_size is: " << quiver_size << endl;
    
// set low and high values quiver values
    high_quiver_fitness=bin_quiver.at(0).fit_rating;
    //low_quiver_fitness=bin_quiver.at(quiver_size).fit_rating;   // ERROR!!!! In bin_quiver both elements are the same... copied in arrange
    low_quiver_fitness=bin_quiver.end()->fit_rating;
}
// --------------------------------------------------
// Compare New Individual
// If an Individual is placed into an already occupied Map_space, currently worse fitness rated Individual is deleted.
// Can be used to sort Individuals if more than one is allowed.
// LATER - can have multiple individuals in vector and sorts based on best fit_rating, possible diversity in Map_space?

// name of function - bin_quiver - vector of individuals that are in a map_space
// control of how many individuals are in quiver.
// Q - need to add first instance of quiver in place_individual function.
// Q - do only the new_individuals need to be added to quiver?
// --------- only add functionality, don't restructure entire code ------------
// I - in addition to quiver function, have low and high int's that correspond to individuals fitness, don't need to access and count every individual in vector.
// --carrying capacity of each bin needs to be dynamic, each map_spaces bin-quiver will size based on how many total individuals can be held by the map.
// carrying capacity number set in: set_map_params. implemented in wrapper
//
// 1) each bin will have a int # correlated with how many individuals it has in bin, (map_space function)
//      - (possibly, give preference to areas of map which system designer deems interesting, allocating more memory to certain bins)
// 2) vector will update at every compare_new_individual(), vector will be ordered greatest to least in number of individuals of the bins. (map_space function)
// 3) when max capacity is achieved, new individuals that will be placed will take the "memory" from the worst fit individual in the most occupied bin. (map_elites function)
//      - accessed as 1st element of combined bin vector, last element of bin_quiver
//      - dynamic allocation of memory will go as follows.
//          - bin accessed will take in individual, as long as memory is not full, individual will be stored into bin.
//          - if memory is full, individual will only be placed if fitness is greater than a current individual in the bin.
// 4)

// ***** Need to have function in Map_Elites that will "oversee" and set the capacity for the quiver_size.
// *** May want to remove quiver function from being directly encoded into compare_new_individual().



void Map_space::compare_new_individual(Individual challenger){
    
    //cout << "COMPARE IN" << endl;
    
    // TODO - If statment for yes or no quiver.
    if(current_individual.size()>1){
        counter++;
        // compare individual fit_ratings and delete worst
        double f1,f2;
        f1=current_individual.at(0).get_fit_rating();
        f2=current_individual.at(1).get_fit_rating();
        
        if (f1>=f2){             //new Individual erased
            new_counter++;
    // CCME function---------------------------------------------------------------- Start
            // check if bin_quiver is full
            if (bin_quiver.size()<quiver_size){
                // if not full, add individual to quiver
                bin_quiver.push_back(current_individual.at(1));
                arrange_quiver(); // rearrange high to low fitness
            }
            else{
                // if full, check to see if new individual has a better fitness, if so replace worst in quiver and rearrange high to low in vector.
                if (current_individual.at(1).fit_rating>low_quiver_fitness) {
                    bin_quiver.erase(bin_quiver.end());
                    bin_quiver.push_back(current_individual.at(1));
                    arrange_quiver(); // rearrange high to low fitness
                }
                else{
                    cout << endl << "individual less than low quiver" << endl;
                }
            }

            // quiver_addition() - will take in new individual, compare agaist current members and take in if condintions are met. // may not use..
            // compare new individual to members of quiver - compare_bin_quiver
//            bin_quiver.push_back()
            
            // CCME function-------------------------- End
            current_individual.erase(current_individual.begin()+1);
            //cout << endl << "new Individual erased" << endl;
        }
        else if (f1<f2){        //old Individual erased
            old_counter++;
            // ------------------ place previous occupant into quiver
            if (bin_quiver.size()<quiver_size){
                // if not full, add individual to quiver
                bin_quiver.push_back(current_individual.at(0));
                arrange_quiver(); // rearrange high to low fitness
            }
            else{
                // if full, check to see if new individual has a better fitness, if so replace worst in quiver and rearrange high to low in vector.
                if (current_individual.at(0).fit_rating>low_quiver_fitness) {
                    bin_quiver.erase(bin_quiver.end());
                    bin_quiver.push_back(current_individual.at(0));
                    arrange_quiver(); // rearrange high to low fitness
                }
                else{
                    cout << endl << "error 265897635: individual less than low quiver" << endl; // This should never display.
                }
            }
            
            // ------------------ VVVVVVVV --------- After quiver placment do -------- VVVVVVVV
            // push_back previous occupant of bin to vector
            previous_genome1.push_back(current_individual.at(0).get_individual1());
            previous_genome2.push_back(current_individual.at(0).get_individual2());
            previous_fit_rating.push_back(current_individual.at(0).fit_rating);
            /// record lineage.
            current_individual.at(1).parents_fitness.push_back(challenger.fit_rating);
            current_individual.at(1).parents_id.push_back(challenger.home_id);
            current_individual.erase(current_individual.begin());
            //cout << endl << "old Individual erased" << endl;
            // ------------------ ^^^^^^^^ --------- After quiver placment do -------- ^^^^^^^^
            
        }
        else {
            cout << endl << "error 567321: compare new individual error.\t" << current_individual.size() << endl;
        }
    }       // LYJF
    else{
       // add to bin_quiver, quiver might be down 1.
        bin_quiver.push_back(current_individual.at(0));     // JFJF - check to see if better way.
        current_individual.at(0).parents_fitness.push_back(challenger.fit_rating);
        current_individual.at(0).parents_id.push_back(challenger.home_id);
    }
    set_best_fit();
}
// --------------------------------------------------
// Counter
// Counts how many times a Map_space has been accessed to place an Individual.
// TODO - 2 seperate counters, one that counts how many times the one being placed is better, possible tuning of mutation.
int Map_space::get_counter(){
    return counter;
}
void Map_space::display_counter(){
    cout << endl << "Number of times bin has been accessed: " << counter << endl;
}
int Map_space::get_new_deleted_counter(){
    return new_counter;
}
void Map_space::display_new_deleted_counter(){
    cout << endl << "Number of times bin deleted new counter: " << new_counter << endl;
}
int Map_space::get_old_deleted_counter(){
    return old_counter;
}
void Map_space::display_old_deleted_counter(){
    cout << endl << "Number of times bin has deleted old counter: " << old_counter << endl;
}
int Map_space::get_mutation_counter(){      // amount of times bin has been choosen for mutation
    return mutate_counter;
}
void Map_space::display_mutation_counter(){
    cout << endl << "Number of times bin has been mutated: " << mutate_counter << endl;
}
// --------------------------------------------------
// Stat_run
void Map_space::set_stat_run(int s){
    stat_run=s;
}
int Map_space::get_stat_run(){
    return stat_run;
}
// --------------------------------------------------
// ID number, associated with bin location
void Map_space::set_id(int i){
    id=i;
}
int Map_space::get_id(){
    return id;
}
// --------------------------------------------------
// Full Bin Check
int Map_space::full_bin_check(){
    if (current_individual.empty()){
        return 1;
    }
    else {
        //cout << endl << "Full bin check successful!" << endl;
        return 0;
    }
}
// --------------------------------------------------



// --------------------------------------------------


#endif /* Map_space_hpp */
