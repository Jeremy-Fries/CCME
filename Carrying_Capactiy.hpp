//
//  Carrying_Capactiy.hpp
//  Map_Elites_Testing_v1
//
//  Created by Jeremy Fries on 7/8/16.
//  Copyright © 2016 Jeremy Fries. All rights reserved.
//

class Carrying_Capacity;

#ifndef Carrying_Capactiy_hpp
#define Carrying_Capactiy_hpp

#include <stdio.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <ostream>

using namespace std;

class Carrying_Capacity{
    friend class Map_Elites;
protected:
    
public:
    int bin;
    int bin_quiver_size=0;
    void decrease_quiver_size();
    void increase_quiver_size();
    
    
    
// --------------------------------------------------
    
    
    
private:
    
// --------------------------------------------------


};
// ------------------------------------------------------------------------------------------------ ^^ Declarations
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ------------------------------------------------------------------------------------------------ vv Definitions

void Carrying_Capacity::decrease_quiver_size(){
    --bin_quiver_size;
}
void Carrying_Capacity::increase_quiver_size(){
    ++bin_quiver_size;
}

#endif /* Carrying_Capactiy_hpp */
