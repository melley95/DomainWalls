#ifndef FLAT_HPP_
#define FLAT_HPP_

#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"



// Class which computes flat (Minkoswki) initial conditions

class Flat
{

    template <class data_t>
    using Vars = ADMConformalVars::VarsWithGauge<data_t>;
    
    
    public:

    template <class data_t> void compute(Cell<data_t> current_cell) const
{
        
            Vars<data_t> vars;

         
            VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below

            vars.chi = 1.0;
            vars.lapse = 1.0;

   

            FOR(i, j)
        {
            if (i==j)
            {
            
                vars.h[i][j] = 1.0;

            }
        
        }

        current_cell.store_vars(vars);


    }
    
       
    



   



};



    



#endif /* FLAT_HPP_ */