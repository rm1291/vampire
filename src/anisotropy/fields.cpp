//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>
#include <sstream>

// Vampire headers
#include "anisotropy.hpp"
#include "errors.hpp"
#include "units.hpp"
#include "vio.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //---------------------------------------------------------------------------
   // Function to calculate magnetic fields from anisotropy tensors
   //
   // Hx = -dE/dSx = - 2      (  Txx Sx  +  Txy Sy  +  Txz Sz  )
   //                - 4 Sx   ( Txx Sx^2 + Txy Sy^2 + Txz Sz^2 )
   //                - 6 Sx^2 ( Txx Sx^3 + Txy Sy^3 + Txz Sz^3 )
   //
   // Hy = -dE/dSy = - 2      (  Tyx Sx  +  Tyy Sy  +  Tyz Sz  )
   //                - 4 Sy   ( Tyx Sx^2 + Tyy Sy^2 + Tyz Sz^2 )
   //                - 6 Sy^2 ( Tyx Sx^3 + Tyy Sy^3 + Tyz Sz^3 )
   //
   // Hz = -dE/dSz = - 2      (  Tzx Sx  +  Tzy Sy  +  Tzz Sz  )
   //                - 4 Sz   ( Tzx Sx^2 + Tzy Sy^2 + Tzz Sz^2 )
   //                - 6 Sz^2 ( Tzx Sx^3 + Tzy Sy^3 + Tzz Sz^3 )
   //
   // Note - tensors are defined with units of Tesla and a double calcellation
   //        of minus (-) signs, one from the Hamiltonian and the other from the
   //        derivative in the field calculation.
   //
   //---------------------------------------------------------------------------
   void fields(std::vector<double>& spin_array_x,
               std::vector<double>& spin_array_y,
               std::vector<double>& spin_array_z,
               std::vector<int>&    type_array,
               std::vector<double>& field_array_x,
               std::vector<double>& field_array_y,
               std::vector<double>& field_array_z,
               const int start_index,
               const int end_index,
               const double temperature){

      // Lop over all atoms between start and end index
      for(int atom = start_index; atom<end_index; atom++){

         const double sx = spin_array_x[atom]; // store spin direction in temporary variables
         const double sy = spin_array_y[atom];
         const double sz = spin_array_z[atom];


           // get atom material(Roberto)
            const unsigned int mat = type_array[atom];

            // Vectors defining the easy axis in cubic anisotropy(Roberto)          
            double e1[3] = { internal::mp[mat].kc_vector1[0],
                                   internal::mp[mat].kc_vector1[1],
                                   internal::mp[mat].kc_vector1[2] };

            double e2[3] = { internal::mp[mat].kc_vector2[0],
                                   internal::mp[mat].kc_vector2[1],
                                   internal::mp[mat].kc_vector2[2] };

            double e3[3] = { (internal::mp[mat].kc_vector1[1]*internal::mp[mat].kc_vector2[2] - internal::mp[mat].kc_vector1[2]*internal::mp[mat].kc_vector2[1]),
                                   (internal::mp[mat].kc_vector1[2]*internal::mp[mat].kc_vector2[0] - internal::mp[mat].kc_vector1[0]*internal::mp[mat].kc_vector2[2]),
                                   (internal::mp[mat].kc_vector1[0]*internal::mp[mat].kc_vector2[1] - internal::mp[mat].kc_vector1[1]*internal::mp[mat].kc_vector2[0]) };

             double mod_e1;           // Normalization (Roberto)
             double mod_e2;           
             double mod_e3; 



             // Strore constant including prefactor and conversion to Tesla (-dE/dS)  (Roberto)
            const double kc4 = 0.5 * internal::mp[mat].kc4;
            double lambda;
             
//           if ( kc4 >=  0.0) {           // (Roberto)
//                 lambda = 1.0;
//           }
//            else if ( kc4 < 0.0){       //  (Roberto)
//                 lambda = 1.0/3.0;
//            } 
 
           if ( kc4 >=  0.0) {           // (Roberto)
                 lambda = 0.0;
           }    
            else if ( kc4 < 0.0){       //  (Roberto)
                 lambda = 0.0/3.0;
            }



                // Normalise Spin Length (Roberto). Anisotropy directions are good enough being orthogonal
                mod_e1 = 1.0/sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
                mod_e2 = 1.0/sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
                mod_e3 = 1.0/sqrt(e3[0]*e3[0] + e3[1]*e3[1] + e3[2]*e3[2]);


                e1[0]=e1[0]*mod_e1;
                e1[1]=e1[1]*mod_e1;
                e1[2]=e1[2]*mod_e1;

                e2[0]=e2[0]*mod_e2;
                e2[1]=e2[1]*mod_e2;
                e2[2]=e2[2]*mod_e2;

                e3[0]=e3[0]*mod_e3;
                e3[1]=e3[1]*mod_e3;
                e3[2]=e3[2]*mod_e3;



         const unsigned int index = 9*atom; // get atom index in tensor array

         // Second order
         double hx = 2.0 * ( internal::second_order_tensor[index + 0] * sx +
                             internal::second_order_tensor[index + 1] * sy +
                             internal::second_order_tensor[index + 2] * sz );

         double hy = 2.0 * ( internal::second_order_tensor[index + 3] * sx +
                             internal::second_order_tensor[index + 4] * sy +
                             internal::second_order_tensor[index + 5] * sz );

         double hz = 2.0 * ( internal::second_order_tensor[index + 6] * sx +
                             internal::second_order_tensor[index + 7] * sy +
                             internal::second_order_tensor[index + 8] * sz );

         // New version choosing uniaxial anisotropy axis (Roberto)
                 
            // Changing the basis
            // sx' = (sx * e1[0] + sy * e1[1]  + sz * e1[2])
            // sy' = (sx * e2[0] + sy * e2[1]  + sz * e2[2])
            // sz' = (sx * e3[0] + sy * e3[1]  + sz * e3[2])
    
         hx += 4.0 * (  internal::fourth_order_tensor[index + 0] * (pow((sx * e1[0] + sy * e1[1]  + sz * e1[2]),3)*e1[0] + pow((sx * e2[0] + sy * e2[1]  + sz * e2[2]),3)*e2[0] + pow((sx * e3[0] + sy * e3[1]  + sz * e3[2]),3)*e3[0] ));
         hy += 4.0 * (  internal::fourth_order_tensor[index + 4] * (pow((sx * e1[0] + sy * e1[1]  + sz * e1[2]),3)*e1[1] + pow((sx * e2[0] + sy * e2[1]  + sz * e2[2]),3)*e2[1] + pow((sx * e3[0] + sy * e3[1]  + sz * e3[2]),3)*e3[1] ));
         hz += 4.0 * (  internal::fourth_order_tensor[index + 8] * (pow((sx * e1[0] + sy * e1[1]  + sz * e1[2]),3)*e1[2] + pow((sx * e2[0] + sy * e2[1]  + sz * e2[2]),3)*e2[2] + pow((sx * e3[0] + sy * e3[1]  + sz * e3[2]),3)*e3[2] ));

         // store net field
         field_array_x[atom] += hx;
         field_array_y[atom] += hy;
         field_array_z[atom] += hz;

      }

      // optionally calclulate lattice anisotropy fields
      if(internal::enable_lattice_anisotropy){
         internal::calculate_lattice_anisotropy_fields(spin_array_x,  spin_array_y,  spin_array_z, type_array,
                                                       field_array_x, field_array_y, field_array_z,
                                                       start_index, end_index, temperature);
      }

      return;

   }

} // end of anisotropy namespace
