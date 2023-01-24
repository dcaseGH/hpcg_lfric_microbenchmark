// read dinodump.dat and save some stuff
// recreating microbenchmark

#include "GenerateGeometry.hpp"
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>
using namespace std;

void read_part_vector_from_string(string line, double** local_array, int &array_counter, string space_delimiter){
    vector<string> words{};
    size_t pos = 0;
    while ((pos = line.find(space_delimiter)) != string::npos) {
        words.push_back(line.substr(0, pos));
           line.erase(0, pos + space_delimiter.length());
        }
//    std::cout << "last bit " << line << std::endl; exit(9);
    // add the last bit additionally
    words.push_back(line);
    for (const auto &str : words) {
       if(!str.empty()){
//              temp_stencil[array_counter] = stod(str);
              (*local_array)[array_counter] = stod(str);
          array_counter++;
       }
    }
}

void apply_helmholtz_operator_code(local_int_t&nlayers, double** map_w3, double** yvec, double** xvec,
                                   double** op1, double** op2, double** op3, double** op4, double** op5, double** op6,
                                   double** op7, double** op8, double** op9, double** ans, local_int_t& undf_w3)
{
  std::cout << "Inside the subroutine " << std::endl;
//  int temp_int = 0;
//  for (int i = 0; i < undf_w3; i++)
//  {
//      std::cout << (*y_vec)[i] << std::endl;
//  }
  //std::cout << "ttemp int " << temp_int << std::endl;
  std::cout << "map " << (*map_w3)[0] << " " << (*map_w3)[1] << " " << (*map_w3)[3456 - 1] << std::endl;
  std::cout << "map " << (*map_w3)[4] << " " << (*map_w3)[5] << " " << (*map_w3)[6] << std::endl;
  std::cout << "yvec " << (*yvec)[0] << " " << (*yvec)[1] << " " << (*yvec)[undf_w3 - 1] << std::endl;
  std::cout << "xvec " << (*xvec)[0] << " " << (*xvec)[1] << " " << (*xvec)[undf_w3 - 1] << std::endl;

}




void read_dinodump(local_int_t& loop0_start, local_int_t& loop0_stop, local_int_t& nlayers, local_int_t& undf_w3, local_int_t& x_vec_max_branch_length,
                   double** map_w3, double** yvec, double** xvec, double** op1, double** op2, double** op3, double** op4, double** op5, double** op6,
                   double** op7, double** op8, double** op9, double** ans)
{
  // read dinodump.dat
  // remember from fortran
  // temp objects made as vectors for higher ranks, then rearranged at bottom
  const int numbers_per_line = 6;
  const int wide_numbers_per_line = 3;
  bool constants_set = false;
  string line;
  int counter = 0;
  ifstream myfile ("dinodump.dat");
//  double tempvec[];
/*  call steggy%input_array(x_vec_stencil_size, 4, undf_w3 loop0stop)
  call steggy%input_array(x_vec_stencil_dofmap,&
       ndf_w3,x_vec_max_branch_length, 4, loop0)
  call steggy%input_array(map_w3,ndf_w3,loop0)
*/
  size_t pos = 0;
  const string space_delimiter = " ";
  double* temp_stencil;
  double* temp_dofmap;
  int array_counter = 0; // use to add to arrays
  int scalars_stop  = 4;
  int stencil_start;
  int stencil_stop;
  int dofmap_start;
  int dofmap_stop;
  int map_start;
  int map_stop;
  int yvec_start;
  int yvec_stop;
  int xvec_start;
  int xvec_stop;
  int op1_start;
  int op1_stop;
  int op2_start;
  int op2_stop;
  int op3_start;
  int op3_stop;
  int op4_start;
  int op4_stop;
  int op5_start;
  int op5_stop;
  int op6_start;
  int op6_stop;
  int op7_start;
  int op7_stop;
  int op8_start;
  int op8_stop;
  int op9_start;
  int op9_stop;
  int ans_start;
  int ans_stop;
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      if (!constants_set){
        if (counter == 0) loop0_start = (local_int_t) stoi(line);
        if (counter == 1) loop0_stop  = (local_int_t) stoi(line);
        if (counter == 2) nlayers     = (local_int_t) stoi(line);
        if (counter == 3) undf_w3     = (local_int_t) stoi(line);
        if (counter == 4) {
                x_vec_max_branch_length = (local_int_t) stoi(line);
                constants_set = true;
//              set bounds
                stencil_start = scalars_stop + 1;
                stencil_stop  = stencil_start + 4*loop0_stop/numbers_per_line -1;
                dofmap_start  = stencil_stop + 1;
                dofmap_stop   = dofmap_start + 4*loop0_stop*x_vec_max_branch_length/numbers_per_line -1;
                map_start     = dofmap_stop + 1;
                map_stop      = map_start + loop0_stop/numbers_per_line -1;
                yvec_start    = map_stop + 1;
                yvec_stop     = yvec_start + undf_w3/wide_numbers_per_line -1;
                xvec_start    = yvec_stop + 1;
                xvec_stop     = xvec_start + undf_w3/wide_numbers_per_line -1;
                op1_start     = xvec_stop + 1;
                op1_stop      = op1_start + undf_w3/wide_numbers_per_line -1;
                op2_start     = op1_stop + 1;
                op2_stop      = op2_start + undf_w3/wide_numbers_per_line -1;
                op3_start     = op2_stop + 1;
                op3_stop      = op3_start + undf_w3/wide_numbers_per_line -1;
                op4_start     = op3_stop + 1;
                op4_stop      = op4_start + undf_w3/wide_numbers_per_line -1;
                op5_start     = op4_stop + 1;
                op5_stop      = op5_start + undf_w3/wide_numbers_per_line -1;
                op6_start     = op5_stop + 1;
                op6_stop      = op6_start + undf_w3/wide_numbers_per_line -1;
                op7_start     = op6_stop + 1;
                op7_stop      = op7_start + undf_w3/wide_numbers_per_line -1;
                op8_start     = op7_stop + 1;
                op8_stop      = op8_start + undf_w3/wide_numbers_per_line -1;
                op9_start     = op8_stop + 1;
                op9_stop      = op9_start + undf_w3/wide_numbers_per_line -1;
                ans_start     = op9_stop + 1;
                ans_stop      = ans_start + undf_w3/wide_numbers_per_line -1;

                // allocate memory
                temp_stencil  = (double*) malloc(4*loop0_stop * sizeof(double));
                temp_dofmap   = (double*) malloc(4*loop0_stop*x_vec_max_branch_length * sizeof(double));
                *map_w3       = (double*) malloc(loop0_stop * sizeof(double));
                *yvec         = (double*) malloc(undf_w3 * sizeof(double));
                *xvec         = (double*) malloc(undf_w3 * sizeof(double));
                *op1          = (double*) malloc(undf_w3 * sizeof(double));
                *op2          = (double*) malloc(undf_w3 * sizeof(double));
                *op3          = (double*) malloc(undf_w3 * sizeof(double));
                *op4          = (double*) malloc(undf_w3 * sizeof(double));
                *op5          = (double*) malloc(undf_w3 * sizeof(double));
                *op6          = (double*) malloc(undf_w3 * sizeof(double));
                *op7          = (double*) malloc(undf_w3 * sizeof(double));
                *op8          = (double*) malloc(undf_w3 * sizeof(double));
                *op9          = (double*) malloc(undf_w3 * sizeof(double));
                *ans          = (double*) malloc(undf_w3 * sizeof(double));

                std::cout << "Scalars " << loop0_start << " "<< loop0_stop << " "<< nlayers << " "<< undf_w3 << " "<< x_vec_max_branch_length << std::endl;
                std::cout << "Points in file " << stencil_start << " " << stencil_stop << " " << dofmap_start << " " << dofmap_stop << " "
                << map_start << " " << map_stop << std::endl;
         }
      }
      else{
//     stencil
        if (counter >= stencil_start and counter <= stencil_stop)
        {
            vector<string> words{};
            while ((pos = line.find(space_delimiter)) != string::npos) {
              words.push_back(line.substr(0, pos));
              line.erase(0, pos + space_delimiter.length());
            }
            for (const auto &str : words) {
              if(!str.empty()){
                 temp_stencil[array_counter] = stod(str);
                 array_counter++;
                 }
              }
        }
//     dofmap       (ndf_w3,x_vec_max_branch_length, 4, loop0)
        if (counter == dofmap_start) array_counter = 0;
        if (counter >= dofmap_start and counter <= dofmap_stop)
        {
            vector<string> words{};
            while ((pos = line.find(space_delimiter)) != string::npos) {
              words.push_back(line.substr(0, pos));
              line.erase(0, pos + space_delimiter.length());
            }
            for (const auto &str : words) {
              if(!str.empty()){
                 temp_dofmap[array_counter] = stod(str);
                 array_counter++;
                 }
              }
        }//dofmap */
        if (counter == map_start) array_counter = 0;
        if (counter >= map_start and counter <= map_stop)
        {
            /*vector<string> words{};
            while ((pos = line.find(space_delimiter)) != string::npos) {
              words.push_back(line.substr(0, pos));
              line.erase(0, pos + space_delimiter.length());
            }
            for (const auto &str : words) {
              if(!str.empty()){
                 (*map_w3)[array_counter] = stod(str);
                 array_counter++;
                 }
              }*/
        read_part_vector_from_string(line, map_w3, array_counter, space_delimiter);
        }//map_w3
        //yvec
        if (counter == yvec_start) array_counter = 0;
        if (counter >= yvec_start and counter <= yvec_stop)
        {
           read_part_vector_from_string(line, yvec, array_counter, space_delimiter);
        }
        if (counter == xvec_start) array_counter = 0;
        if (counter >= xvec_start and counter <= xvec_stop)
        {
           read_part_vector_from_string(line, xvec, array_counter, space_delimiter);
        }
        if (counter == op1_start) array_counter = 0;
        if (counter >= op1_start and counter <= op1_stop)
        {
           read_part_vector_from_string(line, op1, array_counter, space_delimiter);
        }
        if (counter == op2_start) array_counter = 0;
        if (counter >= op2_start and counter <= op2_stop)
        {
           read_part_vector_from_string(line, op2, array_counter, space_delimiter);
        }
        if (counter == op3_start) array_counter = 0;
        if (counter >= op3_start and counter <= op3_stop)
        {
           read_part_vector_from_string(line, op3, array_counter, space_delimiter);
        }
        if (counter == op4_start) array_counter = 0;
        if (counter >= op4_start and counter <= op4_stop)
        {
           read_part_vector_from_string(line, op4, array_counter, space_delimiter);
        }
        if (counter == op5_start) array_counter = 0;
        if (counter >= op5_start and counter <= op5_stop)
        {
           read_part_vector_from_string(line, op5, array_counter, space_delimiter);
        }
        if (counter == op6_start) array_counter = 0;
        if (counter >= op6_start and counter <= op6_stop)
        {
           read_part_vector_from_string(line, op6, array_counter, space_delimiter);
        }
        if (counter == op7_start) array_counter = 0;
        if (counter >= op7_start and counter <= op7_stop)
        {
           read_part_vector_from_string(line, op7, array_counter, space_delimiter);
        }
        if (counter == op8_start) array_counter = 0;
        if (counter >= op8_start and counter <= op8_stop)
        {
           read_part_vector_from_string(line, op8, array_counter, space_delimiter);
        }
        if (counter == op9_start) array_counter = 0;
        if (counter >= op9_start and counter <= op9_stop)
        {
           read_part_vector_from_string(line, op9, array_counter, space_delimiter);
        }
        if (counter == ans_start) array_counter = 0;
        if (counter >= ans_start and counter <= ans_stop)
        {
           read_part_vector_from_string(line, ans, array_counter, space_delimiter);
        }

      }//constants set
      counter++;
    }//while
    if (counter >= 4){
//    for (int i = 0;i< 4*loop0_stop; i++)
//    {
//       std::cout << "temp ve " << temp_stencil[i] << "\n";
 //   }
    free(temp_stencil);
    free(temp_dofmap);
    }//if counter large
    myfile.close();
  }//file open
};
