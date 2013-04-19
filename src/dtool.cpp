/* ----------------------------------------------------------------------
   dtool - A DT Tool
                                 
   Written by Po-Jen Hsu, clusterga@gmail.com
   Copyright (c),2012-, Po-Jen Hsu. 

   This software is distributed under the GNU General Public License. 
-------------------------------------------------------------------------
   | File Name : dtool.cpp
   | Creation Time : 2013-01-09 16:27:54
   | Last Modified : 2013-04-02 10:38:02
------------------------------------------------------------------------- */
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#define PAUSE std::cout << "Press any key to continue...\n" ;std::cin.get();// fgetc(stdin);
// ========= global ===================
const int atom_num_max=800,name_max=101,atom_name_max=31,res_name_max=4,mat_dim=3; //Remember there is an additional \0 for string
const float pi=3.1415927;
int check_pdb_atom_num(char*);
int main(int argc, char * argv[]) {
// ========== general =================
   int atom_num,I0,I1,I2,I3,sel,help_sel,atom_num_dim,file_exist;
   int res_id[atom_num_max];
   float frame_num;
   char read_file_name[name_max],write_file_name[name_max];
   char flag[name_max],atom_name[atom_num_max][atom_name_max],res_name[atom_num_max][res_name_max];
   char dummy1[atom_name_max],dummy2[atom_name_max],dummy3[atom_name_max],dummy4[atom_name_max];
   char dummy5[atom_name_max],dummy6[atom_name_max],dummy7[atom_name_max];
// =========== ASA ====================
   int asa_num,group_sel;
   float asa_eta,read_asa_value[atom_num_max],r_ave;
   const float l_ave=4;
   int output_asa_num,output_asa_id[atom_num_max],output_asa_res_id[atom_num_max];
   float output_asa_value[atom_num_max],output_r_value[atom_num_max],output_xi_value[atom_num_max];
   char output_asa_atom_name[atom_num_max][atom_name_max],output_asa_res_name[atom_num_max][res_name_max];
// =========== weighting ==============
   int weight_num,weight_sel;
   char weight_name[atom_num_max][atom_name_max];
   float weight_value[atom_num_max], weight_numerator,weight_denominator; 
// =========== grouping ===============
   int read_group_target,read_group_element;
   char read_group_weight_name[atom_name_max];
   int group_target_num,previous_group_target_id;
   int group_target_id[atom_num_max],group_target_element_num[atom_num_max],group_target_element_id[atom_num_max][atom_num_max];
   float group_target_element_weight[atom_num_max][atom_num_max];
// =========== bead pdb ================
   int bead_sel,bead_num,bond_num,output_bond_num,output_bond_target_id,output_bond_end_id;
   int bond_target_res_id[atom_num_max];
   int bond_end_res_id[atom_num_max];
   char bond_target_atom_name[atom_num_max][atom_name_max],bond_target_res_name[atom_num_max][res_name_max];
   char bond_end_atom_name[atom_num_max][atom_name_max],bond_end_res_name[atom_num_max][res_name_max];
   float bead_mat_x,bead_mat_y,bead_mat_z,bond_dist;
   enum selection {zero,one,two,three,four,five,six,seven,eight,nine,ten};
   enum bool_enum {is_false,is_true};

   std::cout << "============= Welcome to dtool ==================\n\n";
   std::cout << "Welcome to DT Tool V1.1 by Po-Jen Hsu 2013\n\n";
   std::cout << "Bugs report: clusterga@gmail.com\n\n";
   std::cout << "This software is distributed under the\n";
   std::cout << "GNU General Public License.\n\n"; 
   std::cout << "============= Main Menu =========================\n\n";
   std::cout << "1. Convert or group a pdb trajectory to binary his format for ASA and DT (group and weight table required)\n";
   std::cout << "2. Convert or group ASA to friction-related parameters (ASA.out and group_table.dat required)\n";
   std::cout << "3. Convert or group a pdb trajectory from atomic to bead representation (bond_list.dat required)\n";
   std::cout << "4. How to use this program?\n";
   std::cout << "5. Quit\n";
   std::cin >> sel;

   if(sel == four){
     std::cout <<"============= Help Menu =========================\n\n";
     std::cout <<"dtool can convert parameters and trajectory for usage of ASA and DT.\n";
     std::cout <<"All the files are self-explained. There are mainly three manually type-in files.\n";
     std::cout <<"1. Usage of dtool\n";
     std::cout <<"2. Explain group_table.dat\n";
     std::cout <<"3. Explain weight_table.dat\n";
     std::cout <<"4. Explain bond_list.dat\n";
     std::cout <<"5. What are other required input and corresponding output files?\n";
     std::cout <<"6. Quit\n";
     std::cin >> help_sel;
     if(help_sel == one) {
       std::cout << "============= Usage of dtool ====================\n\n";
       std::cout << "First, convert pdb trajectory to binary his file without any weighting procedure (use: *.pdb1->1->*.his)\n";
       std::cout << "Second, calculate ASA from binary his file (*.his->ASA_ARNALDO->ASA.out)\n";
       std::cout << "Third, convert ASA to frictions (use: ASA.out->group_table.dat->2->ASA.dat)\n";
       std::cout << "Fourth, use ASA to construct the weighting and grouping table(weight_table.dat and group_table.dat)\n";
       std::cout << "Fifth, use weighting and grouping table to generate new binary his file for DT (use: 1-> 2 or 3->*.his)\n";
       std::cout << "Sixth, prepare bond_list.dat for DT. Output bonds and beads to a new pdb file (use: 3)\n";
       std::cout << "Seventh, prepare another ASA group table that match to the beads you defined. Calculate frictions again.(use: 2)\n\n\n";
       std::cout << "The contents of bond_list.dat and ASA.dat (frictions) of beads and bonds as well as the new grouped binary trajectory in his format will be imported to DT\n\n";
       std::cout << "Note!! There might be two kinds of group_table.dat.\n";
       std::cout << "One for generating the weighting factors for atomic grouping.\n";
       std::cout << "Another for generating the corresponding frictions to the beads in DT (final stage).\n\n";

     } else if(help_sel == two) {
       std::cout <<"========= The content of group_table.dat ========\n\n";
       std::cout <<"# target_id element_id weight_name (The first line must be preserved. See example)\n";
       std::cout <<"The elements will be grouped into a target.\n";
       std::cout <<"The weight_name determines the weighting factor for each element.\n";
       std::cout <<"The weight_name only affects grouping of pdb trajectory.\n";
       std::cout <<"No effect on ASA grouping so all weight_names are none in ASA_group_table.\n";
       std::cout <<"Example:\n\n";
       std::cout <<"# target_id element_id weight_name\n";
       std::cout <<"3  3  ACE1H1\n";
       std::cout <<"3  4  ACE1H2\n";
       std::cout <<"3  5  ACE1H3\n\n";
       std::cout <<"3,4,5 will be grouped together into 3 using the weighting factor H1, H2, and H3.\n";
       std::cout <<"The corresponding weighting factors are recorded in weight_table.dat file.\n";
       std::cout <<"Currently, dtool supports two weighting functions(square-root and non-square-root).\n";
       std::cout <<"New functions can be added in the future.\n";
       std::cout <<"Note: This file can be used for both pdb to his conversion or ASA grouping calculation.\n";
       std::cout <<"One may prepare more than one group_table.dat for different purposes.\n";
     } else if(help_sel == three) {
       std::cout <<"======== The content of weight_table.dat ========\n\n";
       std::cout <<"# weight_name weight_value  (The first line must be preserved. See example)\n";
       std::cout <<"It simply records the weighting name and its corresponding weighting factor.\n";
       std::cout <<"The weighting factors (values) can be obtained via ASA calculation.\n";
       std::cout <<"After calculating frictions from ASA, the program will generate weight_name and weight_value in ASA.dat.\n";
       std::cout <<"Currently, this file only works in pdb to his conversion. It has no effect on ASA grouping.\n";
       std::cout <<"Example:\n\n";
       std::cout <<"# weight_name weight_value\n";
       std::cout <<"ACE1H1  3.2423\n";
       std::cout <<"ACE1H2  4.0421\n";
       std::cout <<"ACE1H3  2.9401\n\n";
     } else if(help_sel == four) {
       std::cout <<"========= The content of bond_list.dat ==========\n\n";
       std::cout <<"# bead_a res_a res_a_id bead_b res_b res_b_id  (The first line must be preserved. See example)\n";
       std::cout <<"It determines how beads connect with each other in DT calculation.\n";
       std::cout <<"Each line represents a bond between bead_a and bead_b.\n";
       std::cout <<"The content can be directly imported to DT.\n";
       std::cout <<"One can also use it to generate the trajectory of beads with bond.\n";     
       std::cout <<"Example:\n\n";
       std::cout <<"# bead_a  res_a  res_a_id  bead_b  res_b  res_b_id\n";
       std::cout <<"H1  ACE  1   H3  ACE  1\n";
       std::cout <<"H3  ACE  1   H2  ACE  1\n\n";
       std::cout <<"Three beads are defined as follow: H1-H3-H2.\n";
     } else if(help_sel == five) {
       std::cout <<"=============== Other files =====================\n\n";
       std::cout <<"dtool can convert a trajectory from pdb format to his binary format for both ASA and DT.\n";
       std::cout <<"--> *.pdb trajectory required.\n";
       std::cout <<"--> Output binary trajectory *.his\n\n";
       std::cout <<"dtool can calculate frictions and weighting factors from ASA.\n";
       std::cout <<"--> ASA.out by ASA main program required.\n";
       std::cout <<"--> Output ASA.dat, friction-related parameters as well as the weighting parameters for group table.\n\n";
       std::cout <<"dtool can convert a trajectory from atomic representation to bead representation.\n\n";
       std::cout <<"--> *.pdb trajectory required.\n";
       std::cout <<"--> Output *_bead.pdb with positions of beads and the bond distances.\n\n";
     }
     return(0);
   }


   if(sel == five){
     return(0);
   }
// =========== ASA reading ============
   if(sel == two){
     std::cout <<"============= Convert ASA to frictions ==========\n\n";
     std::cout << "Please input the file name of the source ASA(ex:ASA.out)\n";
     std::cin >> read_file_name;
     std::ifstream asa_input_file;
     file_exist=is_false;
     while (file_exist == is_false) {
       asa_input_file.open(read_file_name,std::ios::in);
       if (!asa_input_file.good()) {
         std::cout << "Cannot open file= "<< read_file_name<<std::endl;
         std::cout << "Please type again or Ctrl+C to quit\n";
         std::cin >>read_file_name;
       } else {
         file_exist=is_true;
       }
     }

// ===== group ASA and Frictions =======
     std::cout << "Please input the output file name for grouping/non-grouping ASA and frictions(ex:ASA.dat)\n";
     std::cin >> write_file_name;
     std::cout << "Please input the value of eta for frictions (def=0.001 Pa.sec).\n";
     std::cin >> asa_eta;
     while (asa_input_file >> flag) {
       if (!strncmp(flag,"[ASA_OUTPUT]",12)) {
         asa_input_file >> asa_num;
         for (I0=0;I0<asa_num;I0++){
           asa_input_file >> dummy1>>atom_name[I0]>>res_name[I0]>>res_id[I0]>>dummy2>>read_asa_value[I0]>>dummy3>>dummy4>>dummy5>>dummy6;
         } 
         asa_input_file.close();
         std::cout << "Number of ASA ="<<asa_num<<std::endl; 
         std::cout << "1. Do not group ASA (frictions of all atoms)\n";
         std::cout << "2. Group ASA (For calculating weights of grouped atom)\n";
         std::cin >> group_sel;
         if (group_sel == one){
           output_asa_num=asa_num-1;
           r_ave=0;
           for (I0=0;I0<asa_num;I0++){
              output_asa_id[I0]=I0+1;
              strcpy(output_asa_atom_name[I0],atom_name[I0]);
              strcpy(output_asa_res_name[I0],res_name[I0]);
              output_asa_res_id[I0]=res_id[I0];
              output_asa_value[I0]=read_asa_value[I0];
              output_r_value[I0]=sqrt(output_asa_value[I0]/(4.*pi));
              output_xi_value[I0]=6.*pi*asa_eta*output_r_value[I0]; //The real friction is xi*10^-13
              r_ave=r_ave+sqrt(output_asa_value[I0]/(4*pi));
           }
           r_ave=r_ave/((float)output_asa_num+1);
         } else if(group_sel == two){
           std::cout << "Please input the file name of ASA grouping table(ex:group_table.dat)\n";
           std::cin >> read_file_name;
           std::ifstream asa_group_file;
           file_exist = is_false;
           while (file_exist == is_false) {
             asa_group_file.open(read_file_name,std::ios::in);
             if (!asa_group_file.good()) {
               std::cout << "Cannot open file= "<< read_file_name<<std::endl;
               std::cout << "Please type again or Ctrl+C to quit\n";
               std::cin >> read_file_name;
             } else {
               asa_group_file >> dummy1>>dummy2 >>dummy3>>dummy4;
               file_exist=is_true;
             }
           }
// =========== group table ============
           std::cout << "Constructing ASA group table...\n";
           group_target_num=-1;
           previous_group_target_id=0;
           while (asa_group_file >> read_group_target >> read_group_element>>read_group_weight_name){
             if(previous_group_target_id!=read_group_target){  //The only limitation of group table is that all elements of a target should not be separated. They should be put together.
               group_target_num++;    
               group_target_id[group_target_num]=read_group_target;
               previous_group_target_id=read_group_target;
               group_target_element_num[group_target_num]=0;
               group_target_element_id[group_target_num][group_target_element_num[group_target_num]]=read_group_element;
             } else {
               group_target_element_num[group_target_num]++;
               group_target_element_id[group_target_num][group_target_element_num[group_target_num]]=read_group_element;
             }
           }
           asa_group_file.close();
// ==== calculate frictions ===========
//           cout << group_target_num <<std::endl;
           output_asa_num=group_target_num;
           r_ave=0;
           for (I0=0;I0<=group_target_num;I0++){
             output_asa_id[I0]=group_target_id[I0];
             strcpy(output_asa_atom_name[I0],atom_name[group_target_id[I0]-1]);
             strcpy(output_asa_res_name[I0],res_name[group_target_id[I0]-1]);
             output_asa_res_id[I0]=res_id[group_target_id[I0]-1];
             output_asa_value[I0]=0;
             for (I1=0;I1<=group_target_element_num[I0];I1++){
//               cout<<"target="<<I0<<",target_id="<<group_target_id[I0]<<",ele_num="<<group_target_element_num[I0]<<",ele_id="<<group_target_element_id[I0][I1]<<std::endl;
               output_asa_value[I0]=output_asa_value[I0]+read_asa_value[group_target_element_id[I0][I1]-1];
             }
             r_ave=r_ave+sqrt(output_asa_value[I0]/(4*pi));
             output_r_value[I0]=sqrt(output_asa_value[I0]/(4.*pi));
             output_xi_value[I0]=6.*pi*asa_eta*output_r_value[I0]; //The real friction is xi*10^-13
              
           }
           r_ave=r_ave/((float)output_asa_num+1);
         }
// === output new ASA and frictions ===
         std::ofstream asa_output_file;
         asa_output_file.open(write_file_name,std::ios::out);
         asa_output_file << "ZETA_ERRE " << r_ave/l_ave<<"   # ZETA_ERRE*ELLE=<r_h>= "<<r_ave<<std::endl;
         asa_output_file << "ELLE " << l_ave<<"   # bead_num= "<<output_asa_num+1<<" ;bond_num= "<<output_asa_num<<std::endl; 
         asa_output_file << "FRICTIONS    "<<" # Frictions in unit of Pa.sec.m. These are calculated from the ASA of each beads. The repeated terms must be discarded according to BONDS (top of this file)#\n";      
         for (I0=0; I0<=output_asa_num;I0++){
           asa_output_file << output_xi_value[I0]*0.0000000001<<std::endl;
         }
         asa_output_file << "END_FRICTIONS\n\n";
         asa_output_file << "ETA  "<<asa_eta<<"    # Viscosity of solvent in unita' Pa.sec. 1.e-3 is for water #\n\n";
         asa_output_file << "# ID  atom  res/res_id  ASA(Group)  Xi  Friction(Xi*10^-13)  Bead_num= "<<output_asa_num+1<<std::endl;
         for (I0=0;I0<=output_asa_num;I0++){
           asa_output_file <<output_asa_id[I0]<<" "<<output_asa_atom_name[I0]<<" "<<output_asa_res_name[I0]<<" "<<output_asa_res_id[I0]<<" "<<output_asa_value[I0]<<" "<<output_xi_value[I0]<<" "<<output_xi_value[I0]*0.0000000001<<std::endl;
         }
         asa_output_file <<std::endl<<std::endl<<"# Below is for weight_table.dat\n\n"<<"# weight_name  weight_value\n";
         for (I0=0;I0<=output_asa_num;I0++){
           asa_output_file <<output_asa_res_name[I0]<<output_asa_res_id[I0]<<output_asa_atom_name[I0]<<" "<<output_asa_value[I0]<<std::endl;
         }
         asa_output_file <<std::endl<<std::endl<<"# Below is for group_table.dat\n\n"<<"# target_id element_id weight_name\n";
         for (I0=0;I0<=output_asa_num;I0++){
           asa_output_file << I0+1 <<"  "<<I0+1<<"  "<<output_asa_res_name[I0]<<output_asa_res_id[I0]<<output_asa_atom_name[I0]<<std::endl;   
         }
         asa_output_file.close();
         return(0);
       }
     }


// ======== pdb to his ================

   }
   else if(sel == one){
     std::cout <<"============= Convert pdb to his ================\n\n";
     std::cout << "Please input the source pdb file name(ex:*.pdb)\n";
     std::cin >> read_file_name;
     atom_num_dim=check_pdb_atom_num(read_file_name);
     std::cout <<"Detect atom_num="<<atom_num_dim<<std::endl;
     std::ifstream pdb_file;
     file_exist = is_false;
     while (file_exist == is_false) {
       pdb_file.open(read_file_name,std::ios::in);
       if (!pdb_file.good()) {
         std::cout << "Cannot open file= "<< read_file_name<<std::endl;
         std::cout << "Please type again or Ctrl+C to quit\n";
         std::cin >> read_file_name;
       } else {
         file_exist= is_true;   
       }
     }
     std::cout << "Please input file name of binary his file(ex:*.his)\n";
     std::cin >> write_file_name;
     std::ofstream his_file;
     his_file.open(write_file_name,std::ios::out|std::ios::binary);
     std::cout << "1. Do not group atoms (Just convert to binary format)\n";
     std::cout << "2. Use original weightings (weight_table.dat and group_table.dat required)\n";
     std::cout << "3. Use square root of weightings (weight_table.dat and group_table.dat required)\n";
     std::cin >> weight_sel;
// ====== dynamical float array mat====
     if (weight_sel == one){
       atom_num=0;
       frame_num=0;
       float *mat0= new float [atom_num_dim*mat_dim];
       while (pdb_file >> flag) {
         if (!strncmp(flag,"ATOM",3)) {
           pdb_file >> dummy1 >> atom_name[atom_num] >> res_name[atom_num] >> dummy2 >>mat0[atom_num*3] >>mat0[atom_num*3+1] >>mat0[atom_num*3+2] >>dummy3>>dummy4;
//       cout << atom_num  << " "<<atom_name[atom_num] <<" "<< res_name[atom_num] <<" "<<mat0[atom_num*3]<<" "<<mat0[atom_num*3+1]<<" "<<mat0[atom_num*3+2] <<std::endl;
           atom_num++;
         }
         if (!strncmp(flag,"ENDMDL",3)) {
           frame_num++;
           std::cout << "Frame=" << frame_num <<" complete!\n";
           his_file.write((char *) &frame_num, sizeof(frame_num)).write((char *) mat0, sizeof(*mat0)*atom_num_dim*mat_dim); 
           atom_num=0;
         }
       }
       delete [] mat0;
       pdb_file.close();
       his_file.close();
       return(0);
     } 

     std::cout << "Please input file name of weighting parameters(ex:weight_table.dat)\n";
     std::cin >> read_file_name;
     std::ifstream weight_file;
     file_exist=is_false;
     while (file_exist == is_false) {
       weight_file.open(read_file_name,std::ios::in);
       if (!weight_file.good()) {
         std::cout << "Cannot open file= "<< read_file_name<<std::endl;
         std::cout << "Please type again or Ctrl+C to quit\n";
         std::cin >> read_file_name;
       } else {
         file_exist=is_true;
         weight_file >> dummy1 >> dummy2 >> dummy3;
       }
     }
     std::cout << "Please input file name of atomic grouping list(ex:group_table.dat)\n";
     std::cin >> read_file_name;
     std::ifstream group_file;
     file_exist=is_false;
     while (file_exist==is_false) {
       group_file.open(read_file_name,std::ios::in);
       if (!group_file.good()) {
         std::cout << "Cannot open file= "<< read_file_name<<std::endl;
         std::cout << "Please type again or Ctrl+C to quit\n";
         std::cin >> read_file_name;
       } else {
         file_exist=is_true;
         group_file >> dummy1 >> dummy2 >> dummy3 >> dummy4;
       }
     }
     std::cout << "Checking atomic group and weightings\n"; 
     weight_num=0;
     while (weight_file >> weight_name[weight_num] >> weight_value[weight_num]){
       if (weight_sel == three){
         weight_value[weight_num]=sqrt(weight_value[weight_num]);
       }
       std::cout << weight_num <<" "<<weight_name[weight_num] <<" "<<weight_value[weight_num]<<std::endl;
       weight_num++;
     }
     std::cout << "There are "<<weight_num <<" weightings\n";  //Note it is from 0 to weight_num-1.
     weight_file.close();
// =========== group table ============
     std::cout << "Constructing atomic group table...\n";
     group_target_num=-1;
     previous_group_target_id=0;
     while (group_file >> read_group_target >> read_group_element >> read_group_weight_name){
       if(previous_group_target_id!=read_group_target){
         group_target_num++;
         group_target_id[group_target_num]=read_group_target;
         previous_group_target_id=read_group_target;
         group_target_element_num[group_target_num]=0;
         group_target_element_id[group_target_num][group_target_element_num[group_target_num]]=read_group_element;
       } else {
         group_target_element_num[group_target_num]++;
         group_target_element_id[group_target_num][group_target_element_num[group_target_num]]=read_group_element;
       }
       I0=0;
       while(strncmp(read_group_weight_name,weight_name[I0],atom_name_max)) {
       I0++;    
       }
       group_target_element_weight[group_target_num][group_target_element_num[group_target_num]]=weight_value[I0];
     }
     group_file.close();
// ==== weighting and grouping ========
     frame_num=0;
     atom_num=0;
     float *mat1= new float [atom_num_dim*mat_dim];
     while (pdb_file >> flag) {
       if (!strncmp(flag,"ATOM",3)) {
         pdb_file >> dummy1 >> atom_name[atom_num] >> res_name[atom_num] >> dummy2 >>mat1[atom_num*3] >>mat1[atom_num*3+1] >>mat1[atom_num*3+2] >>dummy3>>dummy4;
//        cout << atom_num  << " "<<atom_name[atom_num] <<" "<< res_name[atom_num] <<" "<<mat1[atom_num*3]<<" "<<mat1[atom_num*3+1]<<" "<<mat1[atom_num*3+2] <<std::endl;
         atom_num++;
       }
       if (!strncmp(flag,"ENDMDL",3)) {    //Completly read a frame
         frame_num++;
//         cout << "Frame=" << frame_num <<" complete!\n";  
         for (I0=0;I0<=group_target_num;I0++){
            weight_denominator=0;
            for (I1=0;I1<=group_target_element_num[I0];I1++){
              weight_denominator=weight_denominator+group_target_element_weight[I0][I1];
            }
            for (I2=1;I2<=3;I2++){
              weight_numerator=0;
              for (I1=0;I1<=group_target_element_num[I0];I1++){
                weight_numerator=weight_numerator+group_target_element_weight[I0][I1]*mat1[group_target_element_id[I0][I1]*3-I2];
              }
              mat1[group_target_id[I0]*3-I2]=weight_numerator/weight_denominator;
            }
         }
         his_file.write((char *) &frame_num, sizeof(frame_num)).write((char *) mat1, sizeof(*mat1)*atom_num_dim*mat_dim);
//         std::cout << frame_num;
//         for (I1=0;I1<=atom_num;I1++){
//           std::cout << mat1[I1*3]<<" "<<mat1[I1*3+1]<<" "<<mat1[I1*3+2]<<std::endl;
//         }
         atom_num=0;
       }
     }
// ====================================
     delete [] mat1;
     pdb_file.close();
     his_file.close();
     return(0);
   }
   else if(sel == three) {
     std::cout <<"============= Trajectory of beads ===============\n\n";
     std::cout << "Please input the source pdb file name(ex:*.pdb)\n";
     std::cin >> read_file_name;
     atom_num_dim=check_pdb_atom_num(read_file_name);
     std::cout <<"Detect atom_num="<<atom_num_dim<<std::endl;
     std::ifstream atom_pdb_file;
     file_exist=is_false;
     while (file_exist== is_false) {
       atom_pdb_file.open(read_file_name,std::ios::in);
       if (!atom_pdb_file.good()) {
         std::cout << "Cannot open file= "<< read_file_name<<std::endl;
         std::cout << "Please type again or Ctrl+C to quit\n";
         std::cin >> read_file_name;
       } else {
         file_exist=is_true;
       }
     }
     std::cout << "Please input the output pdb file name(ex:*.pdb)\n";
     std::cin >> write_file_name;
     std::cout << "Please input the file name of bond list(ex: bond_list.dat)\n";
     std::cin >> read_file_name;
     std::ifstream bead_bond_file;
     file_exist=is_false;
     while (file_exist == is_false) {
       bead_bond_file.open(read_file_name,std::ios::in);
       if (!bead_bond_file.good()) {
         std::cout << "Cannot open file= "<< read_file_name<<std::endl;
         std::cout << "Please type again or Ctrl+C to quit\n";
         std::cin >> read_file_name;
       } else {
         file_exist=is_true;
         bead_bond_file >>dummy1>>dummy2>>dummy3>>dummy4>>dummy5>>dummy6>>dummy7;
       }
     }
     std::cout << "1. Use original coordinates\n";
     std::cout << "2. Use geometric center of grouped atom(group_table.dat required)\n";
     std::cout << "3. Use weighting position (group_table.dat and weight_table.dat required)\n";
     std::cin >> bead_sel;
// =========== bead bond list =========
     bond_num=0;
     while(bead_bond_file>>bond_target_atom_name[bond_num]>>bond_target_res_name[bond_num]>>bond_target_res_id[bond_num]>>bond_end_atom_name[bond_num]>>bond_end_res_name[bond_num]>>bond_end_res_id[bond_num]) {
      bond_num++;
     }
// =========== read weight table ======
     if (bead_sel == three) {
       std::cout << "1. Use original weightings\n";
       std::cout << "2. Use square root of weightings\n";
       std::cin >> weight_sel;
       std::cout << "Please input file name of weighting parameters(ex:weight_table.dat)\n";
       std::cin >> read_file_name;
       std::ifstream weight_file;
       file_exist=is_false;
       while (file_exist == is_false) {
         weight_file.open(read_file_name,std::ios::in);
         if (!weight_file.good()) {
           std::cout << "Cannot open file= "<< read_file_name<<std::endl;
           std::cout << "Please type again or Ctrl+C to quit\n";
           std::cin >> read_file_name;
         } else {
           file_exist=is_true;
           weight_file >> dummy1 >> dummy2 >> dummy3;
         }
       }
       std::cout << "Please input file name of atomic grouping list(ex:group_table.dat)\n";
       std::cin >> read_file_name;
       std::ifstream group_file;
       file_exist=is_false;
       while (file_exist==is_false) {
         group_file.open(read_file_name,std::ios::in);
         if (!group_file.good()) {
           std::cout << "Cannot open file= "<< read_file_name<<std::endl;
           std::cout << "Please type again or Ctrl+C to quit\n";
           std::cin >> read_file_name;
         } else {
           file_exist=is_true;
           group_file >> dummy1 >> dummy2 >> dummy3 >> dummy4;
         }
       }
       std::cout << "Checking atomic group and weightings\n";
       weight_num=0;
       while (weight_file >> weight_name[weight_num] >> weight_value[weight_num]){
         if (weight_sel == two){
           weight_value[weight_num]=sqrt(weight_value[weight_num]);
         }
         std::cout << weight_num <<" "<<weight_name[weight_num] <<" "<<weight_value[weight_num]<<std::endl;
         weight_num++;
       }
       std::cout << "There are "<<weight_num <<" weightings\n";  //Note it is from 0 to weight_num-1.
       weight_file.close();
// =========== group table ============
       std::cout << "Constructing atomic group table...\n";
       group_target_num=-1;
       previous_group_target_id=0;
       while (group_file >> read_group_target >> read_group_element >> read_group_weight_name){
         if(previous_group_target_id!=read_group_target){
           group_target_num++;
           group_target_id[group_target_num]=read_group_target;
           previous_group_target_id=read_group_target;
           group_target_element_num[group_target_num]=0;
           group_target_element_id[group_target_num][group_target_element_num[group_target_num]]=read_group_element;
         } else {
           group_target_element_num[group_target_num]++;
           group_target_element_id[group_target_num][group_target_element_num[group_target_num]]=read_group_element;
         }
         I0=0;
         while(strncmp(read_group_weight_name,weight_name[I0],atom_name_max)) {
         I0++;
         }
         group_target_element_weight[group_target_num][group_target_element_num[group_target_num]]=weight_value[I0];
       }
       group_file.close();



     }
// =========== read group table =======
     if (bead_sel == two) {
       std::cout << "Please input the file name of pdb grouping(ex: bead_group_table.dat)\n";
       std::cin >> read_file_name;
       std::ifstream group_file;
       file_exist=is_false;
       while (file_exist== is_false) {
         group_file.open(read_file_name,std::ios::in);
         if (!group_file.good()) {
           std::cout << "Cannot open file= "<< read_file_name<<std::endl;
           std::cout << "Please type again or Ctrl+C to quit\n";
           std::cin >> read_file_name;
         } else {
           file_exist=is_true;
           group_file >> dummy1 >> dummy2 >> dummy3 >> dummy4;
         }
       }
       std::cout << "Constructing atomic group table...\n";
       group_target_num=-1;
       previous_group_target_id=0;
       while (group_file >> read_group_target >> read_group_element >> read_group_weight_name){
         if (previous_group_target_id!=read_group_target){
           group_target_num++;
           group_target_id[group_target_num]=read_group_target;
           previous_group_target_id=read_group_target;
           group_target_element_num[group_target_num]=0;
           group_target_element_id[group_target_num][group_target_element_num[group_target_num]]=read_group_element;
         } else {
           group_target_element_num[group_target_num]++;
           group_target_element_id[group_target_num][group_target_element_num[group_target_num]]=read_group_element;
         }
// ====================================
       }
       group_file.close();
     }
// =========== read source pdb ========
     std::ofstream bead_pdb_file;
     bead_pdb_file.open(write_file_name,std::ios::out);
     atom_num=0;
     frame_num=0;
     float *mat2= new float [atom_num_dim*mat_dim];
     while (atom_pdb_file >> flag) {
       if (!strncmp(flag,"ATOM",3)) {
         atom_pdb_file >> dummy1 >> atom_name[atom_num] >> res_name[atom_num] >> res_id[atom_num] >>mat2[atom_num*3] >>mat2[atom_num*3+1] >>mat2[atom_num*3+2] >>dummy3>>dummy4;
//       cout << atom_num  << " "<<atom_name[atom_num] <<" "<< res_name[atom_num] <<" "<<mat2[atom_num*3]<<" "<<mat2[atom_num*3+1]<<" "<<mat2[atom_num*3+2] <<std::endl;
         atom_num++;
       }
       if (!strncmp(flag,"ENDMDL",3)) {
// === use geometric center of ASA ====
         frame_num++;
         if(bead_sel == two) {
           for (I0=0;I0<=group_target_num;I0++){
             bead_mat_x=0;
             bead_mat_y=0;
             bead_mat_z=0;
             for (I1=0;I1<=group_target_element_num[I0];I1++){
               bead_mat_x=bead_mat_x+mat2[group_target_element_id[I0][I1]*3-3];
               bead_mat_y=bead_mat_y+mat2[group_target_element_id[I0][I1]*3-2];
               bead_mat_z=bead_mat_z+mat2[group_target_element_id[I0][I1]*3-1];
             }
             mat2[group_target_id[I0]*3-3]=bead_mat_x/(group_target_element_num[I0]+1);
             mat2[group_target_id[I0]*3-2]=bead_mat_y/(group_target_element_num[I0]+1);
             mat2[group_target_id[I0]*3-1]=bead_mat_z/(group_target_element_num[I0]+1);
           }

         } else if (bead_sel == three) {
           for (I0=0;I0<=group_target_num;I0++){
             weight_denominator=0;
             for (I1=0;I1<=group_target_element_num[I0];I1++){
                weight_denominator=weight_denominator+group_target_element_weight[I0][I1];
             }
             for (I2=1;I2<=3;I2++){
               weight_numerator=0;
               for (I1=0;I1<=group_target_element_num[I0];I1++){
                 weight_numerator=weight_numerator+group_target_element_weight[I0][I1]*mat2[group_target_element_id[I0][I1]*3-I2];
               }
               mat2[group_target_id[I0]*3-I2]=weight_numerator/weight_denominator;
             } 
           }      
         }
// ======== output bead pdb ===========
         bead_pdb_file << "TITLE      Bead connection\n"<<"REMARK     Bond number= "<<bond_num<<std::endl<<"REMARK     Bead number= "<<bond_num+1<<std::endl<<"REMARK     Frame ID= "<<frame_num<<std::endl<<"REMARK     Value right after coordinates is the bond distance\n"<<"MODEL        0\n";
 
         for (I0=0;I0<bond_num;I0++) {
           for (I1=0;I1<atom_num;I1++) {
             if(!strncmp(atom_name[I1],bond_target_atom_name[I0],5) && !strncmp(res_name[I1],bond_target_res_name[I0],3) && res_id[I1] == bond_target_res_id[I0]) {
               output_bond_target_id=I1;  //starting from 0
             }
             if(!strncmp(atom_name[I1],bond_end_atom_name[I0],5) && !strncmp(res_name[I1],bond_end_res_name[I0],3) && res_id[I1] == bond_end_res_id[I0]) {
               output_bond_end_id=I1; //starting from 0
             }
           }
           bond_dist=sqrt(pow((mat2[output_bond_target_id*3]-mat2[output_bond_end_id]),2)+pow((mat2[output_bond_target_id*3+1]-mat2[output_bond_end_id*3+1]),2)+pow((mat2[output_bond_target_id*3+2]-mat2[output_bond_end_id*3+2]),2));
           bead_pdb_file.setf(std::ios::left, std::ios::adjustfield);
           bead_pdb_file.width(6);
           bead_pdb_file <<"ATOM  ";
           bead_pdb_file.setf(std::ios::right, std::ios::adjustfield);
           bead_pdb_file.width(5);
           bead_pdb_file <<I0*2+1;
           bead_pdb_file.width(1);
           bead_pdb_file <<" ";
           bead_pdb_file.setf(std::ios::internal, std::ios::adjustfield);
           bead_pdb_file.width(4);
           bead_pdb_file <<atom_name[output_bond_target_id];
           bead_pdb_file.width(1);
           bead_pdb_file <<" ";
           bead_pdb_file.width(3);
           bead_pdb_file <<res_name[output_bond_target_id];
           bead_pdb_file.width(2);
           bead_pdb_file <<"  ";
           bead_pdb_file.width(4);
           bead_pdb_file <<res_id[output_bond_target_id];
           bead_pdb_file.width(4);
           bead_pdb_file <<"    ";
           bead_pdb_file.width(8);
           bead_pdb_file <<mat2[output_bond_target_id*3];
           bead_pdb_file.width(8);
           bead_pdb_file <<mat2[output_bond_target_id*3+1];
           bead_pdb_file.width(8);
           bead_pdb_file <<mat2[output_bond_target_id*3+2]<<std::endl;       
           bead_pdb_file.setf(std::ios::left, std::ios::adjustfield);
           bead_pdb_file.width(6);
           bead_pdb_file <<"ATOM  ";
           bead_pdb_file.setf(std::ios::right, std::ios::adjustfield);
           bead_pdb_file.width(5);
           bead_pdb_file <<I0*2+2;
           bead_pdb_file.width(1);
           bead_pdb_file <<" ";
           bead_pdb_file.setf(std::ios::internal, std::ios::adjustfield);
           bead_pdb_file.width(4);
           bead_pdb_file <<atom_name[output_bond_end_id];
           bead_pdb_file.width(1);
           bead_pdb_file <<" ";
           bead_pdb_file.width(3);
           bead_pdb_file <<res_name[output_bond_end_id];
           bead_pdb_file.width(2);
           bead_pdb_file <<"  ";
           bead_pdb_file.width(4);
           bead_pdb_file <<res_id[output_bond_end_id];
           bead_pdb_file.width(4);
           bead_pdb_file <<"    ";
           bead_pdb_file.width(8);
           bead_pdb_file <<mat2[output_bond_end_id*3];
           bead_pdb_file.width(8);
           bead_pdb_file <<mat2[output_bond_end_id*3+1]; 
           bead_pdb_file.width(8);
           bead_pdb_file <<mat2[output_bond_end_id*3+2];
           bead_pdb_file.width(6);
           bead_pdb_file <<"      ";
           bead_pdb_file.width(6);
           bead_pdb_file <<bond_dist<<std::endl;

           bead_pdb_file.width(6);
           bead_pdb_file<<"CONECT";
           bead_pdb_file.width(5);
           bead_pdb_file<<I0*2+1;
           bead_pdb_file.width(5);
           bead_pdb_file<<I0*2+2<<std::endl;
         }
         bead_pdb_file <<"ENDMDL";
         atom_num=0;
         std::cout << "Frame=" << frame_num <<" complete!\n";
       }
     }
     delete [] mat2;
     atom_pdb_file.close();
     bead_pdb_file.close();
     return(0);
   }
}


/* ---------------------------------------------------------------------- */

int check_pdb_atom_num(char input_pdb_file_name[name_max])
{
int atom_num=0;
char flag[name_max];
std::ifstream pdb_file;
pdb_file.open(input_pdb_file_name,std::ios::in);
if (!pdb_file.good()){
  std::cerr << "Cannot open file! " <<std::endl;
  return(0);
}
while (pdb_file >> flag){
  if(!strncmp(flag,"ATOM",3)) {
    atom_num++;
  }
  if(!strncmp(flag,"ENDMDL",6)) {
    pdb_file.close();
    return atom_num;
  }
}
pdb_file.close();
std::cerr << "Cannot find a single frame in this pdb file, please check if it includes ENDMDL"<<std::endl;
return atom_num;
}

/* ---------------------------------------------------------------------- */

