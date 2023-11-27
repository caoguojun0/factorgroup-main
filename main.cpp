#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include <set>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/matrix.h>
#include <optional>
#include "create_factorgroup.hpp"
#include <sstream>
//g++ -std=c++11 -I/opt/homebrew/Cellar/ntl/11.5.1/include main.cpp create_factorgroup.cpp -o main -L/opt/homebrew/Cellar/ntl/11.5.1/lib -lntl -lgmp -lpthread
int main() 
{
    int n = 2;                                  // prime number GF(n), field characteristic
    int max_degree = 2;                         // Maximum degree of a polynomial
    vector<int> input_polynomial = {1, 1, 1};   // Setting the coefficients of a polynomial manually

    
    bool print_factorgroup = true;              // whether to take out the Factor Group
    bool print_debag_poly = true;               //   derive intermediate polynomials
    bool find_primitive = true;                 // primitive element
    bool print_all_list = true;                 // output the list of all polynomials in the factor field
    bool create_table = true;                   // create a table


    FactorGroupCreator FactorGroupCreatorObj;
    FactorGroupCreatorObj.print_all_list = print_all_list;

    FactorGroupCreatorObj.FactorGroupCreate(n, max_degree, input_polynomial,print_all_list,print_factorgroup,  print_debag_poly,find_primitive,create_table);

    return 0;
}



