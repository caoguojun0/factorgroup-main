#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include <set>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/matrix.h>
#include "create_factorgroup.hpp"
#include <sstream>
using namespace std;
using namespace NTL;

class FactorGroupCreator {
public:
    FactorGroupCreator();

    std::vector<long> computeDivisors(long number);
    void printCayleyTable(const vector<vector<ZZ_pX>>& table, const ZZ_pX& modulus);
    void createCayleyTable(const vector<ZZ_pX>& primitive_elements, const ZZ_pX& modulus, int p);
    void printDivisors(const std::vector<long>& divisors);
    static bool compare_ntl_polynomials(const NTL::ZZ_pX& a, const NTL::ZZ_pX& b);
    bool checkPrimitive(const ZZ_pX& candidate, const ZZ_pX& modulus, long p, const std::vector<long>& divisors, bool print_debag_poly);
    void FactorGroupCreate(int n, int max_degree,vector<int> input_polynomial,bool print_all_list,bool print_factorgroup,bool print_debag_poly,bool find_primitive,bool create_table);
    std::vector<int> read_polynomial(int degree);
    void generate_all_polynomials(int n, int degree, std::vector<int>& p, int current_degree, std::vector<std::vector<int> >& all_polynomials);
    std::set<vector<int> > build_factor_group(const std::vector<int>& f, int n, int max_degree);
    void print_polynomial(ostream& out, const vector<int>& p);
    bool is_zero(const vector<int>& p);
    vector<int> polynomial_evaluation(const std::vector<int>& polynomial, int x, int n);
    vector<int> polynomial_multiplication(const std::vector<int>& a, const std::vector<int>& b, int n);
    ZZ_pX convert_to_ntl_polynomial(const std::vector<int>& p, int n);
    static string ntl_to_string(const ZZ_pX& p);
    bool is_zero(const std::vector<int>& p, int n);
    vector<ZZ_pX> find_primitive_elements(const std::vector<ZZ_pX>& ntl_polynomials, const ZZ_pX& modulus, long p, long max_degree, bool print_debag_poly);
    std::vector<int> divide_polynomials(const std::vector<int>& dividend, const std::vector<int>& divisor, int n);
    void printPrimitiveElementsList(const std::vector<ZZ_pX>& primitive_elements);
    void printRandomPrimitiveElement(const std::vector<ZZ_pX>& primitive_elements);
    void printFactorGroup(const std::vector<ZZ_pX>& ntl_polynomials, const std::vector<int>& input_polynomial, int n);
    int n;
    int max_degree;
    vector<int> input_polynomial;
    std::set<vector<int>> factor_group;
    vector<ZZ_pX> ntl_polynomials;
    vector<ZZ_pX> primitive_elements;
    bool print_all_list;
    bool print_factorgroup;
    bool print_debag_poly;
    bool find_primitive;
    bool create_table;
private:
    
};
