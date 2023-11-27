#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include <set>
#include <iomanip>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/matrix.h>
#include <optional>
#include "create_factorgroup.hpp"
#include <sstream>
#include <algorithm>
#include <random>

using namespace std;
using namespace NTL;


FactorGroupCreator::FactorGroupCreator()
{
    
}

std::vector<int> FactorGroupCreator::divide_polynomials(const std::vector<int>& dividend, const std::vector<int>& divisor, int n) {
    int dividend_degree = dividend.size() - 1;
    int divisor_degree = divisor.size() - 1;

    std::vector<int> quotient(dividend_degree - divisor_degree + 1, 0);
    std::vector<int> remainder = dividend;

    while (dividend_degree >= divisor_degree) {
        int coef = (remainder[dividend_degree] * ((n - divisor[divisor_degree]) % n)) % n;
        int degree_diff = dividend_degree - divisor_degree;

        quotient[degree_diff] = (quotient[degree_diff] + coef) % n;

        for (int i = 0; i <= divisor_degree; ++i) {
            remainder[degree_diff + i] = (remainder[degree_diff + i] + coef * divisor[i]) % n;
        }

        while (dividend_degree >= 0 && remainder[dividend_degree] == 0) {
            --dividend_degree;
        }
    }

    return remainder;
}

void FactorGroupCreator::generate_all_polynomials(int n, int degree, std::vector<int>& p, int current_degree, std::vector<std::vector<int> >& all_polynomials) {
    if (current_degree > degree) {
        all_polynomials.push_back(p);
        return;
    }

    for (int j = 0; j < n; ++j) {
        p[current_degree] = j;
        generate_all_polynomials(n, degree, p, current_degree + 1, all_polynomials);
    }
}

std::set<std::vector<int> > FactorGroupCreator::build_factor_group(const std::vector<int>& f, int n, int max_degree) {
    std::vector<int> p(max_degree + 1, 0);
    std::vector<std::vector<int> > all_polynomials;
    generate_all_polynomials(n, max_degree, p, 0, all_polynomials);

    std::set<std::vector<int> > factor_group;

    for (size_t i = 0; i < all_polynomials.size(); ++i) {
        std::vector<int> remainder = divide_polynomials(all_polynomials[i], f, n);
        factor_group.insert(remainder);
    }

    return factor_group;
}

void FactorGroupCreator::print_polynomial(std::ostream& out, const std::vector<int>& p) {
    for (int i = p.size() -1; i >= 0; --i) {
        if (p[i] != 0) {
            if (i != p.size() - 1) {
                out << " + ";
            }
            if (i > 0) {
                out << p[i] << "x^" << i;
            } else {
                out << p[i];
            }
        }
    }
    out << std::endl;
}

std::vector<int> FactorGroupCreator::read_polynomial(int degree) {
    std::vector<int> polynomial(degree + 1, 0);
    std::cout << "Enter the coefficients of the polynomial (from the higher degree to the lower degree): ";
    for (int i = degree; i >= 0; --i) {
        std::cin >> polynomial[i];
    }
    return polynomial;
}

bool FactorGroupCreator::is_zero(const std::vector<int>& p, int n) {
    for (int coef : p) {
        if (coef % n != 0) {
            return false;
        }
    }
    return true;
}

std::vector<int> FactorGroupCreator::polynomial_evaluation(const std::vector<int>& polynomial, int x, int n) {
    std::vector<int> evaluated(polynomial.size(), 0);
    int power_of_x = 1;
    for (int i = 0; i < polynomial.size(); ++i) {
        evaluated[i] = (polynomial[i] * power_of_x) % n;
        power_of_x = (power_of_x * x) % n;
    }
    return evaluated;
}

std::vector<int> FactorGroupCreator::polynomial_multiplication(const std::vector<int>& a, const std::vector<int>& b, int n) {
    std::vector<int> product(a.size() + b.size() - 1, 0);

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            product[i + j] = (product[i + j] + a[i] * b[j]) % n;
        }
    }

    return product;
}

ZZ_pX FactorGroupCreator::convert_to_ntl_polynomial(const std::vector<int>& p, int n) {
    ZZ_p::init(conv<ZZ>(n));
    ZZ_pX ntl_poly;

    for (int i = 0; i < p.size(); ++i) {
        SetCoeff(ntl_poly, i, conv<ZZ_p>(p[i]));
    }

    return ntl_poly;
}

std::string FactorGroupCreator::ntl_to_string(const ZZ_pX& p) {
    std::ostringstream oss;
    bool first_term = true;
    for (int i = deg(p); i >= 0; --i) {
        ZZ_p coef = coeff(p, i);
        if (coef != 0) {
            if (!first_term) {
                oss << " + ";
            }
            first_term = false;
            if (i > 0) {
                oss << coef << "x^" << i;
            } else {
                oss << coef;
            }
        }
    }
    return oss.str();
}


std::vector<long> FactorGroupCreator::computeDivisors(long number) {
    std::vector<long> divisors;
    for (long i = 1; i <= number; ++i) {
        if (number % i == 0) {
            divisors.push_back(i);
        }
    }
    return divisors;
}

void FactorGroupCreator::printDivisors(const std::vector<long>& divisors) {
    std::cout << "Simple divisors of a number " << divisors.back() << ": ";
    for (long divisor : divisors) {
        std::cout << divisor << " ";
    }
    std::cout << std::endl;
}

bool FactorGroupCreator::checkPrimitive(const ZZ_pX& candidate, const ZZ_pX& modulus, long p, const std::vector<long>& divisors, bool print_debag_poly) {
    std::string candidate_str = ntl_to_string(candidate);
    std::string modulus_str = ntl_to_string(modulus);
    for (long divisor : divisors) {
        ZZ_pX power_mod_result;
        PowerMod(power_mod_result, candidate, divisor, modulus);
        if (power_mod_result == 1) {
            return false;
        }
        if (print_debag_poly)
        {
            std::cout <<" [ [ (" <<candidate_str <<")^" <<divisor << " ] mod "<<modulus_str<<" ] mod "<<p<<" = ";
            std::string  result= ntl_to_string(power_mod_result);
            std::cout << result <<endl;
             

        }
       
    }
    return true;
}

std::vector<ZZ_pX> FactorGroupCreator::find_primitive_elements(const std::vector<ZZ_pX>& ntl_polynomials, const ZZ_pX& modulus, long p, long max_degree, bool print_debag_poly) {
    long field_size = p * max_degree;
    long phi_field_size = field_size - 1;
    std::vector<long> divisors = computeDivisors(phi_field_size);
    if(print_factorgroup || find_primitive)
    {
        std::cout << "( p * max_degree ) - 1 = ( "<<p<<" * "<<max_degree<<" )"<<" - 1 = " << phi_field_size << std::endl;
        printDivisors(divisors);
        printf("\n");
    }

    std::vector<ZZ_pX> primitive_elements;
    for (const auto& candidate : ntl_polynomials) {
        if (IsZero(candidate) || candidate==1) {
            continue;
        }
        
        if(print_debag_poly)
        {
            printf("\n");
            std::cout << "For a polynomial ";
            std::string candidate_str = ntl_to_string(candidate);
            std::cout << candidate_str<<":"<<endl;
        }
   
        if (checkPrimitive(candidate, modulus, p, divisors, print_debag_poly)) {
            primitive_elements.push_back(candidate);
        }
    }

    return primitive_elements;
}

void FactorGroupCreator::printFactorGroup(const std::vector<ZZ_pX>& ntl_polynomials, const std::vector<int>& input_polynomial, int n) {
    NTL::ZZ_pX temp = convert_to_ntl_polynomial(input_polynomial, n);
    std::string temp_ = ntl_to_string(temp);

    std::cout << "Z_p[x]/f(x) = Z_" << n << "[x]/(" << temp_ << ") \n";
    bool is_irreducible = NTL::IterIrredTest(temp);
    if (is_irreducible) {
        std::cout << "The factor ring is a field since the polynomial f(x) is irreducible and consists of the following polynomials:" << std::endl;
    } else {
        std::cout << "The factor ring is not a field since the polynomial f(x) is reducible and consists of the following polynomials:" << std::endl;
    }
    for (const auto& poly : ntl_polynomials) {
        std::string answer = ntl_to_string(poly);
        std::cout << answer << std::endl;
    }

    
}

void FactorGroupCreator::printRandomPrimitiveElement(const std::vector<ZZ_pX>& primitive_elements) {
    int size = primitive_elements.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, size - 1);
    int random_number = dist(gen);
    printf("\n");
    std::cout << "A primitive element has been found:\n";
    std::string temp2_ = ntl_to_string(primitive_elements[random_number]);
    std::cout << temp2_ << std::endl;
}

void FactorGroupCreator::printPrimitiveElementsList(const std::vector<ZZ_pX>& primitive_elements) {
    std::cout << "List of primitive elements:\n";
    for (const auto& poly : primitive_elements) {
        std::string answer = ntl_to_string(poly);
        std::cout << answer << std::endl;
    }
}

void FactorGroupCreator::printCayleyTable(const vector<vector<ZZ_pX>>& table, const ZZ_pX& modulus) {
    int max_length = 0;
    for (const auto& row : table) {
        for (const auto& elem : row) {
            int length = ntl_to_string(elem).length();
            if (length > max_length) {
                max_length = length;
            }
        }
    }

    std::string separator = std::string((max_length + 3) * table.size() + 1, '-');

    for (size_t i = 0; i < table.size(); ++i) {
        std::cout << separator << std::endl;
        for (size_t j = 0; j < table[i].size(); ++j) {
            std::string elem_str;
            if (table[i][j] == modulus) {
                elem_str = "*";
            } else {
                elem_str = ntl_to_string(table[i][j]);
            }
            if (table[i][j] == 0)
            {
                elem_str = "0";
            }
            std::cout << "| " << std::setw(max_length) << std::left << elem_str << " ";
        }
        std::cout << "|" << std::endl;
    }
    std::cout << separator << std::endl;
}

void FactorGroupCreator::createCayleyTable(const vector<ZZ_pX>& primitive_elements, const ZZ_pX& modulus, int p) {

    size_t n = primitive_elements.size();
    vector<vector<ZZ_pX>> cayley_table(n + 1, vector<ZZ_pX>(n + 1));

    cayley_table[0][0] = modulus;
    for (size_t i = 0; i < n; ++i) {
        cayley_table[0][i + 1] = primitive_elements[i];
        cayley_table[i + 1][0] = primitive_elements[i];
    }

    ZZ_pX temp_result;
    for (size_t i = 1; i <= n ; ++i) {
        for (size_t j = 1; j <= n ; ++j) {
            mul(temp_result, cayley_table[i][0], cayley_table[0][j]);
            rem(temp_result, temp_result, modulus);
            if (temp_result == 0)
            {
                cayley_table[i][j] = 0;
            }
            cayley_table[i][j] = temp_result;
        }
    }


    printf("\n");
    printCayleyTable(cayley_table, modulus);
    printf("\n");
}

void FactorGroupCreator::FactorGroupCreate(int n, int max_degree, vector<int> input_polynomial, bool print_all_list, bool print_factorgroup, bool print_debag_poly, bool find_primitive,bool create_table) {
    std::set<std::vector<int>> factor_group = build_factor_group(input_polynomial, n, max_degree);
    std::vector<ZZ_pX> ntl_polynomials;

    for (const auto &polynomial : factor_group) {
        ZZ_pX ntl_poly = convert_to_ntl_polynomial(polynomial, n);
        ntl_polynomials.push_back(ntl_poly);
        ntl_to_string(ntl_poly);
    }

    vector<ZZ_pX> primitive_elements = find_primitive_elements(ntl_polynomials, convert_to_ntl_polynomial(input_polynomial, n), n, max_degree, print_debag_poly);

    if (print_factorgroup) {
        printf("\n");
        printFactorGroup(ntl_polynomials, input_polynomial, n);
    }

    if (create_table)
    {
        createCayleyTable(ntl_polynomials, convert_to_ntl_polynomial(input_polynomial, n), n);
    }
    
    if (find_primitive) {
        printRandomPrimitiveElement(primitive_elements);
    }

    if (print_all_list) {
        printPrimitiveElementsList(primitive_elements);
    }
}


bool FactorGroupCreator::compare_ntl_polynomials(const NTL::ZZ_pX& a, const NTL::ZZ_pX& b) {
    int a_degree = NTL::deg(a);
    int b_degree = NTL::deg(b);

    if (a_degree != b_degree) {
        return a_degree < b_degree;
    }

    for (int i = a_degree; i >= 0; --i) {
        if (a[i] != b[i]) {
//             std::string a_str = ntl_to_string(a[i]);
//             std::string b_str = ntl_to_string(b[i]);
            return 0;
        }
    }

    return false;
}





