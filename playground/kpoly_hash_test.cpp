#include <iostream>

#include <krowkee/hash/hash.hpp>

typedef boost::multiprecision::uint128_t uint128_t;

typedef krowkee::hash::kPolynomialMersenne_uint64 kpoly_uint64_t;
typedef krowkee::hash::kPolynomialMersenne kpoly_t;

int main(int argc, char** argv) {

    std::uint64_t range = 2;
    std::uint64_t seed = krowkee::hash::default_seed;

    std::uint16_t k = 4;
    std::uint16_t mersenne_power = 31;
    uint128_t base = 2;
    uint128_t prime = boost::multiprecision::pow(base, mersenne_power) - 1;
    std::cout << "p = 2^" << mersenne_power << " - 1 = " << prime << std::endl << std::endl;

    //Try implementation in hash.hpp

    kpoly_t hash{range, seed, k, mersenne_power};

    std::cout << "[" << hash.state() << "]" << std::endl;

    for (int x = 0; x < 10; x++){
        std::cout << "hash(" << x << ") = " << hash(x) << std::endl;
    }
    std::cout << std::endl;

    //Try seeing what's going on from scratch

    boost::random::mt19937_64                             rnd_gen(krowkee::hash::wang64(seed));
    boost::random::uniform_int_distribution<uint128_t>    udist(0, prime);

    std::vector<uint128_t> coefficients (k, 0);
    for (uint16_t i = 0; i < k; i++){
      coefficients[i] = udist(rnd_gen);
    }

    // std::cout << "From scratch. Coefficients ((k-1)th power to constant): "; 
    // for (int i = k - 1; i >= 0; i--){
    //   std::cout << coefficients[i] << " ";
    // }
    // std::cout << std::endl;


    // for (int x = 0; x < 3; x++){
    //     std::cout << std::endl << "Evaluating x = " << x << std::endl;

    //     //Evaluate the polynomial mod p
    //     uint128_t y = coefficients[k-1];
    //     std::cout << "k = " << k-1 << ", y = " << y << std::endl;
    //     for (int i = k - 2; i >= 0; i--){
    //         y = y * x + coefficients[i];
    //         std::cout << "k = " << i << ", y = " << y << std::endl;
    //         y = (y & prime) + (y >> mersenne_power); //calculate y mod prime
    //         std::cout << "y = " << y << " mod p" << std::endl;
    //     }

    //     //Make sure y is in [0, prime - 1]
    //     while (y >= prime){
    //         y = y - prime;
    //     }
    //     std::cout << "After subtracting, y = " << y << std::endl;

    //     //Return y mod range, where we assume range = 2^\ell
    //     y = y & (range - 1);
    //     std::cout << "After bitwise mod range, y = " << y << std::endl;

    //     uint64_t y_final = static_cast<uint64_t>(y);
    //     std::cout << "Final hash(x) = " << y_final << std::endl;

    // }


    return 0;
}