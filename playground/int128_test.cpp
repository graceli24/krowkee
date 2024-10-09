#include <iostream>
#include <string>
#include <typeinfo>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>

int main(int argc, char** argv) {

    //Test uint128_t
    typedef boost::multiprecision::uint128_t uint128_t;

    uint128_t base = 2;
    uint128_t prime = boost::multiprecision::pow(base,89) - 1;

    std::cout << typeid(prime).name() << std::endl;
    std::cout << "Size of prime: " << sizeof(prime) << std::endl;
    std::cout << prime << std::endl << std::endl;

    //uint64_t for comparison
    std::uint64_t small_int = std::pow(2,31) - 1;

    std::cout << typeid(small_int).name() << std::endl;
    std::cout << "Size of small_int: " << sizeof(small_int) << std::endl;
    std::cout << small_int << std::endl << std::endl;

    //Test converting from uint64_t to uint128_t
    uint128_t big_int2 = small_int;

    std::cout << typeid(big_int2).name() << std::endl;
    std::cout << "Size of big_int2: " << sizeof(big_int2) << std::endl;
    std::cout << big_int2 << std::endl << std::endl;

    //Test converting from uint64_t to uint128_t
    std::uint64_t small_int2 = static_cast<uint64_t>(big_int2);

    std::cout << typeid(small_int2).name() << std::endl;
    std::cout << "Size of small_int2: " << sizeof(small_int2) << std::endl;
    std::cout << small_int2 << std::endl << std::endl;

    //Test random number generation
    std::uint64_t seed = 42;
    boost::random::mt19937_64 rnd_gen(seed);
    boost::random::uniform_int_distribution<uint128_t> udist(0, prime);
    // std::uniform_int_distribution<std::uint64_t> udist(0, prime);
    std::cout << "Generating uniformly random numbers" << std::endl;
    std::uint16_t k = 16;
    std::vector<uint128_t> coefficients (k, 0);
    for(uint16_t i = 0; i < k; ++i){
        uint128_t coefficient = udist(rnd_gen);
        coefficients[i] = coefficient;
        std::cout << coefficients[i] << std::endl;
    }

    return 0;
}