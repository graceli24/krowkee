// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other
// krowkee Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

#include <krowkee/hash/util.hpp>
// #include "xxhash.h"

#if __has_include(<cereal/types/base_class.hpp>)
#include <cereal/types/base_class.hpp>
#include <cereal/access.hpp>
#endif

#include <cstdint>
#include <random>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>


//Serialization for boost uint128
#ifndef CEREAL_TYPES_BOOST_UINT128_
#define CEREAL_TYPES_BOOST_UINT128_

namespace cereal{

using uint128_t = boost::multiprecision::uint128_t;

//Documentation on custom cereal serialization functions: https://uscilab.github.io/cereal/serialization_functions.html

template<class Archive> inline
void save(Archive & archive, uint128_t const & x){ 

  uint128_t leading_bits128 = x >> 64;
  uint128_t limit64 = 18446744073709551615U;
  uint128_t ending_bits128 = x & limit64;

  std::uint64_t leading_bits = static_cast<std::uint64_t> (leading_bits128);
  std::uint64_t ending_bits =  static_cast<std::uint64_t> (ending_bits128);

  archive(leading_bits, ending_bits); 

}

template<class Archive> inline
void load(Archive & archive, uint128_t & x){

  std::uint64_t leading_bits;
  std::uint64_t ending_bits;
  archive(leading_bits, ending_bits); 
  x = (static_cast<uint128_t> (leading_bits) << 64) + (static_cast<uint128_t> (ending_bits));

}

//cereal no longer has any ambiguity about which functions to use for serializing boost::multiprecision::uint128_t
//See: https://github.com/USCiLab/cereal/blob/master/include/cereal/specialize.hpp

template <class Archive> 
struct specialize<Archive, boost::multiprecision::uint128_t, cereal::specialization::non_member_load_save> {};
//cereal no longer has any ambiguity about which functions to use for serializing boost::multiprecision::uint128_t

}
// }
#endif //CEREAL_TYPES_BOOST_UINT128_


namespace krowkee {
namespace hash {


using uint128_t = boost::multiprecision::uint128_t;

/**
 * Hash functor base class
 *
 * Right now only supports hashing to a power of 2 range.
 */
struct Base {
  /**
   * @param range the desired range for the hash function. Currently we round up
   * to the nearest power of 2 for `_m'.
   * @param seed the random seed controlling any randomness. Might be ignored.
   */
  Base(const std::uint64_t range, const std::uint64_t seed = default_seed)
      : _m(64 - ceil_log2_64(range)), _seed(seed) {}

  Base() {}

  /**
   * Truncate a hashed 64-bit value to the desired range specified by `_m`.
   *
   * @param val the hashed value to be truncated.
   */
  constexpr std::uint64_t truncate(const std::uint64_t val) const {
    return val >> _m;
  }

  constexpr std::uint64_t get_m() const { return _m; }
  constexpr std::size_t   size() const { return _1u64 << (64 - _m); }
  constexpr std::size_t   seed() const { return _seed; }
  inline std::string      state() const {
    std::stringstream ss;
    ss << "size: " << size() << ", seed: " << seed();
    return ss.str();
  }

  void swap(Base &rhs) {
    std::swap(_m, rhs._m);
    std::swap(_seed, rhs._seed);
  }

#if __has_include(<cereal/types/base_class.hpp>)
  template <class Archive>
  void serialize(Archive &archive) {
    archive(_m, _seed);
  }
#endif

  friend void swap(Base &lhs, Base &rhs) { lhs.swap(rhs); }

 protected:
  std::uint64_t              _m;
  std::uint64_t              _seed;
  const static std::uint64_t _1u64 = 1;
};

/**
 * Thomas Wang hash functor
 *
 * https://naml.us/blog/tag/thomas-wang
 */
struct WangHash : public Base {
  /**
   * @tparam ARGS types of additional (ignored) hash functor arguments.
   * @param range the power of 2 range for the hash function.
   */
  template <typename... ARGS>
  WangHash(std::uint64_t range, ARGS &...) : Base(range) {}

  WangHash() : Base() {}

  template <typename OBJ>
  constexpr std::uint64_t operator()(const OBJ &x) const {
    return truncate(wang64(std::uint64_t(x)));
  }

  /**
   * Print functor name.
   */
  static inline std::string name() { return "WangHash"; }

  inline std::string state() const {
    std::stringstream ss;
    ss << "size: " << size();
    return ss.str();
  }

  friend constexpr bool operator==(const WangHash &lhs, const WangHash &rhs) {
    return lhs._m == rhs._m && lhs._seed == rhs._seed;
  }

  friend constexpr bool operator!=(const WangHash &lhs, const WangHash &rhs) {
    return !operator==(lhs, rhs);
  }
};

// /**
//  * XXHash functor
//  *
//  * https://github.com/Cyan4973/xxHash
//  */
// struct XXHash : public Base {
//   template <typename... ARGS>
//   XXHash(const std::uint64_t m, const std::uint64_t seed = default_seed,
//          const ARGS &&...)
//       : Base(m, std::uint32_t(seed)) {}

//   template <typename OBJ>
//   constexpr std::uint64_t operator()(const OBJ &x) const {
//     std::cout << "gets here" << std::endl;
//     return xxh::xxhash<64>(&x, 1, _seed);
//     // return std::uint64_t(XXH64(&x, sizeof(x), _seed));
//   }

//   inline std::string name() { return "XXHash"; }
// };

/**
 * multiply-shift hash
 *
 * https://en.wikipedia.org/wiki/Universal_hashing
 */
struct MulShift : public Base {
  /**
   * Initialize `_a` given the specified random `seed`.
   *
   * `_a` is randomly sampled from [0, 2^64 - 1].
   *
   * @tparam ARGS types of additional (ignored) hash functor arguments.
   *
   * @param range the range for the hash function.
   * @param seed the random seed.
   */
  template <typename... ARGS>
  MulShift(const std::uint64_t range, const std::uint64_t seed = default_seed,
           const ARGS &...)
      : Base(range, seed) {
    std::mt19937_64                              rnd_gen(wang64(_seed));
    std::uniform_int_distribution<std::uint64_t> udist(
        1, std::numeric_limits<std::uint64_t>::max());
    std::uint64_t a;
    do {
      a = udist(rnd_gen);
    } while (is_even(a));
    _a = a;
  }

  MulShift() : Base() {}

  /**
   * Compute `_a * x mod _m`.
   *
   * @tparam OBJ the object to be hashed. Presently must fit into a 64-bit
   *     register.
   *
   * @param x the obejct to be hashed.
   */
  template <typename OBJ>
  constexpr std::uint64_t operator()(const OBJ &x) const {
    return truncate(_a * x);
  }

  /**
   * Print functor name.
   */
  static inline std::string name() { return "MulShift"; }

  inline std::string state() const {
    std::stringstream ss;
    ss << Base::state() << ", multiplicand: " << _a;
    return ss.str();
  }

  friend void swap(MulShift &lhs, MulShift &rhs) {
    std::swap(lhs._a, rhs._a);
    lhs.swap(rhs);
  }

  friend constexpr bool operator==(const MulShift &lhs, const MulShift &rhs) {
    return lhs._m == rhs._m && lhs._seed == rhs._seed && lhs._a == rhs._a;
  }

  friend constexpr bool operator!=(const MulShift &lhs, const MulShift &rhs) {
    return !operator==(lhs, rhs);
  }

#if __has_include(<cereal/types/base_class.hpp>)
  template <class Archive>
  void serialize(Archive &archive) {
    archive(cereal::base_class<Base>(this), _a);
  }
#endif

 private:
  std::uint64_t _a;
};

/**
 * multiply-add-shift hash
 *
 * https://en.wikipedia.org/wiki/Universal_hashing
 */
struct MulAddShift : public Base {
  /**
   * Initialize `_a` and `_b` given the specified random `seed`.
   *
   * `_a` and `_b` are randomly sampled from [0, 2^64 - 1] and
   * [0, 2^size() - 1], respectively.
   *
   * @tparam ARGS types of additional (ignored) hash functor arguments.
   *
   * @param range the power of 2 range for the hash function.
   * @param seed the random seed.
   */
  template <typename... ARGS>
  MulAddShift(const std::uint64_t range,
              const std::uint64_t seed = default_seed, const ARGS &...)
      : Base(range, seed) {
    std::mt19937_64                              rnd_gen(wang64(_seed));
    std::uniform_int_distribution<std::uint64_t> udist_a(
        0, std::numeric_limits<std::uint64_t>::max());
    std::uniform_int_distribution<std::uint64_t> udist_b(0, size());
    std::uint64_t                                a;
    do {
      a = udist_a(rnd_gen);
    } while (is_even(a));
    _a = a;
    _b = udist_b(rnd_gen);
  }

  MulAddShift() : Base() {}

  /**
   * Compute `_a * x + _b mod _m`.
   *
   * @tparam OBJ the object to be hashed. Presently must fit into a 64-bit
   *     register.
   *
   * @param x the obejct to be hashed.
   */
  template <typename OBJ>
  constexpr std::uint64_t operator()(const OBJ &x) const {
    return truncate(_a * x + _b);
  }

  /**
   * Print functor name.
   */
  static inline std::string name() { return "MulAddShift"; }

  inline std::string state() const {
    std::stringstream ss;
    ss << Base::state() << ", multiplicand: " << _a << ",  summand: " << _b;
    return ss.str();
  }

  friend void swap(MulAddShift &lhs, MulAddShift &rhs) {
    std::swap(lhs._a, rhs._a);
    std::swap(lhs._b, rhs._b);
    lhs.swap(rhs);
  }

  friend constexpr bool operator==(const MulAddShift &lhs,
                                   const MulAddShift &rhs) {
    return lhs._m == rhs._m && lhs._seed == rhs._seed && lhs._a == rhs._a &&
           lhs._b == rhs._b;
  }

  friend constexpr bool operator!=(const MulAddShift &lhs,
                                   const MulAddShift &rhs) {
    return !operator==(lhs, rhs);
  }

#if __has_include(<cereal/types/base_class.hpp>)
  template <class Archive>
  void serialize(Archive &archive) {
    archive(cereal::base_class<Base>(this), _a, _b);
  }
#endif

 private:
  std::uint64_t _a;
  std::uint64_t _b;
};



/**
 * Polynomial hash with a Mersenne prime to generate k-universal hash functions
 *
 * Algorithm 1 from: https://arxiv.org/abs/2008.08654
 */
struct kPolynomialMersenne : public Base {
  /**
   * 
   * Calculate h(x) = ([\sum_{i=0}^{k-1} a_i x^i] mod prime) mod range
   *
   * Initialize the polynomial coefficients a_i = `_coefficients[i]` 
   * given the specified random `seed`.
   * Coefficients are randomly sampled from [0, prime].
   *
   * @tparam ARGS types of additional (ignored) hash functor arguments.
   *
   * @param range: the power of 2 range for the hash function.
   * @param seed: the random seed.
   * @param k: the k universality of the family of hash functions. If unspecified,
               the default is k = 4
   * @param mersenne_power: the power for the mersenne prime (2^mersenne_power - 1) 
            we use for the prime mod.
   */
  template <typename... ARGS>
  kPolynomialMersenne(const std::uint64_t range,
                       const std::uint64_t seed = default_seed, 
                       const std::uint16_t k = 4,
                       const std::uint16_t mersenne_power = 89,
                       const ARGS &...)
      : Base(range, seed) {

    _k = k;
    _mersenne_power = mersenne_power;
    uint128_t base = 2;
    _prime = boost::multiprecision::pow(base, mersenne_power) - 1;

    boost::random::mt19937_64                             rnd_gen(wang64(_seed));
    boost::random::uniform_int_distribution<uint128_t>    udist(0, _prime);
    
    std::vector<uint128_t> coefficients (k, 0);
    for (uint16_t i = 0; i < k; i++){
      coefficients[i] = udist(rnd_gen);
    }
    _coefficients = coefficients;
  }

  kPolynomialMersenne() : Base() {}

  /**
   * Compute ([\sum_{i=0}^{k-1} a[i] x^i] mod prime) mod range
   * Algorithm 1 in https://arxiv.org/abs/2008.08654 
   *
   * @tparam OBJ the object to be hashed. Presently must fit into a 64-bit
   *     register.
   *
   * @param x the object to be hashed.
   */
  template <typename OBJ>
  constexpr std::uint64_t operator()(const OBJ &x) const {

    //Evaluate the polynomial mod p
    uint128_t y = _coefficients[_k-1];
    for (int i = _k - 2; i >= 0; i--){
      y = y * x + _coefficients[i];
      y = (y & _prime) + (y >> _mersenne_power); //calculate y mod prime
    }

    //Make sure y is in [0, prime - 1]
    while (y >= _prime){
      y = y - _prime;
    }

    //Return y mod range, where we assume range = 2^\ell
    uint64_t range = std::pow(2,(64 - _m)); //from base class
    y = y & (range - 1);

    return static_cast<uint64_t>(y);

  }

  /**
   * Print functor name.
   */
  static inline std::string name() { return "kPolynomialMersenne"; }

  inline std::string state() const {
    std::stringstream ss;
    ss << Base::state() << ", coefficients ((k-1)th power to constant): "; 
    for (int i = _k - 1; i >= 0; i--){
      ss << _coefficients[i] << " ";
    }
    return ss.str();
  }

  friend void swap(kPolynomialMersenne &lhs, kPolynomialMersenne &rhs) {
    std::swap(lhs._k, rhs._k);
    std::swap(lhs._mersenne_power, rhs._mersenne_power);
    std::swap(lhs._prime, rhs._prime);
    std::swap(lhs._coefficients, rhs._coefficients);
    lhs.swap(rhs);
  }

  friend constexpr bool operator==(const kPolynomialMersenne &lhs,
                                   const kPolynomialMersenne &rhs) {
    return lhs._m == rhs._m && lhs._seed == rhs._seed && lhs._k == rhs._k &&
           lhs._mersenne_power == rhs._mersenne_power &&
           lhs._prime == rhs._prime &&
           lhs._coefficients == rhs._coefficients;
  }

  friend constexpr bool operator!=(const kPolynomialMersenne &lhs,
                                   const kPolynomialMersenne &rhs) {
    return !operator==(lhs, rhs);
  }

#if __has_include(<cereal/types/base_class.hpp>)
  template <class Archive>
  void serialize(Archive &archive) {
    archive(cereal::base_class<Base>(this), _k, _mersenne_power, _prime, _coefficients);
  }
#endif

 private:
  std::uint16_t _k;
  std::uint16_t _mersenne_power;
  uint128_t _prime;
  std::vector<uint128_t> _coefficients;
};


/**
 * Polynomial hash with a Mersenne prime to generate k-universal hash functions
 * Uses uint64's only and supports inputs in the range [0, 2^31 - 1]
 * Algorithm 1 from: https://arxiv.org/abs/2008.08654
 */
struct kPolynomialMersenne_uint64 : public Base {
  /**
   * 
   * Calculate h(x) = ([\sum_{i=0}^{k-1} a_i x^i] mod prime) mod range
   *
   * Initialize the polynomial coefficients a_i = `_coefficients[i]` 
   * given the specified random `seed`.
   * Coefficients are randomly sampled from [0, prime].
   *
   * @tparam ARGS types of additional (ignored) hash functor arguments.
   *
   * @param range: the power of 2 range for the hash function.
   * @param seed: the random seed.
   * @param k: the k universality of the family of hash functions. If unspecified,
               the default is k = 4
   * @param mersenne_power: the power for the mersenne prime (2^mersenne_power - 1) 
            we use for the prime mod.
   */
  template <typename... ARGS>
  kPolynomialMersenne_uint64(const std::uint64_t range,
                             const std::uint64_t seed = default_seed, 
                             const std::uint16_t k = 4,
                             const std::uint16_t mersenne_power = 31,
                             const ARGS &...)
      : Base(range, seed) {

    _k = k;
    _mersenne_power = mersenne_power;
    _prime = std::pow(2,mersenne_power) - 1;
    
    std::mt19937_64                              rnd_gen(wang64(_seed));
    std::uniform_int_distribution<std::uint64_t> udist(0, _prime);

    std::vector<uint64_t> coefficients (k, 0);
    for (uint16_t i = 0; i < k; i++){
      coefficients[i] = udist(rnd_gen);
    }
    _coefficients = coefficients;

  }

  kPolynomialMersenne_uint64() : Base() {}

  /**
   * Compute ([\sum_{i=0}^{k-1} a[i] x^i] mod prime) mod range
   * Algorithm 1 in https://arxiv.org/abs/2008.08654 
   *
   * @tparam OBJ the object to be hashed. Presently must fit into a 64-bit
   *     register.
   *
   * @param x the object to be hashed.
   */
  template <typename OBJ>
  constexpr std::uint64_t operator()(const OBJ &x) const {

    //Evaluate the polynomial mod p
    std::uint64_t y = _coefficients[_k-1];
    for (int i = _k - 2; i >= 0; i--){
      y = y * x + _coefficients[i];
      y = (y & _prime) + (y >> _mersenne_power); //calculate y mod prime
    }

    //Make sure y is in [0, prime - 1]
    while (y >= _prime){
      y = y - _prime;
    }

    //Return y mod range, where we assume range = 2^\ell
    uint64_t range = std::pow(2,(64 - _m)); //from base class
    y = y & (range - 1);

    return y;

    // return truncate(_a * x + _b);
  }

  /**
   * Print functor name.
   */
  static inline std::string name() { return "kPolynomialMersenne_uint64"; }

  inline std::string state() const {
    std::stringstream ss;
    ss << Base::state() << ", coefficients ((k-1)th power to constant): "; 
    for (int i = _k - 1; i >= 0; i--){
      ss << _coefficients[i] << " ";
    }
    return ss.str();
  }

  friend void swap(kPolynomialMersenne_uint64 &lhs, kPolynomialMersenne_uint64 &rhs) {
    std::swap(lhs._k, rhs._k);
    std::swap(lhs._mersenne_power, rhs._mersenne_power);
    std::swap(lhs._prime, rhs._prime);
    std::swap(lhs._coefficients, rhs._coefficients);
    lhs.swap(rhs);
  }

  friend constexpr bool operator==(const kPolynomialMersenne_uint64 &lhs,
                                   const kPolynomialMersenne_uint64 &rhs) {
    return lhs._m == rhs._m && lhs._seed == rhs._seed && lhs._k == rhs._k &&
           lhs._mersenne_power == rhs._mersenne_power &&
           lhs._prime == rhs._prime &&
           lhs._coefficients == rhs._coefficients;
  }

  friend constexpr bool operator!=(const kPolynomialMersenne_uint64 &lhs,
                                   const kPolynomialMersenne_uint64 &rhs) {
    return !operator==(lhs, rhs);
  }

#if __has_include(<cereal/types/base_class.hpp>)
  template <class Archive>
  void serialize(Archive &archive) {
    archive(cereal::base_class<Base>(this), _k, _mersenne_power, _prime, _coefficients);
  }
#endif

 private:
  std::uint16_t _k;
  std::uint64_t _mersenne_power;
  std::uint64_t _prime;
  std::vector<uint64_t> _coefficients;
};


std::ostream &operator<<(std::ostream &os, const Base &func) {
  os << func.state();
  return os;
}

}  // namespace hash
}  // namespace krowkee