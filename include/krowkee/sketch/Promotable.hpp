// Copyright 2021-2022 Lawrence Livermore National Security, LLC and other
// krowkee Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: MIT

#ifndef _KROWKEE_SKETCH_PROMOTABLE_HPP
#define _KROWKEE_SKETCH_PROMOTABLE_HPP

#include <krowkee/sketch/Dense.hpp>
#include <krowkee/sketch/Sparse.hpp>

#if __has_include(<cereal/types/memory.hpp>)
#include <cereal/types/memory.hpp>
#endif

#include <memory>

namespace krowkee {
namespace sketch {

enum class promotable_mode_t : std::uint8_t { sparse, dense };

/**
 * @brief Sparse to Dense data structure for Sketch objects.
 *
 * Provides a unified interface to the Sparse and Dense containers, allowing for
 * seamless "promotion" of a Sparse representation to a Dense one.
 *
 * @tparam RegType The type held by each register.
 * @tparam MergeOp An element-wise merge operator to combine two sketches,
 * templated on RegType.
 * @tparam MapType The underlying map class to use for buffered updates to the
 * underlying compacting_map.
 * @tparam KeyType The key type to use for this underlying compacting map.
 */
template <typename RegType, typename MergeOp,
          template <typename, typename> class MapType, typename KeyType>
class Promotable {
 public:
  /** Alias for the fully-templated Dense container. */
  typedef krowkee::sketch::Dense<RegType, MergeOp> dense_t;
  typedef std::unique_ptr<dense_t>                 dense_ptr_t;
  /** Alias for the fully-templated Sparse container. */
  typedef krowkee::sketch::Sparse<RegType, MergeOp, MapType, KeyType> sparse_t;
  typedef typename sparse_t::map_t                                    map_t;
  typedef std::unique_ptr<sparse_t> sparse_ptr_t;
  /** Alias for the fully-templated Promotable container. */
  typedef Promotable<RegType, MergeOp, MapType, KeyType> promotable_t;

 private:
  dense_ptr_t       _dense_ptr;
  sparse_ptr_t      _sparse_ptr;
  std::size_t       _promotion_threshold;
  promotable_mode_t _mode;

 public:
  /**
   * @brief Initialize the size parameters and create the base container.
   *
   * @param range_size Size argument for constructing dense container.
   * @param compaction_threshold Size argument for constructing sparse
   *     container.
   * @param promotion_threshold Threshold from promoting from Sparse
   * representation to Dense representation.
   */
  Promotable(const std::size_t range_size,
             const std::size_t compaction_threshold,
             const std::size_t promotion_threshold)
      : _promotion_threshold(promotion_threshold),
        _mode(promotable_mode_t::sparse),
        _sparse_ptr(
            std::make_unique<sparse_t>(range_size, compaction_threshold)) {}

  /**
   * @brief Copy constructor.
   *
   * @note Must explicitly copy contents of rhs's container, not just copy the
   * pointer.
   *
   * @param rhs container to be copied.
   */
  Promotable(const promotable_t &rhs)
      : _promotion_threshold(rhs._promotion_threshold), _mode(rhs._mode) {
    if (_mode == promotable_mode_t::sparse) {
      sparse_ptr_t sparse_ptr(std::make_unique<sparse_t>(*rhs._sparse_ptr));
      std::swap(sparse_ptr, _sparse_ptr);
    } else {
      dense_ptr_t dense_ptr(std::make_unique<dense_t>(*rhs._dense_ptr));
      std::swap(dense_ptr, _dense_ptr);
    }
  }

  /**
   * @brief Default constructor
   *
   * @note Only use for move constructor.
   */
  Promotable()
      : _mode(promotable_mode_t::sparse), _sparse_ptr(), _dense_ptr() {}

  // /**
  //  * move constructor
  //  *
  //  * @param rhs r-value promotable_t to be destructively copied.
  //  */
  // Promotable(promotable_t &&rhs) noexcept : Promotable() {
  //   std::swap(*this, rhs);
  // }

  //////////////////////////////////////////////////////////////////////////////
  // Swaps
  //////////////////////////////////////////////////////////////////////////////
  /**
   * @brief Swap two Promotable objects.
   *
   * @param lhs The left-hand container.
   * @param rhs The right-hand containr.
   */
  friend void swap(promotable_t &lhs, promotable_t &rhs) {
    std::swap(lhs._promotion_threshold, rhs._promotion_threshold);
    std::swap(lhs._mode, rhs._mode);
    if (lhs._mode == promotable_mode_t::sparse) {
      std::swap(lhs._sparse_ptr, rhs._sparse_ptr);
    } else {
      std::swap(lhs._dense_ptr, rhs._dense_ptr);
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Cereal Archives
  //////////////////////////////////////////////////////////////////////////////

#if __has_include(<cereal/types/memory.hpp>)
  /**
   * @brief Save Promotable object to `cereal` archive.
   *
   * @tparam Archive `cereal` archive type.
   * @param archive The `cereal` archive to which to serialize the sketch.
   */
  template <class Archive>
  void save(Archive &oarchive) const {
    oarchive(_promotion_threshold, _mode);
    if (_mode == promotable_mode_t::sparse) {
      oarchive(_sparse_ptr);
    } else {
      oarchive(_dense_ptr);
    }
  }

  /**
   * @brief Load Promotable object to `cereal` archive.
   *
   * @tparam Archive `cereal` archive type.
   * @param archive The `cereal` archive from which to serialize the sketch.
   */
  template <class Archive>
  void load(Archive &iarchive) {
    iarchive(_promotion_threshold, _mode);
    if (_mode == promotable_mode_t::sparse) {
      iarchive(_sparse_ptr);
    } else {
      iarchive(_dense_ptr);
    }
  }
#endif

  //////////////////////////////////////////////////////////////////////////////
  // Assignment
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Copy-and-swap assignment operator
   *
   * @note
   * https://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
   *
   * @param rhs The other container.
   * @return sparse_t& `this`, having been swapped with `rhs`.
   */
  promotable_t &operator=(promotable_t rhs) {
    std::swap(*this, rhs);
    return *this;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Getters
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Returns a description of the type of container.
   *
   * @return std::string "Promotable"
   */
  static inline std::string name() { return "Promotable"; }

  /**
   * @brief Returns a description of the fully-qualified type of container
   *
   * @return std::string "Promotable" plus the full name of the fully-templated
   * compacting_map type.
   */
  static inline std::string full_name() {
    std::stringstream ss;
    ss << name() << " using " << krowkee::hash::type_name<map_t>();
    return ss.str();
  }

  /** The size of the underlying Sparse or Dense container, depending on the
   * mode. */
  constexpr std::size_t size() const {
    return (_mode == promotable_mode_t::sparse) ? _sparse_ptr->size()
                                                : _dense_ptr->size();
  }

  /** True unless we are in Sparse mode and uncompacted. */
  constexpr bool is_compact() const {
    return (_mode == promotable_mode_t::sparse) ? _sparse_ptr->is_compact()
                                                : true;
  }

  /** Return `true` if in sparse mode. */
  constexpr bool is_sparse() const {
    return (_mode == promotable_mode_t::sparse) ? true : false;
  }

  /** The number of bytes used by each register. */
  constexpr std::size_t reg_size() const { return sizeof(RegType); }

  /** The size at which we switch from Sparse mode to Dense mode. */
  constexpr std::size_t promotion_threshold() const {
    return _promotion_threshold;
  }

  /** Get the current mode. */
  constexpr promotable_mode_t mode() const { return _mode; }

  /** The size at which the compaction buffer will flush. */
  constexpr std::size_t compaction_threshold() const {
    return (_mode == promotable_mode_t::sparse)
               ? _sparse_ptr->compaction_threshold()
               : _dense_ptr->compaction_threshold();
  }

  /**
   * @brief Get a copy of the raw register vector.
   *
   * The vector will be sparse, with zeros for unset indices, if we are in
   * Sparse mode.
   *
   * @return std::vector<RegType> Copy of the register vector.
   */
  std::vector<RegType> register_vector() const {
    return (_mode == promotable_mode_t::sparse) ? _sparse_ptr->register_vector()
                                                : _dense_ptr->register_vector();
  }

  //////////////////////////////////////////////////////////////////////////////
  // Equality operators
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Checks whether another Promotable container has the same promotion
   * threshold.
   *
   * @param rhs The other container.
   * @return true The thresholds agree.
   * @return false The thresholds disagree.
   */
  constexpr bool same_parameters(const promotable_t &rhs) const {
    return _promotion_threshold == rhs._promotion_threshold;
  }

  /**
   * @brief Checks equality between two Promotable containers
   *
   * @note[bwp] do we want to be able to compare sparse to dense containers
   * directly? no for now but will want to revisit.
   *
   * @param lhs The left-hand container.
   * @param rhs The right-hand container.
   * @return true The registers agree and promotion thresholds are the same.
   * @return false At least one register disagrees or the promotion thresholds
   * disagree.
   */
  friend constexpr bool operator==(const promotable_t &lhs,
                                   const promotable_t &rhs) {
    if (lhs.same_parameters(rhs) == false || lhs._mode != rhs._mode) {
      return false;
    }
    if (lhs._mode == promotable_mode_t::sparse) {
      return *lhs._sparse_ptr == *rhs._sparse_ptr;
    } else {
      return *lhs._dense_ptr == *rhs._dense_ptr;
    }
  }

  /**
   * @brief Checks inequality between two Promotable containers
   *
   * @param lhs The left-hand container.
   * @param rhs The right-hand container.
   * @return true At least one register disagrees or the promotion thresholds
   * disagree.
   * @return false The registers agree and promotion thresholds are the same.
   */
  friend constexpr bool operator!=(const promotable_t &lhs,
                                   const promotable_t &rhs) {
    return !operator==(lhs, rhs);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Compaction
  //////////////////////////////////////////////////////////////////////////////

  /** If in Sparse mode, flush the update buffer. */
  void compactify() {
    if (_mode == promotable_mode_t::sparse) {
      _sparse_ptr->compactify();
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Erase
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Trigger the underlying container to erase an index.
   *
   * @param index The index to erase.
   */
  inline void erase(const std::uint64_t index) {
    if (_mode == promotable_mode_t::sparse) {
      _sparse_ptr->erase(index);
    } else {
      _dense_ptr->erase(index);
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Register iterators
  //////////////////////////////////////////////////////////////////////////////

  /** Mutable begin iterator. */
  constexpr auto begin() {
    if (_mode == promotable_mode_t::sparse) {
      return std::begin(*_sparse_ptr);
    } else {
      return std::begin(*_dense_ptr);
    }
  }
  /** Const begin iterator. */
  constexpr auto begin() const {
    if (_mode == promotable_mode_t::sparse) {
      return std::cbegin(*_sparse_ptr);
    } else {
      return std::cbegin(*_dense_ptr);
    }
  }
  /** Const begin iterator. */
  constexpr auto cbegin() const {
    if (_mode == promotable_mode_t::sparse) {
      return std::cbegin(*_sparse_ptr);
    } else {
      return std::cbegin(*_dense_ptr);
    }
  }
  /** Mutable end iterator. */
  constexpr auto end() {
    if (_mode == promotable_mode_t::sparse) {
      return std::end(*_sparse_ptr);
    } else {
      return std::end(*_dense_ptr);
    }
  }
  /** Const end iterator. */
  constexpr auto end() const {
    if (_mode == promotable_mode_t::sparse) {
      return std::cend(*_sparse_ptr);
    } else {
      return std::cend(*_dense_ptr);
    }
  }
  /** Const end iterator. */
  constexpr auto cend() const {
    if (_mode == promotable_mode_t::sparse) {
      return std::cend(*_sparse_ptr);
    } else {
      return std::cend(*_dense_ptr);
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Accessors
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Access reference to the specified element.
   *
   * @param index The index to access.
   * @return RegType& The accessed register.
   */
  RegType &operator[](const std::uint64_t index) {
    if (_mode == promotable_mode_t::sparse) {
      if (size() == _promotion_threshold) {
        compactify();
        promote();
        return (*_dense_ptr)[index];
      }
      return (*_sparse_ptr)[index];
    } else {
      return (*_dense_ptr)[index];
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Merge Operators
  //////////////////////////////////////////////////////////////////////////////
  /**
   * @brief Merge other Promotable into `this`.
   *
   * @param rhs The other promotable. Must have the same parameters.
   * @return promotable_t& `this` after merging rhs.
   */
  promotable_t &operator+=(const promotable_t &rhs) {
    if (same_parameters(rhs) == false) {
      throw std::invalid_argument(
          "containers do not have congruent parameters!");
    }
    if (_mode == rhs._mode) {
      if (_mode == promotable_mode_t::sparse) {
        // sparse_t tmp(*rhs._sparse_ptr);
        // *_sparse_ptr += tmp;
        *_sparse_ptr += *rhs._sparse_ptr;
        if (size() >= _promotion_threshold) {
          promote();
        }
      } else {
        *_dense_ptr += *rhs._dense_ptr;
      }
    } else {
      if (_mode == promotable_mode_t::sparse) {
        // we are sparse; promote to dense and add rhs._dense_ptr
        // This will be slow; try to avoid it.
        promote();
        *_dense_ptr += *rhs._dense_ptr;
      } else {
        // we are dense; add rhs._sparse_ptr
        merge_from_sparse(rhs);
      }
    }
    return *this;
  }

  /**
   * @brief Merget two Promotable objects together.
   *
   * lhs and rhs must have the same parameters.
   *
   * @param lhs The left-hand container.
   * @param rhs The right-hand container.
   * @return promotable_t The merge of the two containers.
   */
  inline friend promotable_t operator+(const promotable_t &lhs,
                                       const promotable_t &rhs) {
    if (rhs._mode == promotable_mode_t::dense) {
      promotable_t ret(rhs);
      ret += lhs;
      return ret;
    }
    promotable_t ret(lhs);
    ret += rhs;
    return ret;
  }

  /**
   * Promote sparse container into a dense container.
   *
   * Initializes dense pointer while releasing sparse pointer.
   *
   * @throws std::logic_error if the container is not in sparse mode.
   * @throws std::logic_error if the sparse container is not compacted.
   */

  /**
   * @brief Promote sparse container into a dense container.
   *
   * Initializes dense pointer while releasing sparse pointer.
   *
   * @throws std::logic_error if the container is not in sparse mode.
   * @throws std::logic_error if the sparse container is not compacted.
   */
  void promote() {
    if (_mode != promotable_mode_t::sparse) {
      throw std::logic_error("Attempt to promote non-sparse container!");
    }
    if (_sparse_ptr->is_compact() == false) {
      throw std::logic_error("Attempt to promote uncompacted container!");
    }

    _dense_ptr = std::make_unique<dense_t>(_sparse_ptr->range_size());
    merge_from_sparse(*this);
    _sparse_ptr.reset();
    _mode = promotable_mode_t::dense;
  }

  /**
   * @brief Incorporate the information in a sparse container into a dense
   * container.
   *
   * Initializes dense pointer while releasing sparse pointer.
   *
   * @param rhs The sparse container.
   * @throws std::logic_error if the rhs container is not in sparse mode.
   */
  void merge_from_sparse(const promotable_t &rhs) {
    if (rhs._mode != promotable_mode_t::sparse) {
      throw std::logic_error("Attempt to dense merge a non-sparse rhs!");
    }
    if (rhs.is_compact() == false) {
      throw std::logic_error("Attempt to dense merge a non-compact rhs!");
    }
    for_each(*rhs._sparse_ptr, [&](const auto &p) {
      RegType &val((*_dense_ptr)[p.first]);
      val = MergeOp()(val, p.second);
    });
  }

  //////////////////////////////////////////////////////////////////////////////
  // I/O Operators
  //////////////////////////////////////////////////////////////////////////////
  /**
   * @brief Serialize a Sparse container to human-readable output stream.
   *
   * Output format is a list of (key, value) pairs, skipping empty registers if
   * in Sparse mode.
   *
   * @note Intended for debugging only.
   *
   * @param os The output stream.
   * @param con The Promotable to print.
   * @return std::ostream& The new stream state.
   * @throw std::logic_error if invoked on uncompacted sketch.
   */
  friend std::ostream &operator<<(std::ostream &os, const promotable_t &con) {
    if (con._mode == promotable_mode_t::sparse) {
      os << *con._sparse_ptr;
    } else {
      os << *con._dense_ptr;
    }
    return os;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Accumulation
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Accumulate sum of register values.
   *
   * @tparam RetType The value type to return.
   * @param con The Promotable to accumulate.
   * @param init Initial value of return.
   * @return RetType The sum over all register values + `init`.
   */
  template <typename RetType>
  friend RetType accumulate(const promotable_t &con, const RetType init) {
    if (con._mode == promotable_mode_t::sparse) {
      return accumulate(*con._sparse_ptr, init);
    } else {
      return accumulate(*con._dense_ptr, init);
    }
  }
};

}  // namespace sketch
}  // namespace krowkee

#endif