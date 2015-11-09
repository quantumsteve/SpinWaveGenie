#ifndef __ThreeVectors__
#define __ThreeVectors__

#include <vector>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>
#include <array>

//! Structure of Arrays used for storing vectors with three components.
/*!
 Vectors containing three elements are not ideally coellesced in memory.
 Therefore, we use an alternative where where each component is stored
 in a separate vector and accessed using an Iterator.
 */

namespace SpinWaveGenie
{

template <typename T>
class ThreeVectorsIterator
    : public boost::iterator_facade<ThreeVectorsIterator<T>, std::array<T, 3> const, boost::random_access_traversal_tag>
{
public:
  ThreeVectorsIterator(typename std::vector<T>::iterator _it0, typename std::vector<T>::iterator _it1,
                       typename std::vector<T>::iterator _it2)
      : iteratorArray{{_it0, _it1, _it2}}, values{{*_it0, *_it1, *_it2}}
  {
  }

private:
  friend class boost::iterator_core_access;
  mutable std::array<typename std::vector<T>::iterator, 3> iteratorArray;
  mutable std::array<T, 3> values;
  void increment()
  {
    for (int i = 0; i < 3; ++i)
    {
      iteratorArray[i] += 1;
      values[i] = *iteratorArray[i];
    }
  }

  void decrement()
  {
    for (int i = 0; i < 3; ++i)
    {
      iteratorArray[i] -= 1;
      values[i] = *iteratorArray[i];
    }
  }

  void advance(std::ptrdiff_t n)
  {
    for (int i = 0; i < 3; ++i)
    {
      iteratorArray[i] += n;
      values[i] = *iteratorArray[i];
    }
  }

  std::ptrdiff_t distance_to(ThreeVectorsIterator j) { return std::distance(iteratorArray(0), j.iteratorArray(0)); }

  bool equal(ThreeVectorsIterator const &other) const
  {
    for (int i = 0; i < 3; ++i)
    {
      if (iteratorArray[0] != other.iteratorArray[0])
        return false;
    }
    return true;
  }

  std::array<T, 3> &dereference() const { return values; }
};

template <typename T>
class ConstThreeVectorsIterator : public boost::iterator_facade<ConstThreeVectorsIterator<T>, std::array<T, 3> const,
                                                                boost::random_access_traversal_tag>
{
public:
  ConstThreeVectorsIterator(typename std::vector<T>::iterator _it0, typename std::vector<T>::iterator _it1,
                            typename std::vector<T>::iterator _it2)
      : iteratorArray{{_it0, _it1, _it2}}, values{{*_it0, *_it1, *_it2}}
  {
  }

private:
  friend class boost::iterator_core_access;
  mutable std::array<typename std::vector<T>::iterator, 3> iteratorArray;
  std::array<T, 3> const values;
  void increment()
  {
    for (int i = 0; i < 3; ++i)
    {
      iteratorArray[i] += 1;
      values[i] = *iteratorArray[i];
    }
  }

  void decrement()
  {
    for (int i = 0; i < 3; ++i)
    {
      iteratorArray[i] -= 1;
      values[i] = *iteratorArray[i];
    }
  }

  void advance(std::ptrdiff_t n)
  {
    for (int i = 0; i < 3; ++i)
    {
      iteratorArray[i] += n;
      values[i] = *iteratorArray[i];
    }
  }

  std::ptrdiff_t distance_to(ConstThreeVectorsIterator j)
  {
    return std::distance(iteratorArray(0), j.iteratorArray(0));
  }

  bool equal(ConstThreeVectorsIterator const &other) const
  {
    for (int i = 0; i < 3; ++i)
    {
      if (iteratorArray[0] != other.iteratorArray[0])
        return false;
    }
    return true;
  }

  std::array<T, 3> const &dereference() const { return values; }
};

template <typename T> class ThreeVectors
{
protected:
  // typedef typename std::vector<T>::iterator ValueIterator;
  // typedef typename std::vector<T>::const_iterator ConstValueIterator;

public:
  bool empty();
  //! insert three elements x,y,z
  //! \param x zeroth element of type T
  //! \param y first element of type T
  //! \param z second element of type T
  void insert(T x, T y, T z);
  // typedef boost::zip_iterator<boost::tuple<ValueIterator, ValueIterator, ValueIterator>> Iterator;
  // typedef boost::zip_iterator<boost::tuple<ConstValueIterator, ConstValueIterator, ConstValueIterator>>
  // ConstIterator;
  typedef ThreeVectorsIterator<T> Iterator;
  typedef ConstThreeVectorsIterator<T> ConstIterator;

  // typedef boost::zip_iterator<boost::tuple<ConstValueIterator, ConstValueIterator, ConstValueIterator>>
  // ConstIterator;
  //! \return number of elements in the ThreeVector
  size_t size();
  //! Clears all data stored in the ThreeVector.
  void clear();
  //! \return Returns an Iterator pointing to the first element
  Iterator begin();
  //! \return Returns an Iterator pointing to the end of the vector
  Iterator end();
  //! \return Returns a ConstIterator pointing to the first element
  ConstIterator cbegin() const;
  //! \return Returns an ConstIterator pointing to the end of the vector
  ConstIterator cend() const;

protected:
  std::vector<T> valuesX;
  std::vector<T> valuesY;
  std::vector<T> valuesZ;
};

template <typename T> bool ThreeVectors<T>::empty() { return valuesX.empty(); }

template <typename T> void ThreeVectors<T>::insert(T x, T y, T z)
{
  valuesX.push_back(x);
  valuesY.push_back(y);
  valuesZ.push_back(z);
}

template <typename T>

size_t ThreeVectors<T>::size()
{
  return valuesX.size();
}

template <typename T> typename ThreeVectors<T>::Iterator ThreeVectors<T>::begin()
{
  return boost::make_zip_iterator(boost::make_tuple(valuesX.begin(), valuesY.begin(), valuesZ.begin()));
}

template <typename T> typename ThreeVectors<T>::Iterator ThreeVectors<T>::end()
{
  return boost::make_zip_iterator(boost::make_tuple(valuesX.end(), valuesY.end(), valuesZ.end()));
}

template <typename T> typename ThreeVectors<T>::ConstIterator ThreeVectors<T>::cbegin() const
{
  return boost::make_zip_iterator(boost::make_tuple(valuesX.cbegin(), valuesY.cbegin(), valuesZ.cbegin()));
}

template <typename T> typename ThreeVectors<T>::ConstIterator ThreeVectors<T>::cend() const
{
  return boost::make_zip_iterator(boost::make_tuple(valuesX.cend(), valuesY.cend(), valuesZ.cend()));
}

template <typename T> void ThreeVectors<T>::clear()
{
  valuesX.clear();
  valuesY.clear();
  valuesZ.clear();
}
}
#endif /* defined(__ThreeVectors__) */
