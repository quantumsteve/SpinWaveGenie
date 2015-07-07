#ifndef spin_wave_genie_PointsOnAPlane_h
#define spin_wave_genie_PointsOnAPlane_h

#include <iostream>
#include <array>
#include "SpinWaveGenie/Containers/ThreeVectors.h"

namespace SpinWaveGenie
{

//! generates k-points used by the SpinWavePlot routines.
/*!
 This class generates a grid of evenly spaced k-points on a plane.
 The result is stored in a ThreeVectors<double> container.
 */

struct Extent
{
  Extent() : minimumValue(0.0), maximumValue(0.0), numberOfBins(1) {}
  Extent(double min, double max, std::size_t bins) : minimumValue(min), maximumValue(max), numberOfBins(bins) {}
  double minimumValue, maximumValue;
  std::size_t numberOfBins;
};

class PointsOnAPlane
{
public:
  void setProjections(const std::array<double, 3> &proj1, const std::array<double, 3> &proj2);
  void setExtents(const Extent &extent1, const Extent &extent2);
  const std::array<Extent, 2> &getExtents();
  ThreeVectors<double> getPoints();

private:
  std::array<Extent, 2> m_extents;
  std::array<double, 3> m_proj1, m_proj2;
};
}

#endif
