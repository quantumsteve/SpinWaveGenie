//
//  PointsOnAPlane.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 7/1/15.
//
//

#include <stdio.h>
#include "SpinWaveGenie/Containers/PointsOnAPlane.h"

namespace SpinWaveGenie
{

void PointsOnAPlane::setProjections(const std::array<double, 3> &proj1, const std::array<double, 3> &proj2)
{
  m_proj1 = proj1;
  m_proj2 = proj2;
}

void PointsOnAPlane::setExtents(const SpinWaveGenie::Extent &extent1, const SpinWaveGenie::Extent &extent2)
{
  m_extents[0] = extent1;
  m_extents[1] = extent2;
}

const std::array<Extent, 2> &PointsOnAPlane::getExtents() { return m_extents; }

ThreeVectors<double> PointsOnAPlane::getPoints()
{
  double numberPoints1 = m_extents[0].numberOfBins;
  double numberPoints2 = m_extents[1].numberOfBins;

  double H0 = m_extents[0].minimumValue * m_proj1[0] + m_extents[1].minimumValue * m_proj2[0];
  double K0 = m_extents[0].minimumValue * m_proj1[1] + m_extents[1].minimumValue * m_proj2[1];
  double L0 = m_extents[0].minimumValue * m_proj1[2] + m_extents[1].minimumValue * m_proj2[2];

  double H1 = m_extents[0].maximumValue * m_proj1[0] + m_extents[1].minimumValue * m_proj2[0];
  double K1 = m_extents[0].maximumValue * m_proj1[1] + m_extents[1].minimumValue * m_proj2[1];
  double L1 = m_extents[0].maximumValue * m_proj1[2] + m_extents[1].minimumValue * m_proj2[2];

  double H2 = m_extents[0].minimumValue * m_proj1[0] + m_extents[1].maximumValue * m_proj2[0];
  double K2 = m_extents[0].minimumValue * m_proj1[1] + m_extents[1].maximumValue * m_proj2[1];
  double L2 = m_extents[0].minimumValue * m_proj1[2] + m_extents[1].maximumValue * m_proj2[2];

  ThreeVectors<double> Kpoints;
  for (auto m = 0; m < numberPoints2; ++m)
  {
    double Hmi = H0 + (H2 - H0) * static_cast<double>(m) / static_cast<double>(numberPoints2 - 1);
    double Kmi = K0 + (K2 - K0) * static_cast<double>(m) / static_cast<double>(numberPoints2 - 1);
    double Lmi = L0 + (L2 - L0) * static_cast<double>(m) / static_cast<double>(numberPoints2 - 1);
    for (auto n = 0; n < numberPoints1; ++n)
    {
      double Hni = H0 + (H1 - H0) * static_cast<double>(n) / static_cast<double>(numberPoints1 - 1);
      double Kni = K0 + (K1 - K0) * static_cast<double>(n) / static_cast<double>(numberPoints1 - 1);
      double Lni = L0 + (L1 - L0) * static_cast<double>(n) / static_cast<double>(numberPoints1 - 1);
      Kpoints.insert(Hmi + Hni, Kmi + Kni, Lmi + Lni);
    }
  }
  return Kpoints;
}
}