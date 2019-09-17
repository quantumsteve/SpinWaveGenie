//
//  OneDimensionalLorentzian.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//

#ifndef __OneDimensionalLorentzian__
#define __OneDimensionalLorentzian__

#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"
#include <iostream>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT OneDimensionalLorentzian : public OneDimensionalShapes
{
public:
  void setFWHM(double InFWHM);
  void setTolerance(double InTolerance) override;
  double getMinimumEnergy() override;
  double getMaximumEnergy() override;
  void setFrequency(double frequency) override;
  double getFunction(double energy) override;
  std::unique_ptr<OneDimensionalShapes> clone() const override;

private:
  double getExponentialFactor();
  double FWHM, Tolerance;
  double m_frequency;
};
} // namespace SpinWaveGenie
#endif /* defined(__OneDimensionalLorentzian__) */
