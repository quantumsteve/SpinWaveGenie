//
//  ThreeDimensionalCut.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 7/2/15.
//
//

#ifndef spin_wave_genie_ThreeDimensionalCut_h
#define spin_wave_genie_ThreeDimensionalCut_h

#include <memory>
#include "SpinWaveGenie/Containers/PointsOnAPlane.h"
#include "SpinWaveGenie/Plot/SpinWavePlot.h"

namespace SpinWaveGenie
{
    
    class ThreeDimensionalCut
    {
    public:
        ThreeDimensionalCut();
        ThreeDimensionalCut(const ThreeDimensionalCut &other);
        ThreeDimensionalCut &operator=(const ThreeDimensionalCut &other);
        ThreeDimensionalCut(ThreeDimensionalCut &&other);
        ThreeDimensionalCut &operator=(ThreeDimensionalCut &&other);
        void setFilename(const std::string& name);
        void setPoints(const PointsOnAPlane& pos);
        void setEnergyPoints(double min, double max, std::size_t numberpoints);
        void setPlotObject(std::unique_ptr<SpinWavePlot> object);
        void save();
        Eigen::MatrixXd getMatrix();
        ~ThreeDimensionalCut();
        
    private:
        class CutImpl;
        std::unique_ptr<CutImpl> m_p;
    };
}

#endif
