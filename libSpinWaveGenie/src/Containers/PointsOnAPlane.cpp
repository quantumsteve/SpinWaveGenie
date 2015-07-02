//
//  PointsOnAPlane.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 7/1/15.
//
//

#include <stdio.h>
#include "SpinWaveGenie/Containers/PointsOnAPlane.h"
#include "SpinWaveGenie/Containers/Cell.h"

namespace SpinWaveGenie
{
    
    void PointsOnAPlane::setOriginPoint(double H0, double K0, double L0)
    {
        m_H0 = H0;
        m_K0 = K0;
        m_L0 = L0;
    }
    
    void PointsOnAPlane::setFinalPointFirstDirection(double H1, double K1, double L1)
    {
        m_H1 = H1;
        m_K1 = K1;
        m_L1 = L1;
    }
    
    void PointsOnAPlane::setFinalPointSecondDirection(double H2, double K2, double L2)
    {
        m_H2 = H2;
        m_K2 = K2;
        m_L2 = L2;
    }
    
    void PointsOnAPlane::setNumberPoints(std::size_t points1, std::size_t points2)
    {
        m_numberPoints1 = points1;
        m_numberPoints2 = points2;
    }
    
    std::tuple<std::size_t,std::size_t> PointsOnAPlane::getNumberPoints()
    {
        return std::tie(m_numberPoints1,m_numberPoints2);
    }
       
    ThreeVectors<double> PointsOnAPlane::getPoints()
    {
        ThreeVectors<double> Kpoints;
        if (m_numberPoints1 == 1 && m_numberPoints2 == 1)
        {
            Kpoints.insert(m_H0, m_K0, m_L0);
        }
        else
        {
            for (auto m = 0; m < m_numberPoints2; ++m)
            {
                double Hmi = m_H0 + (m_H2 - m_H0) * static_cast<double>(m) / static_cast<double>(m_numberPoints2 - 1);
                double Kmi = m_K0 + (m_K2 - m_K0) * static_cast<double>(m) / static_cast<double>(m_numberPoints2 - 1);
                double Lmi = m_L0 + (m_L2 - m_L0) * static_cast<double>(m) / static_cast<double>(m_numberPoints2 - 1);
                for (auto n = 0; n < m_numberPoints1; ++n)
                {
                    double Hni = m_H0 + (m_H1 - m_H0) * static_cast<double>(m) / static_cast<double>(m_numberPoints1 - 1);
                    double Kni = m_K0 + (m_K1 - m_K0) * static_cast<double>(m) / static_cast<double>(m_numberPoints1 - 1);
                    double Lni = m_L0 + (m_L1 - m_L0) * static_cast<double>(m) / static_cast<double>(m_numberPoints1 - 1);
                    Kpoints.insert(Hmi+Hni,Kmi+Kni,Lmi+Lni);
                }
            }
        }
        return Kpoints;
    }
}