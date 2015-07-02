#ifndef spin_wave_genie_PointsOnAPlane_h
#define spin_wave_genie_PointsOnAPlane_h

#include <iostream>
#include <tuple>
#include "SpinWaveGenie/Containers/ThreeVectors.h"

namespace SpinWaveGenie
{
   
    //! generates k-points used by the SpinWavePlot routines.
    /*!
     This class generates a grid of evenly spaced k-points on a plane.
     The result is stored in a ThreeVectors<double> container.
     */
    
    class PointsOnAPlane
    {
    public:
        //! Set starting point of line. (inclusive)
        //! \param H0 H component of k-point,in rlu.
        //! \param K0 K component of k-point,in rlu.
        //! \param L0 L component of k-point,in rlu.
        void setOriginPoint(double H0, double K0, double L0);
        //! Set final point of line. (inclusive)
        //! \param H1 H component of k-point,in rlu.
        //! \param K1 K component of k-point,in rlu.
        //! \param L1 L component of k-point,in rlu.
        void setFinalPointFirstDirection(double H1, double K1, double L1);
        //! Set number of k-points along line in the first direction.
        //! \param points number of k-points.
        void setFinalPointSecondDirection(double H2, double K2, double L2);
        //! Set number of k-points along line.
        //! \param points number of k-points.
        void setNumberPoints(std::size_t points1, std::size_t points2);
        //! Get number of points along the first and second directions.
        //! \return points along two dimensions.
        std::tuple<std::size_t,std::size_t> getNumberPoints();
        //! Get resulting k-points along line.
        //! \return k-points along line.
        ThreeVectors<double> getPoints();
        
    private:
        double m_H0,m_K0,m_L0,m_H1,m_K1,m_L1,m_H2,m_K2,m_L2;
        long m_numberPoints1,m_numberPoints2;
    };
}

#endif
