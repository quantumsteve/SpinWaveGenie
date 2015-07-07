//
//  ThreeDimensionalCut.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/16/14.
//
//
#include <fstream>
#include <Eigen/Dense>
#include "SpinWaveGenie/Plot/ThreeDimensionalCut.h"
#include "SpinWaveGenie/Plot/EnergyResolutionFunction.h"
#include "SpinWaveGenie/Containers/Energies.h"
#include "SpinWaveGenie/Containers/PointsOnAPlane.h"
#include "External/ezRateProgressBar.hpp"
#include <thread>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif
#ifdef USE_THREADS
#include "tbb/tbb.h"
using namespace tbb;
#endif

// fix for gcc 4.4 only having cstdatomic
#ifdef HAVE_ATOMIC_H
#include <atomic>
#else
#include <cstdatomic>
#endif

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkNew.h>
#include <vtkDoubleArray.h>
#include "vtkImagePointsArray.h"

using namespace std;

namespace SpinWaveGenie
{
class ThreeDimensionalCut::CutImpl
{
public:
  std::string filename;
  atomic_size_t counter;
  Eigen::MatrixXd mat;
  unique_ptr<SpinWaveGenie::SpinWavePlot> cut;
  SpinWaveGenie::PointsOnAPlane GeneratePoints;
  SpinWaveGenie::ThreeVectors<double> points;
  CutImpl() { counter = 0; };
  CutImpl(unique_ptr<SpinWaveGenie::SpinWavePlot> inCut, SpinWaveGenie::PointsOnAPlane inPoints)
      : cut(move(inCut)), GeneratePoints(inPoints), points(inPoints.getPoints())
  {
    counter = 0;
  };
  std::unique_ptr<CutImpl> clone()
  {
    std::unique_ptr<CutImpl> newCut(new CutImpl(cut->clone(), GeneratePoints));
    newCut->filename = filename;
    newCut->mat = mat;
    return std::move(newCut);
  };

  void progressBar(std::size_t numberPoints)
  {
    ez::ezRateProgressBar<std::size_t> p(numberPoints);
    p.units = "Q-points";
    p.start();
    while (counter < numberPoints)
    {
      p.update(counter);
#ifdef _WIN32
      Sleep(1);
#else
      sleep(1);
#endif
    }
    p.update(numberPoints);
  }

  void partialCut(size_t begin, size_t end)
  {
    unique_ptr<SpinWaveGenie::SpinWavePlot> cutclone = cut->clone();
    for (size_t m = begin; m < end; m++)
    {
      auto it = points.begin() + m;
      vector<double> val = cutclone->getCut(it->get<0>(), it->get<1>(), it->get<2>());
      for (size_t n = 0; n < val.size(); n++)
      {
        mat(n, m) = val[n];
      }
      counter++;
    }
  }
#ifdef USE_THREADS
  Eigen::MatrixXd generateMatrix()
  {
    mat.resize(cut->getEnergies().size(), points.size());
    TbbExecutor tbbExec(this);
    thread myThread(&CutImpl::progressBar, this, points.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, points.size()), tbbExec);
    myThread.join();
    return mat;
  }
  struct TbbExecutor
  {
  public:
    TbbExecutor(CutImpl *w) : w_(w) {}
    void operator()(const tbb::blocked_range<size_t> r) const { w_->partialCut(r.begin(), r.end()); }

  private:
    CutImpl *w_;
  };
#else
  Eigen::MatrixXd generateMatrix()
  {
    mat.resize(cut->getEnergies().size(), points.size());
    thread pbar(&CutImpl::progressBar, this, points.size());
    partialCut(0, points.size());
    pbar.join();
    return mat;
  }
#endif
};

ThreeDimensionalCut::ThreeDimensionalCut() : m_p{new CutImpl{}} {};
ThreeDimensionalCut::ThreeDimensionalCut(const ThreeDimensionalCut &other) : m_p(other.m_p->clone()) {}
ThreeDimensionalCut &ThreeDimensionalCut::operator=(const ThreeDimensionalCut &other)
{
  m_p = move(other.m_p->clone());
  return *this;
}
ThreeDimensionalCut::ThreeDimensionalCut(ThreeDimensionalCut &&other)
{
  m_p = move(other.m_p);
  other.m_p = NULL;
}
ThreeDimensionalCut &ThreeDimensionalCut::operator=(ThreeDimensionalCut &&other)
{
  if (m_p != other.m_p)
  {
    m_p = move(other.m_p);
    other.m_p = NULL;
  }
  return *this;
}

ThreeDimensionalCut::~ThreeDimensionalCut(){};

void ThreeDimensionalCut::setFilename(const string& name) { m_p->filename = name; }

void ThreeDimensionalCut::setPlotObject(unique_ptr<SpinWavePlot> object) { m_p->cut = move(object); }

void ThreeDimensionalCut::setPoints(const PointsOnAPlane& pts)
{
  m_p->GeneratePoints = pts;
  m_p->points = m_p->GeneratePoints.getPoints();
}

void ThreeDimensionalCut::setEnergyPoints(double min, double max, size_t points)
{
  m_p->cut->setEnergies(Energies(min, max, points));
}

Eigen::MatrixXd ThreeDimensionalCut::getMatrix() { return m_p->generateMatrix(); }

void ThreeDimensionalCut::save()
{
    Eigen::MatrixXd figure = getMatrix();
    
    // Create a grid
    vtkNew<vtkStructuredGrid> structuredGrid;
    vtkNew<vtkPoints> points;
    
    std::string m_scalarName("scalarData");
    vtkSmartPointer<vtkDoubleArray> signal = vtkSmartPointer<vtkDoubleArray>::New();
    signal->SetName(m_scalarName.c_str());
    signal->SetNumberOfComponents(1);
    signal->SetArray(figure.data(),figure.size(),1);
    
    const Matrix3& recip = m_p->cut->getCell().getReciprocalVectors();
    const Energies& energies = m_p->cut->getEnergies();
    for(auto it = m_p->points.begin(); it!= m_p->points.end();++it)
    {
      Vector3 K;
      K << it->get<0>(),it->get<1>(),it->get<2>();
      K = K.transpose() * recip;
      for(auto it2 = energies.cbegin(); it2 != energies.cend(); ++it2)
      {
        points->InsertNextPoint(*it2,K[0],K[1]);
      }
    }
    
    std::size_t xdim,ydim;
    std::tie(xdim,ydim) = m_p->GeneratePoints.getNumberPoints();
    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(energies.size(),xdim,ydim);
    structuredGrid->SetPoints(points.GetPointer());
    structuredGrid->GetPointData()->SetScalars(signal.GetPointer());
    
    std::vector<double> u = {1.0,0.0,0.0};
    std::vector<double> v = {recip(0,2),recip(0,0),recip(0,1)};
    std::vector<double> w = {recip(1,2),recip(1,0),recip(1,1)};
    
    vtkSmartPointer<vtkMatrix4x4> cobMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
    cobMatrix->Identity();
    std::copy(u.begin(), u.end(), cobMatrix->Element[0]);
    std::copy(v.begin(), v.end(), cobMatrix->Element[1]);
    std::copy(w.begin(), w.end(), cobMatrix->Element[2]);
    
    cobMatrix->Transpose();
    
    vtkNew<vtkDoubleArray> cobArray;
    cobArray->SetName("ChangeOfBasisMatrix");
    cobArray->SetNumberOfComponents(16);
    cobArray->SetNumberOfTuples(1);
    std::copy(&cobMatrix->Element[0][0], (&cobMatrix->Element[0][0]) + 16,cobArray->GetPointer(0));
    structuredGrid->GetFieldData()->AddArray(cobArray.GetPointer());
    
    std::vector<double> bbox = {0.0,10.0,0.0,3.0,0.0,3.0};
    vtkNew<vtkDoubleArray> bounds;
    bounds->SetName("BoundingBoxInModelCoordinates");
    bounds->SetNumberOfComponents(6);
    bounds->SetNumberOfTuples(1);
    std::copy(bbox.begin(), bbox.end(), bounds->GetPointer(0));
    structuredGrid->GetFieldData()->AddArray(bounds.GetPointer());
    
    // Write file
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName("output.vts");
    writer->SetInputData(structuredGrid.GetPointer());
    writer->Write();

}
    
}

