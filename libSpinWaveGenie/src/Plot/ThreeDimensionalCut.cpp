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

#include <vtkOpenGLExtensionManager.h>


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

/*void ThreeDimensionalCut::save()
{
  Eigen::MatrixXd figure = getMatrix();
  std::ofstream file(m_p->filename + ".mat");
  if (file.is_open())
  {
    file << figure << endl;
  }
  file.close();

  file.open(m_p->filename + ".x");
  if (file.is_open())
  {
    ThreeVectors<double> pts = m_p->points;
    for (auto it = pts.begin(); it != pts.end(); ++it)
      file << it->get<0>() << "\t" << it->get<1>() << "\t" << it->get<2>() << endl;
  }

  file.close();
  file.open(m_p->filename + ".y");
  if (file.is_open())
  {
    Energies energies = m_p->cut->getEnergies();
    for (auto it = energies.begin(); it != energies.end(); ++it)
      file << (*it) << endl;
  }
  file.close();
}
*/
void ThreeDimensionalCut::save()
{
    Eigen::MatrixXd figure = getMatrix();
    
    // Create a grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid =
    vtkSmartPointer<vtkStructuredGrid>::New();
    
    vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
    double x, y, z;
    
    x = 0.0;
    y = 0.0;
    z = 0.0;
    
    std::string m_scalarName("scalarData");
    vtkSmartPointer<vtkDoubleArray> signal = vtkSmartPointer<vtkDoubleArray>::New();
    signal->SetName(m_scalarName.c_str());
    signal->SetNumberOfComponents(1);
    
    const Matrix3& recip = m_p->cut->getCell().getReciprocalVectors();
    const Energies& energies = m_p->cut->getEnergies();
    for(auto it = m_p->points.begin(); it!= m_p->points.end();++it)
    {
      Vector3 K;
      K << it->get<0>(),it->get<1>(),it->get<2>();
      K = K.transpose() * recip;
      auto pos1 = std::distance(m_p->points.begin(),it);
      for(auto it2 = energies.cbegin(); it2 != energies.cend(); ++it2)
      {
        auto pos2 = std::distance(energies.cbegin(),it2);
        points->InsertNextPoint(K[0],K[1],*it2);
        signal->InsertNextValue(figure(pos2,pos1));
      }
    }
    
    std::size_t xdim,ydim;
    std::tie(xdim,ydim) = m_p->GeneratePoints.getNumberPoints();
    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(xdim,ydim,energies.size());
    structuredGrid->SetPoints(points);
    structuredGrid->GetPointData()->SetScalars(signal.GetPointer());
    
    std::vector<double> u = {recip(0,0),recip(0,1),recip(0,2)};
    std::vector<double> v = {recip(1,0),recip(1,1),recip(1,2)};
    std::vector<double> w = {recip(2,0),recip(2,1),recip(2,2)};
    //u[0] = 6.28319; u[1]= v[1]=1.0; w[0]=1.0;
    
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
    
    std::vector<double> bbox = {0.0,recip(0,0)*3.0,0.0,recip(1,1)*3.0,0.0,5.0};
    vtkNew<vtkDoubleArray> bounds;
    bounds->SetName("BoundingBoxInModelCoordinates");
    bounds->SetNumberOfComponents(6);
    bounds->SetNumberOfTuples(1);
    std::copy(bbox.begin(), bbox.end(), bounds->GetPointer(0));
    structuredGrid->GetFieldData()->AddArray(bounds.GetPointer());
    
    // Write file
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName("output.vts");
    writer->SetInputData(structuredGrid);
    writer->Write();

}

/*
    int main(int, char *[])
    {
 
        // Create a grid
        vtkSmartPointer<vtkStructuredGrid> structuredGrid =
        vtkSmartPointer<vtkStructuredGrid>::New();
 
        vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
        double x, y, z;
 
        x = 0.0;
        y = 0.0;
        z = 0.0;
 
        std::string m_scalarName("scalarData");
        vtkSmartPointer<vtkDoubleArray> signal = vtkSmartPointer<vtkDoubleArray>::New();
        signal->SetName(m_scalarName.c_str());
        signal->SetNumberOfComponents(1);
        
        for(unsigned int k = 0; k < 11; k++)
        {
            y = 0.0;
            for(unsigned int j = 0; j < 11; j++)
            {
                x = 0.0;
                for(unsigned int i = 0; i < 11; i++)
                {
                    //std::cout << x << " " << y << " " << z <<std::endl;
                    points->InsertNextPoint(.628319*x, .36276*x+.7252*y, .418879*z);
                    signal->InsertNextValue(x);
                    x += 0.1;
                }
                y += 0.1;
            }
            z += 0.1;
        }
        
        // Specify the dimensions of the grid
        structuredGrid->SetDimensions(11,11,11);
        structuredGrid->SetPoints(points);
        structuredGrid->GetPointData()->SetScalars(signal.GetPointer());
        
        int* dims = structuredGrid->GetDimensions();
        
        // Retrieve the entries from the grid and print them to the screen
        unsigned int counter = 0;
        
        for (int k = 0; k < dims[2]; k++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                for (int i = 0; i < dims[0]; i++)
                {
                    double p[3];
                    structuredGrid->GetPoint(counter, p);
                    
                    double pNew[3];
                    structuredGrid->GetPoint(i, j, k, pNew);
                    
                    //std::cout << "P   : "
                    //          << p[0] << " "
                    //          << p[1] << " "
                    //          << p[2] << std::endl;
                    //std::cout << "PNew: "
                    //          << pNew[0] << " "
                    //          << pNew[1] << " "
                    //          << pNew[2] << std::endl;
                    
                    counter++;
                }
            }
        }
        
        std::vector<double> u = {.628319,.36276,0.0};
        std::vector<double> v = {0.0,.7252,0.0};
        std::vector<double> w = {0.0,0.0,.418879};
        //u[0] = 6.28319; u[1]= v[1]=1.0; w[0]=1.0;
        
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
        
        std::vector<double> bbox = {0.0,1.0,0.0,1.0,0.0,1.0};
        vtkNew<vtkDoubleArray> bounds;
        bounds->SetName("BoundingBoxInModelCoordinates");
        bounds->SetNumberOfComponents(6);
        bounds->SetNumberOfTuples(1);
        std::copy(bbox.begin(), bbox.end(), bounds->GetPointer(0));
        structuredGrid->GetFieldData()->AddArray(bounds.GetPointer());
        
        // Write file
        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName("output.vts");
        writer->SetInputData(structuredGrid);
        writer->Write();
        
        // Create a mapper and actor
        vtkSmartPointer<vtkDataSetMapper> mapper =
        vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
        mapper->SetInputConnection(structuredGrid->GetProducerPort());
#else
        mapper->SetInputData(structuredGrid);
#endif
        
        vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        
        // Create a renderer, render window, and interactor
        vtkSmartPointer<vtkRenderer> renderer =
        vtkSmartPointer<vtkRenderer>::New();
        vtkSmartPointer<vtkRenderWindow> renderWindow =
        vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
        vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);
        
        vtkSmartPointer<vtkOpenGLExtensionManager> asdf = vtkOpenGLExtensionManager::New();
        asdf->SetRenderWindow(renderWindow.GetPointer());
        asdf->Update();
        std::cout << asdf->DriverIsNvidia() << std::endl;
        std::cout << asdf->DriverIsATI() << std::endl;
        std::cout << asdf->DriverIsMesa() << std::endl;
        std::cout << asdf->DriverIsIntel() << std::endl;
        const char* vendor = asdf->GetDriverGLVendor();
        std::cout << vendor << std::endl;
        std::cout << asdf->GetDriverGLRenderer() << std::endl;
        std::cout << asdf->GetDriverGLVersionMajor() << std::endl;
        std::cout << asdf->GetDriverGLVersionMinor() << std::endl;
        std::cout << asdf->GetDriverGLVersionPatch() << std::endl;
        std::cout << asdf->GetDriverVersionMajor() << std::endl;
        std::cout << asdf->GetDriverVersionMinor() << std::endl;
        std::cout << asdf->GetDriverVersionPatch() << std::endl;
        
        // Add the actor to the scene
        renderer->AddActor(actor);
        renderer->SetBackground(1.,1.,1.); // Background color white
        
        // Render and interact
        renderWindow->Render();
        renderWindowInteractor->Start();
        
        return EXIT_SUCCESS;
    }

*/
    
}

