#include <cmath>
#include <nlopt.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

/*#include "log4cplus/logger.h"
#include "log4cplus/loglevel.h"
#include "log4cplus/layout.h"
#include "log4cplus/loggingmacros.h"
#include "log4cplus/ndc.h"

#include "log4cplus/consoleappender.h"
#include "log4cplus/fileappender.h"
*/
#include "Genie/SpinWave.h"
#include "Initializer.h"
#include "Cell/Cell.h"
#include "Cell/Neighbors.h"
#include "Interactions/Exch_Interaction.h"
#include "Interactions/AnisotropyInteraction.h"


using namespace std;
using namespace Eigen;
//using namespace log4cplus;
//using namespace log4cplus::helpers;

double fitting(SpinWave &SW, double KX, double KY, double KZ, double freq, long pos, double error)
{
    
    double diff_sq = 0.0;
    SW.createMatrix(KX,KY,KZ);
    SW.Calc();
    vector<point> points1 = SW.getPoints();
    
    SW.createMatrix(KZ,KX,KY);
    SW.Calc();
    vector<point> points2 = SW.getPoints();
    
    SW.createMatrix(KY,KZ,KX);
    SW.Calc();
    vector<point> points3 = SW.getPoints();
    
    double average = (points1[pos].frequency*points1[pos].intensity + points2[pos].frequency*points2[pos].intensity  + points3[pos].frequency*points3[pos].intensity);
    average = average/(points1[pos].intensity+points2[pos].intensity+points3[pos].intensity);
    diff_sq += pow((average-freq)/error,2);
    
    /*double lower_limit = (points1[2].frequency*points1[2].intensity + points2[2].frequency*points2[2].intensity  + points3[2].frequency*points3[2].intensity)/(points1[2].intensity+points2[2].intensity+points3[2].intensity);
    if (lower_limit < 15.0)
    {
        cout << "missing spin gap" << endl;
        diff_sq += pow(lower_limit-15.0,2);
    }*/
    
    return diff_sq;
}

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    double *angle = (double *) my_func_data;
    if (!grad.empty())
    {
        cout << "error: no gradient available" << endl;
    }
    
    cout << "x[i] = ";
    for (size_t i = 0; i<x.size();i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;
    
    double SA = 2.3;
    double SB = 0.9;
    double theta = M_PI - (*angle)*M_PI/180.0;
    
    Cell cell;
    cell.setBasisVectors(8.5,8.5,8.5,90.0,90.0,90.0);
    
    Sublattice Mn0;
    string name = "Mn0";
    Mn0.setName(name);
    Mn0.setType("MN2");
    Mn0.setMoment(SA,0.0,0.0);
    cell.addSublattice(Mn0);
    cell.addAtom(name,0.0,0.0,0.0);
    cell.addAtom(name,0.0,0.5,0.5);
    cell.addAtom(name,0.5,0.0,0.5);
    cell.addAtom(name,0.5,0.5,0.0);
    
    Sublattice Mn1;
    name = "Mn1";
    Mn1.setName(name);
    Mn1.setType("MN2");
    Mn1.setMoment(SA,0.0,0.0);
    cell.addSublattice(Mn1);
    cell.addAtom(name,0.75,0.25,0.75);
    cell.addAtom(name,0.75,0.75,0.25);
    cell.addAtom(name,0.25,0.25,0.25);
    cell.addAtom(name,0.25,0.75,0.75);
    
    Sublattice V0;
    name = "V0";
    V0.setName(name);
    V0.setType("V3");
    V0.setMoment(SB,theta,3.0*M_PI/4.0);
    cell.addSublattice(V0);
    cell.addAtom(name,0.875,0.125,0.375);
    cell.addAtom(name,0.875,0.625,0.875);
    cell.addAtom(name,0.375,0.125,0.875);
    cell.addAtom(name,0.375,0.625,0.375);
    
    Sublattice V1;
    name = "V1";
    V1.setName(name);
    V1.setType("V3");
    V1.setMoment(SB,theta,7.0*M_PI/4.0);
    cell.addSublattice(V1);
    cell.addAtom(name,0.125,0.375,0.875);
    cell.addAtom(name,0.125,0.875,0.375);
    cell.addAtom(name,0.625,0.375,0.375);
    cell.addAtom(name,0.625,0.875,0.875);
    
    Sublattice V2;
    name = "V2";
    V2.setName(name);
    V2.setType("V3");
    V2.setMoment(SB,theta,M_PI/4.0);
    cell.addSublattice(V2);
    cell.addAtom(name,0.375,0.875,0.125);
    cell.addAtom(name,0.375,0.375,0.625);
    cell.addAtom(name,0.875,0.875,0.625);
    cell.addAtom(name,0.875,0.375,0.125);
    
    Sublattice V3;
    name = "V3";
    V3.setName(name);
    V3.setType("V3");
    V3.setMoment(SB,theta,5.0*M_PI/4.0);
    cell.addSublattice(V3);
    cell.addAtom(name,0.625,0.625,0.625);
    cell.addAtom(name,0.625,0.125,0.125);
    cell.addAtom(name,0.125,0.625,0.125);
    cell.addAtom(name,0.125,0.125,0.625);
    
    SW_Builder builder(cell);

    builder.Add_Interaction(new Exch_Interaction("Jaa",x[0],"Mn0","Mn1",3.0,5.0));
    
    builder.Add_Interaction(new Exch_Interaction("Jbb",x[1],"V0","V1",2.975,3.06));
    builder.Add_Interaction(new Exch_Interaction("Jbb",x[1],"V2","V3",2.975,3.06));

    builder.Add_Interaction(new Exch_Interaction("Jbbp",x[2],"V0","V2",2.975,3.06));
    builder.Add_Interaction(new Exch_Interaction("Jbbp",x[2],"V0","V3",2.975,3.06));
    builder.Add_Interaction(new Exch_Interaction("Jbbp",x[2],"V1","V2",2.975,3.06));
    builder.Add_Interaction(new Exch_Interaction("Jbbp",x[2],"V1","V3",2.975,3.06));
    
    Vector3 zhat(0.0,0.0,1.0);

    builder.Add_Interaction(new AnisotropyInteraction("Daz",x[3],zhat,"Mn0"));
    builder.Add_Interaction(new AnisotropyInteraction("Daz",x[3],zhat,"Mn1"));
    
    builder.Add_Interaction(new AnisotropyInteraction("Dbz",x[4],zhat,"V0"));
    builder.Add_Interaction(new AnisotropyInteraction("Dbz",x[4],zhat,"V1"));
    builder.Add_Interaction(new AnisotropyInteraction("Dbz",x[4],zhat,"V2"));
    builder.Add_Interaction(new AnisotropyInteraction("Dbz",x[4],zhat,"V3"));
    
    Vector3 direction(-1.0,1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dby",x[5],direction,"V0"));
    direction = Vector3(1.0,-1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dbx",x[5],direction,"V1"));
    direction = Vector3(1.0,1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dby",x[5],direction,"V2"));
    direction = Vector3(-1.0,-1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dbx",x[5],direction,"V3"));

    
    double JBB = x[1];
    double JBBP = x[2];
    double DB = x[5];
    double DBz = x[4];
    double JAB = SB*((6.0*JBB+6.0*JBBP+DB-3.0*DBz)*cos(theta)*sin(theta) -sqrt(2.0)*DB*(2.0*pow(cos(theta),2)-1.0))/(-9.0*SA*sin(theta));
    
    cout << "JAB= " << JAB << endl;
    
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn0","V0",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn0","V1",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn0","V2",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn0","V3",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn1","V0",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn1","V1",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn1","V2",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn1","V3",3.48,3.57));
    
    SpinWave SW = builder.Create_Element();
    
    double diff_sq = 0.0;
    
    
    //MnV2O4
    diff_sq += fitting(SW,2.0,0.0,0.0898357,9.67548,1,1.0);
    diff_sq += fitting(SW,2.0,0.0,0.300308,9.59135,1,1.0);
    diff_sq += fitting(SW,2.0,0.0,0.490246,9.59135,1,1.0);
    diff_sq += fitting(SW,2.0,0.0,0.703285,9.63341,1,1.0);
    diff_sq += fitting(SW,2.0,0.0,0.898357,9.54928,1,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.10370,8.96034,1,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.29877,7.86659,0,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.48871,6.35216,0,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.68121,4.33293,0,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.85062,2.52404,0,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.94559,1.93510,0,1.0);
    diff_sq += fitting(SW,2.0,0.0,2.14836,2.52404,0,1.0);
    diff_sq += fitting(SW,2.0,0.0,2.31519,4.16466,0,1.0);
    diff_sq += fitting(SW,2.0,0.0,2.49743,6.05769,0,1.0);
    
    diff_sq += fitting(SW,4.0-1.09151,0.0,1.09151,9.16866,0,1.0);
    diff_sq += fitting(SW,4.0-1.29285,0.0,1.29285,7.78708,0,1.0);
    diff_sq += fitting(SW,4.0-1.49251,0.0,1.49251,6.57297,0,1.0);
    diff_sq += fitting(SW,4.0-1.65225,0.0,1.65225,5.35885,0,1.0);
    diff_sq += fitting(SW,4.0-1.75541,0.0,1.75541,4.31220,0,1.0);
    diff_sq += fitting(SW,4.0-1.84859,0.0,1.84859,3.01435,0,1.0);
    diff_sq += fitting(SW,4.0-1.94842,0.0,1.94842,1.92584,0,1.0);
    
    //MnV2O4
    //diff_sq += fitting(SW,2.0,0.0,0.0,9.67548,1,1.0);
    //diff_sq += fitting(SW,2.0,0.0,1.94559,1.93510,0,1.0);
    
    //diff_sq += fitting(SW,4.0-1.0,0.0,1.0,9.16866,0,1.0);
    //diff_sq += fitting(SW,4.0-2.0,0.0,2.0,2.0,0,1.0);
    
    
    //Mn0.6Co0.4V2O4
    //diff_sq += fitting(SW,1.0,1.0,0.0,10.06597219,0,1.0);;
    //diff_sq += fitting(SW,2.0,2.0,0.0,2.0,0,0.1);
    
    
    //Mn0.4=6Co0.4V2O4
    //diff_sq += fitting(SW,0.0,0.0,4.00,2.0,0,1.0);
    //diff_sq += fitting(SW,0.0,0.0,3.00,10.66,0,1.0);

    //diff_sq += fitting(SW,2.00,2.00,0.0,2.0,0,0.5);
    //diff_sq += fitting(SW,1.00,1.00,0.0,10.49,0,1.0);
    
    
    cout << "diff_sq= " << diff_sq << endl;
    return diff_sq;
}

int main()
{
    
    /*SharedAppenderPtr myConsoleAppender(new ConsoleAppender());
    myConsoleAppender->setName("myConsoleAppender");
    
    SharedAppenderPtr myFileAppender(new FileAppender("myLogFile.log"));
    myFileAppender->setName("myFileAppender");
    
    std::auto_ptr<Layout> myLayout = std::auto_ptr<Layout>(new log4cplus::TTCCLayout());
    Logger myLogger= Logger::getInstance("myLoggerName");
    myLogger.addAppender(myConsoleAppender);
    
    LOG4CPLUS_FATAL(myLogger, "logEvent");//for a fatal priority event</pre>
    LOG4CPLUS_ERROR(myLogger, "logEvent");//for a error priority event</pre>
    LOG4CPLUS_WARN(myLogger, "logEvent") ;//for a warn priority event</pre>
    LOG4CPLUS_INFO(myLogger, "logEvent");      //for a info priority event</pre>
    LOG4CPLUS_DEBUG(myLogger, "logEvent");//for a debug priority event</pre>
    */
    nlopt::opt opt(nlopt::LN_COBYLA,6);
    std::vector<double> ub(6);
    ub[0] =  0.0;
    ub[1] = -9.0;
    ub[2] =  5.0;
    ub[3] =  1.0;
    ub[4] =  4.0;
    ub[5] = 0.0;
    opt.set_upper_bounds(ub);
    std::vector<double> lb(6);
    lb[0] = 0.0;
    lb[1] = -11.0;
    lb[2] =  0.0;
    lb[3] = -1.0;
    lb[4] = -4.0;
    lb[5] = -10.0;
    opt.set_lower_bounds(lb);
    opt.set_xtol_rel(1.0e-4);
    std::vector<double> x(6);
    
    double angle = 36;
       
    
    for (int i = 36 ; i <= 36 ; i++)
    {
        x[0] = 0.0;
        x[1] = -10.0;
        x[2] = 3.6;
        x[3] = 0.37;
        x[4] = -2.6;
        x[5] = -2.9;
        
        angle = (double)i;
        cout << "angle = " << angle << endl;
        opt.set_min_objective(myfunc,&angle);
 
        double minf = 0.0;
    
        nlopt::result result = opt.optimize(x, minf);
    
        //cout << result << endl;
    
        //cout << minf << endl;
    
        cout << "results:" << endl;
        for (vector<double>::iterator it = x.begin(); it!=x.end();it++)
            cout << (*it) << " ";
        cout << endl;
        cout << "chi_squared = " << minf << endl;
        //cout << "reduced chi_squared = " << minf/21.0 << endl;
    }
    
    /*double minf = 2.55511;
    vector<double> gradient;
    

    //cout << endl;
    long i = 1;
    vector<double> xtemp = x;
    for (int j=-200;j!=201;j++)
    {
        xtemp[i] = x[i] + x[i]*(double)j/1000.0;
        double diffsq = myfunc(xtemp,gradient,(void*)&angle);
        cout << "x[" << i << "] = " << xtemp[i];
        cout << ", reduced chi^2= " << (diffsq - minf)/4.0 -0.0 << endl;
    }*/

    return 0;
}
