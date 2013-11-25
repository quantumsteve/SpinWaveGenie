#include "Anis_X_Interaction.h"
#include <iostream>

using namespace std;
using namespace Eigen;

Anis_X_Interaction::Anis_X_Interaction(double value_in, string sl_r_in)
{
    this->Update_Interaction(value_in, sl_r_in);
}

Interaction* Anis_X_Interaction::do_clone() const
{
    return new Anis_X_Interaction(*this);

}

void Anis_X_Interaction::Update_Interaction(double value_in, string sl_r_in)
{
    value = value_in;
    sl_r = sl_r_in;
}

vector<string> Anis_X_Interaction::sublattices() const
{
    vector<string> sl;
    sl.push_back(sl_r);
    return sl;
}

void Anis_X_Interaction::calcConstantValues(Cell& cell)
{
    complex<double> XI (0.0,1.0);
    //find location of r,s
    r= -1;
    M = 0;
    for (Cell::Iterator sl=cell.begin(); sl!=cell.end(); ++sl)
    {
        //cout << sl_r << " " << sl->getName() << endl;
        if ( sl_r == sl->getName())
            r = M;
        M++;
    }
    assert(r!=-1);
    
    double S = cell.getSublattice(sl_r).getMoment();
    double theta = cell.getSublattice(sl_r).getTheta();
    double phi = cell.getSublattice(sl_r).getPhi();
    double X = value;
    LNrr =   -0.5*X*S*(pow(cos(theta),2)*pow(cos(phi),2)+pow(sin(phi),2)-2.0*pow(sin(theta),2)*pow(cos(phi),2));
    LNrrM =  -0.5*X*S*pow(cos(theta)*cos(phi)-XI*sin(phi),2);
    LNrMr =  -0.5*X*S*pow(cos(theta)*cos(phi)+XI*sin(phi),2);
    LNrMrM = -0.5*X*S*(pow(cos(theta),2)*pow(cos(phi),2)+pow(sin(phi),2)-2.0*pow(sin(theta),2)*pow(cos(phi),2));
    
    cout << "working X implementation" << endl;
    cout << LNrr << " "<< LNrMr << " " << LNrrM << " " << endl;
}

void Anis_X_Interaction::checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements)
{
        complex<double> XI (0.0,1.0);
        double S = cell.getSublattice(sl_r).getMoment();
        double theta = cell.getSublattice(sl_r).getTheta();
        double phi = cell.getSublattice(sl_r).getPhi();
        elements[r] -= sqrt(pow(S,3)*2.0)*value*sin(theta)*cos(phi)*(cos(theta)*cos(phi)+XI*sin(phi));
        elements[r+M] -= sqrt(pow(S,3)*2.0)*value*sin(theta)*cos(phi)*(cos(theta)*cos(phi)-XI*sin(phi));
}

void Anis_X_Interaction::Update_Matrix(Vector3d K, MatrixXcd &LN, int quadrant)
{
    switch (quadrant)
    {
        case 0:
            LN(r,r) += LNrr;
            break;
        case 1:
            LN(r,r+M) += LNrrM;
            break;
        case 2:
            LN(r+M,r) += LNrMr;
            break;
        case 3:
            LN(r+M,r+M) += LNrMrM;
            break;
        default:
            //cout << "error: case must be between 0 and 3" << endl;
            break;
    }
}
