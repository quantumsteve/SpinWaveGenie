#ifndef __DMY_Interaction_H__
#define __DMY_Interaction_H__ 1

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "Cell.h"
#include "Interaction.h"

class DM_Y_Interaction: public Interaction
{
public:
    DM_Y_Interaction(double value_in, std::string sl_r_in,std::string sl_s_in, double min_in, double max_in);
    void Update_Interaction(double value_in, std::string sl_r_in,std::string sl_s_in, double min_in, double max_in);
    void calcConstantValues(boost::shared_ptr<Cell> cell);
    void calcChangingValues(boost::shared_ptr<Cell> cell, Eigen::Vector3d K);
    void Update_Matrix(Eigen::Vector3d K, boost::shared_ptr<Cell> cell, Eigen::MatrixXcd &LN,int quadrant);
    std::vector<std::string> sublattices() const;
private:
    std::string sl_r,sl_s;
    int r,s,M;
    double value,min,max;
    double value0,value1,value2,value3;
    double z_rs;
    std::complex<double> gamma_rs;
};

#endif