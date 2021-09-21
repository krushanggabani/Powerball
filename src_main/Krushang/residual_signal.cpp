#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <chrono>
#include <ostream> // included for color output to the terminal

// boost headers
#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>


// schunk powerball headers
#include "powerball/schunk_powerball.h"
#include "vrep/v_repClass.h"
#include "powerball/schunk_kinematics.h"

// Toon headers
#include <TooN/TooN.h>
#include <TooN/LU.h>
#include <TooN/SVD.h>

// dynamixel headers
#include <utils.h>
#include <USB2Dynamixel.h>
#include <commonOptions.h>

// phidget headers
#include <phidget21.h>


// plotting and linear algebra libraries
#include "matplotlibcpp.h"
// #include "sigpack.h"
// #include <armadillo>

//the following are UBUNTU/LINUX ONLY terminal color codes.
#define RESET       "\033[0m"
#define RED         "\033[31m"              /* Red */
#define GREEN       "\033[32m"              /* Green */
#define BOLDRED     "\033[1m\033[31m"       /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"       /* Bold Green */


using namespace::std;
using namespace::TooN;
namespace plt = matplotlibcpp;
using boost::asio::ip::tcp;

const double dt = 0.01; // sampling time
float integral_sum      = 0;
float residual          = 0; 
float mass              = 0.2260;                //  (100+2*45+2*18)/1000;   // mass of the gripper
float Fv                = 0.45;                   // viscous coeeficient
Vector<3,float> xe      = makeVector(1,0,0);      // coordinate of X axis of end effector in end effector frame
Vector<3,float> gravity = makeVector(0, 0,-9.81); // gravity vector in base frame 
float KI                = 40;                     // detection gain


// Magnet force model
double magForce(double s){
  double a = 5.5e-6;
  double b = 6.85e-5;
  double c = 2.22e-7;

  if (s<0){
      s=abs(s);
  }
  return a/(pow(s,3)+b*s+c);
}

// Skew Symmetry Matrix 
void skew_symmetry(Vector<3,float> w, Matrix<3,3,float> *S){
    float w1= w[0];
    float w2= w[1];
    float w3= w[2];


    *S = Data(0 , -w3, w2 ,
            w3 , 0 , -w1,
            -w2, w1, 0  );
       
}


int main(int argc, char **argv){

    Kin kin;

    cout<< "Reading file...."<<endl;
    std::vector< std::vector<double> > collisiondata;
    kin.inputFile("../data/Gripper_collision_data/MyFile.txt",&collisiondata);

    int  ncols, nrows=0;
    int traj_col = 0;
    for (std::vector< std::vector<double> >::const_iterator it = collisiondata.begin(); it != collisiondata.end(); ++ it)
    {
        nrows++;
        ncols = 0;
        for (std::vector<double>::const_iterator itit = it->begin(); itit != it->end(); ++ itit)
        {
            ncols++;
        }
    }
    cout << "size of imported matrix = " << nrows << "*" << ncols << endl;

    TooN::Matrix<Dynamic,Dynamic,double> imported_data(nrows, ncols);

    if (ncols == 33){

        // Put the matrix into TooN matrix
        nrows = 0;  ncols = 0;
        for (std::vector< std::vector<double> >::const_iterator it = collisiondata.begin(); it != collisiondata.end(); ++ it)
        {
            ncols = 0;
            for (std::vector<double>::const_iterator itit = it->begin(); itit != it->end(); ++ itit)
            {
                imported_data(nrows,ncols) = *itit;
                ncols++;
            }
            nrows++;
            traj_col = ncols;
        }
    }else{
        cout<< BOLDRED <<"Inconsistent matrix size for trajectory generation" << RESET << endl;
        return 0;
    }
    cout << "File read succesful" <<endl;
    /*------------------------------------*/

    // file values
    // 0-5   Q
    // 6-11  Qdot
    // 12-17 joint torque
    // 18-20 force
    // 21-23 torque
    // 24-27 airgap
    // 28-31 mag force
    // 32    outermag pose

    nrows=1700;
    Matrix<1700,3,float> Omega =Zeros;
    std::vector<float> end_effector_frame (nrows,0);

    std::vector<float> x1_dot (nrows,0);
    std::vector<float> Fm(nrows,0);
    std::vector<float> residual_history(nrows,0);
    std::vector<float> pe_x (nrows,0);
    std::vector<float> pe_y (nrows,0);
    std::vector<float> pe_z (nrows,0);
    std::vector<float> timestamp(nrows,0);



    // nrows=300;
    float previos_gripper_finger1 = imported_data[1][32]-imported_data[1][24];
    
   
    for(int i=1; i<=nrows; i++){

        
        Vector<33,float> Row_data = imported_data[i];
        Vector<6,float> Q_present = makeVector(Row_data[0],Row_data[1],Row_data[2],Row_data[3],Row_data[4],Row_data[5]) ;
        Vector<6,float> Qdot      = makeVector(Row_data[6],Row_data[7],Row_data[8],Row_data[9],Row_data[10],Row_data[11]);

        float timestep = 0.01; 
        timestamp[i] = timestep * i;
       
        Matrix<3,6,float> J_lin = Zeros;     // 3 x 6 linear part of jacobian Matrix
        Matrix<3,6,float> J_ang = Zeros;     // 3 x 6 angular part of jacobian Matrix 
        Matrix<3,3,float> R     = Zeros;     // 3 x 3 forward rotation matrix 
        
        kin.JacobPos(Q_present,&J_lin);
        kin.JacobRot(Q_present,&J_ang);
        kin.FK_R(Q_present,&R);
        

        // Matrix<3,3,float> R_trans = R.T(); 
        Matrix<3,3,float> R_trans = R;
        Matrix<3,3,float> Sw =Zeros; 
        Omega[i]=R_trans * J_ang * Qdot;
        skew_symmetry(Omega[i],&Sw);


        float Gripper_finger1 = Row_data[32]-Row_data[24];
        float Gripper_finger2 = Row_data[32]-Row_data[26];
        
        Vector<3,float> r_end_effector_pos = makeVector(Gripper_finger1-Gripper_finger2,0,0.17);
        x1_dot[i]       = (Gripper_finger1-previos_gripper_finger1)/timestep;  


        end_effector_frame[i] = r_end_effector_pos[0];
        previos_gripper_finger1 = Gripper_finger1;

        float Fm1 = magForce(Row_data[24]) + magForce(Row_data[25]);
        float Fm2 = magForce(Row_data[26]) + magForce(Row_data[27]);
        float Fe  = Fm1-Fm2;
        Fm[i] = magForce(Row_data[24]) - magForce(Row_data[25])- magForce(Row_data[26]) + magForce(Row_data[27]);
        
        //calculating absolute translation velocity of grasped object in end-effector frame
        Vector<3,float> pe_dot = R_trans*J_lin*Qdot + Sw * r_end_effector_pos + makeVector(x1_dot[i],0,0);
        pe_x[i] = pe_dot[0];
        pe_y[i] = pe_dot[1];
        pe_z[i] = pe_dot[2];
       


        // calculating linear momentum and its derivative
        float p_linear_momentum = mass* (xe * pe_dot);
        float p_linear_mom_dot  = mass * (xe *(R_trans*gravity - Sw *pe_dot) ) + Fm[i] + Fv * x1_dot[i]; 
    
        float integral_part = (p_linear_mom_dot+residual_history[i-1])*timestep;
        integral_sum = integral_sum+ integral_part;
        
        // *pre_sum = integral_sum; 
        residual_history[i] = KI * (p_linear_momentum-integral_sum);

        
                
    }

    
    plt::figure(); 
    plt::plot(residual_history);
    // plt::plot(pe_x,"r--");
    // plt::plot(pe_y,"g--");
    // plt::plot(pe_z,"b--");
    plt::show();



}
