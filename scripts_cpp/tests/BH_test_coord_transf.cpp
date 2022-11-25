#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <stack>
#include <stdio.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>
#include <unordered_set>

#include "../causets_cpp/spacetimes.h"
#include "../causets_cpp/functions.h"
#include "../causets_cpp/vecfunctions.h"

#include "boost/math/tools/roots.hpp"
#include "boost/numeric/odeint.hpp"

using namespace boost::numeric::odeint;
using std::cout;
using std::endl;
using std::vector;

int main()
{
    cout<<"\n================= TESTING BH COORD TRANSF ==================\n";
    double mass = 1;

    // SOME COMPLETELY RANDOM COORDINATES IN SCHWARZ
    vector<vector<double>> coords_2D = {  { 0,   1  },
                                          { 0.5, 1  },
                                          { 1,   1  },
                                          {-1,   1  },
                                          { 1,   3  },
                                          { 1,   0.5},
                                          {-1,   3  },
                                          { 1,  10  },
                                          {-5,   7  }};
    
    vector<vector<double>> coords_3D = {{ 0,   1  ,1   },
                                        { 0.5, 1  ,2   },
                                        { 1,   1  ,3   },
                                        {-1,   1  ,1   },
                                        { 1,   3  ,0.80},
                                        { 1,   0.5,0.2 },
                                        {-1,   3  ,1   },
                                        { 1,  10  ,1.5 },
                                        {-5,   7  ,1.2 }};
    
    vector<vector<double>> coords_4D = {  { 0,   1  ,0.3 ,1  },
                                          { 0.5, 1  ,1   ,2  },
                                          { 1,   1  ,0.5 ,3  },
                                          {-1,   1  ,1   ,1  },
                                          { 1,   3  ,0.80,0.80},
                                          { 1,   0.5,1   ,0.2},
                                          {-1,   3  ,1   ,1},
                                          { 1,  10  ,1.5 ,1.5},
                                          {-5,   7  ,1.2 ,1.2}};

    // A COPY OF THE COORDINATES IN SCHWARZ
    vector<vector<double>> S_coords_2D = {  { 0,   1  },
                                          { 0.5, 1  },
                                          { 1,   1  },
                                          {-1,   1  },
                                          { 1,   3  },
                                          { 1,   0.5},
                                          {-1,   3  },
                                          { 1,  10  },
                                          {-5,   7  }};
    
    vector<vector<double>> S_coords_3D = {{ 0,   1  ,1   },
                                        { 0.5, 1  ,2   },
                                        { 1,   1  ,3   },
                                        {-1,   1  ,1   },
                                        { 1,   3  ,0.80},
                                        { 1,   0.5,0.2 },
                                        {-1,   3  ,1   },
                                        { 1,  10  ,1.5 },
                                        {-5,   7  ,1.2 }};
    
    vector<vector<double>> S_coords_4D = {  { 0,   1  ,0.3 ,1  },
                                          { 0.5, 1  ,1   ,2  },
                                          { 1,   1  ,0.5 ,3  },
                                          {-1,   1  ,1   ,1  },
                                          { 1,   3  ,0.80,0.80},
                                          { 1,   0.5,1   ,0.2},
                                          {-1,   3  ,1   ,1},
                                          { 1,  10  ,1.5 ,1.5},
                                          {-5,   7  ,1.2 ,1.2}};
    
    // THEIR EF CORRESPONDENT VALUES
    vector<vector<double>> EForig_coords_2D = { {-1.386294, 1},
                                                {-0.886294, 1},
                                                {-0.386294, 1},
                                                {-2.386294, 1},
                                                {-0.386294, 3},
                                                {0.424636, 0.5},
                                                {-2.386294, 3},
                                                {3.772588 , 10},
                                                {-3.167418, 7} };
    
    vector<vector<double>> EForig_coords_3D = { {-1.386294, 1, 1},
                                                {-0.886294, 1, 2},
                                                {-0.386294, 1, 3},
                                                {-2.386294, 1, 1},
                                                {-0.386294, 3, 0.8},
                                                { 0.424636, 0.5, 0.2},
                                                {-2.386294, 3, 1},
                                                {3.772588 , 10, 1.5},
                                                {-3.167418, 7, 1.2}  };
    
    vector<vector<double>> EForig_coords_4D = { {-1.386294, 1, 0.3, 1},
                                                {-0.886294, 1, 1  , 2},
                                                {-0.386294, 1, 0.5, 3},
                                                {-2.386294, 1, 1  , 1},
                                                {-0.386294, 3, 0.8, 0.8},
                                                {0.424636, 0.5, 1  , 0.2},
                                                {-2.386294, 3 , 1  , 1},
                                                {3.772588 , 10, 1.5,1.5},
                                                {-3.167418, 7 ,1.2,1.2}  };
    

    vector<vector<double>> EFuv_coords_2D = { {-0.386294, 1},
                                              {-0.886294+1, 1},
                                              {-0.386294+1, 1},
                                              {-1.386294, 1},
                                              {-0.386294+3, 3},
                                              {0.924636, 0.5},
                                              {-2.386294+3, 3},
                                              {13.772588 , 10},
                                              {-3.167418+7, 7} };
    
    vector<vector<double>> EFuv_coords_3D = { {-1.386294+1, 1, 1},
                                              {-0.886294+1, 1, 2},
                                              {-0.386294+1, 1, 3},
                                              {-2.386294+1, 1, 1},
                                              {-0.386294+3, 3, 0.8},
                                              { 0.424636+0.5, 0.5, 0.2},
                                              {-2.386294+3, 3, 1},
                                              {3.772588 +10, 10, 1.5},
                                              {-3.167418+7, 7, 1.2}  };
    
    vector<vector<double>> EFuv_coords_4D = { {-1.386294+1, 1, 0.3, 1},
                                              {-0.886294+1, 1, 1  , 2},
                                              {-0.386294+1, 1, 0.5, 3},
                                              {-2.386294+1, 1, 1  , 1},
                                              {-0.386294+3, 3, 0.8, 0.8},
                                              {0.424636+0.5, 0.5, 1  , 0.2},
                                              {-2.386294+3, 3 , 1  , 1},
                                              {3.772588 +10, 10, 1.5,1.5},
                                              {-3.167418+7, 7 ,1.2,1.2}  };

    Spacetime::StoInEF(coords_2D);
    Spacetime::StoInEF(coords_3D);
    Spacetime::StoInEF(coords_4D);

    for (int i = 0; i<coords_2D.size(); i++)
    {   
        for (int j = 0; j<4; j++)
        {
            if (j<2){
            if (std::abs(coords_2D[i][j] - EForig_coords_2D[i][j])>1e-4)
            printf("\n2D: From S to EForig failed at (%i, %i)", i, j);}

            if (j<3){
            if (std::abs(coords_3D[i][j] - EForig_coords_3D[i][j])>1e-4)
            printf("\n3D: From S to EForig failed at (%i, %i)", i, j);}

            if (std::abs(coords_4D[i][j] - EForig_coords_4D[i][j])>1e-4)
            printf("\n4D: From S to EForig failed at (%i, %i)", i, j);
        }
    }

    Spacetime::switchInEF(coords_2D);
    Spacetime::switchInEF(coords_3D);
    Spacetime::switchInEF(coords_4D);
    for (int i = 0; i<coords_2D.size(); i++)
    {   
        for (int j = 0; j<4; j++)
        {
            if (j<2){
            if (std::abs(coords_2D[i][j] - EFuv_coords_2D[i][j])>1e-4)
            printf("\n2D: From orig to uv failed at (%i, %i)", i, j);}

            if (j<3){
            if (std::abs(coords_3D[i][j] - EFuv_coords_3D[i][j])>1e-4)
            printf("\n3D: From orig to uv failed at (%i, %i)", i, j);}

            if (std::abs(coords_4D[i][j] - EFuv_coords_4D[i][j])>1e-4)
            printf("\n4D: From orig to uv failed at (%i, %i)", i, j);
        }
    }

    Spacetime::InEFtoS(coords_2D, 1, "uv");
    Spacetime::InEFtoS(coords_3D, 1, "uv");
    Spacetime::InEFtoS(coords_4D, 1, "uv");
    for (int i = 0; i<coords_2D.size(); i++)
    {   
        for (int j = 0; j<4; j++)
        {
            if (j<2){
            if (std::abs(coords_2D[i][j] - S_coords_2D[i][j])>1e-4)
            printf("\n2D: From uv to S failed at (%i, %i)", i, j);}

            if (j<3){
            if (std::abs(coords_3D[i][j] - S_coords_3D[i][j])>1e-4)
            printf("\n3D: From uv to S failed at (%i, %i)", i, j);}

            if (std::abs(coords_4D[i][j] - S_coords_4D[i][j])>1e-4)
            printf("\n4D: From uv to S failed at (%i, %i)", i, j);
        }
    }


    Spacetime::StoInEF(coords_2D, 1, "uv");
    Spacetime::StoInEF(coords_3D, 1, "uv");
    Spacetime::StoInEF(coords_4D, 1, "uv");
    for (int i = 0; i<coords_2D.size(); i++)
    {   
        for (int j = 0; j<4; j++)
        {
            if (j<2){
            if (std::abs(coords_2D[i][j] - EFuv_coords_2D[i][j])>1e-4)
            printf("\n2D: From S to uv failed at (%i, %i)", i, j);}

            if (j<3){
            if (std::abs(coords_3D[i][j] - EFuv_coords_3D[i][j])>1e-4)
            printf("\n3D: From S to uv failed at (%i, %i)", i, j);}

            if (std::abs(coords_4D[i][j] - EFuv_coords_4D[i][j])>1e-4)
            printf("\n4D: From S to uv failed at (%i, %i)", i, j);
        }
    }


    Spacetime::switchInEF(coords_2D, "uv");
    Spacetime::switchInEF(coords_3D, "uv");
    Spacetime::switchInEF(coords_4D, "uv");
    for (int i = 0; i<coords_2D.size(); i++)
    {   
        for (int j = 0; j<4; j++)
        {
            if (j<2){
            if (std::abs(coords_2D[i][j] - EForig_coords_2D[i][j])>1e-4)
            printf("\n2D: From uv to EForig failed at (%i, %i)", i, j);}

            if (j<3){
            if (std::abs(coords_3D[i][j] - EForig_coords_3D[i][j])>1e-4)
            printf("\n3D: From uv to EForig failed at (%i, %i)", i, j);}

            if (std::abs(coords_4D[i][j] - EForig_coords_4D[i][j])>1e-4)
            printf("\n4D: From uv to EForig failed at (%i, %i)", i, j);
        }
    }


    Spacetime::InEFtoS(coords_2D);
    Spacetime::InEFtoS(coords_3D);
    Spacetime::InEFtoS(coords_4D);
    for (int i = 0; i<coords_2D.size(); i++)
    {   
        for (int j = 0; j<4; j++)
        {
            if (j<2){
            if (std::abs(coords_2D[i][j] - S_coords_2D[i][j])>1e-4)
            printf("\n2D: From orig to S failed at (%i, %i)", i, j);}

            if (j<3){
            if (std::abs(coords_3D[i][j] - S_coords_3D[i][j])>1e-4)
            printf("\n3D: From orig to S failed at (%i, %i)", i, j);}

            if (std::abs(coords_4D[i][j] - S_coords_4D[i][j])>1e-4)
            printf("\n4D: From orig to S failed at (%i, %i)", i, j);
        }
    }

    cout<<"\n======= IF NOTHING PRINTED, EVERYTHING IS FINE =========\n";
    cout<<"(In conversion between S, EForig and EFuv)\n";


    return 0;
}
