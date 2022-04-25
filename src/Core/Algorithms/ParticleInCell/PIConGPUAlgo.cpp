/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2020 Scientific Computing and Imaging Institute,
   University of Utah.

   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/


//
#include <openPMD/openPMD.hpp>
//#include <algorithm>
#include <array>
#include <iostream>
#include <memory>
//


#include<Core/Algorithms/ParticleInCell/PIConGPUAlgo.h>
#include<Core/Datatypes/MatrixTypeConversions.h>
#include <chrono>

using namespace SCIRun;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Algorithms::ParticleInCell;

//ALGORITHM_PARAMETER_DEF(ParticleInCell, NumTimeSteps);

const AlgorithmOutputName PIConGPUAlgo::x_coordinates("x_coordinates");
const AlgorithmOutputName PIConGPUAlgo::y_coordinates("y_coordinates");
const AlgorithmOutputName PIConGPUAlgo::z_coordinates("z_coordinates");

PIConGPUAlgo::PIConGPUAlgo()
    {
    addParameter(Variables::Method,0);;
//    addParameter(Parameters::NumTimeSteps,5000);
    }

AlgorithmOutput PIConGPUAlgo::run(const AlgorithmInput&) const
    {
    AlgorithmOutput output;

/*
************************************************ Variable Setup
*/

#include <stdlib.h>
//#include <string>
using namespace std;

    const int buffer_size   = 100;
    auto buffer_pos_x       = new double[buffer_size];
    auto buffer_pos_y       = new double[buffer_size];
    auto buffer_pos_z       = new double[buffer_size];

/*
************************************************ Output Data Setup
*/

    DenseMatrixHandle output_mat_0(new DenseMatrix(buffer_size, 1));
    DenseMatrixHandle output_mat_1(new DenseMatrix(buffer_size, 1));
    DenseMatrixHandle output_mat_2(new DenseMatrix(buffer_size, 1));
                          
/*
************************************************ Simulation Code

The LaserWakefield.py python program is:

import os
import sys
os.system('bash LaserWakefield_compile_run')
*/


    string str="cd $HOME && python3 LaserWakefield_run.py";
	const char *command=str.c_str();
	system(command);

/* Doesn't work:

    string str="cd $HOME && source LaserWakefield_compile_run";
	const char *command=str.c_str();
	system(command);
*/


/*
************************************************End of the simulation code
*/

/*
************************************************Wrap Up
                              includes passing buffered data to SCIRun


    double *data0=output_mat_0->data();
    double *data1=output_mat_1->data();
    double *data2=output_mat_2->data();

    std::copy(buffer_pos_x, buffer_pos_x+buffer_size, data0);
    std::copy(buffer_pos_y, buffer_pos_y+buffer_size, data1);
    std::copy(buffer_pos_z, buffer_pos_z+buffer_size, data2);

    output[x_coordinates] = output_mat_0;
    output[y_coordinates] = output_mat_1;
    output[z_coordinates] = output_mat_2;
*/
    return output;
    }
