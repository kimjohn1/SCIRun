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

#include <filesystem>
#include <stdlib.h>

#include <Modules/ParticleInCell/OpenVisusViewer.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/DenseColumnMatrix.h>
#include <Core/Datatypes/MatrixTypeConversions.h>

#include <Core/GeometryPrimitives/Vector.h>
#include <Core/GeometryPrimitives/Point.h>
#include <Core/Datatypes/Legacy/Field/Field.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>

#include <Visus/IdxDataset.h>
#include <Visus/Encoder.h>

using namespace SCIRun;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Dataflow::Networks; 
using namespace SCIRun::Modules::ParticleInCell;
using namespace SCIRun::Core::Geometry;
using namespace SCIRun::Core::Thread;

using std::cout;
using namespace std;
using namespace Visus;

MODULE_INFO_DEF(OpenVisusViewer,ParticleInCell,SCIRun);

OpenVisusViewer::OpenVisusViewer() : Module(staticInfo_) {}

void OpenVisusViewer::setStateDefaults() {}

void OpenVisusViewer::execute()
    {

    while (!std::filesystem::exists(v_idxdataout_dir)) std::this_thread::sleep_for(std::chrono::seconds(1));

    string V_file;
    //V_file = "PYTHONPATH=~/OpenVisus/build/Release python3 -m OpenVisus viewer /dev/shm/idx_data.idx &";
    V_file = v_VisFile;
    const char *command_v = V_file.c_str();
    system(command_v);
    //std::cout << "\n\t" << "Debug 1 First Pass\n";

/*
    //Track run time required to execute using the OpenVisus Viewer
    v_t2 = std::chrono::high_resolution_clock::now();
    float v_duration     = std::chrono::duration_cast<std::chrono::milliseconds>( v_t2 - v_t1 ).count();
    float v_big_duration = std::chrono::duration_cast<std::chrono::milliseconds>( v_t2 - v_big_time ).count();
    //std::cout << "\n\t" << "Debug 3, v_duration is " << v_duration << "\n";
    std::cout << "Viewer visualization time for iteration " << v_data_counter << " is ";
    std::cout << "\t" << v_duration/1000.0 << " seconds\n";
    std::cout << "Total Viewer visualization time is\t\t" << v_big_duration/1000.0 << " seconds\n\n";

    v_out.open(v_out_dir, ios::app);
    v_out << "\nViewer visualization time for iteration " << v_data_counter << " is " << "\t" << v_duration/1000.0 << " seconds\n";
    v_out << "Total Viewer visualization time is\t\t" << v_big_duration/1000.0 << " seconds\n";
    v_out.close();

    v_data_counter++;
    v_t1 = v_t2;
*/
    } //end of OpenVisusViewer::execute()

