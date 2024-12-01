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

#include <Modules/ParticleInCell/GodotVR.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/DenseColumnMatrix.h>
#include <Core/Datatypes/MatrixTypeConversions.h>

#include <Core/GeometryPrimitives/Vector.h>
#include <Core/GeometryPrimitives/Point.h>
#include <Core/Datatypes/Legacy/Field/Field.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>

using namespace SCIRun;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Dataflow::Networks; 
using namespace SCIRun::Modules::ParticleInCell;
using namespace SCIRun::Core::Geometry;
using namespace SCIRun::Core::Thread;

using std::cout;
using namespace std;

MODULE_INFO_DEF(GodotVR,ParticleInCell,SCIRun);

GodotVR::GodotVR() : Module(staticInfo_,false) {}

void GodotVR::setStateDefaults() {}

void GodotVR::execute()
    {
    //while (!std::filesystem::exists(vr_out_dir)) std::this_thread::sleep_for(std::chrono::seconds(1));
    while (!std::filesystem::exists(vr_SST_dir)) std::this_thread::sleep_for(std::chrono::seconds(1));

    if(vr_data_count < 1)
        {

        //Direct execution
        //const char *command_dr_vr = vr_GodotLaunch_dir.c_str();
        const char *command_dr_vr = vr_Run_Godot.c_str();
        system(command_dr_vr);
/*
        Script launch exeecution
        const char *command_vr = vr_GodotLaunch.c_str();
        system(command_vr);
        vr_data_count++;
*/
        cout <<"Debug 01"<< "\n";

        }

/*
    else
        {
        const char *command_run = vr_Run_Godot.c_str();
        system(command_run);
        }

        cout <<"Debug 02"<< "\n";
*/
    } //end of GodotVR::execute()

