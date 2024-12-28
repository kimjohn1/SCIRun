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

#ifndef MODULES_PARTICLEINCELL_GodotVR_H
#define MODULES_PARTICLEINCELL_GodotVR_H

#include <Modules/Fields/share.h>
#include <Dataflow/Network/Module.h>
#include <Modules/Basic/share.h>

#include <iostream>
#include <fstream>

#include <Core/Algorithms/Base/AlgorithmVariableNames.h>

namespace SCIRun         {
namespace Modules        {
namespace ParticleInCell {

int vr_data_count                  = 0;
const std::string& vr_home_        = std::getenv("HOME");
const std::string& vr_out_dir      = vr_home_+"/scratch/runs/SST/iteration/00000.png";  // This was changed to include the file name -kj 21 Nov 2024
const std::string& vr_GodotLaunch_dir  = "godot --path "+ vr_home_ + "/Documents/Godot/Projects/simple_wave -r & exit";
const std::string& vr_GodotLaunch  = vr_home_+"/launch_godot.sh";
const std::string& vr_Run_Godot    = vr_home_+"/launch_run_godot.sh";
const std::string& vr_SST_dir      = vr_home_+"/scratch/runs/SST/simOutput/openPMD/simData.sst";

class SCISHARE GodotVR : public SCIRun::Dataflow::Networks::Module,
    public HasNoInputPorts,
    public HasNoOutputPorts
        {
        public:
            GodotVR();
            virtual void execute();
            virtual void setStateDefaults();

            MODULE_TRAITS_AND_INFO(SCIRun::Modules::ModuleFlags::NoAlgoOrUI);
        };
}}}
#endif