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

#ifndef MODULES_PARTICLEINCELL_PIConGPUReader_H
#define MODULES_PARTICLEINCELL_PIConGPUReader_H

#include <openPMD/openPMD.hpp>
#include <Modules/Fields/share.h>
#include <Dataflow/Network/Module.h>
#include <Modules/Basic/share.h>

#include <iostream>                                              //here, out
#include <fstream>                                               //here, out

#include <Core/Algorithms/Base/AlgorithmVariableNames.h>

namespace SCIRun         {
namespace Modules        {
namespace ParticleInCell {

using namespace openPMD;
#define openPMDIsAvailable 1

#if openPMDIsAvailable
Series series;
SeriesIterator it, end;
#endif

int data_counter           = 0;                                  //here
int iteration_filter_i     = 1;
int iteration_filter_j     = 1;
int iteration_filter_k     = 1;
bool setup_                = false;
bool particlesPresent      = false;
bool vectorFieldPresent    = false;
bool scalarFieldPresent    = false;
const std::string& home_   = std::getenv("HOME");
const std::string& SST_dir = home_+"/scratch/runs/SST/simOutput/openPMD/simData.sst";
const std::string& visout_dir = home_+"/scratch/runs/SST/simOutput/visout.txt";
const std::string& rawdataout_dir = home_+"/scratch/runs/SST/simOutput/raw_data_out.bin";     //the raw data output task 1 Nov
const std::string& idxdataout_dir = "/dev/shm/idx_data.idx";                                  //the idx data output task 9 Nov
const std::string& VisFile = home_+"/launch_visus.sh";

std::string stringDirCreate;
std::string stringDirCreate2;
std::string stringDir;
std::string stringDirRemove;
std::string stringDirRemoveZip;                                  //the SF2PNG task 23 Nov 2024
std::string stringDirZip;                                        //

bool        DataSet1;
bool        DataSet2;
bool        DataSet3;
int         Dim_i_max;
int         Dim_j_max;
int         Dim_k_max;
int         SampleRate;
std::string ParticleType;
std::string ScalarFieldComp;
std::string VectorFieldType;


std::ofstream vis_out;                                           //here, out
//std::ofstream v_out;


auto t1       = std::chrono::high_resolution_clock::now();       //here
auto t2       = std::chrono::high_resolution_clock::now();       //here
auto big_time = std::chrono::high_resolution_clock::now();       //here


class SCISHARE PIConGPUReader : public SCIRun::Dataflow::Networks::Module,
    public HasNoInputPorts,
    public Has3OutputPorts<FieldPortTag, FieldPortTag, FieldPortTag>
        {
        public:
            PIConGPUReader();
            void setupStream();
            void showDataSet();
            void shutdownStream();
            virtual void execute();
            virtual void setStateDefaults();

            OUTPUT_PORT(0, Particles,   Field);
            OUTPUT_PORT(1, ScalarField, Field);
            OUTPUT_PORT(2, VectorField, Field);

            MODULE_TRAITS_AND_INFO(SCIRun::Modules::ModuleFlags::ModuleHasUI);

        private:
            std::unique_ptr<class SimulationStreamingReaderBaseImpl> impl_;
        };
}}}
#endif
