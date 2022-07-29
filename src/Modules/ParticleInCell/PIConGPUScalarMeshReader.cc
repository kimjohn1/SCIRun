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

#include <openPMD/openPMD.hpp>
#include <filesystem>

#include <Modules/ParticleInCell/PIConGPUScalarMeshReader.h>
#include <Core/Algorithms/ParticleInCell/PIConGPUScalarMeshReaderAlgo.h>


using namespace SCIRun;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Modules::ParticleInCell;
using namespace SCIRun::Core::Algorithms::ParticleInCell;

using std::cout;
using namespace openPMD;

MODULE_INFO_DEF(PIConGPUScalarMeshReader,ParticleInCell,SCIRun);

const AlgorithmOutputName PIConGPUScalarMeshReaderAlgo::ScalarMesh_value("ScalarMesh_value");
const AlgorithmOutputName PIConGPUScalarMeshReaderAlgo::ScalarMesh_unused1("ScalarMesh_unused1");
const AlgorithmOutputName PIConGPUScalarMeshReaderAlgo::ScalarMesh_unused2("ScalarMesh_unused2");

PIConGPUScalarMeshReader::PIConGPUScalarMeshReader() : Module(staticInfo_)
    {
    INITIALIZE_PORT(ScalarMesh_value);
    INITIALIZE_PORT(ScalarMesh_unused1);
    INITIALIZE_PORT(ScalarMesh_unused2);
    }

void PIConGPUScalarMeshReader::setStateDefaults() {}

void PIConGPUScalarMeshReader::execute()
    {
    AlgorithmInput input;
    if(needToExecute())
        {
        auto state = get_state();
        auto output=algo().run(input);

/*
************************************************Start the openPMD Reader function and loop
*/
                                                        //Wait for simulation output data to be generated and posted via SST
                                                        // TODO: figure out how to use a general reference for the home directory in these two lines of code
        while(!std::filesystem::exists("/home/kj/scratch/runs/SST/simOutput/openPMD/simData.sst")) sleep(1);
//        while(!std::filesystem::exists("scratch/runs/SST/simOutput/openPMD/simData.sst")) sleep(1);

        Series series = Series("/home/kj/scratch/runs/SST/simOutput/openPMD/simData.sst", Access::READ_ONLY);
//        Series series = Series("scratch/runs/SST/simOutput/openPMD/simData.sst", Access::READ_ONLY);

        for (IndexedIteration iteration : series.readIterations())
            {
            cout << "\nFrom PIConGPUScalarMeshReader: Current iteration is: " << iteration.iterationIndex << std::endl;

//*********************Preamble: Output information about the Series content to the terminal
/*
                                                        //From https://openpmd-api.readthedocs.io/en/latest/usage/serial.html#c
            Iteration iter = series.iterations[iteration.iterationIndex];
            cout << "Iteration " << iteration.iterationIndex << " contains "
                 << iter.meshes.size()    << " meshes " << "and "
                 << iter.particles.size() << " particle species\n";
            cout << "The Series contains " << series.iterations.size() << " iterations\n";

                                                        //Output data about particles
            cout << "\nParticle data \n";
            for (auto const &ps : iter.particles)
                {
                cout << "\n\t" << ps.first;
                cout << "\n";
                for (auto const &r : ps.second) cout << "\n\t" << r.first;
                }
            cout << '\n';

                                                        //Output data about meshes
            cout << "\nMesh data \n";

            for (auto const &pm : iter.meshes) cout << "\n\t" << pm.first;
            cout << "\n";

            MeshRecordComponent B_x = iter.meshes["B"]["x"];
            Extent extent_B = B_x.getExtent();
            cout << "\nField B_x has shape (";
            for (auto const &dim : extent_B) cout << dim << ',';
            cout << ") and has datatype " << B_x.getDatatype() << '\n';

            MeshRecordComponent E_charge_density = iter.meshes["e_all_chargeDensity"][MeshRecordComponent::SCALAR];
            Extent extent_cd = E_charge_density.getExtent();
            cout << "\nField E_charge_density has shape (";
            for (auto const &dim : extent_cd) cout << dim << ',';
            cout  << ") and has datatype " << E_charge_density.getDatatype() << '\n';
            cout << "\n----------" << std::endl;

//*********************End of Preamble

                                                 //Load particles xyz position (back to https://openpmd-api.readthedocs.io/en/latest/usage/streaming.html#c)
            Record electronPositions = iteration.particles["e"]["position"];
            std::array<std::shared_ptr<position_t>, 3> loadedChunks;
            std::array<Extent, 3> extents;
            std::array<std::string, 3> const dimensions{{"x", "y", "z"}};

            for (size_t i_dim = 0; i_dim < 3; ++i_dim)
                {
                std::string dim_str = dimensions[i_dim];
                RecordComponent rc = electronPositions[dim_str];
                loadedChunks[i_dim] = rc.loadChunk<position_t>(Offset(rc.getDimensionality(), 0), rc.getExtent());
                extents[i_dim] = rc.getExtent();
                }

                                                        //Load mesh data; ijk values at xyz node points (from Franz Poschel email, 17 May 2022)
            auto mesh = iteration.meshes["E"];
            auto E_x = mesh["x"].loadChunk<float>();
            auto E_y = mesh["y"].loadChunk<float>();
            auto E_z = mesh["z"].loadChunk<float>();
            iteration.seriesFlush();                    //Data is now available


                                                        //Add code to extract scalar values from e_all_chargeDensity mesh here


                                                        //Output some useful particle information to the terminal

            cout << "\nAfter loading particle position data\n";
            Extent const &extent_0 = extents[0];
            int num_particles = extent_0[0];
            cout << "\nNumber of particles is " << num_particles;

                                                        //Set the particle sample rate (this may eventually be a config file setting)
//            int particle_sample_rate = 1000000;         //Number of samples in the final frame is 4
//            int particle_sample_rate = 100000;          //Number of samples in the final frame is 38
            int particle_sample_rate = 100;             //Number of samples in the final frame is 37154
            cout << "\nParticle sample_rate is " << particle_sample_rate << "\n";
            cout << "The number of particles sampled is " << 1+(num_particles/particle_sample_rate) << "\n";

            iteration.close();

                                                      //testing and debug: Output the single 3 dim vector from E field data at a single point
            cout << "\nAfter loading Mesh data\n";

            auto extent_x = mesh["x"].getExtent();
            size_t i = 6;
            size_t j = 512;
            size_t k = 86;
            size_t flat_index = i * extent_x[1] * extent_x[2] + j * extent_x[2] + k;
            cout << "\nxyz values at mesh E node point (" << i << ", " << j << ", " << k << ") are:\n";
            cout << "\t x: " << E_x.get()[flat_index] << "\t y: " << E_y.get()[flat_index] << "\t z: " << E_z.get()[flat_index] << "\n----------\n";


//    ***************************************************** Set up and load the module output buffers


            const int buffer_size   = 1+(num_particles/particle_sample_rate);
            auto buffer_pos_x       = new double[buffer_size];
            auto buffer_pos_y       = new double[buffer_size];
            auto buffer_pos_z       = new double[buffer_size];

            for (size_t i_pos = 0; i_pos < 3; ++i_pos)
                {
                std::string dim = dimensions[i_pos];
                auto chunk = loadedChunks[i_pos];

                cout <<"\nThe sampled values for particle position in dimension " << dim << " are\n";
                cout <<"not printed\n";
//                for (size_t j = 0; j<num_particles ; j+=particle_sample_rate) cout << "\t" << chunk.get()[j] << ", ";
//                cout << "\n----------" << std::endl;

                if(i_pos==0) for (size_t k = 0; k<num_particles ; k+=particle_sample_rate) buffer_pos_x[k/particle_sample_rate]=chunk.get()[k];
                if(i_pos==1) for (size_t i = 0; i<num_particles ; i+=particle_sample_rate) buffer_pos_y[i/particle_sample_rate]=chunk.get()[i];
                if(i_pos==2) for (size_t m = 0; m<num_particles ; m+=particle_sample_rate) buffer_pos_z[m/particle_sample_rate]=chunk.get()[m];
                }


//    *****************************************************  Set up the output data


            DenseMatrixHandle output_mat_0(new DenseMatrix(buffer_size, 1));
            DenseMatrixHandle output_mat_1(new DenseMatrix(buffer_size, 1));
            DenseMatrixHandle output_mat_2(new DenseMatrix(buffer_size, 1));

            double *data0=output_mat_0->data();
            double *data1=output_mat_1->data();
            double *data2=output_mat_2->data();

            std::copy(buffer_pos_x, buffer_pos_x+buffer_size, data0);
            std::copy(buffer_pos_y, buffer_pos_y+buffer_size, data1);
            std::copy(buffer_pos_z, buffer_pos_z+buffer_size, data2);


//    *****************************************************  Send data to the output ports
    

            sendOutput(ScalarMesh_value, output_mat_0);
            sendOutput(ScalarMesh_unused1, output_mat_1);
            sendOutput(ScalarMesh_unused2, output_mat_2);
*/                
            }  //end of the openPMD reader loop
        }  //end of the "needToExecute" block
    }  //end of the "PIConGPU::execute()" function

