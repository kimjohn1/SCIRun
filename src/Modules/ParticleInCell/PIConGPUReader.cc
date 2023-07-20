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
#include <stdlib.h>

#include <Modules/ParticleInCell/PIConGPUReader.h>
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
using namespace openPMD;

#define openPMDIsAvailable 1

MODULE_INFO_DEF(PIConGPUReader,ParticleInCell,SCIRun);

PIConGPUReader::PIConGPUReader() : Module(staticInfo_)
    {
    INITIALIZE_PORT(Particles);
    INITIALIZE_PORT(ScalarField);
    INITIALIZE_PORT(VectorField);
    }

void PIConGPUReader::setStateDefaults()
    {
    auto state = get_state();
    state->setValue(Variables::Method, 1);
    state->setValue(Variables::Dim_i_max, 0);
    state->setValue(Variables::Dim_j_max, 0);
    state->setValue(Variables::Dim_k_max, 0);
    state->setValue(Variables::SampleRate, 100);
    state->setValue(Variables::ParticleType, std::string("e"));
    state->setValue(Variables::VectorFieldType, std::string("E"));
    state->setValue(Variables::ScalarFieldComp, std::string("e_all_chargeDensity"));
    }

namespace SCIRun::Modules::ParticleInCell
{

class SimulationStreamingReaderBaseImpl
    {
    public:
#if openPMDIsAvailable
    FieldHandle particleData(int buffer_size, float component_x[], float component_y[], float component_z[])
        {
        FieldInformation pcfi("PointCloudMesh",0,"int");
        FieldHandle ofh = CreateField(pcfi);
        VMesh* omesh    = ofh->vmesh();
        
        for(VMesh::Node::index_type p=0; p < buffer_size; p++) omesh->add_point(Point(component_z[p],component_y[p],component_x[p]));

        delete [] component_x;
        delete [] component_y;
        delete [] component_z;

        return ofh;
        }

    FieldHandle scalarField(std::shared_ptr<float> scalarFieldData_buffer, std::vector<long unsigned int> extent_sFD)
        {
        int sFD_0 = extent_sFD[0];
        int sFD_1 = extent_sFD[1];
        int sFD_2 = extent_sFD[2];

        if(extent_sFD[0] > Dim_i_max && Dim_i_max > 0)
            {
            sFD_0 = Dim_i_max;
            iteration_filter_i = (extent_sFD[0]) / sFD_0;
            }

        if(extent_sFD[1] > Dim_j_max && Dim_j_max > 0)
            {
            sFD_1 = Dim_j_max;
            iteration_filter_j = (extent_sFD[1]) / sFD_1;
            }

        if(extent_sFD[2] > Dim_k_max && Dim_k_max > 0)
            {
            sFD_2 = Dim_k_max;
            iteration_filter_k = (extent_sFD[2]) / sFD_2;
            }

        const int buffer_size_sFD = (sFD_0 * sFD_1 * sFD_2)+1;                    //added a plus 1 here that might not be needed: 6 March - kj
        FieldInformation lfi("LatVolMesh",1,"float");
        std::vector<float> values(buffer_size_sFD);                                               //look at values created here, and maybe in mesh created in the statment below, for the data leak: 6 March - kj

        MeshHandle mesh = CreateMesh(lfi, sFD_0, sFD_1, sFD_2, Point(0.0,0.0,0.0), Point(sFD_0,sFD_1,sFD_2));
        FieldHandle ofh = CreateField(lfi,mesh);                                                  //mesh is used here: 6 March - kj

        for(int i=0; i < sFD_0; i++) for(int j=0; j < sFD_1; j++) for(int k=0; k < sFD_2; k++)
            {
            int flat_index    = (i * iteration_filter_i)*sFD_1*sFD_2+j*sFD_2+k;
            int c_m_index     = k*sFD_0*sFD_1+j*sFD_0+i;
            values[c_m_index] = scalarFieldData_buffer.get()[flat_index];                         //values is used here: 6 March - kj
            }

        VField* ofield = ofh->vfield();
        ofield->set_values(values);

        return ofh;
        }

    FieldHandle vectorField(std::vector<long unsigned int> extent_vFD, std::shared_ptr<float> vFD_component_x, std::shared_ptr<float> vFD_component_y, std::shared_ptr<float> vFD_component_z)
        {

        int vFD_0 = extent_vFD[0];
        int vFD_1 = extent_vFD[1];
        int vFD_2 = extent_vFD[2];

        if(extent_vFD[0] > Dim_i_max && Dim_i_max > 0)
            {
            vFD_0 = Dim_i_max;
            iteration_filter_i = (extent_vFD[0]) / vFD_0;
            }

        if(extent_vFD[1] > Dim_j_max && Dim_j_max > 0)
            {
            vFD_1 = Dim_j_max;
            iteration_filter_j = (extent_vFD[1]) / vFD_1;
            }

        if(extent_vFD[2] > Dim_k_max && Dim_k_max > 0)
            {
            vFD_2 = Dim_k_max;
            iteration_filter_k = (extent_vFD[2]) / vFD_2;
            }

        FieldInformation lfi("LatVolMesh",1,"float");
        lfi.make_vector();
        MeshHandle mesh = CreateMesh(lfi, vFD_0, vFD_1, vFD_2, Point(0.0,0.0,0.0), Point(vFD_0,vFD_1,vFD_2));
        FieldHandle ofh = CreateField(lfi,mesh);
        VField* ofield  = ofh->vfield();

        for(int i=0; i<vFD_0; i++) for(int j=0; j<vFD_1; j++) for(int k=0; k<vFD_2; k++)
            {
            int flat_index = (i * iteration_filter_i)*vFD_1*vFD_2+j*vFD_2+k;
            int c_m_index  = k*vFD_0*vFD_1+j*vFD_0+i;

            Vector v;
            v[0] = vFD_component_x.get()[flat_index];
            v[1] = vFD_component_y.get()[flat_index];
            v[2] = vFD_component_z.get()[flat_index];
            ofield->set_value(v, c_m_index);
            }

        return ofh;
        }

    FieldHandle makeParticleOutput(openPMD::IndexedIteration iteration, int particle_sample_rate, std::string particle_type)
        {

                                                                 //Read particle data
        Record particlePositions       = iteration.particles[particle_type]["position"];
        Record particlePositionOffsets = iteration.particles[particle_type]["positionOffset"];

        std::array<std::shared_ptr<float>, 3> loadedChunks;
        std::array<std::shared_ptr<int>,   3> loadedChunks1;
        std::array<Extent,                 3> extents;
        std::array<std::string,            3> const dimensions{{"x", "y", "z"}};

        for (size_t i_dim = 0; i_dim < 3; ++i_dim)
            {
            std::string dim_str  = dimensions[i_dim];
            RecordComponent rc   = particlePositions[dim_str];
            RecordComponent rc1  = particlePositionOffsets[dim_str];

            loadedChunks[i_dim]  = rc.loadChunk<float>(Offset(rc.getDimensionality(), 0), rc.getExtent());
            loadedChunks1[i_dim] = rc1.loadChunk<int>(Offset(rc1.getDimensionality(), 0), rc1.getExtent());
            extents[i_dim]       = rc.getExtent();
            }

        iteration.seriesFlush();                                 //Data is now available

        Extent const &extent_0 = extents[0];
        int num_particles      = extent_0[0];

        const int buffer_size  = 1+(num_particles/particle_sample_rate);
        auto component_x       = new float[buffer_size];
        auto component_y       = new float[buffer_size];
        auto component_z       = new float[buffer_size];

        for (size_t i_pos = 0; i_pos < 3; ++i_pos)
            {
            std::string dim_str = dimensions[i_pos];
            auto chunk          = loadedChunks[i_pos];
            auto chunk1         = loadedChunks1[i_pos];
                                                                 //Calculate (dimensionless) particle xyz position
            if(i_pos==0) for (size_t k = 0; k<num_particles; k+=particle_sample_rate) component_x[k/particle_sample_rate] = chunk1.get()[k] + chunk.get()[k];
            if(i_pos==1) for (size_t i = 0; i<num_particles; i+=particle_sample_rate) component_y[i/particle_sample_rate] = chunk1.get()[i] + chunk.get()[i];
            if(i_pos==2) for (size_t m = 0; m<num_particles; m+=particle_sample_rate) component_z[m/particle_sample_rate] = chunk1.get()[m] + chunk.get()[m];
            }
                                                                 //Send data to output
        return particleData(buffer_size, component_x, component_y, component_z);
        }

    FieldHandle makeScalarOutput(openPMD::IndexedIteration iteration, std::string scalar_field_component)
        {
                                                                 //Read scalar field data
        auto scalarFieldData        = iteration.meshes[scalar_field_component][MeshRecordComponent::SCALAR];
        auto scalarFieldData_buffer = scalarFieldData.loadChunk<float>();

        iteration.seriesFlush();                                 //Data is now available

        auto extent_sFD             = scalarFieldData.getExtent();

                                                                 //Send data to output
        return scalarField(scalarFieldData_buffer, extent_sFD);
        }

    FieldHandle makeVectorOutput(openPMD::IndexedIteration iteration, std::string vector_field_type)
        {

                                                                 //Read Vector field data
        auto vectorFieldData      = iteration.meshes[vector_field_type];
        auto vFD_component_x      = vectorFieldData["x"].loadChunk<float>();
        auto vFD_component_y      = vectorFieldData["y"].loadChunk<float>();
        auto vFD_component_z      = vectorFieldData["z"].loadChunk<float>();

        iteration.seriesFlush();                                 //Data is now available

        auto extent_vFD           = vectorFieldData["x"].getExtent();

                                                                 //Send data to output
        return vectorField(extent_vFD, vFD_component_x, vFD_component_y, vFD_component_z);
#endif
        }
    }; //end of class SimulationStreamingReaderBaseImpl
} //end of namespace SCIRun::Modules::ParticleInCell

void PIConGPUReader::execute()
    {
    AlgorithmInput input;
    SimulationStreamingReaderBaseImpl P;
    if (!setup_) setupStream();


    t2 = std::chrono::high_resolution_clock::now();                                                      //here
    float duration     = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();       //here
    float big_duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - big_time ).count();
    std::cout << "Visualization time for iteration " << data_counter << " is ";                          //here
    std::cout << "\t" << duration/1000.0 << " seconds\n";                                                //here
    std::cout << "Total visualization time is\t\t" << big_duration/1000.0 << " seconds\n\n";

    //vis_out.open("visout.txt", ios::app);                                                               //here, out
    //vis_out << "\nTotal visualization time is\t\t" << big_duration/1000.0 << " seconds\n";              //here, out
    //vis_out.close();                                                                                    //here, out

    t1 = std::chrono::high_resolution_clock::now();                                                      //here
    data_counter++;                                                                                      //here


#if openPMDIsAvailable
    IndexedIteration iteration = *it;

    if(particlesPresent  ) sendOutput(Particles,   P.makeParticleOutput(iteration, SampleRate, ParticleType));
    if(scalarFieldPresent) sendOutput(ScalarField, P.makeScalarOutput(iteration,   ScalarFieldComp));
    if(vectorFieldPresent) sendOutput(VectorField, P.makeVectorOutput(iteration,   VectorFieldType));
    iteration.close();
#endif
    ++it;
#if openPMDIsAvailable
    if(it != end) enqueueExecuteAgain(false);
    else shutdownStream();
#endif

    }

void PIConGPUReader::setupStream()
    {
    auto state      = get_state();
    DataSet         = state->getValue(Variables::Method).toInt();
    Dim_i_max       = state->getValue(Variables::Dim_i_max).toInt();
    Dim_j_max       = state->getValue(Variables::Dim_j_max).toInt();
    Dim_k_max       = state->getValue(Variables::Dim_k_max).toInt();
    SampleRate      = state->getValue(Variables::SampleRate).toInt();
    ParticleType    = state->getValue(Variables::ParticleType).toString();
    ScalarFieldComp = state->getValue(Variables::ScalarFieldComp).toString();
    VectorFieldType = state->getValue(Variables::VectorFieldType).toString();

    while (!std::filesystem::exists(SST_dir)) std::this_thread::sleep_for(std::chrono::seconds(1));


    t1       = std::chrono::high_resolution_clock::now();        //here
    big_time = t1;                                               //here


#if openPMDIsAvailable
    series = Series(SST_dir, Access::READ_ONLY);
    end    = series.readIterations().end();
    it     = series.readIterations().begin();
    setup_ = true;
    IndexedIteration iter_ss = *it;

    if(iter_ss.particles.size())
        for (auto const &ps : iter_ss.particles)
            if(ps.first == ParticleType)      particlesPresent = true;

    if(iter_ss.meshes.size())
        for (auto const &pm : iter_ss.meshes)
            {
            if(pm.first == ScalarFieldComp) scalarFieldPresent = true;
            if(pm.first == VectorFieldType) vectorFieldPresent = true;
            }

    if(DataSet==0) showDataSet();
#endif
    }

void PIConGPUReader::shutdownStream()
    {
    string text_file;
    text_file = "rm ~/picongpu.profile ~/picongpu_reRun.profile ~/Sim.py ~/Sim_run";
    const char *command_shutDown = text_file.c_str();
    system(command_shutDown);

    setup_             = false;
    particlesPresent   = false;
    vectorFieldPresent = false;
    scalarFieldPresent = false;
    }

void PIConGPUReader::showDataSet()
    {
#if openPMDIsAvailable
    IndexedIteration iteration_00 = *it;
    cout << "\nShowing the data set content for iteration " << iteration_00.iterationIndex << "\n";

                                                //Output data about the Series

    Iteration iter = series.iterations[iteration_00.iterationIndex];
    cout << "Iteration " << iteration_00.iterationIndex << " contains "
         << iter.meshes.size()    << " meshes " << "and "
         << iter.particles.size() << " particle species\n";

                                                //Output data about particles
    if(iter.particles.size())
        {
        cout << "\nParticle data \n";
        for (auto const &ps : iter.particles)
            {
            cout << "\nSpecies\t" << ps.first;
            for (auto const &r : ps.second) cout << "\n\t" << r.first;
            }
        cout << "\nParticle visualization will use a sampled subset of particles using a sampling rate of: " << SampleRate << "\n";
        }
    else cout << "\nThere is no particle data in this data set\n";

                                                //Output data about meshes
    if(iter.meshes.size())
        {
        cout << "Mesh data \n";
        for (auto const &pm : iter.meshes) cout << "\n\t" << pm.first;
        cout << "\n";

        if(vectorFieldPresent)
            {
            MeshRecordComponent FieldType_x = iter.meshes[VectorFieldType]["x"];
            Extent extent_FieldType = FieldType_x.getExtent();
            cout << "\nField " << VectorFieldType << " is vector valued, has shape (";
            for (auto const &dim : extent_FieldType) cout << dim << ',';
            cout << ") and datatype " << FieldType_x.getDatatype() << '\n';

            cout << "Vector field data output to the vectorField output port has shape ";
            if(extent_FieldType[0] > Dim_i_max) cout << "(" << Dim_i_max << ", " << extent_FieldType[1] << ", " << extent_FieldType[2] << ")\n";
            else cout << "Same as listed above\n";
            }

        if(scalarFieldPresent)
            {
            MeshRecordComponent E_density = iter.meshes[ScalarFieldComp][MeshRecordComponent::SCALAR];
            Extent extent_d = E_density.getExtent();
            cout << "\nField " << ScalarFieldComp << " is scalar valued, has shape (";
            for (auto const &dim : extent_d) cout << dim << ',';
            cout  << ") and datatype " << E_density.getDatatype() << '\n';

            cout << "Scalar field data output to the scalarField output port has shape ";
            if(extent_d[0] > Dim_i_max) cout << "(" << Dim_i_max << ", " << extent_d[1] << ", " << extent_d[2] << ")\n";
            else cout << "Same as listed above\n";
            }
        }
    else cout << "\nThere is no mesh data in this data set\n";
    cout << "\n";
#endif
    }

