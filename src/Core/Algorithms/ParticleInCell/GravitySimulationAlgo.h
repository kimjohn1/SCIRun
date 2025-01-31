
#ifndef CORE_ALGORITHMS_ParticleInCell_GravitySimulationAlgo_H
#define CORE_ALGORITHMS_ParticleInCell_GravitySimulationAlgo_H

#include <sstream>
#include <vector>
#include <algorithm>

#include <Core/Algorithms/Base/AlgorithmVariableNames.h>
#include <Core/Algorithms/Base/AlgorithmBase.h>
#include <Core/Algorithms/ParticleInCell/share.h>

namespace SCIRun         {
namespace Core           {
namespace Algorithms     {
namespace ParticleInCell {

    class SCISHARE GravitySimulationAlgo : public AlgorithmBase
        {
        public:
            GravitySimulationAlgo();
            AlgorithmOutput run(const AlgorithmInput& input) const;
            static const AlgorithmOutputName x_coordinates;
            static const AlgorithmOutputName y_coordinates;
            static const AlgorithmOutputName z_coordinates;
        };

/*
************************ Global variables used by the GravitySim program ************************
*/

using namespace std;

int thread_count        = 8;
const int Dim           = 100;
const int num_particles = 1000000;
const int iterations    = 400;
const int sample_size_p = 100000;
const int sample_size_i = 2;
const int buffer_size   = (iterations/sample_size_i)*(num_particles/sample_size_p);
int j                   = 0;
float delta_t           = 0.01;

//float alpha             = 0.3;

float initial_vel_x     = 1.0;
float initial_vel_y     = 20.0;
float initial_vel_z     = 0.0;

float x_min             = 0.0;
float y_min             = 0.0;
float z_min             = 0.0;

float x_max             = 100.0;
float y_max             = 100.0;
float z_max             = 100.0;

string output_file      = "P2_GravitySim.txt";
//bool save               = 1;                        //here; save is a UI entered value  save = 0 means do not save output to file P2_GravitySim.txt.

                                                      //The Clementine (H1) computer processes 1,600,000 particles in a 10 second run at 100 steps per second
                                                      //in 16 seconds run time using 8 threads, 11 seconds using 16 threads and 8 seconds using 32 threads.
                                                      //H1 processed 4.8Mil particles in 21 seconds, no problems.

                                                      //The HP laptop takes 63 seconds to process 1.6 Mil particles using 1 thread, 26 seconds using 4 threads,
                                                      //and 18 seconds using 8 threads.  Processing 4.8 Mil particles using 8 threads took 70 seconds.

auto pos_x        = new float[num_particles];
auto pos_y        = new float[num_particles];
auto pos_z        = new float[num_particles];

auto buffer_pos_x = new double[buffer_size];
auto buffer_pos_y = new double[buffer_size];
auto buffer_pos_z = new double[buffer_size];

auto vel_x        = new float[num_particles];
auto vel_y        = new float[num_particles];
auto vel_z        = new float[num_particles];
auto acc_x        = new float[num_particles];
auto acc_y        = new float[num_particles];
auto acc_z        = new float[num_particles];
auto acc_field    = new float[Dim][Dim][Dim][3];
auto p_alive      = new int  [num_particles];
auto t_alive      = new int  [thread_count];
auto t_blk_size   = new int  [thread_count];
auto t_index      = new int  [thread_count];

/*
************************ End of global variables used by the GravitySim program ************************
*/

}}}}
#endif
