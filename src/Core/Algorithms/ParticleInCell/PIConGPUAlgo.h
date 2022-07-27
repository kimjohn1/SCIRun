
#ifndef CORE_ALGORITHMS_ParticleInCell_PIConGPUAlgo_H
#define CORE_ALGORITHMS_ParticleInCell_PIConGPUAlgo_H

#include <string>
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

ALGORITHM_PARAMETER_DECL(SimulationFile);
ALGORITHM_PARAMETER_DECL(ConfigFile);
ALGORITHM_PARAMETER_DECL(CloneDir);
ALGORITHM_PARAMETER_DECL(OutputDir);

class SCISHARE PIConGPUAlgo : public AlgorithmBase
    {
    public:
        PIConGPUAlgo();
        AlgorithmOutput run(const AlgorithmInput& input) const;
//            static const AlgorithmOutputName x_coordinates;
//            static const AlgorithmOutputName y_coordinates;
//            static const AlgorithmOutputName z_coordinates;
    private:
        bool StartPIConGPU(const std::string, const std::string, const std::string, const std::string, const int) const;
 //       bool simulationstarted_;
    };

}}}}
#endif
