/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2015 Scientific Computing and Imaging Institute,
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

#ifndef ALGORITHMS_MATH_CreateImageAlgo_H
#define ALGORITHMS_MATH_CreateImageAlgo_H

#include <Core/Algorithms/Base/AlgorithmBase.h>
#include <Core/Algorithms/Field/share.h>

namespace SCIRun {
namespace Core {
namespace Algorithms {
namespace Fields {
  
  ALGORITHM_PARAMETER_DECL(Width);
  ALGORITHM_PARAMETER_DECL(Height);
  ALGORITHM_PARAMETER_DECL(Depth);
  ALGORITHM_PARAMETER_DECL(PadPercent);
  
  ALGORITHM_PARAMETER_DECL(Mode);
  
  ALGORITHM_PARAMETER_DECL(Axis);
  
  ALGORITHM_PARAMETER_DECL(CenterX);
  ALGORITHM_PARAMETER_DECL(CenterY);
  ALGORITHM_PARAMETER_DECL(CenterZ);
  
  ALGORITHM_PARAMETER_DECL(NormalX);
  ALGORITHM_PARAMETER_DECL(NormalY);
  ALGORITHM_PARAMETER_DECL(NormalZ);
  
  ALGORITHM_PARAMETER_DECL(Position);
  ALGORITHM_PARAMETER_DECL(Index);
  
  ALGORITHM_PARAMETER_DECL(DataLocation);
  
  class SCISHARE CreateImageAlgo : public AlgorithmBase
  {

  public:
    CreateImageAlgo();
    static const AlgorithmInputName OVMatrix;
    static const AlgorithmInputName SizeMatrix;
    virtual AlgorithmOutput run(const AlgorithmInput& input) const override;
    
  private:
    enum DataTypeEnum { SCALAR, VECTOR, TENSOR };
  };

}}}}

#endif
