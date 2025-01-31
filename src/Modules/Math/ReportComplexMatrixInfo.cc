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


/// @todo Documentation Modules/Math/ReportMatrixInfo.cc

#include <Modules/Math/ReportComplexMatrixInfo.h>
#include <Core/Algorithms/Math/ReportComplexMatrixInfo.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/Scalar.h>

using namespace SCIRun::Modules::Math;
using namespace SCIRun::Dataflow::Networks;
using namespace SCIRun::Core::Datatypes;

MODULE_INFO_DEF(ReportComplexMatrixInfo, Math, SCIRun)

ReportComplexMatrixInfo::ReportComplexMatrixInfo() : Module(staticInfo_)
{
  INITIALIZE_PORT(InputMatrix);
  INITIALIZE_PORT(NumRows);
  INITIALIZE_PORT(NumCols);
  INITIALIZE_PORT(NumElements);
}

void ReportComplexMatrixInfo::execute()
{
  auto matrix = getRequiredInput(InputMatrix);

  if (needToExecute())
  {
    auto output = algo().run(withInputData((InputMatrix, matrix)));
    get_state()->setTransientValue("ReportedInfo", output.getTransient());

    auto info = transient_value_cast<SCIRun::Core::Algorithms::Math::ReportComplexMatrixInfoAlgo::Outputs>(output.getTransient());
    /// @todo: requires knowledge of algorithm type
    sendOutput(NumRows, makeShared<Int32>(info.get<1>()));
    sendOutput(NumCols, makeShared<Int32>(info.get<2>()));
    sendOutput(NumElements, makeShared<Int32>(info.get<3>()));
  }
}
