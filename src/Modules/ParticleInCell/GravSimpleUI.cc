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


#include <Modules/ParticleInCell/GravSimpleUI.h>
#include <Core/Datatypes/String.h>

using namespace SCIRun;
using namespace SCIRun::Modules::StringManip;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Dataflow::Networks;

SCIRun::Core::Algorithms::AlgorithmParameterName GravSimpleUI::FormatString("FormatString");

const ModuleLookupInfo GravSimpleUI::staticInfo_("GravSimpleUI","ParticleInCell","SCIRun");

GravSimpleUI::GravSimpleUI() : Module(staticInfo_)
    {
    INITIALIZE_PORT(InputString);
    INITIALIZE_PORT(OutputString);
    }

void GravSimpleUI::setStateDefaults()
    {
    auto state = get_state();
    state->setValue(FormatString,std::string ("[Insert message here]"));
    }

void GravSimpleUI::execute()
    {
    std::string message_string;
    auto stringH = getOptionalInput(InputString);
    auto state = get_state();
    if (stringH && *stringH)
        {
        state -> setValue(FormatString, (*stringH) -> value());
        }

    message_string = state -> getValue(FormatString).toString();
    StringHandle msH(new String(message_string));
    sendOutput(OutputString,msH);
    }
