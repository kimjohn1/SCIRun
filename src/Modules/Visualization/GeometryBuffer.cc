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


#include <Modules/Visualization/GeometryBuffer.h>
#include <Graphics/Datatypes/GeometryImpl.h>
#include <chrono>

using namespace SCIRun;
using namespace Modules::Visualization;
using namespace Core::Datatypes;
using namespace Dataflow::Networks;
using namespace Core::Algorithms;
using namespace Core::Algorithms::Visualization;
using namespace Core::Geometry;
using namespace Graphics::Datatypes;

ALGORITHM_PARAMETER_DEF(Visualization, MaxBufferSize);
ALGORITHM_PARAMETER_DEF(Visualization, FrameDelay);
ALGORITHM_PARAMETER_DEF(Visualization, SendFlag);

MODULE_INFO_DEF(GeometryBuffer, Visualization, SCIRun)

namespace SCIRun::Modules::Visualization
{
  class GeometryBufferImpl
  {
  public:
    std::map<std::string, std::vector<GeometryBaseHandle>> buffer_;
  };
}

GeometryBuffer::GeometryBuffer() : ModuleWithAsyncDynamicPorts(staticInfo_, true),
  impl_(new GeometryBufferImpl)
{
  INITIALIZE_PORT(GeometryInput);
  INITIALIZE_PORT(GeometryOutputSeries0);
  INITIALIZE_PORT(GeometryOutputSeries1);
  INITIALIZE_PORT(GeometryOutputSeries2);
  INITIALIZE_PORT(GeometryOutputSeries3);
  INITIALIZE_PORT(GeometryOutputSeries4);
  INITIALIZE_PORT(GeometryOutputSeries5);
  INITIALIZE_PORT(GeometryOutputSeries6);
  INITIALIZE_PORT(GeometryOutputSeries7);
}

GeometryBuffer::~GeometryBuffer() = default;

void GeometryBuffer::setStateDefaults()
{
  auto state = get_state();
  state->setValue(Parameters::MaxBufferSize, 50);
  state->setValue(Parameters::FrameDelay, 1.0);
  state->setValue(Parameters::SendFlag, false);
  state->connectSpecificStateChanged(Parameters::SendFlag, [this]()
    {
      Core::Thread::Util::launchAsyncThread([this]() { sendAllGeometries(); });
    });
}

void GeometryBuffer::execute()
{

}

void GeometryBuffer::sendAllGeometries()
{
  using namespace std::chrono_literals;
  auto state = get_state();
  if (state->getValue(Parameters::SendFlag).toBool())
  {
    logCritical("Send all geoms module {}", true);

    int i = 0;
    //TODO: need to transpose this loop
    for (const auto& geomsByPort : impl_->buffer_)
    {

      for (const auto& geom : geomsByPort.second)
      {
        logCritical("Outputting geom on port number {}", i);
        sendOutput(get_output_port(i), geom);
      }
      std::this_thread::sleep_for(10ms);
      i++;
    }
  }
  state->setValue(Parameters::SendFlag, false);
}

void GeometryBuffer::asyncExecute(const PortId& pid, DatatypeHandle data)
{
  logCritical("Received object!");

  const auto geom = std::dynamic_pointer_cast<GeometryObject>(data);
  impl_->buffer_[pid.toString()].push_back(geom);
  logCritical("Buffer is size {}", impl_->buffer_.size());
}

void GeometryBuffer::portRemovedSlotImpl(const PortId& pid)
{
  impl_->buffer_.remove(pid.toString());
}
