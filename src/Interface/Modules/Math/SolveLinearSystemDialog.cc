/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2012 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
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

#include <Interface/Modules/Math/SolveLinearSystemDialog.h>
#include <Core/Algorithms/Math/SolveLinearSystemWithEigen.h>
#include <Dataflow/Network/ModuleStateInterface.h>  //TODO: extract into intermediate
#include <QtGui>

using namespace SCIRun::Gui;
using namespace SCIRun::Dataflow::Networks;
using namespace SCIRun::Core::Algorithms::Math;

SolveLinearSystemDialog::SolveLinearSystemDialog(const std::string& name, ModuleStateHandle state,
  QWidget* parent /* = 0 */)
  : ModuleDialogGeneric(state, parent)
{
  setupUi(this);
  setWindowTitle(QString::fromStdString(name));
  fixSize();
  //executeButton_->setEnabled(false);
  
  //connect(executeButton_, SIGNAL(clicked()), this, SIGNAL(executeButtonPressed()));

  //TODO: clean these up...still getting circles of push/pull
  connect(buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()), this, SLOT(pushParametersToState()));
  connect(targetErrorLineEdit_, SIGNAL(editingFinished()), this, SLOT(pushParametersToState()));
  connect(maxIterationsSpinBox_, SIGNAL(valueChanged(int)), this, SLOT(pushParametersToState()));
}

void SolveLinearSystemDialog::pushParametersToState()
{
  //std::cout << "SLS pushing" << std::endl;
  //TODO: need pattern for this, to avoid silly recursion of push/pull.
  int max = maxIterationsSpinBox_->value();
  if (max != state_->getValue(SolveLinearSystemAlgorithm::MaxIterations).getInt())
    state_->setValue(SolveLinearSystemAlgorithm::MaxIterations, max);

  double error = targetErrorLineEdit_->text().toDouble();
  if (error != state_->getValue(SolveLinearSystemAlgorithm::Tolerance).getDouble())
    state_->setValue(SolveLinearSystemAlgorithm::Tolerance, error);
}

void SolveLinearSystemDialog::pull()
{
  //std::cout << "SLS pulling" << std::endl;
  auto iterations = state_->getValue(SolveLinearSystemAlgorithm::MaxIterations).getInt();
  //std::cout << "SLS state value: " << iterState << std::endl;
  //std::cout << "SLS QString value: " << QString::number(iterState).toStdString() << std::endl;
  
  auto tolerance = state_->getValue(SolveLinearSystemAlgorithm::Tolerance).getDouble();
  //std::cout << "SLS state value: " << stateValue << std::endl;
  //std::cout << "SLS QString value: " << QString::number(stateValue).toStdString() << std::endl;

  maxIterationsSpinBox_->setValue(iterations);
  targetErrorLineEdit_->setText(QString::number(tolerance));
}

