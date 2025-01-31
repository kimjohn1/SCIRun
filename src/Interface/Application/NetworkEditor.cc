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


#include <sstream>
#include <Interface/qt_include.h>
#include <QtConcurrent>
#include <Interface/Application/NetworkEditor.h>
#include <Interface/Application/Connection.h>
#include <Interface/Application/ModuleWidget.h>
#include <Interface/Application/ModuleProxyWidget.h>
#include <Interface/Application/Utility.h>
#include <Interface/Application/Port.h>
#include <Interface/Application/PortWidgetManager.h>
#include <Interface/Application/GuiLogger.h>
#include <Interface/Application/Subnetworks.h>
#include <Interface/Application/StateViewer.h>
#include <Interface/Application/NetworkEditorControllerGuiProxy.h>
#include <Interface/Application/ClosestPortFinder.h>
#include <Dataflow/Serialization/Network/NetworkDescriptionSerialization.h>
#include <Dataflow/Engine/Controller/NetworkEditorController.h> //TODO: remove
#include <Dataflow/Network/Module.h> //TODO: remove
#include <Core/Application/Preferences/Preferences.h>
#include <Core/Application/Application.h>
#include <Dataflow/Serialization/Network/XMLSerializer.h>
#include <Interface/Application/MainWindowCollaborators.h>
#ifdef BUILD_WITH_PYTHON
#include <Dataflow/Engine/Python/NetworkEditorPythonAPI.h>
#endif

#include <boost/algorithm/string/find.hpp>

using namespace SCIRun;
using namespace Core;
using namespace Logging;
using namespace Algorithms;
using namespace Gui;
using namespace NetworkBoundaries;
using namespace Dataflow::Networks;
using namespace Dataflow::Engine;

NetworkEditor::NetworkEditor(const NetworkEditorParameters& params, QWidget* parent)
  : QGraphicsView(parent),
  ctorParams_(params),
  tagColor_(params.tagColor),
  tagName_(params.tagName),
  moduleSelectionGetter_(params.moduleSelectionGetter),
  defaultNotePositionGetter_(params.dnpg),
  preexecute_(params.preexecuteFunc),
  highResolutionExpandFactor_(params.highResolutionExpandFactor),
  scene_(new QGraphicsScene(parent)),
  moduleEventProxy_(new ModuleEventProxy),
  zLevelManager_(new ZLevelManager(scene_))
{
  setBackgroundBrush(QPixmap(networkBackgroundImage()));

  scene_->setItemIndexMethod(QGraphicsScene::NoIndex);

  Preferences::Instance().forceGridBackground.connectValueChanged([this](bool value) { updateBackground(value); });
  Preferences::Instance().moduleExecuteDownstreamOnly.connectValueChanged([this](bool value) { updateExecuteButtons(value); });

  setHighResolutionExpandFactor(highResolutionExpandFactor_);

  setScene(scene_);
  setDragMode(RubberBandDrag);
  setAcceptDrops(true);
  setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);

  setSceneRect(QRectF());
  centerOn(100, 100);

  setMouseAsDragMode();

#ifdef BUILD_WITH_PYTHON
  NetworkEditorPythonAPI::setExecutionContext(this);
#endif

  connect(this, &NetworkEditor::moduleMoved, this, &NetworkEditor::redrawTagGroups);

  setObjectName(QString::fromUtf8("networkEditor_"));
  setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
  setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
  verticalScrollBar()->setValue(0);
  horizontalScrollBar()->setValue(0);

#ifdef MODULE_POSITION_LOGGING
  setViewUpdateFunc([](const QString& q) { qDebug() << q; });
#endif

  if (allowModificationSignalConnection())
    connect(this, &NetworkEditor::modified, [this]() { setSceneRect(QRectF()); });
}

void NetworkEditor::setHighResolutionExpandFactor(double factor)
{
  highResolutionExpandFactor_ = factor;
  ModuleWidget::highResolutionExpandFactor_ = highResolutionExpandFactor_;
}

void NetworkEditor::setNetworkEditorController(SharedPointer<NetworkEditorControllerGuiProxy> controller)
{
  if (controller_ == controller)
    return;

  if (controller_)
  {
    disconnect(controller_.get(), &NetworkEditorControllerGuiProxy::moduleAdded,
      this, &NetworkEditor::addModuleWidget);

    disconnect(this, &NetworkEditor::connectionDeleted,
      controller_.get(), &NetworkEditorControllerGuiProxy::removeConnection);

    disconnect(controller_.get(), &NetworkEditorControllerGuiProxy::connectionAdded,
      this, &NetworkEditor::connectionAddedQueued);

    controller_->setExecutableLookup(nullptr);
  }

  controller_ = controller;

  if (controller_)
  {
    connect(controller_.get(), &NetworkEditorControllerGuiProxy::moduleAdded, this, &NetworkEditor::addModuleWidget);

    connect(this, &NetworkEditor::connectionDeleted,
      controller_.get(), &NetworkEditorControllerGuiProxy::removeConnection);

    connect(controller_.get(), &NetworkEditorControllerGuiProxy::connectionAdded,
      this, &NetworkEditor::connectionAddedQueued);

    connect(controller_.get(), &NetworkEditorControllerGuiProxy::connectionStatusChanged,
      this, &NetworkEditor::connectionStatusChanged);

    controller_->setExecutableLookup(this);
  }
}

SharedPointer<NetworkEditorControllerGuiProxy> NetworkEditor::getNetworkEditorController() const
{
  return controller_;
}

NetworkEditor::ViewUpdateFunc NetworkEditor::viewUpdateFunc_;
  // = [](const QString& s) { logCritical("Gui info: {}", s.toStdString()); };
QGraphicsView* NetworkEditor::miniview_ {nullptr};

static const int macModulePositionWorkaroundTimerValue = 5;

void NetworkEditor::addModuleWidget(const std::string& name, ModuleHandle module, const ModuleCounter& count)
{
  if (!fileLoading_ && inEditingContext_ != this)
  {
    //std::cout << "\t\tnew condition: !" << fileLoading_ << " and " << inEditingContext_ << " != " << this << std::endl;
    return;
  }

  latestModuleId_ = module->id().id_;
  auto moduleWidget = new ModuleWidget(this, QString::fromStdString(name), module);
  moduleEventProxy_->trackModule(module);

#ifdef MODULE_POSITION_LOGGING

#endif
  auto proxy = setupModuleWidget(moduleWidget);

#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();

  logCritical("module {} added at pos {},{} proxy pos {}, {} proxy scenePos {},{}",
    name,
    moduleWidget->pos().x(), moduleWidget->pos().y(),
    proxy->pos().x(), proxy->pos().y(),
    proxy->scenePos().x(), proxy->scenePos().y());
#endif

  if (!fileLoading_)
  {
    moduleWidget->postLoadAction();
  }
#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  count.increment();
  Q_EMIT modified();
  Q_EMIT newModule(QString::fromStdString(module->id()), module->hasUI());

#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

#ifdef __APPLE__
  QTimer::singleShot(macModulePositionWorkaroundTimerValue, [proxy]()
    {
      proxy->setSelected(true);
      proxy->setSelected(false);
      proxy->setSelected(true);
      proxy->show();
    }
  );
#endif

#if 0
  if (name.find("Loop") != std::string::npos)
  {
    qDebug() << "Special shape logic for loopers";
    auto diamondBackground = new LoopDiamondPolygon();
    scene_->addItem(diamondBackground);
    ensureVisible(diamondBackground);
    diamondBackground->setPos(proxy->pos());
    diamondBackground->setZValue(-99999);
    proxy->setBackgroundPolygon(diamondBackground);
  }
#endif
}

void NetworkEditor::connectionAddedQueued(const ConnectionDescription&)
{
  //std::cout << "Received queued connection request: " << ConnectionId::create(cd).id_ << std::endl;
}

SharedPointer<DisableDynamicPortSwitch> NetworkEditor::createDynamicPortDisabler()
{
  return controller_->createDynamicPortSwitch();
}

std::optional<ConnectionId> NetworkEditor::requestConnection(const PortDescriptionInterface* from, const PortDescriptionInterface* to)
{
  auto id = controller_->requestConnection(from, to);
  Q_EMIT modified();
  return id;
}

std::optional<ConnectionId> NetworkEditor::requestConnectionWidget(const PortWidget* from, const PortWidget* to)
{
  return requestConnection(from->description(), to->description());
}

namespace
{
  const int TagTextKey = 123;

  ModuleProxyWidget* findById(const QList<QGraphicsItem*>& list, const std::string& id)
  {
    Q_FOREACH(QGraphicsItem* item, list)
    {
      if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
      {
        if (id == w->getModuleWidget()->getModuleId())
          return w;
      }
    }
    return nullptr;
  }

  ModuleProxyWidget* findFirstByName(const QList<QGraphicsItem*>& list, const std::string& name)
  {
    Q_FOREACH(QGraphicsItem* item, list)
    {
      if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
      {
        if (w->getModuleWidget()->getModuleId().find(name) != std::string::npos)
          return w;
      }
    }
    return nullptr;
  }
}

void NetworkEditor::duplicateModule(const ModuleHandle& module)
{
  InEditingContext iec(this);

  auto widget = findById(scene_->items(), module->id());
  modulePlacement_.updateLatestFromDuplicate(widget->scenePos());
  //TODO: need better duplicate placement. hard code it for now.
  controller_->duplicateModule(module);
}

namespace
{
  static const QPointF incr1 {100, 0};
  static const QPointF replaceIncr {-15, -15};
}

QPointF ModuleWidgetPlacementManager::getLastForDoubleClickedItem(const QPointF& sceneCenter) const
{
  static int counter = 0;
  counter = (counter + 1) % 5;
  double coord = 10*counter;
  return sceneCenter + QPointF{coord, coord};
}

QPointF ModuleWidgetPlacementManager::connectNewIncrement(bool isInput)
{
  return {0.0, isInput ? -90.0 : 90.0};
}

void ModuleWidgetPlacementManager::updateLatestFromDuplicate(const QPointF& scenePos)
{
  lastModulePosition_ = scenePos + incr1;
}

void ModuleWidgetPlacementManager::updateLatestFromConnectNew(const QPointF& scenePos, bool isInputPort)
{
  lastModulePosition_ = scenePos + connectNewIncrement(isInputPort);
}

void ModuleWidgetPlacementManager::updateLatestFromReplace(const QPointF& scenePos)
{
  lastModulePosition_ = scenePos + replaceIncr;
}

void NetworkEditor::insertNewModule(const ModuleHandle& moduleToConnectTo, const PortDescriptionInterface* portToConnect, const QMap<QString, std::string>& info)
{
  auto startWidget = findById(scene_->items(), moduleToConnectTo->id());

  if (startWidget)
  {
    InEditingContext iec(this);
    modulePlacement_.updateLatestFromConnectNew(startWidget->scenePos(), portToConnect->isInput());
    controller_->insertNewModule(portToConnect, info);
  }
}

void NetworkEditor::connectNewModule(const ModuleHandle& moduleToConnectTo, const PortDescriptionInterface* portToConnect, const std::string& newModuleName)
{
  connectNewModuleImpl(moduleToConnectTo, portToConnect, newModuleName);
}

void NetworkEditor::connectNewModuleImpl(const ModuleHandle& moduleToConnectTo, const PortDescriptionInterface* portToConnect, const std::string& newModuleName)
{
  auto widget = findById(scene_->items(), moduleToConnectTo->id());

  if (widget)
  {
    InEditingContext iec(this);
    modulePlacement_.updateLatestFromConnectNew(widget->scenePos(), portToConnect->isInput());
    controller_->connectNewModule(portToConnect, newModuleName);
    return;
  }

#if 0
  for (auto& child : childrenNetworks_)
  {
    child.second->get()->connectNewModuleImpl(moduleToConnectTo, portToConnect, newModuleName);
  }
#endif
}

void NetworkEditor::replaceModuleWith(const ModuleHandle& moduleToReplace, const std::string& newModuleName)
{
  InEditingContext iec(this);

  auto oldModule = findById(scene_->items(), moduleToReplace->id());
  modulePlacement_.updateLatestFromReplace(oldModule->scenePos());
  controller_->addModule(newModuleName);

  // connect up same ports
  auto newModule = findById(scene_->items(), latestModuleId_);

  const auto& oldModPorts = oldModule->getModuleWidget()->ports();
  const auto& newModPorts = newModule->getModuleWidget()->ports();

  {
    size_t nextInputIndex = 0;
    for (const auto& iport : oldModPorts.inputs())
    {
      if (iport->isConnected())
      {
        const auto& newInputs = newModPorts.inputs();
        auto toConnect = std::find_if(newInputs.begin(), newInputs.end(),
          [&](const PortWidget* port) { return port->description()->get_typename() == iport->description()->get_typename() && port->description()->getIndex() >= nextInputIndex; });
        if (toConnect == newInputs.end())
        {
          guiLogCritical("Logical error: could not find input port to connect to {}, {}", iport->name().toStdString(), nextInputIndex);
          break;
        }
        nextInputIndex = (*toConnect)->description()->getIndex() + 1;
        requestConnectionWidget(iport->connectedPorts()[0], *toConnect);
      }
    }
  }

  {
    size_t nextOutputIndex = 0;
    auto newOutputs = newModPorts.outputs();
    for (const auto& oport : oldModPorts.outputs())
    {
      if (oport->isConnected())
      {
        auto toConnect = std::find_if(newOutputs.begin(), newOutputs.end(),
          [&](const PortWidget* port) { return port->description()->get_typename() == oport->description()->get_typename() && port->description()->getIndex() >= nextOutputIndex; });
        if (toConnect == newOutputs.end())
        {
          guiLogCritical("Logical error: could not find output port to connect to {}", oport->name().toStdString());
          break;
        }
        auto connectedPorts = oport->connectedPorts();
        std::vector<PortWidget*> dynamicPortsNeedSpecialHandling;
        std::copy_if(connectedPorts.begin(), connectedPorts.end(), std::back_inserter(dynamicPortsNeedSpecialHandling), [](const PortWidget* p) { return p->description()->isDynamic(); });
        connectedPorts.erase(std::remove_if(connectedPorts.begin(), connectedPorts.end(), [](const PortWidget* p) { return p->description()->isDynamic(); }), connectedPorts.end());
        oport->deleteConnections();
        for (const auto& connected : connectedPorts)
        {
          requestConnectionWidget(connected, *toConnect);
        }
        nextOutputIndex = (*toConnect)->description()->getIndex() + 1;
      }
    }
  }

  oldModule->deleteLater();
}

ModuleProxyWidget* NetworkEditor::setupModuleWidget(ModuleWidget* module)
{
  logViewerDims("Scene bounds pre-add:");

  auto proxy = new ModuleProxyWidget(module);

  connect(module, &ModuleWidget::removeModule, controller_.get(), &NetworkEditorControllerGuiProxy::removeModule);
  connect(module, &ModuleWidget::removeModule, this, &NetworkEditor::modified);
  connect(module, &ModuleWidget::noteChanged, this, &NetworkEditor::modified);
  connect(module, &ModuleWidget::executionDisabled, this, &NetworkEditor::modified);
  connect(module, &ModuleWidget::requestConnection, this, &NetworkEditor::requestConnection);
  connect(module, &ModuleWidget::duplicateModule, this, &NetworkEditor::duplicateModule);
  connect(this, &NetworkEditor::networkEditorMouseButtonPressed, module, &ModuleWidget::cancelConnectionsInProgress);
  connect(controller_.get(), &NetworkEditorControllerGuiProxy::connectionAdded, module, &ModuleWidget::connectionAdded);
  connect(module, &ModuleWidget::executedManually, this, &NetworkEditor::executeModule);
  connect(module, &ModuleWidget::connectionDeleted, this, &NetworkEditor::connectionDeleted);
  connect(module, &ModuleWidget::connectionDeleted, this, &NetworkEditor::modified);
  connect(module, &ModuleWidget::connectNewModule, this, &NetworkEditor::connectNewModule);
  connect(module, &ModuleWidget::insertNewModule, this, &NetworkEditor::insertNewModule);
  connect(module, &ModuleWidget::replaceModuleWith, this, &NetworkEditor::replaceModuleWith);
  connect(module, &ModuleWidget::disableWidgetDisabling, this, &NetworkEditor::disableWidgetDisabling);
  connect(module, &ModuleWidget::reenableWidgetDisabling, this, &NetworkEditor::reenableWidgetDisabling);
  connect(this, &NetworkEditor::connectionStatusChanged, module, &ModuleWidget::connectionStatusChanged);

  if (module->hasDynamicPorts())
  {
    connect(controller_.get(), &NetworkEditorControllerGuiProxy::portAdded, module, &ModuleWidget::addDynamicPort);
    connect(controller_.get(), &NetworkEditorControllerGuiProxy::portRemoved, module, &ModuleWidget::removeDynamicPort);
    connect(module, &ModuleWidget::dynamicPortChanged, proxy, &ModuleProxyWidget::createPortPositionProviders);
  }

  LOG_TRACE("NetworkEditor connecting to state.");
  module->getModule()->get_state()->connectStateChanged([this]() { modified(); });

  connect(this, &NetworkEditor::networkExecuted, module, &ModuleWidget::resetLogButtonColor);
  connect(this, &NetworkEditor::networkExecuted, module, &ModuleWidget::resetProgressBar);

#ifdef MODULE_POSITION_LOGGING
  qDebug() << "__NW__" << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  proxy->setZValue(zLevelManager_->get_max());

#ifdef MODULE_POSITION_LOGGING
  qDebug() << "__NW__" << __LINE__ << modulePlacement_.getLast();
#endif

  proxy->setPos(modulePlacement_.getLast());

#ifdef MODULE_POSITION_LOGGING
  qDebug() << "__NW__" << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  proxy->setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemSendsGeometryChanges);

#ifdef MODULE_POSITION_LOGGING
  qDebug() << "__NW__" << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  connect(scene_, &QGraphicsScene::selectionChanged, proxy, &ModuleProxyWidget::highlightIfSelected);
  connect(proxy, &ModuleProxyWidget::selected, this, &NetworkEditor::bringToFront);
  connect(proxy, &ModuleProxyWidget::widgetMoved, this, &NetworkEditor::modified);
  connect(proxy, &ModuleProxyWidget::widgetMoved, this, &NetworkEditor::moduleMoved);
  connect(this, &NetworkEditor::snapToModules, proxy, &ModuleProxyWidget::snapToGrid);
  connect(this, &NetworkEditor::highlightPorts, proxy, &ModuleProxyWidget::highlightPorts);
  connect(this, &NetworkEditor::resetModulesDueToCycle, module, &ModuleWidget::changeExecuteButtonToPlay);
  connect(this, &NetworkEditor::defaultNotePositionChanged, proxy, &ModuleProxyWidget::setDefaultNotePosition);
  connect(this, &NetworkEditor::defaultNoteSizeChanged, proxy, &ModuleProxyWidget::setDefaultNoteSize);
  connect(module, &ModuleWidget::displayChanged, this, &NetworkEditor::updateViewport);
  connect(module, &ModuleWidget::displayChanged, proxy, &ModuleProxyWidget::createPortPositionProviders);
  connect(proxy, &ModuleProxyWidget::tagChanged, this, &NetworkEditor::highlightTaggedItem);

#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  proxy->setDefaultNotePosition(defaultNotePositionGetter_->position());

#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  proxy->setDefaultNoteSize(defaultNotePositionGetter_->size());

#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  proxy->createPortPositionProviders();

#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  proxy->highlightPorts(Preferences::Instance().highlightPorts ? 1 : 0);

#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  if (highResolutionExpandFactor_ > 1)
  {
#ifdef MODULE_POSITION_LOGGING
    qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

    auto multiplier = std::min(highResolutionExpandFactor_, 1.2);
    module->setFixedHeight(proxy->size().height() * multiplier);
    proxy->setMaximumHeight(proxy->size().height() * multiplier);
    module->setFixedWidth(proxy->size().width() * std::max(multiplier*0.9, 1.0));
    proxy->setMaximumWidth(proxy->size().width() * std::max(multiplier*0.9, 1.0));
  }

#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

#ifdef __APPLE__
  proxy->setVisible(false);
#endif
  scene_->addItem(proxy);
  ensureVisible(proxy);

#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  proxy->createStartupNote();

  scene_->clearSelection();
  proxy->setSelected(true);
  bringToFront();
  proxy->setVisible(visibleItems_);

#ifdef MODULE_POSITION_LOGGING
  qDebug() << __LINE__ << "mpw pos" << proxy->pos() << proxy->scenePos();
#endif

  guiLogDebug("Module added: {}", module->getModuleId());

  logViewerDims("Scene bounds post-add:");

  return proxy;
}

void NetworkEditor::logViewerDims(const QString& msg)
{
  if (!viewUpdateFunc_)
    return;

  const auto rect = sceneRect();
  const auto itemBound = scene_->itemsBoundingRect();
  viewUpdateFunc_(msg);
  viewUpdateFunc_(tr(" sceneRect topLeft %1,%2 bottomRight %3,%4")
    .arg(rect.topLeft().x())
    .arg(rect.topLeft().y())
    .arg(rect.bottomRight().x())
    .arg(rect.bottomRight().y())
  );
  viewUpdateFunc_(tr("itemsBoundingRect topLeft %1,%2 bottomRight %3,%4")
    .arg(itemBound.topLeft().x())
    .arg(itemBound.topLeft().y())
    .arg(itemBound.bottomRight().x())
    .arg(itemBound.bottomRight().y())
  );
  const auto visibleRect = mapToScene(viewport()->geometry()).boundingRect();
  viewUpdateFunc_(tr("visibleRect topLeft %1,%2 bottomRight %3,%4")
    .arg(visibleRect.topLeft().x())
    .arg(visibleRect.topLeft().y())
    .arg(visibleRect.bottomRight().x())
    .arg(visibleRect.bottomRight().y())
  );
}

void NetworkEditor::setMouseAsDragMode()
{
  setDragMode(ScrollHandDrag);
  //tailRecurse(&NetworkEditor::setMouseAsDragMode);
}

void NetworkEditor::setMouseAsSelectMode()
{
  setDragMode(RubberBandDrag);
  //tailRecurse(&NetworkEditor::setMouseAsSelectMode);
}

void NetworkEditor::bringToFront()
{
  zLevelManager_->bringToFront();
}

void ZLevelManager::bringToFront()
{
  ++maxZ_;
  setZValue(maxZ_);
}

void NetworkEditor::sendToBack()
{
  zLevelManager_->sendToBack();
}

void ZLevelManager::sendToBack()
{
  --minZ_;
  setZValue(minZ_);
}

void ZLevelManager::setZValue(int z)
{
  auto node = selectedModuleProxy();
  if (node)
  {
    node->setZValue(z);
  }
}

ModuleProxyWidget* getModuleProxy(QGraphicsItem* item)
{
  return dynamic_cast<ModuleProxyWidget*>(item);
}

ModuleWidget* Gui::getModule(QGraphicsItem* item)
{
  auto proxy = getModuleProxy(item);
  if (proxy)
    return static_cast<ModuleWidget*>(proxy->widget());
  return nullptr;
}

void NetworkEditor::setVisibility(bool visible)
{
  visibleItems_ = visible;
  ConnectionFactory::setVisibility(visible);
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto p = getModuleProxy(item))
      p->setVisible(visibleItems_);
    else if (auto c = dynamic_cast<ConnectionLine*>(item))
      c->setVisible(visibleItems_);
  }
}

void NetworkEditor::hidePipesByType(const std::string& type)
{
  if (connectionDimmingTimeLine_)
  {
    qDebug() << "Already dimming certain connections, please wait.";
    return;
  }

  std::vector<ConnectionLine*> conns;
  Q_FOREACH(auto item, scene_->items())
  {
    if (auto c = dynamic_cast<ConnectionLine*>(item))
    {
      if (type == c->connectedPorts().first->description()->get_typename())
      {
        guiLogDebug("dimming {}", c->id().id_);
        conns.push_back(c);
      }
    }
  }
  if (!conns.empty())
  {
    connectionDimmingTimeLine_ = new QTimeLine(ConnectionHideTimeMS_, this);
    connect(connectionDimmingTimeLine_, &QTimeLine::valueChanged,
      [conns](qreal q)
      {
        std::for_each(conns.begin(), conns.end(), [q](ConnectionLine* c) { c->setOpacity(q); });
      });
    connect(connectionDimmingTimeLine_, &QTimeLine::finished, [this]()
      {
        connectionDimmingTimeLine_->deleteLater();
        connectionDimmingTimeLine_ = nullptr;
      });
    connectionDimmingTimeLine_->start();
  }

}

//TODO copy/paste
ModuleWidget* NetworkEditor::selectedModule() const
{
  auto items = scene_->selectedItems();
  if (items.count() == 1)
  {
    return getModule(items.first());
  }
  return nullptr;
}

ModuleProxyWidget* ZLevelManager::selectedModuleProxy() const
{
  auto items = scene_->selectedItems();
  if (items.count() == 1)
  {
    return getModuleProxy(items.first());
  }
  return nullptr;
}

ConnectionLine* NetworkEditor::selectedLink() const
{
  auto items = scene_->selectedItems();
  if (items.count() == 1)
    return dynamic_cast<ConnectionLine*>(items.first());
  return nullptr;
}

NetworkEditor::ModulePair NetworkEditor::selectedModulePair() const
{
  auto items = scene_->selectedItems();
  if (items.count() == 2)
  {
    auto first = getModule(items.first());
    auto second = getModule(items.last());
    if (first && second)
		return ModulePair(first, second);
  }
  return ModulePair();
}

void NetworkEditor::del()
{
  if (!isActiveWindow())
    return;

  deleteImpl(scene_->selectedItems());
  updateViewport();
  Q_EMIT modified();
}

void NetworkEditor::deleteImpl(QList<QGraphicsItem*> items)
{
  QMutableListIterator<QGraphicsItem*> i(items);
  while (i.hasNext())
  {
    auto link = dynamic_cast<QGraphicsPathItem*>(i.next());
    if (link)
    {
      scene()->removeItem(link);
      delete link;
      i.remove();
    }
  }
  qDeleteAll(items);
}

void NetworkEditor::cut()
{
  if (!isActiveWindow())
    return;

  copy();
  del();
}

void NetworkEditor::copy()
{
  if (!isActiveWindow())
    return;

  auto selected = scene_->selectedItems();
  if (selected.empty())
    return;

  auto modSelected = [&selected](ModuleHandle mod)
  {
    for (const auto& item : selected)
    {
      if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
      {
        if (w->getModuleWidget()->getModuleId() == mod->id().id_)
          return true;
      }
    }
    return false;
  };
  auto connSelected = [&selected](const ConnectionDescription& conn)
  {
    for (const auto& item : selected)
    {
      if (auto connLine = dynamic_cast<ConnectionLine*>(item))
      {
        if (connLine->id().describe() == conn)
          return true;
      }
    }
    return false;
  };

  auto file = controller_->serializeNetworkFragment(modSelected, connSelected);

  if (file)
  {
    //TODO encapsulate
    std::ostringstream ostr;
    XMLSerializer::save_xml(*file, ostr, "networkFragment");
    auto xml = QString::fromStdString(ostr.str());

    QApplication::clipboard()->setText(xml);
    Q_EMIT newSubnetworkCopied(xml);
  }
  else
  {
    qDebug() << "null net fragment";
  }
}

void NetworkEditor::paste()
{
  pasteImpl(QApplication::clipboard()->text());
}

void NetworkEditor::pasteImpl(const QString& xml)
{
  if (parentNetwork_)
    return;

  std::istringstream istr(xml.toStdString());
  try
  {
    appendToNetwork(XMLSerializer::load_xml<NetworkFile>(istr));
  }
  catch (...)
  {
    QMessageBox::critical(this, "Paste error", "Invalid clipboard contents: " + xml);
  }
}

void NetworkEditor::contextMenuEvent(QContextMenuEvent *event)
{
  //TODO: this menu needs to check for certain editing conditions. Disabling for now.
  if (false)
  {
    auto items = scene_->items(mapToScene(event->pos()));
    if (items.isEmpty())
    {
      QMenu menu(this);
      menu.addActions(actions());
      menu.exec(event->globalPos());
    }
  }
}

static auto eventPos(QDropEvent* event)
{
#ifdef SCIRUN_QT6_ENABLED
  return event->position().toPoint();
#else
  return event->pos();
#endif
}

void NetworkEditor::dropEvent(QDropEvent* event)
{
  auto data = event->mimeData();
  if (data->hasUrls())
  {
    auto urls = data->urls();
    if (!urls.isEmpty())
    {
      auto file = urls[0].toLocalFile();
      QFileInfo check_file(file);
      if (check_file.exists() && check_file.isFile() && file.endsWith("srn5"))
      {
        Q_EMIT requestLoadNetwork(file);
        return;
      }
    }
  }

  if (moduleSelectionGetter_->isModule())
  {
    addNewModuleAtPosition(mapToScene(eventPos(event)));
  }
  else if (moduleSelectionGetter_->isClipboardXML())
    pasteImpl(moduleSelectionGetter_->clipboardXML());
}

void NetworkEditor::addNewModuleAtPosition(const QPointF& position)
{
  InEditingContext iec(this);

  modulePlacement_.setLastFromAddingNew(position);

#ifdef MODULE_POSITION_LOGGING
  logCritical("{},{}", __FILE__, __LINE__);
#endif

  controller_->addModule(moduleSelectionGetter_->text().toStdString());

#ifdef MODULE_POSITION_LOGGING
  logCritical("{},{}", __FILE__, __LINE__);
#endif

  Q_EMIT modified();
}

void NetworkEditor::addModuleViaDoubleClickedTreeItem()
{
  if (parentNetwork_)
    return;

  if (moduleSelectionGetter_->isModule())
  {
    addNewModuleAtPosition(modulePlacement_.getLastForDoubleClickedItem(scene_->itemsBoundingRect().center()));
  }
  else if (moduleSelectionGetter_->isClipboardXML())
    pasteImpl(moduleSelectionGetter_->clipboardXML());
}

void NetworkEditor::dragEnterEvent(QDragEnterEvent* event)
{
  event->acceptProposedAction();
}

void NetworkEditor::dragMoveEvent(QDragMoveEvent* event)
{
  event->acceptProposedAction();
}

void NetworkEditor::updateViewport()
{
  viewport()->update();
}

void NetworkEditor::mouseMoveEvent(QMouseEvent *event)
{
	if (event->button() != Qt::LeftButton && event->button() != Qt::NoButton)
		Q_EMIT networkEditorMouseButtonPressed();

  if (auto cL = getSingleConnectionSelected())
  {
    if (event->buttons() & Qt::LeftButton)
    {
      if (!(event->modifiers() & Qt::ControlModifier))
      {
        auto selectedPair = cL->getConnectedToModuleIds();
        auto c1 = findById(scene_->items(), selectedPair.first);
        if (c1)
          c1->setSelected(true);
        auto c2 = findById(scene_->items(), selectedPair.second);
        if (c2)
          c2->setSelected(true);

        modulesSelectedByCL_ = true;
      }
    }
  }
  QGraphicsView::mouseMoveEvent(event);
}

void NetworkEditor::mousePressEvent(QMouseEvent *event)
{
  if (event->button() == Qt::MiddleButton)
    Q_EMIT middleMouseClicked();
  QGraphicsView::mousePressEvent(event);
}

void NetworkEditor::mouseReleaseEvent(QMouseEvent *event)
{
  //logViewerDims("mouseReleaseEvent 0:");

  if (modulesSelectedByCL_)
  {
    unselectConnectionGroup();
    Q_EMIT modified();
  }
  modulesSelectedByCL_ = false;

  //logViewerDims("mouseReleaseEvent 1:");

  QGraphicsView::mouseReleaseEvent(event);

  //logViewerDims("mouseReleaseEvent 2:");
}

NetworkSearchWidget::NetworkSearchWidget(NetworkEditor* ned)
{
  setupUi(this);
  connect(searchLineEdit_, &QLineEdit::textChanged, ned, &NetworkEditor::searchTextChanged);
  connect(clearToolButton_, &QPushButton::clicked, searchLineEdit_, &QLineEdit::clear);
}

SearchResultItem::SearchResultItem(const QString& text, const QColor& color, std::function<void()> action, QGraphicsItem* parent)
  : FloatingTextItem(text, action, parent)
{
  setDefaultTextColor(color);
  const auto backgroundGray = QString("background:rgba(%1, %1, %1, 30%)").arg(200);
  setHtml("<div style='" + backgroundGray + ";font: 15px Lucida, sans-serif'>" + toPlainText() + "</div>");
  items_.insert(this);
}

SearchResultItem::~SearchResultItem()
{
  items_.erase(this);
}

void SearchResultItem::removeAll()
{
  auto copyOfItems(items_);
  for (auto& sri : copyOfItems)
    delete sri;
  items_.clear();
}

std::set<SearchResultItem*> SearchResultItem::items_;

enum SearchTupleParts
{
  ItemType,
  ItemName,
  ItemAction,
  ItemColor
};

class NetworkSearchEngine
{
public:
  NetworkSearchEngine(QGraphicsScene* scene, TagColorFunc tagColor) : scene_(scene), tagColor_(tagColor) {}

  using Result = std::tuple<QString, QString, std::function<void()>, QColor>;
  using ResultList = std::vector<Result>;
  ResultList search(const QString& text) const
  {
    ResultList results;
    Q_FOREACH(auto item, scene_->items())
    {
      ResultList subresults;
      if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
      {
        subresults = searchItem(w, text);
      }
      else if (dynamic_cast<FloatingTextItem*>(item))
      {
        // skip--don't search errors or search results
      }
      else if (auto t = dynamic_cast<QGraphicsTextItem*>(item))
      {
        subresults = searchItem(t, text);
      }
      else if (auto s = dynamic_cast<QGraphicsSimpleTextItem*>(item))
      {
        subresults = searchItem(s, text);
      }
      else
      {
      }
      results.insert(results.end(), subresults.begin(), subresults.end());
    }
    std::sort(results.begin(), results.end(), [](const Result& r1, const Result& r2) { return std::get<ItemType>(r1) < std::get<ItemType>(r2); });
    return results;
  }
private:
  ResultList searchItem(ModuleProxyWidget* mod, const QString& text) const
  {
    ResultList results;
    auto id = mod->getModuleWidget()->getModuleId();
    if (boost::ifind_first(id, text.toStdString()))
    {
      auto tag = mod->data(TagDataKey).toInt();
      results.emplace_back("Module",
        QString::fromStdString(id),
        [mod]() { mod->showAndColor(Qt::green); },
        tag != NoTag ? tagColor_(tag) : Qt::white);
    }

    auto metadata = mod->getModuleWidget()->metadataToString();
    if (metadata.contains(text, Qt::CaseInsensitive))
    {
      results.emplace_back("Module metadata match in",
        QString::fromStdString(id),
        [mod]() { mod->showAndColor(Qt::yellow); },
        Qt::yellow);
    }

    auto dialog = mod->getModuleWidget()->dialog();
    if (dialog && text.length() > 5)
    {
      auto widgetMatches = dialog->findChildren<QWidget*>(QRegularExpression(".*" + text + ".*", QRegularExpression::CaseInsensitiveOption));
      Q_FOREACH(auto widget, widgetMatches)
      {
        results.emplace_back("Module UI widget match",
          QString::fromStdString(id) + "::" + widget->objectName(),
          [mod]() { mod->showAndColor("#AA3333"); },
          "#AA3333");
      }
    }
    return results;
  }

  ResultList searchItem(QGraphicsTextItem* note, const QString& text) const
  {
    ResultList results;
    auto cursor = note->document()->find(text);
    if (!cursor.isNull())
    {
      results.emplace_back("Note",
        "..." + note->toPlainText().mid(cursor.position() - 10, 20) + "...",
        [note, text]() { ModuleProxyWidget::ensureItemVisible(note); selectNote(note, text); },
        Qt::white);
    }
    return results;
  }

  ResultList searchItem(QGraphicsSimpleTextItem* tag, const QString& text) const
  {
    ResultList results;
    if (tag->text().contains(text, Qt::CaseInsensitive))
    {
      results.emplace_back("Tag",
        tag->text(),
        [tag]() { ModuleProxyWidget::ensureItemVisible(tag); },
        tagColor_(tag->data(TagTextKey).toInt()));
    }
    return results;
  }

  static void selectNote(QGraphicsTextItem* note, const QString& text)
  {
    auto doc = note->document();
    QTextCursor cur(doc->find(text));
    note->setTextCursor(cur);
  }

  QGraphicsScene* scene_;
  TagColorFunc tagColor_;
};

void scaleTextItem(QGraphicsTextItem* item, double textScale)
{
  item->setTransform(QTransform::fromScale(textScale, textScale));
}

void NetworkEditor::searchTextChanged(const QString& text)
{
  if (text.isEmpty())
  {
    SearchResultItem::removeAll();
    return;
  }
  if (text.length() > 2)
  {
    SearchResultItem::removeAll();

    NetworkSearchEngine engine(scene(), tagColor_);
    auto results = engine.search(text);
    auto textScale = 1.0 / currentScale_;

    if (!results.empty())
    {
      auto title = new SearchResultItem("Search results:", Qt::green, {});
      title->setPos(positionOfFloatingText(title->num(), true, 20, textScale * 22));
      scene()->addItem(title);
      scaleTextItem(title, textScale);
    }
    for (const auto& result : results)
    {
      auto searchItem = new SearchResultItem(std::get<SearchTupleParts::ItemType>(result) + ": " + std::get<SearchTupleParts::ItemName>(result),
        std::get<ItemColor>(result), std::get<ItemAction>(result));
      searchItem->setPos(positionOfFloatingText(searchItem->num(), true, 50, textScale * 22));
      scene()->addItem(searchItem);
      scaleTextItem(searchItem, textScale);
    }
  }
}

QPointF NetworkEditor::positionOfFloatingText(int num, bool top, int horizontalIndent, int verticalSpacing) const
{
  auto visibleRect = mapToScene(viewport()->geometry()).boundingRect();
  auto corner = top ? visibleRect.topLeft() : visibleRect.bottomLeft();
  return (corner + QPointF(horizontalIndent, (top ? 1 : -1)* (verticalSpacing * num + 100)));
}

void NetworkEditor::displayError(const QString& msg, std::function<void()> showModule)
{
  if (Preferences::Instance().showModuleErrorInlineMessages)
  {
    auto errorItem = new ErrorItem(msg, showModule);
    scene()->addItem(errorItem);

    errorItem->setPos(positionOfFloatingText(errorItem->num(), false, 30, 40));
  }
}

ConnectionLine* NetworkEditor::getSingleConnectionSelected()
{
	ConnectionLine* connectionSelected = nullptr;
	auto item = scene_->selectedItems();
	if(item.count() == 1 && (connectionSelected = qgraphicsitem_cast<ConnectionLine*>(item.first())))
		return connectionSelected;
	return connectionSelected;
}

void NetworkEditor::unselectConnectionGroup()
{
  auto items = scene_->selectedItems();
	if (items.count() == 3)
	{
		int hasConnection = 0;
		int	hasWidgets = 0;

		Q_FOREACH(QGraphicsItem* item, items)
		{
			if (auto cL = qgraphicsitem_cast<ConnectionLine*>(item))
			{
				++hasConnection;
				items.push_front(cL);
			}
			if (qgraphicsitem_cast<ModuleProxyWidget*>(item))
				++hasWidgets;
		}
		if(hasConnection == 1 && hasWidgets == 2)
		{
			if (auto cL = qgraphicsitem_cast<ConnectionLine*>(items.first()))
			{
				auto selectedPair = cL->getConnectedToModuleIds();

				cL->setSelected(false);
				findById(scene_->items(),selectedPair.first)->setSelected(false);
				findById(scene_->items(),selectedPair.second)->setSelected(false);
			}
		}
	}
}

ModulePositionsHandle NetworkEditor::dumpModulePositions(ModuleFilter filter) const
{
  auto positions(makeShared<ModulePositions>());
  fillModulePositionMap(*positions, filter);
#if 0
  for (const auto& sub : childrenNetworks_)
  {
    sub.second->get()->fillModulePositionMap(*positions, filter);
  }
#endif
  return positions;
}

void NetworkEditor::fillModulePositionMap(ModulePositions& positions, ModuleFilter filter) const
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      if (filter(w->getModuleWidget()->getModule()))
        positions.modulePositions[w->getModuleWidget()->getModuleId()] = std::make_pair(item->scenePos().x(), item->scenePos().y());
    }
  }
}

void NetworkEditor::centerView()
{
  if (!isActiveWindow())
  {
    //tailRecurse(&NetworkEditor::centerView);
    return;
  }

  ModulePositions positions;
  fillModulePositionMap(positions, [](ModuleHandle) { return true; });
  centerOn(findCenterOfNetwork(positions));
}

ModuleNotesHandle NetworkEditor::dumpModuleNotes(ModuleFilter filter) const
{
  auto notes(makeShared<ModuleNotes>());
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      auto note = w->currentNote();
      if (filter(w->getModuleWidget()->getModule()) &&
        !note.plainText_.isEmpty())
        notes->notes[w->getModuleWidget()->getModuleId()] = NoteXML(note.html_.toStdString(), static_cast<int>(note.position_), note.plainText_.toStdString(), note.fontSize_);
    }
  }
  return notes;
}

void NetworkEditor::copyNote(ModuleHandle from, ModuleHandle to) const
{
  Note note;
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      if (w->getModuleWidget()->getModule()->id() == from->id())
      {
        note = w->currentNote();
        break;
      }
    }
  }

  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      if (w->getModuleWidget()->getModule()->id() == to->id())
      {
        w->updateNote(note);
        w->getModuleWidget()->setCurrentNote(note, true);
        break;
      }
    }
  }
}

namespace
{
  std::string connectionNoteId(const ModuleIdPair& ms)
  {
    return ms.first.id_ + "--" + ms.second.id_;
  }
}

ConnectionNotesHandle NetworkEditor::dumpConnectionNotes(ConnectionFilter filter) const
{
  auto notes(makeShared<ConnectionNotes>());
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto conn = dynamic_cast<ConnectionLine*>(item))
    {
      auto note = conn->currentNote();
      if (filter(conn->id().describe()) &&
        !note.plainText_.isEmpty())
      {
        //TODO hacky
        auto id = connectionNoteId(conn->getConnectedToModuleIds());
        notes->notes[id] = NoteXML(note.html_.toStdString(), static_cast<int>(note.position_), note.plainText_.toStdString(), note.fontSize_);
      }
    }
  }
  return notes;
}

ModuleTagsHandle NetworkEditor::dumpModuleTags(ModuleFilter filter) const
{
  auto tags(makeShared<ModuleTags>());
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto mod = dynamic_cast<ModuleProxyWidget*>(item))
    {
      if (filter(mod->getModuleWidget()->getModule()))
        tags->tags[mod->getModuleWidget()->getModuleId()] = mod->data(TagDataKey).toInt();
    }
  }
  tags->labels = tagLabelOverrides_;
  tags->showTagGroupsOnLoad = showTagGroupsOnFileLoad();
  return tags;
}

DisabledComponentsHandle NetworkEditor::dumpDisabledComponents(ModuleFilter modFilter, ConnectionFilter connFilter) const
{
  auto disabled(makeShared<DisabledComponents>());
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto mod = dynamic_cast<ModuleProxyWidget*>(item))
    {
      if (mod->getModuleWidget()->isExecutionDisabled() && modFilter(mod->getModuleWidget()->getModule()))
        disabled->disabledModules.push_back(mod->getModuleWidget()->getModuleId());
    }
    if (auto conn = dynamic_cast<ConnectionLine*>(item))
    {
      if (conn->disabled() && connFilter(conn->id().describe()))
      {
        disabled->disabledConnections.emplace_back(connectionNoteId(conn->getConnectedToModuleIds()));
      }
    }
  }
  return disabled;
}

#if 0
SubnetworksHandle NetworkEditor::dumpSubnetworks(ModuleFilter modFilter) const
{
  auto subnets(makeShared<Subnetworks>());

  for (const auto& child : childrenNetworks_)
  {
    child.second->get()->dumpSubnetworksImpl(child.first, *subnets, modFilter);
  }
  return subnets;
}
#endif

QPointF NetworkEditor::getModulePositionAdjustment(const ModulePositions& modulePositions)
{
  std::vector<QPointF> positions;
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      auto posIter = modulePositions.modulePositions.find(w->getModuleWidget()->getModuleId());
      if (posIter != modulePositions.modulePositions.end())
      {
        positions.emplace_back(posIter->second.first, posIter->second.second);
      }
    }
  }

  QPointF adjustment;
  if (!positions.empty())
  {
    auto minX = *std::min_element(positions.begin(), positions.end(), [](const QPointF& p1, const QPointF& p2) { return p1.x() < p2.x(); });
    auto minY = *std::min_element(positions.begin(), positions.end(), [](const QPointF& p1, const QPointF& p2) { return p1.y() < p2.y(); });
    if (minX.x() < 0 || minY.y() < 0)
    {
      adjustment = { minX.x(), minY.y() };
    }
  }
  return adjustment;
}

void NetworkEditor::updateModulePositions(const ModulePositions& modulePositions, bool selectAll)
{
#ifdef MODULE_POSITION_LOGGING
  logCritical("updateModulePositions {},{}", __FILE__, __LINE__);
#endif

  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      auto posIter = modulePositions.modulePositions.find(w->getModuleWidget()->getModuleId());
      if (posIter != modulePositions.modulePositions.end())
      {
        QPointF p {posIter->second.first, posIter->second.second};
        w->setPos(p);
        ensureVisible(w);
        if (selectAll)
          w->setSelected(true);
      }
    }
  }
#if 0
  for (const auto& child : childrenNetworks_)
  {
    child.second->get()->updateModulePositions(modulePositions, selectAll);
  }
#endif
}

void NetworkEditor::updateModuleNotes(const ModuleNotes& moduleNotes)
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      auto noteIter = moduleNotes.notes.find(w->getModuleWidget()->getModuleId());
      if (noteIter != moduleNotes.notes.end())
      {
        auto noteXML = noteIter->second;
        Note note(QString::fromStdString(noteXML.noteHTML), QString::fromStdString(noteXML.noteText), noteXML.fontSize, NotePosition(noteXML.position));
        w->getModuleWidget()->updateNoteFromFile(note);
      }
    }
  }
}

void NetworkEditor::updateModuleTags(const ModuleTags& moduleTags)
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      auto tagIter = moduleTags.tags.find(w->getModuleWidget()->getModuleId());
      if (tagIter != moduleTags.tags.end())
      {
        w->setData(TagDataKey, tagIter->second);
      }
    }
  }
  setShowTagGroupsOnFileLoad(moduleTags.showTagGroupsOnLoad);
  tagLabelOverrides_ = moduleTags.labels;
  if (showTagGroupsOnFileLoad())
  {
    tagGroupsActive_ = true;
    drawTagGroups();
  }
}

void NetworkEditor::updateConnectionNotes(const ConnectionNotes& notes)
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto conn = dynamic_cast<ConnectionLine*>(item))
    {
      auto id = connectionNoteId(conn->getConnectedToModuleIds());
      auto noteIter = notes.notes.find(id);
      if (noteIter != notes.notes.end())
      {
        auto noteXML = noteIter->second;
        Note note(QString::fromStdString(noteXML.noteHTML), QString::fromStdString(noteXML.noteText), noteXML.fontSize, NotePosition(noteXML.position));
        conn->updateNoteFromFile(note);
      }
    }
  }
}

void NetworkEditor::updateDisabledComponents(const DisabledComponents& disabled)
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto conn = dynamic_cast<ConnectionLine*>(item))
    {
      auto id = connectionNoteId(conn->getConnectedToModuleIds());
      if (std::find(disabled.disabledConnections.begin(), disabled.disabledConnections.end(), id) != disabled.disabledConnections.end())
      {
        conn->setDisabled(true);
      }
    }

    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      if (std::find(disabled.disabledModules.begin(), disabled.disabledModules.end(), w->getModuleWidget()->getModuleId()) != disabled.disabledModules.end())
      {
        w->getModuleWidget()->setExecutionDisabled(true);
      }
    }
  }
}

void NetworkEditor::executeAll()
{
  preexecute_();
  // explicit type needed for older Qt and/or clang
  std::function<void()> exec = [this]() { controller_->executeAll(); };
  (void)QtConcurrent::run(exec);

  //TODO: not sure about this right now.
  //Q_EMIT modified();
  Q_EMIT networkExecuted();
}

void NetworkEditor::executeModule(const ModuleHandle& module, bool fromButton)
{
  preexecute_();
  // explicit type needed for older Qt and/or clang
  std::function<void()> exec = [this, &module, fromButton]() { controller_->executeModule(module, fromButton); };
  (void)QtConcurrent::run(exec);
  //TODO: not sure about this right now.
  //Q_EMIT modified();
  Q_EMIT networkExecuted();
}

ExecutableObject* NetworkEditor::lookupExecutable(const ModuleId& id) const
{
#if 0
  for (const auto& child : childrenNetworks_)
  {
    auto exec = child.second->get()->lookupExecutable(id);
    if (exec)
      return exec;
  }
#endif

  auto widget = findById(scene_->items(), id.id_);
  return widget ? widget->getModuleWidget() : nullptr;
}

void NetworkEditor::resetNetworkDueToCycle()
{
  Q_EMIT resetModulesDueToCycle();
  //TODO: ??reset module colors--right now they stay yellow
}

void NetworkEditor::removeModuleWidget(const ModuleId& id)
{
  auto widget = findById(scene_->items(), id.id_);
  if (widget)
  {
    widget->getModuleWidget()->setDeletedFromGui(id.id_.find("Subnet") != std::string::npos);
    delete widget;
    Q_EMIT modified();
  }
}

void NetworkEditor::clear()
{
  tagLabelOverrides_.clear();
  ModuleWidget::NetworkClearingScope clearing;

  QList<QGraphicsItem*> deleteTheseFirst;
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
#if 0
    if (dynamic_cast<SubnetWidget*>(getModule(item)))
    {
      deleteTheseFirst.append(item);
    }
    else
#endif
    if (auto vsw = dynamic_cast<ModuleWidget*>(getModule(item)))
    {
      auto vs = vsw->dialog();
      if (vs)
      {
        vs->blockSignals(true);
      }
    }
  }
  deleteImpl(deleteTheseFirst);
  scene_->clear();
  //TODO: this (unwritten) method does not need to be called here.  the dtors of all the module widgets get called when the scene_ is cleared, which triggered removal from the underlying network.
  // we'll need a similar hook when programming the scripting interface (moduleWidgets<->modules).
  //controller_->clear();
  if (!parentNetwork_)
    Q_EMIT modified();
}

NetworkFileHandle NetworkEditor::saveNetwork() const
{
  return controller_->saveNetwork();
}

void NetworkEditor::loadNetwork(const NetworkFileHandle& xml)
{
  fileLoading_ = true;
  controller_->loadNetwork(xml);
  fileLoading_ = false;

  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      w->getModuleWidget()->postLoadAction();
    }
  }

  setSceneRect(QRectF());

#ifdef __APPLE__
  QTimer::singleShot(macModulePositionWorkaroundTimerValue, [this]()
  {
    Q_FOREACH(QGraphicsItem* item, scene_->items())
    {
      item->setSelected(true);
      item->setSelected(false);
    }
  });
#endif
}

void NetworkEditor::deselectAll()
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
    {
      w->setSelected(false);
    }
  }
}

void NetworkEditor::appendToNetwork(const NetworkFileHandle& xml)
{
  auto originalItems = scene_->items();
  fileLoading_ = true;
  controller_->appendToNetwork(xml);
  fileLoading_ = false;

  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (!originalItems.contains(item))
    {
      if (auto w = dynamic_cast<ModuleProxyWidget*>(item))
      {
        w->getModuleWidget()->postLoadAction();
        w->setSelected(true);
      }
    }
  }

  setSceneRect(QRectF());
}

void NetworkEditor::disableViewScenes()
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto c = dynamic_cast<ConnectionLine*>(item))
    {
      if (c->id().id_.find("ViewScene") != std::string::npos)
        c->setDisabled(true);
    }
  }
  //TODO: doesn't work yet.
  //Application::Instance().controller()->connectNetworkExecutionFinished([this](int code){ enableViewScenes(); });
}

void NetworkEditor::enableViewScenes()
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto c = dynamic_cast<ConnectionLine*>(item))
    {
      if (c->id().id_.find("ViewScene") != std::string::npos)
        c->setDisabled(false);
    }
  }
}

size_t NetworkEditor::numModules() const
{
  if (!controller_)
    return 0;
  return controller_->numModules();
}

void NetworkEditor::setConnectionPipelineType(int type)
{
  ConnectionFactory::setType(ConnectionDrawType(type));
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto c = dynamic_cast<ConnectionLine*>(item))
    {
      c->setDrawStrategy(ConnectionFactory::getCurrentDrawer());
    }
  }
}

int NetworkEditor::connectionPipelineType() const
{
  return static_cast<int>(ConnectionFactory::getType());
}

int NetworkEditor::errorCode() const
{
  return controller_->errorCode();
}

ModuleEventProxy::ModuleEventProxy()
{
  qRegisterMetaType<std::string>("std::string");
  qRegisterMetaType<ModuleHandle>("SCIRun::Dataflow::Networks::ModuleHandle");
  qRegisterMetaType<ConnectionDescription>("SCIRun::Dataflow::Networks::ConnectionDescription");
  qRegisterMetaType<ModuleId>("SCIRun::Dataflow::Networks::ModuleId");
  qRegisterMetaType<ConnectionId>("SCIRun::Dataflow::Networks::ConnectionId");
  qRegisterMetaType<ModuleCounter>("SCIRun::Dataflow::Engine::ModuleCounter");
}

void ModuleEventProxy::trackModule(ModuleHandle module)
{
  module->connectExecuteBegins([this](const std::string& id) { moduleExecuteStart(id); });
  module->connectExecuteEnds([this](double t, const std::string& id) { moduleExecuteEnd(t, id); });
}

void NetworkEditor::disableInputWidgets()
{
  //deleteAction_->setDisabled(true);
}

void NetworkEditor::enableInputWidgets()
{
  //deleteAction_->setEnabled(true);
}

void NetworkEditor::setBackground(const QBrush& brush)
{
  scene_->setBackgroundBrush(brush);
  Preferences::Instance().networkBackgroundColor.setValue(brush.color().name().toStdString());
}

QBrush NetworkEditor::background() const
{
  return scene_->backgroundBrush();
}

void NetworkEditor::selectAll()
{
  if (!isActiveWindow())
  {
    //tailRecurse(&NetworkEditor::selectAll);
    return;
  }

  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    item->setSelected(true);
  }
}

void NetworkEditor::pinAllModuleUIs()
{
  if (!isActiveWindow())
  {
    //tailRecurse(&NetworkEditor::pinAllModuleUIs);
    return;
  }

  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    auto module = getModule(item);
    if (module)
      module->pinUI();
  }
}

void NetworkEditor::hideAllModuleUIs()
{
  if (!isActiveWindow())
  {
    //tailRecurse(&NetworkEditor::hideAllModuleUIs);
    return;
  }

  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    auto module = getModule(item);
    if (module)
      module->hideUI();
  }
}

void NetworkEditor::seeThroughAllModuleUIs()
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    auto module = getModule(item);
    if (module)
      module->seeThroughUI();
  }
}

void NetworkEditor::normalOpacityAllModuleUIs()
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    auto module = getModule(item);
    if (module)
      module->normalOpacityUI();
  }
}

void NetworkEditor::restoreAllModuleUIs()
{
  if (!isActiveWindow())
  {
    //tailRecurse(&NetworkEditor::restoreAllModuleUIs);
    return;
  }

  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    auto module = getModule(item);
    if (module)
    {
      module->showUI();
      module->collapsePinnedDialog();
    }
  }
}

namespace
{
  const double minScale = 0.3;
  const double maxScale = 2.0;
  const double scaleFactor = 1.15;
}

//TODO!!!
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

void NetworkEditor::wheelEvent(QWheelEvent* event)
{
  logViewerDims("pre-zoom: ");
  if (event->modifiers() & Qt::ShiftModifier)
  {
    setTransformationAnchor(AnchorUnderMouse);

    if (event->angleDelta().y() > 0)
    {
      zoomIn();
    }
    else
    {
      zoomOut();
    }
    // Don't call superclass handler here
    // as wheel is normally used for moving scrollbars
  }
  else
    QGraphicsView::wheelEvent(event);

  logViewerDims("post-zoom: ");
}

#if 0
void NetworkEditor::resizeSubnetPortHolders(double scaleFactor)
{
  for (auto& item : subnetPortHolders_)
  {
    item->resize(QSize(item->size().width() * scaleFactor, item->size().height()));

    auto isInput = item->data(123).toString() == "Inputs";
    if (isInput)
      item->setPos(topSubnetPortHolderPositioner_(visibleRect()));
    else
      item->setPos(bottomSubnetPortHolderPositioner_(visibleRect()));
    item->updateConnections();
  }
}
#endif

void NetworkEditor::zoomIn()
{
  if (!isActiveWindow())
  {
    //tailRecurse(&NetworkEditor::zoomIn);
    return;
  }

  if (currentScale_ < maxScale)
  {
    double factor = std::min(scaleFactor, 4.0/currentScale_);
    scale(factor, factor);
    currentScale_ *= factor;

#if 0
    resizeSubnetPortHolders(1.0 / factor);
#endif

    Q_EMIT zoomLevelChanged(currentZoomPercentage());
  }
}

void NetworkEditor::zoomOut()
{
  if (!isActiveWindow())
  {
    //tailRecurse(&NetworkEditor::zoomOut);
    return;
  }

  if (currentScale_ > minScale)
  {
    scale(1.0 / scaleFactor, 1.0 / scaleFactor);
    currentScale_ /= scaleFactor;

#if 0
    resizeSubnetPortHolders(scaleFactor);
#endif

    Q_EMIT zoomLevelChanged(currentZoomPercentage());
  }
}

void NetworkEditor::zoomReset()
{
  if (!isActiveWindow())
  {
    //tailRecurse(&NetworkEditor::zoomReset);
    return;
  }

  scale(1.0 / currentScale_, 1.0 / currentScale_);
  currentScale_ = 1;
  Q_EMIT zoomLevelChanged(currentZoomPercentage());
}

void NetworkEditor::zoomBestFit()
{
  //TODO not quite working yet.
  auto oldRect = sceneRect();
  setSceneRect(QRectF());
  fitInView(sceneRect(), Qt::KeepAspectRatio);
  currentScale_ *= sceneRect().x() / oldRect.x();
  //scale(1.0 / currentScale_, 1.0 / currentScale_);
  //currentScale_ = 1;
  Q_EMIT zoomLevelChanged(currentZoomPercentage());
}

int NetworkEditor::currentZoomPercentage() const
{
  return static_cast<int>(currentScale_ * 100);
}

bool NetworkEditor::containsViewScene() const
{
  return findFirstByName(scene_->items(), "ViewScene") != nullptr;
}

void NetworkEditor::moduleWindowAction()
{
  auto action = qobject_cast<QAction*>(sender());
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    auto module = getModule(item);
    if (module && module->getModuleId() == action->text().toStdString())
    {
      if (module->guiVisible())
        module->hideUI();
      else
        module->showUI();
      break;
    }
  }
}

void NetworkEditor::updateBackground(bool forceGrid)
{
  setBackgroundBrush(QPixmap(forceGrid ? standardNetworkBackgroundImage() : networkBackgroundImage()));
}

void NetworkEditor::adjustModuleWidth(int delta)
{
#ifdef MODULE_POSITION_LOGGING
  logCritical("{},{}", __FILE__, __LINE__);
#endif
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    auto proxy = getModuleProxy(item);
    if (proxy)
    {
      proxy->adjustWidth(delta);
    }
  }
}

void NetworkEditor::adjustModuleHeight(int delta)
{
#ifdef MODULE_POSITION_LOGGING
  logCritical("{},{}", __FILE__, __LINE__);
#endif
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    auto proxy = getModuleProxy(item);
    if (proxy)
    {
      proxy->adjustHeight(delta);
    }
  }
}

void NetworkEditor::metadataLayer(bool active)
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    item->setOpacity(active ? 0.4 : 1);
    auto module = getModule(item);
    if (module)
      module->updateMetadata(active);
  }
  // TODO: tailRecurse(&NetworkEditor::metadataLayer, active);
}

void NetworkEditor::adjustExecuteButtonsToDownstream(bool downOnly)
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    auto module = getModule(item);
    if (module)
    {
      module->adjustExecuteButtonToDownstream(downOnly);
    }
  }

  //TODO: tailRecurse(&NetworkEditor::adjustExecuteButtonsToDownstream, downOnly);
}

void NetworkEditor::updateExecuteButtons(bool downstream)
{
  adjustExecuteButtonsToDownstream(downstream);
}

void NetworkEditor::saveImages()
{
  Q_FOREACH(auto item, scene_->items())
  {
    auto module = getModule(item);
    if (module)
    {
      module->saveImagesFromViewScene();
    }
  }
}

QColor Gui::defaultTagColor(int tag)
{
  switch (tag)
  {
  case 0:
    return Qt::blue;
  case 1:
    return Qt::green;
  case 2:
    return Qt::darkYellow;
  case 3:
    return Qt::darkMagenta;
  case 4:
    return Qt::darkCyan;
  case 5:
    return Qt::darkRed;
  case 6:
    return Qt::darkGray;
  case 7:
    return Qt::darkGreen;
  case 8:
    return Qt::darkBlue;
  case 9:
    return Qt::black;
  default:
    return Qt::white;
  }
}

QString Gui::colorToString(const QColor& color)
{
  return QString("rgb(%1, %2, %3)").arg(color.red()).arg(color.green()).arg(color.blue());
}

QColor Gui::stringToColor(const QString& s)
{
  static QRegularExpression re("rgb\\((\\d+), (\\d+), (\\d+)\\)");
  auto match = re.match(s);
  if (match.hasMatch())
  {
    return QColor(match.captured(1).toInt(), match.captured(2).toInt(), match.captured(3).toInt());
  }
  return {};
}

QGraphicsEffect* Gui::blurEffect(double radius)
{
  auto blur = new QGraphicsBlurEffect;
  blur->setBlurRadius(radius);
  return blur;
}

void NetworkEditor::tagLayer(bool active, TagValues tag)
{
  tagLayerActive_ = active;

  if (active)
  {
    auto items = scene_->selectedItems();
    Q_FOREACH(QGraphicsItem* item, items)
    {
      if (item->data(TagDataKey).toInt() == NoTag)
      {
        if (validTag(tag))
          item->setData(TagDataKey, tag);
      }
      else if (ClearTags == tag)
      {
        item->setData(TagDataKey, NoTag);
      }
    }
  }

  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    item->setData(static_cast<int>(TagLayerKey), active);
    item->setData(static_cast<int>(CurrentTagKey), tag);
    if (active)
    {
      const auto itemTag = static_cast<TagValues>(item->data(static_cast<int>(TagDataKey)).toInt());
      if (AllTags == tag || ShowGroups == tag)
      {
        highlightTaggedItemImpl(item, itemTag);
      }
      else if (tag != NoTag && tag != ClearTags)
      {
        if (tag == itemTag)
        {
          highlightTaggedItemImpl(item, itemTag);
        }
        else
          item->setGraphicsEffect(blurEffect());
      }
    }
    else
      item->setGraphicsEffect(nullptr);
  }
  if (ShowGroups == tag)
  {
    tagGroupsActive_ = true;
    redrawTagGroups();
  }
  if (HideGroups == tag)
  {
    tagGroupsActive_ = false;
    removeTagGroups();
  }

  //TODO: tailRecurse(&NetworkEditor::tagLayer, active, tag);
}

namespace
{
  class TagGroupBox : public QGraphicsRectItem
  {
  public:
    explicit TagGroupBox(int tagNum, const QRectF& rect, NetworkEditor* ned) : QGraphicsRectItem(rect),
      tagNumber_(tagNum),
      ned_(ned)
    {
      setAcceptHoverEvents(true);
    }
  protected:
    void hoverEnterEvent(QGraphicsSceneHoverEvent*) override
    {
      setPen(QPen(pen().color(), 5));
    }

    void hoverLeaveEvent(QGraphicsSceneHoverEvent*) override
    {
      setPen(QPen(pen().color(), 3));
    }

    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event) override
    {
      QMenu menu;
      auto autoDisplay = menu.addAction("Display in saved network", ned_, &NetworkEditor::saveTagGroupRectInFile);
      autoDisplay->setCheckable(true);
      autoDisplay->setChecked(ned_->showTagGroupsOnFileLoad());
      auto rename = menu.addAction("Rename in saved network...", ned_, &NetworkEditor::renameTagGroupInFile);
      rename->setProperty("tag", tagNumber_);
      menu.exec(event->screenPos());
      QGraphicsRectItem::mouseDoubleClickEvent(event);
    }
  private:
    int tagNumber_;
    NetworkEditor* ned_;
  };
}

void NetworkEditor::saveTagGroupRectInFile()
{
  auto action = qobject_cast<QAction*>(sender());
  setShowTagGroupsOnFileLoad(action->isChecked());
  Q_EMIT modified();
}

void NetworkEditor::renameTagGroupInFile()
{
  auto action = qobject_cast<QAction*>(sender());
  auto tagNum = action->property("tag").toInt();

  bool ok;
  auto text = QInputDialog::getText(this, tr("Rename tag group"),
    tr("Enter new tag group name for this network file:"), QLineEdit::Normal, checkForOverriddenTagName(tagNum), &ok);
  if (ok && !text.isEmpty())
  {
    bool changed = tagLabelOverrides_[tagNum] != text.toStdString();
    tagLabelOverrides_[tagNum] = text.toStdString();
    if (changed)
      renameTagGroup(tagNum, text);
  }

  Q_EMIT modified();
}

void NetworkEditor::scrollContentsBy(int dx, int dy)
{
#if 0
  for (auto& item : subnetPortHolders_)
  {
    item->setPos(item->pos() + QPointF(-dx / currentScale_, -dy / currentScale_));
    item->updateConnections();
  }
#endif
  QGraphicsView::scrollContentsBy(dx, dy);
}

void NetworkEditor::drawTagGroups()
{
  if (tagGroupsActive_)
  {
    QMap<int, QRectF> tagItemRects;

    Q_FOREACH(QGraphicsItem* item, scene_->items())
    {
      if (dynamic_cast<ModuleProxyWidget*>(item))
      {
        const auto itemTag = item->data(TagDataKey).toInt();

        if (itemTag != NoTag)
        {
          auto r = item->boundingRect();
          r.translate(item->pos());
          if (!tagItemRects.contains(itemTag))
          {
            tagItemRects.insert(itemTag, r);
          }
          else
          {
            tagItemRects[itemTag] = tagItemRects[itemTag].united(r);
          }
        }
      }
    }

    for (auto rectIter = tagItemRects.constBegin(); rectIter != tagItemRects.constEnd(); ++rectIter)
    {
      auto rectBounds = rectIter.value().adjusted(-10, -10, 10, 10);
      auto tagNum = rectIter.key();
      QPen pen(tagColor_(tagNum));
      pen.setWidth(3);
      pen.setCapStyle(Qt::RoundCap);
      pen.setJoinStyle(Qt::RoundJoin);
      auto rect = new TagGroupBox(tagNum, rectBounds, this);
      rect->setPen(pen);
      scene_->addItem(rect);
      rect->setFlags(QGraphicsItem::ItemIsFocusable | QGraphicsItem::ItemIsSelectable);
      rect->setZValue(-100000);

      auto fill = new QGraphicsRectItem(rectBounds);
      auto c = pen.color();
      c.setAlphaF(0.15);
      fill->setBrush(c);
      fill->setZValue(-100000);
      scene_->addItem(fill);

      static const QFont labelFont("Courier", 20, QFont::Bold);

      auto label = scene_->addSimpleText(checkForOverriddenTagName(tagNum), labelFont);
      label->setBrush(pen.color());
      label->setData(TagTextKey, tagNum);
      static const QFontMetrics fm(labelFont);

      auto textWidthInPixels = fm.WIDTH_FUNC(label->text());
      label->setPos((rect->rect().topLeft() + rect->rect().topRight()) / 2 + QPointF(-textWidthInPixels / 2, -30));
    }
  }
}

QString NetworkEditor::checkForOverriddenTagName(int tag) const
{
  auto nameOverrideIter = tagLabelOverrides_.find(tag);
  if (nameOverrideIter != tagLabelOverrides_.end() && !nameOverrideIter->second.empty())
    return QString::fromStdString(nameOverrideIter->second);
  return tagName_(tag);
}

void NetworkEditor::removeTagGroups()
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (dynamic_cast<QGraphicsRectItem*>(item) || dynamic_cast<QGraphicsSimpleTextItem*>(item))
    {
      delete item;
    }
  }
}

void NetworkEditor::renameTagGroup(int tag, const QString& name)
{
  Q_FOREACH(QGraphicsItem* item, scene_->items())
  {
    if (auto rectLabel = dynamic_cast<QGraphicsSimpleTextItem*>(item))
    {
      if (rectLabel->data(TagTextKey).toInt() == tag)
      {
        rectLabel->setText(name);
        return;
      }
    }
  }
}

void NetworkEditor::redrawTagGroups()
{
  removeTagGroups();
  drawTagGroups();
}

void NetworkEditor::highlightTaggedItem(int tagValue)
{
  highlightTaggedItemImpl(qobject_cast<QGraphicsItem*>(sender()), static_cast<TagValues>(tagValue));
  Q_EMIT modified();
}

void NetworkEditor::highlightTaggedItemImpl(QGraphicsItem* item, TagValues tagValue)
{
  if (tagValue == NoTag)
  {
    item->setGraphicsEffect(blurEffect());
  }
  else
  {
    auto colorize = new QGraphicsColorizeEffect;
    auto color = tagColor_(tagValue);
    colorize->setColor(color);
    item->setGraphicsEffect(colorize);
  }
}

void NetworkEditor::cleanUpNetwork()
{
  controller_->cleanUpNetwork();
  centerView();
}

void NetworkEditor::showStateViewer()
{
  //TODO: make non-modal, but needs update slot
  StateViewer viewer(this);
  viewer.exec();
}

ErrorItem::ErrorItem(const QString& text, std::function<void()> action, QGraphicsItem* parent) : FloatingTextItem(text, action, parent)
{
  setDefaultTextColor(Qt::red);
}

std::atomic<int> FloatingTextItem::instanceCounter_(0);

FloatingTextItem::FloatingTextItem(const QString& text, std::function<void()> action, QGraphicsItem* parent) : QGraphicsTextItem(text, parent),
  action_(action), counter_(instanceCounter_), rect_(nullptr)
{
  setFlags(ItemIsMovable | ItemIsSelectable | ItemSendsGeometryChanges);
  setZValue(10000);
  ++instanceCounter_;

  {
    timeLine_ = new QTimeLine(10000, this);
    connect(timeLine_, &QTimeLine::valueChanged, this, &FloatingTextItem::animate);
    connect(timeLine_, &QTimeLine::finished, this, &FloatingTextItem::deleteLater);
  }
  timeLine_->start();
}

FloatingTextItem::~FloatingTextItem()
{
  --instanceCounter_;
}

void FloatingTextItem::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  if (event->buttons() & Qt::LeftButton)
  {
    if (action_)
      action_();
  }
  else if (event->buttons() & Qt::RightButton)
  {
    if (rect_)
    {
      scene()->removeItem(rect_);
      rect_ = nullptr;
    }
    scene()->removeItem(this);
  }
  QGraphicsTextItem::mousePressEvent(event);
}

void FloatingTextItem::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
  timeLine_->setCurrentTime(0);
  QGraphicsTextItem::hoverEnterEvent(event);
}

void FloatingTextItem::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event)
{
  if (!rect_)
  {
    timeLine_->stop();
    timeLine_->setCurrentTime(0);
    auto f = font();
    f.setBold(true);
    setFont(f);
    setFlags(flags() & ~ItemIsMovable);
    rect_ = scene()->addRect(boundingRect(), QPen(defaultTextColor(), 2, Qt::DotLine));
    rect_->setPos(pos());
  }
  else
  {
    auto f = font();
    f.setBold(false);
    setFont(f);
    setFlags(flags() & ItemIsMovable);
    scene()->removeItem(rect_);
    rect_ = nullptr;
    timeLine_->start();
  }

  QGraphicsTextItem::mouseDoubleClickEvent(event);
}

void FloatingTextItem::animate(qreal val)
{
  if (val < 1)
    show();
  else
    hide();
  setOpacity(val < 0.5 ? 1 : 2 - 2*val);
}

ZLevelManager::ZLevelManager(QGraphicsScene* scene)
  : scene_(scene), minZ_(INITIAL_Z), maxZ_(INITIAL_Z)
{

}

bool Gui::allowModificationSignalConnection()
{
  auto cmd = Application::Instance().parameters();
  return !cmd->executeNetwork() &&
    !cmd->executeNetworkAndQuit() &&
    !cmd->isRegressionMode();
}
