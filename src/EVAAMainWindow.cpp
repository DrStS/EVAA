/*
 * Copyright &copy; 2019, Dr. Stefan Sicklinger, Munich \n
 *
 *  All rights reserved.
 *
 *  This file is part of EVAA.
 *
 *  EVAA is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EVAA is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EVAA.  If not, see http://www.gnu.org/licenses/.
 */

#include "EVAAMainWindow.h"

#include <QApplication>
#include <QDebug>
#include <QIcon>
#include <QQmlApplicationEngine>
#include <QQmlContext>
#include <QQuickStyle>

#include "Model.h"
#include "ProcessingEngine.h"
#include "QVTKFramebufferObjectItem.h"
#include "QVTKFramebufferObjectRenderer.h"

namespace EVAA {

EVAAMainWindow::EVAAMainWindow(int argc, char **argv)
{
    QApplication app(argc, argv);
    QQmlApplicationEngine engine;

    app.setApplicationName("EVAA: Efficient Vehicle dynAmics simulAtor");
    app.setWindowIcon(QIcon(":/resources/EVAA.png"));

    // Register QML types
    qmlRegisterType<QVTKFramebufferObjectItem>("QtVTK", 1, 0, "VtkFboItem");

    // Create classes instances
    m_processingEngine = std::shared_ptr<ProcessingEngine>(new ProcessingEngine());

    // Expose C++ classes to QML
    QQmlContext *ctxt = engine.rootContext();

    ctxt->setContextProperty("mainWindowHandler", this);

    QQuickStyle::setStyle("Material");

    // Load main QML file
    engine.load(QUrl("qrc:/resources/EVAAMainWindow.qml"));

    // Get reference to the QVTKFramebufferObjectItem created in QML
    // We cannot use smart pointers because this object must be deleted by QML
    QObject *rootObject = engine.rootObjects().first();
    m_vtkFboItem = rootObject->findChild<QVTKFramebufferObjectItem *>("vtkFboItem");

    // Give the vtkFboItem reference to the EVAAMainWindow
    if (m_vtkFboItem) {
        qDebug() << "EVAAMainWindow::EVAAMainWindow: setting vtkFboItem to EVAAMainWindow";

        m_vtkFboItem->setProcessingEngine(m_processingEngine);

        connect(m_vtkFboItem, &QVTKFramebufferObjectItem::rendererInitialized, this,
                &EVAAMainWindow::startApplication);
        connect(m_vtkFboItem, &QVTKFramebufferObjectItem::isModelSelectedChanged, this,
                &EVAAMainWindow::isModelSelectedChanged);
        connect(m_vtkFboItem, &QVTKFramebufferObjectItem::selectedModelPositionXChanged, this,
                &EVAAMainWindow::selectedModelPositionXChanged);
        connect(m_vtkFboItem, &QVTKFramebufferObjectItem::selectedModelPositionYChanged, this,
                &EVAAMainWindow::selectedModelPositionYChanged);
    }
    else {
        qCritical() << "EVAAMainWindow::EVAAMainWindow: Unable to get vtkFboItem instance";
        return;
    }

    int rc = app.exec();

    qDebug() << "EVAAMainWindow::EVAAMainWindow: Execution finished with return code:" << rc;
}

void EVAAMainWindow::startApplication() const
{
    qDebug() << "EVAAMainWindow::startApplication()";

    disconnect(m_vtkFboItem, &QVTKFramebufferObjectItem::rendererInitialized, this,
               &EVAAMainWindow::startApplication);
}

void EVAAMainWindow::openModel(const QUrl &path) const
{
    qDebug() << "EVAAMainWindow::openModel():" << path;

    QUrl localFilePath;

    if (path.isLocalFile()) {
        // Remove the "file:///" if present
        localFilePath = path.toLocalFile();
    }
    else {
        localFilePath = path;
    }

    m_vtkFboItem->addModelFromFile(localFilePath);
}

bool EVAAMainWindow::isModelExtensionValid(const QUrl &modelPath) const
{
    if (modelPath.toString().toLower().endsWith(".stl") ||
        modelPath.toString().toLower().endsWith(".obj")) {
        return true;
    }

    return false;
}

void EVAAMainWindow::mousePressEvent(const int button, const int screenX, const int screenY) const
{
    qDebug() << "EVAAMainWindow::mousePressEvent()";

    m_vtkFboItem->selectModel(screenX, screenY);
}

void EVAAMainWindow::mouseMoveEvent(const int button, const int screenX, const int screenY)
{
    if (!m_vtkFboItem->isModelSelected()) {
        return;
    }

    if (!m_draggingMouse) {
        m_draggingMouse = true;

        m_previousWorldX = m_vtkFboItem->getSelectedModelPositionX();
        m_previousWorldY = m_vtkFboItem->getSelectedModelPositionY();
    }

    CommandModelTranslate::TranslateParams_t translateParams;

    translateParams.screenX = screenX;
    translateParams.screenY = screenY;

    m_vtkFboItem->translateModel(translateParams, true);
}

void EVAAMainWindow::mouseReleaseEvent(const int button, const int screenX, const int screenY)
{
    qDebug() << "EVAAMainWindow::mouseReleaseEvent()";

    if (!m_vtkFboItem->isModelSelected()) {
        return;
    }

    if (m_draggingMouse) {
        m_draggingMouse = false;

        CommandModelTranslate::TranslateParams_t translateParams;

        translateParams.screenX = screenX;
        translateParams.screenY = screenY;
        translateParams.previousPositionX = m_previousWorldX;
        translateParams.previousPositionY = m_previousWorldY;

        m_vtkFboItem->translateModel(translateParams, false);
    }
}

bool EVAAMainWindow::getIsModelSelected() const
{
    // QVTKFramebufferObjectItem might not be initialized when QML loads
    if (!m_vtkFboItem) {
        return 0;
    }

    return m_vtkFboItem->isModelSelected();
}

double EVAAMainWindow::getSelectedModelPositionX() const
{
    // QVTKFramebufferObjectItem might not be initialized when QML loads
    if (!m_vtkFboItem) {
        return 0;
    }

    return m_vtkFboItem->getSelectedModelPositionX();
}

double EVAAMainWindow::getSelectedModelPositionY() const
{
    // QVTKFramebufferObjectItem might not be initialized when QML loads
    if (!m_vtkFboItem) {
        return 0;
    }

    return m_vtkFboItem->getSelectedModelPositionY();
}

void EVAAMainWindow::setModelsRepresentation(const int representationOption)
{
    m_vtkFboItem->setModelsRepresentation(representationOption);
}

void EVAAMainWindow::setModelsOpacity(const double opacity)
{
    m_vtkFboItem->setModelsOpacity(opacity);
}

void EVAAMainWindow::setGouraudInterpolation(const bool gouraudInterpolation)
{
    m_vtkFboItem->setGouraudInterpolation(gouraudInterpolation);
}

void EVAAMainWindow::setModelColorR(const int colorR) { m_vtkFboItem->setModelColorR(colorR); }

void EVAAMainWindow::setModelColorG(const int colorG) { m_vtkFboItem->setModelColorG(colorG); }

void EVAAMainWindow::setModelColorB(const int colorB) { m_vtkFboItem->setModelColorB(colorB); }

}  // namespace EVAA
