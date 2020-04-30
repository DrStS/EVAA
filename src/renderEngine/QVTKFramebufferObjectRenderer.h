/*
 * Copyright &copy; 2019, Stefan Sicklinger, Munich
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

/**
 * \file QVTKFramebufferObjectRenderer.h
 * This file holds the class of QVTKFramebufferObjectRenderer.
 * \date 6/13/2019
 */

#pragma once

#include <vtkActor.h>
#include <vtkCellPicker.h>
#include <vtkCubeSource.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkGenericRenderWindowInteractor.h>
#include <vtkObject.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <QDir>
#include <QOpenGLFramebufferObject>
#include <QOpenGLFunctions>
#include <QQuickFramebufferObject>
#include <QUndoStack>
#include <memory>
#include <mutex>
#include <vector>

class Model;
class QVTKFramebufferObjectItem;
class ProcessingEngine;

class QVTKFramebufferObjectRenderer : public QObject,
                                      public QQuickFramebufferObject::Renderer,
                                      protected QOpenGLFunctions {
    Q_OBJECT

public:
    QVTKFramebufferObjectRenderer();

    void setProcessingEngine(const std::shared_ptr<ProcessingEngine> processingEngine);

    virtual void synchronize(QQuickFramebufferObject *item);
    virtual void render();
    virtual void openGLInitState();
    QOpenGLFramebufferObject *createFramebufferObject(const QSize &size);

    void addModelActor(const std::shared_ptr<Model> model);

    std::shared_ptr<Model> getSelectedModel() const;
    bool isModelSelected() const;

    void setSelectedModelPositionX(const double positionX);
    void setSelectedModelPositionY(const double positionY);

    double getSelectedModelPositionX() const;
    double getSelectedModelPositionY() const;

    void resetCamera();
    const bool screenToWorld(const int16_t screenX, const int16_t screenY, double worldPos[]);

signals:
    void isModelSelectedChanged();

    void selectedModelPositionXChanged();
    void selectedModelPositionYChanged();

private:
    void initScene();
    void generatePlatform();
    void updatePlatform();

    void selectModel(const int16_t x, const int16_t y);
    void clearSelectedModel();
    void setIsModelSelected(const bool isModelSelected);

    void createLine(const double x1, const double y1, const double z1, const double x2,
                    const double y2, const double z2, vtkSmartPointer<vtkPoints> points,
                    vtkSmartPointer<vtkCellArray> cells);
    std::shared_ptr<Model> getSelectedModelNoLock() const;

    std::shared_ptr<ProcessingEngine> m_processingEngine;
    QVTKFramebufferObjectItem *m_vtkFboItem = nullptr;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> m_vtkRenderWindow;
    vtkSmartPointer<vtkRenderer> m_renderer;
    vtkSmartPointer<vtkGenericRenderWindowInteractor> m_vtkRenderWindowInteractor;

    vtkSmartPointer<vtkCellPicker> m_picker;

    std::shared_ptr<Model> m_selectedModel = nullptr;
    vtkSmartPointer<vtkActor> m_selectedActor = nullptr;
    bool m_isModelSelected = false;

    double m_selectedModelPositionX = 0.0;
    double m_selectedModelPositionY = 0.0;

    std::shared_ptr<QMouseEvent> m_mouseLeftButton = nullptr;
    std::shared_ptr<QMouseEvent> m_mouseEvent = nullptr;
    std::shared_ptr<QMouseEvent> m_moveEvent = nullptr;
    std::shared_ptr<QWheelEvent> m_wheelEvent = nullptr;

    vtkSmartPointer<vtkCubeSource> m_platformModel;
    vtkSmartPointer<vtkPolyData> m_platformGrid;
    vtkSmartPointer<vtkActor> m_platformModelActor;
    vtkSmartPointer<vtkActor> m_platformGridActor;

    double m_platformWidth = 200.0;
    double m_platformDepth = 200.0;
    double m_platformHeight = 200.0;
    double m_platformThickness = 2.0;
    double m_gridBottomHeight = 0.15;
    uint16_t m_gridSize = 10;

    double m_camPositionX;
    double m_camPositionY;
    double m_camPositionZ;

    double m_clickPositionZ = 0.0;

    bool m_firstRender = true;

    int m_modelsRepresentationOption = 0;
    double m_modelsOpacity = 1.0;
    bool m_modelsGouraudInterpolation = false;
};
