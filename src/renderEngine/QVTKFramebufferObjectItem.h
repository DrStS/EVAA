/*  Copyright &copy; 2019, Stefan Sicklinger, Munich
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
/*************************************************************************************************
 * \file QVTKFramebufferObjectItem.h
 * This file holds the class of QVTKFramebufferObjectItem.
 * \date 6/13/2019
 **************************************************************************************************/
#pragma once

#include <QtQuick/QQuickFramebufferObject>
#include <memory>
#include <mutex>
#include <queue>

#include "CommandModelTranslate.h"

class CommandModel;
class Model;
class ProcessingEngine;
class QVTKFramebufferObjectRenderer;

class QVTKFramebufferObjectItem : public QQuickFramebufferObject {
    Q_OBJECT

public:
    QVTKFramebufferObjectItem();

    Renderer *createRenderer() const Q_DECL_OVERRIDE;
    void setVtkFboRenderer(QVTKFramebufferObjectRenderer *);
    bool isInitialized() const;
    void setProcessingEngine(const std::shared_ptr<ProcessingEngine> processingEngine);

    // Model releated functions
    bool isModelSelected() const;

    double getSelectedModelPositionX() const;
    double getSelectedModelPositionY() const;

    void selectModel(const int screenX, const int screenY);
    void resetModelSelection();
    void addModelFromFile(const QUrl &modelPath);

    void translateModel(CommandModelTranslate::TranslateParams_t &translateData,
                        const bool inTransition);

    // Camera related functions
    void wheelEvent(QWheelEvent *e) override;
    void mousePressEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void mouseMoveEvent(QMouseEvent *e) override;

    QMouseEvent *getLastMouseLeftButton();
    QMouseEvent *getLastMouseButton();
    QMouseEvent *getLastMoveEvent();
    QWheelEvent *getLastWheelEvent();

    void resetCamera();

    int getModelsRepresentation() const;
    double getModelsOpacity() const;
    bool getGourauInterpolation() const;
    int getModelColorR() const;
    int getModelColorG() const;
    int getModelColorB() const;

    void setModelsRepresentation(const int representationOption);
    void setModelsOpacity(const double opacity);
    void setGouraudInterpolation(const bool gouraudInterpolation);
    void setModelColorR(const int colorR);
    void setModelColorG(const int colorG);
    void setModelColorB(const int colorB);

    CommandModel *getCommandsQueueFront() const;
    void commandsQueuePop();
    bool isCommandsQueueEmpty() const;
    void lockCommandsQueueMutex();
    void unlockCommandsQueueMutex();

signals:
    void rendererInitialized();

    void isModelSelectedChanged();
    void selectedModelPositionXChanged();
    void selectedModelPositionYChanged();

    void addModelFromFileDone();
    void addModelFromFileError(QString error);

private:
    void addCommand(CommandModel *command);

    QVTKFramebufferObjectRenderer *m_vtkFboRenderer = nullptr;
    std::shared_ptr<ProcessingEngine> m_processingEngine;

    std::queue<CommandModel *> m_commandsQueue;
    std::mutex m_commandsQueueMutex;

    std::shared_ptr<QMouseEvent> m_lastMouseLeftButton;
    std::shared_ptr<QMouseEvent> m_lastMouseButton;
    std::shared_ptr<QMouseEvent> m_lastMouseMove;
    std::shared_ptr<QWheelEvent> m_lastMouseWheel;

    int m_modelsRepresentationOption = 2;
    double m_modelsOpacity = 1.0;
    bool m_gouraudInterpolation = false;
    int m_modelColorR = 3;
    int m_modelColorG = 169;
    int m_modelColorB = 244;
};
