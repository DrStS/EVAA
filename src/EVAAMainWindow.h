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
 * \file EVAAMainWindow.h
 * This file holds the class of StartWindow.
 * \date 6/12/2019
 */

#pragma once

#include <QObject>
#include <QUrl>
#include <memory>

class ProcessingEngine;
class QVTKFramebufferObjectItem;
/**
 * \brief Class EVAAMainWindow the core of the GUI
 */
class EVAAMainWindow : public QObject {
    Q_OBJECT
    Q_PROPERTY(bool showFileDialog MEMBER m_showFileDialog NOTIFY showFileDialogChanged)
    Q_PROPERTY(bool isModelSelected READ getIsModelSelected NOTIFY isModelSelectedChanged)
    Q_PROPERTY(
        double modelPositionX READ getSelectedModelPositionX NOTIFY selectedModelPositionXChanged)
    Q_PROPERTY(
        double modelPositionY READ getSelectedModelPositionY NOTIFY selectedModelPositionYChanged)

public:
    /**
     * \brief Constructor
     * \author Stefan Sicklinger
     */
    EVAAMainWindow(int argc, char **argv);

    Q_INVOKABLE void openModel(const QUrl &path) const;
    Q_INVOKABLE void mousePressEvent(const int button, const int mouseX, const int mouseY) const;
    Q_INVOKABLE void mouseMoveEvent(const int button, const int mouseX, const int mouseY);
    Q_INVOKABLE void mouseReleaseEvent(const int button, const int mouseX, const int mouseY);

    bool getIsModelSelected() const;
    double getSelectedModelPositionX() const;
    double getSelectedModelPositionY() const;

    Q_INVOKABLE void setModelsRepresentation(const int representationOption);
    Q_INVOKABLE void setModelsOpacity(const double opacity);
    Q_INVOKABLE void setGouraudInterpolation(const bool gouraudInterpolation);
    Q_INVOKABLE void setModelColorR(const int colorR);
    Q_INVOKABLE void setModelColorG(const int colorG);
    Q_INVOKABLE void setModelColorB(const int colorB);

public slots:
    void startApplication() const;

signals:
    void showFileDialogChanged();

    void isModelSelectedChanged();
    void selectedModelPositionXChanged();
    void selectedModelPositionYChanged();

private:
    bool isModelExtensionValid(const QUrl &modelPath) const;
    std::shared_ptr<ProcessingEngine> m_processingEngine;
    QVTKFramebufferObjectItem *m_vtkFboItem = nullptr;
    double m_previousWorldX = 0;
    double m_previousWorldY = 0;
    bool m_draggingMouse = false;
    bool m_showFileDialog = false;
};
