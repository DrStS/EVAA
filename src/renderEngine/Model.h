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
 * \file Model.h
 * This file holds the class of Model.
 * \date 6/13/2019
 */

#pragma once

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>
#include <vtkTransformPolyDataFilter.h>

#include <QColor>
#include <QObject>
#include <memory>
#include <mutex>

class Model : public QObject {
    Q_OBJECT

public:
    Model(vtkSmartPointer<vtkPolyData> modelData);

    const vtkSmartPointer<vtkActor> &getModelActor() const;

    double getPositionX();
    double getPositionY();

    void translateToPosition(const double x, const double y);

    void setSelected(const bool selected);
    static void setSelectedModelColor(const QColor &selectedModelColor);

    const double getMouseDeltaX() const;
    const double getMouseDeltaY() const;
    void setMouseDeltaXY(const double deltaX, const double deltaY);

    void updateModelColor();

signals:
    void positionXChanged(const double positionX);
    void positionYChanged(const double positionY);

private:
    void setPositionX(const double positionX);
    void setPositionY(const double positionY);

    void setColor(const QColor &color);

    static QColor m_defaultModelColor;
    static QColor m_selectedModelColor;

    vtkSmartPointer<vtkPolyData> m_modelData;
    vtkSmartPointer<vtkPolyDataMapper> m_modelMapper;
    vtkSmartPointer<vtkActor> m_modelActor;

    vtkSmartPointer<vtkTransformPolyDataFilter> m_modelFilterTranslate;

    std::mutex m_propertiesMutex;

    double m_positionX{0.0};
    double m_positionY{0.0};
    double m_positionZ{0.0};

    bool m_selected = false;

    double m_mouseDeltaX = 0.0;
    double m_mouseDeltaY = 0.0;
};
