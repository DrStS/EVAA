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
 * \file ProcessingEngine.h
 * This file holds the class of ProcessingEngine.
 * \date 6/13/2019
 */

#pragma once

#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <QUrl>
#include <array>
#include <cstdint>
#include <memory>
#include <mutex>
#include <vector>

namespace EVAA {

class Model;

class ProcessingEngine {
public:
    ProcessingEngine();

    const std::shared_ptr<Model> &addModel(const QUrl &modelFilePath);

    void placeModel(Model &model) const;

    void setModelsRepresentation(const int modelsRepresentationOption) const;
    void setModelsOpacity(const double modelsOpacity) const;
    void setModelsGouraudInterpolation(const bool enableGouraudInterpolation) const;
    void updateModelsColor() const;

    std::shared_ptr<Model> getModelFromActor(const vtkSmartPointer<vtkActor> modelActor) const;

private:
    vtkSmartPointer<vtkPolyData> preprocessPolydata(const vtkSmartPointer<vtkPolyData> inputData) const;

    std::vector<std::shared_ptr<Model>> m_models;
};

}  // namespace EVAA
