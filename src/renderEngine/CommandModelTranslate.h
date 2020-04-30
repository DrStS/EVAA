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
 * \file CommandModelTranslate.h
 * This file holds the class of CommandModelTranslate.
 * \date 6/12/2019
 */

#pragma once

#include <memory>

#include "CommandModel.h"

namespace EVAA {

class Model;
class QVTKFramebufferObjectRenderer;

class CommandModelTranslate : public CommandModel {
public:
    typedef struct {
        std::shared_ptr<Model> model;
        int screenX{0};
        int screenY{0};
        double previousPositionX{0};
        double previousPositionY{0};
        double targetPositionX{0};
        double targetPositionY{0};
    } TranslateParams_t;

    CommandModelTranslate(QVTKFramebufferObjectRenderer *vtkFboRenderer,
                          const TranslateParams_t &translateVector, bool inTransition);

    bool isReady() const override;
    void execute() override;

private:
    void transformCoordinates();

    TranslateParams_t m_translateParams;
    bool m_inTransition;
    bool m_needsTransformation = true;
};

}  // namespace EVAA
