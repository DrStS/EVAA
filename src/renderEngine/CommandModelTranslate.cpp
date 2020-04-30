/*  Copyright &copy; 2019, Dr. Stefan Sicklinger, Munich \n
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
#include "CommandModelTranslate.h"

#include <array>

#include "Model.h"
#include "QVTKFramebufferObjectRenderer.h"

CommandModelTranslate::CommandModelTranslate(QVTKFramebufferObjectRenderer *vtkFboRenderer,
                                             const TranslateParams_t &translateData,
                                             bool inTransition) :
    m_translateParams{translateData}, m_inTransition{inTransition}
{
    m_vtkFboRenderer = vtkFboRenderer;
}

bool CommandModelTranslate::isReady() const { return true; }

void CommandModelTranslate::transformCoordinates()
{
    std::array<double, 3> worldCoordinates;

    if (m_vtkFboRenderer->screenToWorld(m_translateParams.screenX, m_translateParams.screenY,
                                        worldCoordinates.data())) {
        m_translateParams.targetPositionX =
            worldCoordinates[0] - m_translateParams.model->getMouseDeltaX();
        m_translateParams.targetPositionY =
            worldCoordinates[1] - m_translateParams.model->getMouseDeltaY();
    }
    else {
        m_translateParams.targetPositionX = m_translateParams.model->getPositionX();
        m_translateParams.targetPositionY = m_translateParams.model->getPositionY();
    }

    m_needsTransformation = false;
}

void CommandModelTranslate::execute()
{
    // Screen to world transformation can only be done within the Renderer thread
    if (m_needsTransformation) {
        this->transformCoordinates();
    }

    m_translateParams.model->translateToPosition(m_translateParams.targetPositionX,
                                                 m_translateParams.targetPositionY);
}
