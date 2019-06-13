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
#include "CommandModelAdd.h"
#include "Model.h"
#include "ProcessingEngine.h"
#include "QVTKFramebufferObjectRenderer.h"


CommandModelAdd::CommandModelAdd(QVTKFramebufferObjectRenderer *vtkFboRenderer, std::shared_ptr<ProcessingEngine> processingEngine, QUrl modelPath)
	: m_processingEngine{processingEngine}
	, m_modelPath{modelPath}
{
	m_vtkFboRenderer = vtkFboRenderer;
}


void CommandModelAdd::run()
{
	qDebug() << "CommandModelAdd::run()";

	m_model = m_processingEngine->addModel(m_modelPath);

	m_processingEngine->placeModel(*m_model);

	m_ready = true;
	emit ready();
}


bool CommandModelAdd::isReady() const
{
	return m_ready;
}

void CommandModelAdd::execute()
{
	qDebug() << "CommandModelAdd::execute()";

	m_vtkFboRenderer->addModelActor(m_model);

	emit done();
}
