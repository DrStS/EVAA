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
* \file CommandModelAdd.h
* This file holds the class of CommandModelAdd.
* \date 6/12/2019
**************************************************************************************************/
#pragma once

#include <memory>
#include <QUrl>
#include <QThread>
#include "CommandModel.h"


class Model;
class ProcessingEngine;
class QVTKFramebufferObjectRenderer;

class CommandModelAdd : public QThread, public CommandModel
{
	Q_OBJECT

public:
	CommandModelAdd(QVTKFramebufferObjectRenderer *vtkFboRenderer, std::shared_ptr<ProcessingEngine> processingEngine, QUrl modelPath);

	void run() Q_DECL_OVERRIDE;

	bool isReady() const override;
	void execute() override;

signals:
	void ready();
	void done();

private:
	std::shared_ptr<ProcessingEngine> m_processingEngine;
	std::shared_ptr<Model> m_model = nullptr;
	QUrl m_modelPath;
	double m_positionX;
	double m_positionY;

	bool m_ready = false;
};
