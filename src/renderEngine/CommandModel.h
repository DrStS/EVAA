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
* \file CommandModel.h
* This file holds the class of CommandModel.
* \date 6/13/2019
**************************************************************************************************/
#pragma once

class QVTKFramebufferObjectRenderer;
class CommandModel
{
public:
	CommandModel(){}
	virtual ~CommandModel(){}

	virtual bool isReady() const = 0;
	virtual void execute() = 0;

protected:
	QVTKFramebufferObjectRenderer *m_vtkFboRenderer;
};
