/*********************************  open_iA 2016 06  ******************************** *
* **********  A tool for scientific visualisation and 3D image processing  ********** *
* *********************************************************************************** *
* Copyright (C) 2016  C. Heinzl, M. Reiter, A. Reh, W. Li, M. Arikan, J. Weissenböck, *
*                     Artem & Alexander Amirkhanov, B. Fröhler                        *
* *********************************************************************************** *
* This program is free software: you can redistribute it and/or modify it under the   *
* terms of the GNU General Public License as published by the Free Software           *
* Foundation, either version 3 of the License, or (at your option) any later version. *
*                                                                                     *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY     *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A     *
* PARTICULAR PURPOSE.  See the GNU General Public License for more details.           *
*                                                                                     *
* You should have received a copy of the GNU General Public License along with this   *
* program.  If not, see http://www.gnu.org/licenses/                                  *
* *********************************************************************************** *
* Contact: FH OÖ Forschungs & Entwicklungs GmbH, Campus Wels, CT-Gruppe,              *
*          Stelzhamerstraße 23, 4600 Wels / Austria, Email:                           *
* ************************************************************************************/
#include "iAConsole.h"

#include <iostream>
#include <fstream>

void iAConsole::Log(std::string const & text)
{
	Log(QString::fromStdString(text));
}

void iAConsole::Log(char const * text)
{
	Log(QString(text));
}

void iAConsole::Log(QString const & text)
{
	emit LogSignal(text);
}

void iAConsole::LogSlot(QString const & text)
{
	std::cout << text.toStdString() << std::endl;
	std::ofstream logfile("debug.log", std::ofstream::out | std::ofstream::app);
	logfile << text.toStdString() << std::endl;
	logfile.flush();
}

iAConsole::iAConsole()
{
	connect(this, SIGNAL(LogSignal(QString const &)), this, SLOT(LogSlot(QString const &)));
}

iAConsole::~iAConsole()
{
}

iAConsole& iAConsole::GetInstance()
{
	static iAConsole instance;
	return instance;
}


void iAConsole::Close()
{
}
