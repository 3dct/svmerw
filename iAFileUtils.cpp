/************************************  SVMERW  ****************************************
* **  Command line tool for segmenting multi-modal 3D data with a Support Vector   ** *
* **  Machine (SVM) and the Extended Random Walker (ERW)                           ** *
***************************************************************************************
* Copyright (C) 2017 Bernhard Fröhler                                                 *
***************************************************************************************
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
***************************************************************************************
* Contact: FH OÖ Forschungs & Entwicklungs GmbH, Campus Wels, CT-Gruppe,              *
*       Stelzhamerstraße 23, 4600 Wels / Austria, Email: bernhard.froehler@fh-wels.at *
**************************************************************************************/
#include "iAFileUtils.h"

#include <QString>

QString MakeAbsolute(QString const & baseDir, QString const & fileName)
{
#if _WIN32
	if (fileName.contains(":"))
#else
	if (fileName.startsWith("/"))
#endif
	{
		return fileName;
	}
	// TODO : use canonicalFilePath here? drawback: empty if file doesn't exist!
	return baseDir + "/" + fileName;
}

QString MakeRelative(QString const & baseDir,  QString const & fileName)
{
	if (fileName.startsWith(baseDir))
	{                                                               // for '/'
		return fileName.right(fileName.length() - baseDir.length() - 1);
	}
	return fileName;
}
