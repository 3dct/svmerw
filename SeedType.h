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
#pragma once

#include "iAImageCoordinate.h"

//! class representing a single seed
//! either through its coordinate value, and optionally through its multiple feature values
class Seed
{
public:
	bool hasValues() const
	{
		return m_featureValue && m_featureValue->size() > 0;
	}
	iAImageCoordinate const & coord() const
	{
		return m_coord;
	}
	int featureCount() const
	{
		return hasValues() ? m_featureValue->size() : -1;
	}
	double featureValue(int idx) const
	{
		return (*m_featureValue)[idx];
	}

	void setCoordinate(int x, int y, int z)
	{
		m_coord.x = x;
		m_coord.y = y;
		m_coord.z = z;
	}
	void setValues(QSharedPointer<QVector<double> > featureValues)
	{
		m_featureValue = featureValues;
	}
private:
	iAImageCoordinate m_coord;
	QSharedPointer<QVector<double> > m_featureValue;
};

// a list (QList) of seeds per label (QVector)
typedef QVector<QList<Seed> > Seeds;
typedef QSharedPointer<Seeds> SeedsPointer;
