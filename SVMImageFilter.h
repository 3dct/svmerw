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

#include <vtkSmartPointer.h>

#include <QList>
#include <QSharedPointer>
#include <QVector>

#include "SeedType.h"

struct iAImageCoordinate;
class iAModalityList;

class vtkImageData;

class SVMImageFilter
{
public:
	typedef QVector<vtkSmartPointer<vtkImageData> >	ProbabilityImagesType;
	typedef QSharedPointer<ProbabilityImagesType> ProbabilityImagesPointer;
	typedef QSharedPointer<iAModalityList const> ModalitiesPointer;

	SVMImageFilter(double c, double gamma,
		ModalitiesPointer modalities,
		SeedsPointer seeds,
		int channelCount); // TODO: replace by bool mask to allow selective enabling/disabling of channels
	void Run();
	ProbabilityImagesPointer GetResult();
private:
	//! @{
	//! Input
	double m_c;
	double m_gamma;
	SeedsPointer m_seeds;
	ModalitiesPointer m_modalities;
	int m_channelCount;
	//! @}
	//! @{
	//! Output
	ProbabilityImagesPointer m_probabilities;
	QString m_error;
	//! @}
};
