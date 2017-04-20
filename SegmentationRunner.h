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

#include <QMap>
#include <QMutex>
#include <QSharedPointer>
#include <QThread>

#include "iAPerformanceHelper.h"
#include "iARandomWalker.h"
#include "iASpectrumType.h"


class iAModalityList;
class iAMMSegParameter;

enum RWMode {
	Standard,
	Extended,
	Normalize
};

class SegmentationRunner
{
public:
	SegmentationRunner(
		QSharedPointer<iAModalityList const> modalities,
		SeedsPointer seeds,
		iAMMSegParameter const & parameter,
		RWMode mode
	);
	QSharedPointer<iARWResult> GetResult();
	void run();
	//void Abort();
	//bool IsAborted();
private:
	//! @{
	//! input
	QSharedPointer<iAModalityList const> m_modalities;
	SeedsPointer m_seeds;
	iAMMSegParameter const & m_param;
	RWMode m_mode;
	//! @}

	bool m_aborted;


	//! @{
	//! Performance Measurement
	iAPerformanceTimer m_erwPerfTimer;
	iAPerformanceTimer::DurationType m_calcERWDuration;
	//! @}

	iAExtendedRandomWalker* m_runningERW;
	QMutex m_mutex;
	QSharedPointer<iARWResult> m_result;
};
