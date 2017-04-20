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

#include "iAImageGraphTypes.h"
#include "SeedType.h"

#include <QSharedPointer>
#include <QThread>
#include <QVector>

#include "iAImageTypes.h"

class iANormalizer;
class iALabelMapper;
class iAImageGraph;
class iASpectraDistance;

struct iARWInputChannel
{
	QSharedPointer<iASpectralVoxelData const> image;
	QSharedPointer<iASpectraDistance> distanceFunc;
	QSharedPointer<iANormalizer> normalizeFunc;
	double             weight;
};

struct iARWResult
{
	LabelImagePointer labelledImage;
	QVector<ProbabilityImagePointer> erwProbImages;
	QVector<PriorModelImagePointer> svmProbImages;
};

class iARandomWalker//: public QThread
{
public:
	//! initialize "standard" random walker
	iARandomWalker(iAVoxelIndexType width,
		iAVoxelIndexType height,
		iAVoxelIndexType depth,
		double const spacing[3],
		QSharedPointer<QVector<iARWInputChannel> > inputChannels,
		SeedsPointer seeds,
		int maxIter
	);
	QSharedPointer<iARWResult> GetResult();
	virtual void run();
private:
	int m_vertexCount;
	QSharedPointer<QVector<iARWInputChannel> > m_inputChannels;
	QSharedPointer<iAImageGraph> m_imageGraph;
	SeedsPointer m_seeds;
	int m_minLabel, m_maxLabel;
	double m_spacing[3];
	int m_maxIterations;
	QSharedPointer<iARWResult>  m_result;

};

class vtkImageData;
class vtkPolyData;
class QObject;

class iAExtendedRandomWalker//: public QThread
{
public:
	//! initialize "extended" random walker
	iAExtendedRandomWalker(
		iAVoxelIndexType width,
		iAVoxelIndexType height,
		iAVoxelIndexType depth,
		double const spacing[3],
		QSharedPointer<QVector<iARWInputChannel> > inputChannels,
		QSharedPointer<QVector<PriorModelImagePointer> > priorModel,
		double priorModelWeight,
		int maxIterations
	);
	iAExtendedRandomWalker(
		iAVoxelIndexType width,
		iAVoxelIndexType height,
		iAVoxelIndexType depth,
		double const spacing[3],
		QSharedPointer<QVector<iARWInputChannel> > inputChannels,
		QSharedPointer<QVector<PriorModelImagePointer> > priorModel,
		double priorModelWeight,
		int maxIterations,
		vtkImageData* i,
		vtkPolyData* p,
		QObject* parent
	);
	QSharedPointer<iARWResult> GetResult();

	virtual void run();
private:
	int m_vertexCount;
	QSharedPointer<QVector<iARWInputChannel> > m_inputChannels;
	QSharedPointer<iAImageGraph> m_imageGraph;
	QSharedPointer<QVector<PriorModelImagePointer> > m_priorModel;
	double m_priorModelWeight;
	int m_labelCount;
	double m_spacing[3];
	QSharedPointer<iARWResult>  m_result;
	int m_maxIterations;
};
