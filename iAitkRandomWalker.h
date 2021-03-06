/************************************  SVMERW  ****************************************
* **  Command line tool for segmenting multi-modal 3D data with a Support Vector   ** *
* **  Machine (SVM) and the Extended Random Walker (ERW)                           ** *
***************************************************************************************
* Copyright (C) 2017 Bernhard Fr�hler                                                 *
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
* Contact: FH O� Forschungs & Entwicklungs GmbH, Campus Wels, CT-Gruppe,              *
*       Stelzhamerstra�e 23, 4600 Wels / Austria, Email: bernhard.froehler@fh-wels.at *
**************************************************************************************/
#pragma once

#include "iARandomWalker.h"

// wrapper classes to make the (extended) random walker applicable for arbitrary itk imges

template <class TInputImage>
class iAitkRandomWalker
{
public:
	typedef iAitkRandomWalker				  Self;
	iAitkRandomWalker();
	~iAitkRandomWalker();
	void SetInput(TInputImage* image, SeedVector, double beta);
	bool Success() const;
	LabelImagePointer GetLabelImage();
	QVector<ProbabilityImagePointer> GetProbabilityImages();
	void Calculate();
private:
	// not implemented on purpose:
	iAitkRandomWalker(const Self &);
	void operator=(const Self &);
	
	QSharedPointer<iARandomWalker> m_randomWalker;
	
	TInputImage* m_input;
	SeedVector m_seeds;
	double m_beta;
	QSharedPointer<iARWResult> m_result;
};

template <class TInputImage>
class iAitkExtendedRandomWalker
{
public:
	typedef iAitkExtendedRandomWalker				  Self;
	iAitkExtendedRandomWalker();
	~iAitkExtendedRandomWalker();
	void SetInput(TInputImage* image);
	void AddPriorModel(PriorModelImagePointer priorModel);
	bool Success() const;
	LabelImagePointer GetLabelImage();
	QVector<ProbabilityImagePointer> GetProbabilityImages();
	void Calculate(); 
private:
	// not implemented on purpose:
	iAitkExtendedRandomWalker(const Self &);
	void operator=(const Self &);
	
	QSharedPointer<iAExtendedRandomWalker> m_extendedRandomWalker;
	
	TInputImage* m_input;
	QSharedPointer<QVector<PriorModelImagePointer> > m_priorModel;
	QSharedPointer<iARWResult> m_result;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "iAitkRandomWalker.hxx"
#endif
