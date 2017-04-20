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
#include "iAMMSegParameter.h"

#include "iASpectraDistanceImpl.h"
#include "iAConsole.h"

#include <QStringList>

iAMMSegParameter::iAMMSegParameter(double erw_beta, double erw_gamma, int erw_maxIter,
	double svm_c, double svm_gamma, int svm_channels, bool do_gad, int gad_iter, double gad_step, double gad_cond)
:
	m_id(-1),
	m_erw_beta(erw_beta),
	m_erw_gamma(erw_gamma),
	m_erw_maxIter(erw_maxIter),
	m_svm_c(svm_c),
	m_svm_gamma(svm_gamma),
	m_svm_channels(svm_channels),
	m_do_gad(do_gad),
	m_gad_iter(gad_iter),
	m_gad_step(gad_step),
	m_gad_cond(gad_cond)
{}

void iAMMSegParameter::addModalityParam(int pca, double weight, int distanceFuncIdx)
{
	m_modalityParams.push_back(iAMMSegModalityParameter(pca, weight, distanceFuncIdx));
}

int iAMMSegParameter::modalityCount() const
{
	return m_modalityParams.size();
}

void iAMMSegParameter::setID(int id){
	m_id = id;
}

int iAMMSegParameter::id() const
{
	return m_id;
}
double iAMMSegParameter::erw_beta() const
{
	return m_erw_beta;
}
double iAMMSegParameter::erw_gamma() const
{
	return m_erw_gamma;
}
int iAMMSegParameter::erw_maxIter() const
{
	return m_erw_maxIter;
}
double iAMMSegParameter::svm_c() const
{
	return m_svm_c;
}
double iAMMSegParameter::svm_gamma() const
{
	return m_svm_gamma;
}
int iAMMSegParameter::svm_channels() const
{
	return m_svm_channels;
}


bool iAMMSegParameter::doGAD() const
{
	return m_do_gad;
}

int iAMMSegParameter::gad_iter() const
{
	return m_gad_iter;
}

double iAMMSegParameter::gad_step() const
{
	return m_gad_step;
}

double iAMMSegParameter::gad_cond() const
{
	return m_gad_cond;
}


int iAMMSegParameter::pcaDim(int modalityIdx) const
{
	return m_modalityParams[modalityIdx].pcaDim;
}

QSharedPointer<iASpectraDistance> iAMMSegParameter::distanceFunc(int modalityIdx) const
{
	return GetDistanceMeasure(m_modalityParams[modalityIdx].distanceFuncIdx);
}

int iAMMSegParameter::distanceFuncIdx(int modalityIdx) const
{
	return m_modalityParams[modalityIdx].distanceFuncIdx;
}


double iAMMSegParameter::modalityWeight(int modalityIdx) const
{
	return m_modalityParams[modalityIdx].weight;
}
