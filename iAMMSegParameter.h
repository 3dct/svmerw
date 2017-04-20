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

#include <QSharedPointer>
#include <QVector>

class iASpectraDistance;
class iAMMSegParameterRange;

class iAMMSegModalityParameter
{
public:
	iAMMSegModalityParameter() : weight(0), distanceFuncIdx(0), pcaDim(0) {}
	iAMMSegModalityParameter(int p, double w, int d) :
		weight(w),
		distanceFuncIdx(d),
		pcaDim(p)
	{}
	double weight;
	int distanceFuncIdx;
	int pcaDim;
};

class iAMMSegParameter
{
public:
	//iAMMSegParameter(int modalityCount);

	iAMMSegParameter(double erw_beta, double erw_gamma, int erw_maxIter,
		double svm_c, double svm_gamma, int svm_channels,
		bool do_gad, int gad_iter, double gad_step, double gad_cond);

	//static QSharedPointer<iAMMSegParameter> Create(QString const & Descriptor, int modalityCount);
	//QString GetDescriptor() const;

	void setID(int id);
	int id() const;
	double erw_beta() const;
	double erw_gamma() const;
	int erw_maxIter() const;
	double svm_c() const;
	double svm_gamma() const;
	int svm_channels() const;

	bool doGAD() const;
	int gad_iter() const;
	double gad_step() const;
	double gad_cond() const;

	void addModalityParam(int pca, double weight, int distanceFuncIdx);

	int pcaDim(int modalityIdx) const;
	int distanceFuncIdx(int modalityIdx) const;
	QSharedPointer<iASpectraDistance> distanceFunc(int modalityIdx) const;
	double modalityWeight(int modalityIdx) const;
	int modalityCount() const;
private:
	int m_id;
	// for smoothing with Gradient Anisotropic Diffusion:
	bool m_do_gad;
	int m_gad_iter;
	double m_gad_step;
	double m_gad_cond;
	// for ERW:
	double m_erw_gamma;
	double m_erw_beta;
	int m_erw_maxIter;
	// for SVM:
	double m_svm_c;
	double m_svm_gamma;
	double m_svm_channels;

	int m_modalityCount;
	// for each channel:
	QVector<iAMMSegModalityParameter> m_modalityParams;
};
