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

#include "iASpectraDistance.h"

#include <QString>

int GetDistanceMeasureCount();
char const * const * const GetDistanceMeasureNames();
char const * const * const GetShortMeasureNames();
int GetShortNameIdx(QString const & shortName);
QSharedPointer<iASpectraDistance> GetDistanceMeasure(int idx);
//QSharedPointer<iASpectraDistance> GetDistanceMeasure(QString const & distFuncName);
//QSharedPointer<iASpectraDistance> GetDistanceMeasureFromShortName(QString const & distFuncName);

class iASpectralAngularDistance: public iASpectraDistance
{
public:
	virtual char const * GetName() const;
	virtual char const * GetShortName() const;
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2) const;
};

class iAL1NormDistance: public iASpectraDistance
{
public:
	virtual char const * GetName() const;
	virtual char const * GetShortName() const;
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2) const;
};

class iAL2NormDistance: public iASpectraDistance
{
public:
	virtual char const * GetName() const;
	virtual char const * GetShortName() const;
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2) const;
};

class iALInfNormDistance: public iASpectraDistance
{
public:
	virtual char const * GetName() const;
	virtual char const * GetShortName() const;
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2) const;
};

class iAJensenShannonDistance : public iASpectraDistance
{
public:
	virtual char const * GetName() const;
	virtual char const * GetShortName() const;
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2) const;
};

class iAKullbackLeiblerDivergence : public iASpectraDistance
{
public:
	virtual char const * GetName() const;
	virtual char const * GetShortName() const;
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2) const;
	virtual bool isSymmetric() const {return false; }
};

class iAChiSquareDistance: public iASpectraDistance
{
public:
	virtual char const * GetName() const;
	virtual char const * GetShortName() const;
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2) const;
};

class iAEarthMoversDistance: public iASpectraDistance
{
public:
	virtual char const * GetName() const;
	virtual char const * GetShortName() const;
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2) const;
};

class iASquaredDistance: public iASpectraDistance
{
public:
	virtual char const * GetName() const;
	virtual char const * GetShortName() const;
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2) const;
};

class iANullDistance: public iASpectraDistance
{
public:
	virtual char const * GetName() const { return "Null Dist."; }
	virtual char const * GetShortName() const { return "null"; }
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2) const { return 0; }
};

/*
class iAMutualInformation : public iASpectraDistance
{
public:
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2);
	virtual bool isSymmetric() {return true; }
};

class iACrossCorrelation : public iASpectraDistance
{
public:
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2);
	virtual bool isSymmetric() {return true; }
};

class iAMahalanobisDistance : public iASpectraDistance
{
public:
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2);
	virtual bool isSymmetric() {return true; }
};

class iAHausdorffDistance : public iASpectraDistance
{
public:
	virtual double GetDistance(QSharedPointer<iASpectrumType const> spec1, QSharedPointer<iASpectrumType const> spec2);
	virtual bool isSymmetric() {return true; }
};
*/
