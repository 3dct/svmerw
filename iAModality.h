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

#include <QSharedPointer>
#include <QString>
#include <QVector>

class iAImageCoordConverter;
class iASpectralVoxelData;

class vtkImageData;

class iAModality
{
public:
	//! create uninitialized modality
	iAModality();
	//! create modality from name and file
	iAModality(QString const & name, QString const & filename);
	//! create modality from name, file and image data
	iAModality(QString const & name, QString const & filename, vtkSmartPointer<vtkImageData> imgData);
	//! returns name of the modality
	QString GetName() const;
	//! returns file holding the modality data
	QString GetFileName() const;
	//! set name of the modality
	void SetName(QString const & name);
	//! set file containing the modality data
	void SetFileName(QString const & filename);
	//! set flag indicating location where to render
	void SetRenderFlag(int renderFlag);

	int GetWidth() const;
	int GetHeight() const;
	int GetDepth() const;
	double const * GetSpacing() const;
	iAImageCoordConverter const & GetConverter() const;
	
	// TODO: retrieve modality image/data:
	QSharedPointer<iASpectralVoxelData const> GetData() const;

	vtkSmartPointer<vtkImageData> GetImage() const;
	void SetImage(vtkSmartPointer<vtkImageData>) const;

	bool LoadData();

	// TODO: Refactor
	QString positionSettings;
	QString orientationSettings;
	QString tfFileName;
private:
	QString m_name;
	QString m_filename;

	QSharedPointer<iASpectralVoxelData> m_data;
	QSharedPointer<iAImageCoordConverter> m_converter;
	double m_spacing[3];

	// IO-related methods:
	void SetData(vtkSmartPointer<vtkImageData> imgData);
	
	vtkSmartPointer<vtkImageData> m_imgData;
};


typedef QVector<QSharedPointer<iAModality> > ModalityCollection;


class iAModalityList: public QObject
{
	Q_OBJECT
public:
	iAModalityList();
	bool Load(QString const & filename);
	
	int size() const;
	QSharedPointer<iAModality> Get(int idx);
	QSharedPointer<iAModality const> Get(int idx) const;
	void Add(QSharedPointer<iAModality> mod);
	void Remove(int idx);
	QString const & GetFileName() const;
signals:
	void Added(QSharedPointer<iAModality> mod);
private:
	bool ModalityExists(QString const & filename) const;

	ModalityCollection m_modalities;
	double m_spacing[3];
	QString m_fileName;
};
