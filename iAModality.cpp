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
#include "iAModality.h"

#include "iAConsole.h"
#include "iAImageCoordinate.h"
#include "iAitkImagesMultiChannelAdapter.h"

#include <QFileInfo>
#include <QSettings>
#include <QTextStream>

#include <cassert>

iAModality::iAModality():
	m_name(""),
	m_filename("")
{
}

iAModality::iAModality(QString const & name, QString const & filename):
	m_name(name),
	m_filename(filename)
{
}

iAModality::iAModality(QString const & name, QString const & filename, vtkSmartPointer<vtkImageData> imgData) :
	m_name(name),
	m_filename(filename)
{
	SetData(imgData);
}

QString iAModality::GetName() const
{
	return m_name;
}

QString iAModality::GetFileName() const
{
	return m_filename;
}

void iAModality::SetName(QString const & name)
{
	m_name = name;
}

void iAModality::SetFileName(QString const & filename)
{
	m_filename = filename;
}
	
QSharedPointer<iASpectralVoxelData const> iAModality::GetData() const
{
	assert(m_data);
	return m_data;
}

int iAModality::GetWidth() const
{
	assert(m_converter);
	return m_converter->GetWidth();
}

int iAModality::GetHeight() const
{
	assert(m_converter);
	return m_converter->GetHeight();
}

int iAModality::GetDepth() const
{
	assert(m_converter);
	return m_converter->GetDepth();
}

double const * iAModality::GetSpacing() const
{
	assert(m_data);
	return m_spacing;
}

iAImageCoordConverter const & iAModality::GetConverter() const
{
	assert(m_converter);
	return *m_converter;
}

vtkSmartPointer<vtkImageData> iAModality::GetImage() const
{
	return m_imgData;
}

void iAModality::SetImage(vtkSmartPointer<vtkImageData> modImage) const
{
	m_imgData->DeepCopy(modImage);
}


#include <vtkMetaImageReader.h>

vtkSmartPointer<vtkImageData> loadMetaImageFile(QString const & cfn)
{
	QString fileName(cfn);
	fileName = fileName.replace(QString("/"), QString("\\"));
	fileName = fileName.replace(QString("\""), QString(""));
	vtkSmartPointer<vtkMetaImageReader> reader = vtkSmartPointer<vtkMetaImageReader>::New();
	reader->SetFileName(fileName.toStdString().c_str());
	reader->Update();
	return reader->GetOutput();
}

QString getParameterValues(QString fileName, QString parameter, int index, QString section = "", QString sep = ":")
{
	if (index<0 || index>3)
		return 0;
	QString values[4];
	QFile file(fileName);

	if (file.open(QIODevice::ReadOnly))
	{
		QTextStream textStream(&file);
		QString currentLine, currentSection;
		while (!textStream.atEnd())
		{
			currentLine = textStream.readLine();
			if (!currentLine.isEmpty())
			{
				QString currentParameter;
				if ((currentLine.indexOf("[") == 0) && (currentLine.indexOf("[") < currentLine.indexOf("]")))
					currentSection = currentLine.section(" ", 0, 0);
				else
					currentParameter = currentLine.section(" ", 0, 0);

				if ((section != "") && (!currentSection.startsWith(section))) continue;
				else
				{
					if (currentParameter.startsWith(parameter))
					{
						QString temp = currentLine.remove(0, currentLine.indexOf(sep) + 1);
						temp = temp.simplified();
						if (currentParameter == "name")
						{
							values[index] = temp;
							return values[index];
						}

						values[0] = temp.section(" ", 0, 0).trimmed();
						values[1] = temp.section(" ", 1, 1).trimmed();
						values[2] = temp.section(" ", 2, 2).trimmed();

						if (currentLine.indexOf("%") >= 0)
							values[3] = temp.section("%", 1, 1);
					}
				}
			}
		}
		file.close();
	}
	else return "";

	return values[index];
}

bool iAModality::LoadData()
{
	if (m_filename.endsWith(".volstack"))
	{
		std::vector<vtkSmartPointer<vtkImageData> > volumes;

		QFileInfo fi(m_filename);
		QString fileNamesBase = fi.absolutePath() + "/" + getParameterValues(m_filename, "file_names_base:", 0);
		QString extension = getParameterValues(m_filename, "extension:", 0);
		int digitsInIndex = getParameterValues(m_filename, "number_of_digits_in_index:", 0).toInt();
		int indexRange[2] = {
			getParameterValues(m_filename, "minimum_index:", 0).toInt(),
			getParameterValues(m_filename, "maximum_index:", 0).toInt() };

		std::vector<QString> fileNames;
		for (int i = indexRange[0]; i <= indexRange[1]; i++)
		{
			QString temp = fileNamesBase + QString("%1").arg(i, digitsInIndex, 10, QChar('0')) + extension;
			fileNames.push_back(temp);
		}

		if (fileNames.size() == 0)
		{
			DEBUG_LOG("  No files matched the given criteria!");
			return false;
		}
		for (int m = 0; m < fileNames.size(); m++)
		{
			QString fileName = fileNames[m];
			vtkSmartPointer<vtkImageData> img = loadMetaImageFile(fileName);
			if (!img)
			{
				return false;
			}
			volumes.push_back(img);
		}
		assert(volumes.size() > 0);
		if (volumes.size() == 0)
		{
			DEBUG_LOG("No volume found in stack!");
			return false;
		}
		int channels = volumes.size();
		int extent[6];
		volumes[0]->GetExtent(extent);
		volumes[0]->GetSpacing(m_spacing);

		m_converter = QSharedPointer<iAImageCoordConverter>(new iAImageCoordConverter(
			extent[1]-extent[0]+1, extent[3]-extent[2]+1, extent[5]-extent[4]+1));
		
		QSharedPointer<iAvtkImagesMultiChannelAdapter> data(new iAvtkImagesMultiChannelAdapter(
			m_converter->GetWidth(), m_converter->GetHeight(), m_converter->GetDepth()));

		for (int i=0; i<volumes.size(); ++i)
		{
			data->AddImage(volumes[i]);
		}
		m_data = data;
		// TODO: make all channel images available somehow?
		m_imgData = volumes[0];
	}
	else
	{
		// TODO: use ITK image loading?
		vtkSmartPointer<vtkImageData> img = loadMetaImageFile(m_filename);
		SetData(img);
	}
	return true;
}

bool Str2Vec3D(QString const & str, double vec[3])
{
	QStringList list = str.split(" ");
	if (list.size() != 3)
	{
		return false;
	}
	for (int i = 0; i < 3; ++i)
	{
		bool ok;
		vec[i] = list[i].toDouble(&ok);
		if (!ok)
			return false;
	}
	return true;
}

void iAModality::SetData(vtkSmartPointer<vtkImageData> imgData)
{
	assert(imgData);
	m_imgData = imgData;
	int extent[6];
	imgData->GetExtent(extent);
	imgData->GetSpacing(m_spacing);
	m_converter = QSharedPointer<iAImageCoordConverter>(new iAImageCoordConverter(
		extent[1]-extent[0]+1, extent[3]-extent[2]+1, extent[5]-extent[4]+1));
	QSharedPointer<iAvtkImagesMultiChannelAdapter> data(new iAvtkImagesMultiChannelAdapter(
		m_converter->GetWidth(), m_converter->GetHeight(), m_converter->GetDepth()));
	data->AddImage(imgData);
	m_data = data;
}


// iAModalityList
#include "iAFileUtils.h"

namespace
{
	QString GetModalityKey(int idx, QString key)
	{
		return QString("Modality")+QString::number(idx)+"/"+key;
	}

	static const QString FileVersionKey("FileVersion");
	static const QString ModFileVersion("1.0");
	static const QString SetSpacingToOneKey("SetSpacingToOne");

	static const QString CameraPositionKey("CameraPosition");
	static const QString CameraFocalPointKey("CameraFocalPoint");
	static const QString CameraViewUpKey("CameraViewUp");
}


iAModalityList::iAModalityList()
{
	m_spacing[0] = m_spacing[1] = m_spacing[2] = 1.0;
}

bool iAModalityList::ModalityExists(QString const & filename) const
{
	foreach (QSharedPointer<iAModality> mod, m_modalities)
	{
		if (mod->GetFileName() == filename)
		{
			return true;
		}
	}
	return false;
}

QString const & iAModalityList::GetFileName() const
{
	return m_fileName;
}

QString Vec3D2String(double* vec)
{
	return QString("%1 %2 %3").arg(vec[0]).arg(vec[1]).arg(vec[2]);
}


bool iAModalityList::Load(QString const & filename)
{
	if (filename.isEmpty())
	{
		DEBUG_LOG("No modality file given.");
		return false;
	}
	QFileInfo fi(filename);
	if (!fi.exists())
	{
		DEBUG_LOG(QString("Given modality file '%1' does not exist").arg(filename));
		return false;
	}
	QSettings settings(filename, QSettings::IniFormat );
	
	if (!settings.contains(FileVersionKey) ||
		settings.value(FileVersionKey).toString() != ModFileVersion)
	{
		DEBUG_LOG(QString("Invalid modality file version (was %1, expected %2! Trying to parse anyway, but expect failures.")
			.arg(settings.contains(FileVersionKey) ? settings.value(FileVersionKey).toString() : "not set")
			.arg(ModFileVersion));
		return false;
	}

	bool setSpacingToOne = settings.contains(SetSpacingToOneKey) && settings.value(SetSpacingToOneKey).toBool();
	int currIdx = 0;
	
	double spacingFactor[3] = { 1.0, 1.0, 1.0 };
	while (settings.contains(GetModalityKey(currIdx, "Name")))
	{
		QString modalityName = settings.value(GetModalityKey(currIdx, "Name")).toString();
		QString modalityFile = settings.value(GetModalityKey(currIdx, "File")).toString();
		modalityFile = MakeAbsolute(fi.absolutePath(), modalityFile);
		QString orientationSettings = settings.value(GetModalityKey(currIdx, "Orientation")).toString();
		QString positionSettings = settings.value(GetModalityKey(currIdx, "Position")).toString();
		QString tfFileName = settings.value(GetModalityKey(currIdx, "TransferFunction")).toString();
		tfFileName = MakeAbsolute(fi.absolutePath(), tfFileName);
		if (ModalityExists(modalityFile))
		{
			//DebugOut () << "Modality (name="<<modalityName<<", filename="<<modalityFile<<") already exists!" << std::endl;
		}
		else
		{
			QSharedPointer<iAModality> mod(new iAModality(modalityName, modalityFile));
			if (!mod->LoadData())
			{
				return false;
			}

			// fake a spacing of 1 1 1 for main dataset (to improve transparency renderings)
			if (setSpacingToOne)
			{
				if (currIdx == 0)
				{
					mod->GetImage()->GetSpacing(spacingFactor);
				}
				double spacing[3];
				mod->GetImage()->GetSpacing(spacing);
				mod->GetImage()->SetSpacing(spacing[0] / spacingFactor[0], spacing[1] / spacingFactor[1], spacing[2] / spacingFactor[2]);
			}

			mod->orientationSettings = orientationSettings;
			mod->positionSettings = positionSettings;
			mod->tfFileName = tfFileName;

			m_modalities.push_back(mod);
			emit Added(mod);
		}
		currIdx++;
	}
	m_fileName = filename;
	return true;
}

namespace
{
QString GetMeasurementString(QSharedPointer<iAModality> mod)
{
	return QString("") + QString::number(mod->GetWidth()) + "x" + 
		QString::number(mod->GetHeight()) + "x" +
		QString::number(mod->GetDepth()) + " (" +
		QString::number(mod->GetSpacing()[0]) + ", " +
		QString::number(mod->GetSpacing()[1]) + ", " +
		QString::number(mod->GetSpacing()[2]) + ")";
}
}

void iAModalityList::Add(QSharedPointer<iAModality> mod)
{
	if (m_modalities.size() > 0)
	{
		// make sure that size & spacing fit:
		/*
		if (m_modalities[0]->GetWidth() != mod->GetWidth() ||
			m_modalities[0]->GetHeight() != mod->GetHeight() ||
			m_modalities[0]->GetDepth() != mod->GetDepth() ||
			m_modalities[0]->GetSpacing()[0] != mod->GetSpacing()[0] ||
			m_modalities[0]->GetSpacing()[1] != mod->GetSpacing()[1] ||
			m_modalities[0]->GetSpacing()[2] != mod->GetSpacing()[2])
		{
			DebugOut() << "Measurements of new modality " <<
				GetMeasurementString(mod) << " don't fit measurements of existing one: " <<
				GetMeasurementString(m_modalities[0]) << std::endl;
			return;
		}
		*/
	}
	m_modalities.push_back(mod);
	emit Added(mod);
}

void iAModalityList::Remove(int idx)
{
	m_modalities.remove(idx);
}

QSharedPointer<iAModality> iAModalityList::Get(int idx)
{
	return m_modalities[idx];
}

QSharedPointer<iAModality const> iAModalityList::Get(int idx) const
{
	return m_modalities[idx];
}

int iAModalityList::size() const
{
	return m_modalities.size();
}