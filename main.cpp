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
#include "iAConsole.h"
#include "ParameterHelper.h"
#include "iAMMSegParameter.h"
#include "SeedType.h"
#include "SegmentationRunner.h"
#include "iAModality.h"
#include "iAToolsITK.h"
#include "iAToolsVTK.h"
#include "iASpectraDistanceImpl.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMaskImageFilter.h>

#include <vtkObjectFactory.h>
#include <vtkOutputWindow.h>
#include <vtkSmartPointer.h>
#ifdef WIN32
#include <vtkWindows.h>
#endif

#include <QFileInfo>
#include <QXmlStreamReader>

QSharedPointer<iAModalityList> LoadModalities(char* argv[], int argc, int & currentParam)
{
	QSharedPointer<iAModalityList> modalities(new iAModalityList());

	if (currentParam >= argc)
	{
		throw std::exception(QString("Expected more arguments, only %1 given").arg(argc).toStdString().c_str());
	}
	QString curStr(argv[currentParam]);

	while (!isDouble(curStr))
	{
		QSharedPointer<iAModality> mod(new iAModality(QString("mod%1").arg(currentParam - 1), curStr));
		if (!mod->LoadData())
		{
			DEBUG_LOG(QString("ERROR while loading file '%1'.").arg(curStr));
			throw std::exception(QString("Could not load file %1!").arg(curStr).toStdString().c_str());
		}
		modalities->Add(mod);
		
		currentParam++;
		if (currentParam >= argc)
		{
			throw std::exception(QString("Expected more arguments, only %1 given").arg(argc).toStdString().c_str());
		}
		curStr = argv[currentParam];
	}
	return modalities;
}

RWMode LoadMode(char* argv[], int argc, int & currentParam)
{
	QString modeStr(argv[currentParam++]);
	if (modeStr == "RW")
		return Standard;
	else if (modeStr == "ERW")
		return Extended;
	else if (modeStr == "Norm")
	{
		return Normalize;
	}
	else
	{
		QString msg(QString("Unknown mode %1!").arg(modeStr));
		DEBUG_LOG(msg);
		throw std::exception(msg.toStdString().c_str());
	}
}

iAMMSegParameter LoadParameters(char* argv[], int argc, RWMode mode, bool do_gad, int & currentParam)
{
	double erw_beta = toDoubleAdv(argv, argc, currentParam, "ERW Beta");
	double erw_gamma = 0,
		erw_maxIter = 0,
		svm_c = 0,
		svm_gamma = 0,
		svm_channels = 0;
	int gad_iter = 0;
	double gad_step = 0,
		gad_cond = 0;
	if (mode == Extended)
	{
		erw_gamma = toDoubleAdv(argv, argc, currentParam, "ERW Gamma");
		erw_maxIter = toDoubleAdv(argv, argc, currentParam, "ERW Max Iter");
		svm_c = toDoubleAdv(argv, argc, currentParam, "SVM C");
		svm_gamma = toDoubleAdv(argv, argc, currentParam, "SVM Gamma");
		svm_channels = toDoubleAdv(argv, argc, currentParam, "SVM Channels");
	}
	if (do_gad)
	{
		gad_iter = toIntAdv(argv, argc, currentParam, "GAD Iterations");
		gad_step = toDoubleAdv(argv, argc, currentParam, "GAD Step");
		gad_cond = toDoubleAdv(argv, argc, currentParam, "GAD Conductance");
	}
	return iAMMSegParameter(erw_beta, erw_gamma, erw_maxIter, svm_c, svm_gamma, svm_channels, do_gad, gad_iter, gad_step, gad_cond);
}

void LoadModalityParams(char *argv[], int argc, iAMMSegParameter & param, int modalityCount, int & currentParam)
{
	int * pca = new int[modalityCount];
	double * weight = new double[modalityCount];
	for (int i = 0; i < modalityCount; ++i)
	{
		pca[i] = toIntAdv(argv, argc, currentParam, QString("PCA Mod %1").arg(i).toStdString().c_str());
	}
	for (int i = 0; i < modalityCount; ++i)
	{
		weight[i] = toDoubleAdv(argv, argc, currentParam, QString("Weight Mod %1").arg(i).toStdString().c_str());
	}
	for (int i = 0; i < modalityCount; ++i)
	{
		if (currentParam >= argc)
		{
			throw std::exception(QString("Expected more arguments, only %1 given, param Distance Func Mod %2\n").arg(argc).arg(i).toStdString().c_str());
		}
		QString distFuncShortName(argv[currentParam++]);
		int distFuncIdx = GetShortNameIdx(distFuncShortName);
		if (distFuncIdx == -1)
		{
			throw std::exception(QString("Unrecognized distance function name '%1' given, param Distance Func Mod %2\n").arg(distFuncShortName).arg(i).toStdString().c_str());
		}
		param.addModalityParam(pca[i], weight[i], distFuncIdx);
	}
	delete[] pca;
	delete[] weight;
}

SeedsPointer LoadSeeds(char *argv[], int argc, int & currentParam)
{
	if (currentParam >= argc)
	{
		throw std::exception(QString("Expected more arguments, only %1 given").arg(argc).toStdString().c_str());
	}
	QString filename(argv[currentParam++]);
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		DEBUG_LOG(QString("Seed file loading: Failed to open file '%1'!").arg(filename));
		return SeedsPointer();
	}
	QXmlStreamReader stream(&file);
	stream.readNext();

	SeedsPointer allSeeds(new Seeds);

	QList<Seed> labelSeeds;
	Seed currentSeed;
	QSharedPointer<QVector<double> > seedValues;
	while (!stream.atEnd())
	{
		if (stream.isStartElement())
		{
			if (stream.name() == "Label")
			{
				int id = stream.attributes().value("id").toInt();
				if (allSeeds->size() != id)
				{
					DEBUG_LOG(QString("Invalid seed id %1; note: currently only sequential, consecutive IDs are supported!").arg(id));
				}
			}
			else if (stream.name() == "Seed")
			{
				currentSeed.setCoordinate(
					stream.attributes().value("x").toInt(),
					stream.attributes().value("y").toInt(),
					stream.attributes().value("z").toInt());
				seedValues = QSharedPointer<QVector<double> >(new QVector<double>);
			}
			else if (stream.name() == "Value")
			{
				double value = stream.attributes().value("value").toDouble();
				seedValues->push_back(value);
			}
		}
		if (stream.isEndElement())
		{
			if (stream.name() == "Label")
			{
				if (labelSeeds.empty())
				{
					DEBUG_LOG("Seed file contains a label with no seeds defined!");
				}
				allSeeds->push_back(labelSeeds);
				labelSeeds.clear();
			}
			else if (stream.name() == "Seed")
			{
				if (seedValues->size() > 0)
				{
					currentSeed.setValues(seedValues);
				}
				labelSeeds.push_back(currentSeed);
			}
		}
		stream.readNext();
	}

	if (!labelSeeds.empty())
	{
		DEBUG_LOG("Invalid seed file format: Missing closing label tag?");
	}

	file.close();
	if (stream.hasError())
	{
		DEBUG_LOG(QString("Error: Failed to parse seed xml file '%1': %2")
			.arg(filename)
			.arg(stream.errorString()));
		return SeedsPointer();
	}
	else if (file.error() != QFile::NoError)
	{
		DEBUG_LOG(QString("Error: Cannot read file '%1': %2")
			.arg(filename)
			.arg(file.errorString()));
		return SeedsPointer();
	}
	return allSeeds;
}

class StdOutNoOutputWindow : public vtkOutputWindow
 {
public:
	vtkTypeMacro(StdOutNoOutputWindow, vtkOutputWindow);
	void PrintSelf(ostream& os, vtkIndent indent);
	static StdOutNoOutputWindow * New();
	virtual void DisplayText(const char*);
private:
	StdOutNoOutputWindow();
	StdOutNoOutputWindow(const StdOutNoOutputWindow &) = delete;
	void operator=(const StdOutNoOutputWindow &) =delete;
};

vtkStandardNewMacro(StdOutNoOutputWindow);

StdOutNoOutputWindow::StdOutNoOutputWindow() {}

void StdOutNoOutputWindow::DisplayText(const char* someText)
{
	std::cout << someText;
}

//----------------------------------------------------------------------------
void StdOutNoOutputWindow::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}


bool GetBoolParam(const char * param)
{
	if (QString("true") != param &&  QString("false") != param)
	{
		throw std::exception(QString("Invalid parameter value, expected true or false, got '%1'").arg(param).toStdString().c_str());
	}
	return param == QString("true");
}

template <class T>
class statistics 
{
private:
	std::vector<T> values;		// is there a better way?
	double sum;
public:
	statistics(): sum(0) {}
	~statistics() {}
	void add(T value)
	{
		sum += value;
		values.push_back(value);
	}
	double mean() {
		return sum / values.size();
	}
	double variance() {
		double diff_sum = 0;
		double m = mean();
		for (auto it = values.begin(); it != values.end(); ++it)
		{
			diff_sum += ((*it - m)*(*it - m));
		}
		return diff_sum / values.size();
	}
	double stddev() {
		return sqrt(variance());
	}
};


void NormalizeProbabilities(QVector<ProbabilityImagePointer> probabilities)
{
	int labels = probabilities.size();
	typedef itk::ImageRegionIterator<ProbabilityImageType> IteratorType;
	IteratorType * its = new IteratorType[labels];
	bool atEnd = false;
	for (int l = 0; l < labels; ++l)
	{
		its[l] = IteratorType(probabilities[l], probabilities[l]->GetLargestPossibleRegion());
		its[l].GoToBegin();
		atEnd = its[l].IsAtEnd();
	}
	bool zeroSum = false;
	double factorSum = 0;
	statistics<double> factorStats;
	while (!atEnd)
	{
		double sum = 0;
		for (int l = 0; l < labels; ++l)
		{
			sum += its[l].Get();
		}
		if (sum == 0)
		{
			zeroSum = true;
		}
		else
		{
			double factor = 1.0 / sum;
			factorStats.add(factor);
			for (int l = 0; l < labels; ++l)
			{
				its[l].Set(its[l].Get() * factor);
			}
		}
		for (int l = 0; l < labels; ++l)
		{
			++(its[l]);
			if ((atEnd && !its[l].IsAtEnd()) ||
				(l > 0 && !atEnd && its[l].IsAtEnd()))
			{
				DEBUG_LOG("Image sizes are not equal!");
			}
			atEnd |= its[l].IsAtEnd();
		}
	}
	if (zeroSum)
	{
		DEBUG_LOG("One or more pixels had probability sum of 0!");
	}
	std::cout << "Factor: mean=" << factorStats.mean() << ", stddev=" << factorStats.stddev() << std::endl;
	delete[] its;
}

#include <QDir>

int main(int argc, char* argv[])
{
	std::cout << "svmerw v1.0" << std::endl;

	vtkSmartPointer<StdOutNoOutputWindow> myOutputWindow = vtkSmartPointer<StdOutNoOutputWindow>::New();
	vtkOutputWindow::SetInstance(myOutputWindow);

	int currentParam = 1;
	
	try
	{
		// retrieve arguments

		// "fixed":
		RWMode mode = LoadMode(argv, argc, currentParam);
		if (mode == Normalize)
		{
			QString inputDirectory = argv[currentParam++];
			QDir dir(inputDirectory);
			std::cout << "Normalizing files in directory " << QFileInfo(inputDirectory).absoluteFilePath().toStdString() << std::endl;
			dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
			QFileInfoList list = dir.entryInfoList();
			QVector<ProbabilityImagePointer> probImages;
			QVector<QString> filenames;
			for (int i = 0; i < list.size(); ++i) {
				QFileInfo fileInfo = list.at(i);
				if (fileInfo.fileName().startsWith("prob") && (fileInfo.fileName().endsWith(".mhd")))
				{
					//std::cout << "  Reading image " << fileInfo.fileName().toStdString() << std::endl;
					ProbabilityImagePointer probImg = ReadImage<ProbabilityImageType>(fileInfo.absoluteFilePath());
					probImages.push_back(probImg);
					filenames.push_back(fileInfo.absoluteFilePath());
				}
			}
			NormalizeProbabilities(probImages);
			for (int i = 0; i<probImages.size(); ++i)
			{
				//std::cout << "  Writing image " << filenames[i].toStdString() << std::endl;
				StoreImage<ProbabilityImageType>(probImages[i], filenames[i]);
			}
			return 0;
		}
		SeedsPointer seeds = LoadSeeds(argv, argc, currentParam);
		if (!seeds)
		{
			return 1;
		}
		bool storeSVMProbabilities = false;
		if (mode == Extended)
		{
			storeSVMProbabilities = GetBoolParam(argv[currentParam++]);
		}
		bool storeERWProbabilities = GetBoolParam(argv[currentParam++]);
		bool do_gad = GetBoolParam(argv[currentParam++]);
		bool use_mask = GetBoolParam(argv[currentParam++]);

		LabelImagePointer maskImage;
		if (use_mask)
		{
			maskImage = ReadImage<LabelImageType>(argv[currentParam++]);
		}

		// "output-dir"
		if (currentParam >= argc)
		{
			throw std::exception(QString("Expected more arguments, only %1 given").arg(argc).toStdString().c_str());
		}
		QString outputFile(argv[currentParam++]);

		// modalities
		QSharedPointer<iAModalityList> modalities = LoadModalities(argv, argc, currentParam);
		int modalityCount = modalities->size();

		// parameters
		iAMMSegParameter param = LoadParameters(argv, argc, mode, do_gad, currentParam);
		LoadModalityParams(argv, argc, param, modalityCount, currentParam);

		int labelCount = seeds->size();

		// run segmentation
		DEBUG_LOG(QString("Starting segmentation for %1 modalities into %2 labels").arg(modalityCount).arg(labelCount));
		SegmentationRunner segmRunner(
			modalities,
			seeds,
			param,
			mode
		);
		segmRunner.run();

		QSharedPointer<iARWResult> result = segmRunner.GetResult();

		QString labelFileName(outputFile);
		LabelImagePointer img = result->labelledImage;
		if (use_mask)
		{
			typedef itk::MaskImageFilter< LabelImageType, LabelImageType > MaskFilterType;
			MaskFilterType::Pointer maskFilter = MaskFilterType::New();
			maskFilter->SetInput(img);
			maskFilter->SetMaskImage(maskImage);
			maskFilter->Update();
			img = maskFilter->GetOutput();
		}
		StoreImage<LabelImageType>(img, labelFileName);
		QFileInfo outputFileInfo(outputFile);
		QString outputDir(outputFileInfo.absolutePath());
		for (int i = 0; storeSVMProbabilities && i<result->svmProbImages.size(); ++i)
		{
			QString probFileName(outputDir + "\\svm-prob" + QString::number(i) + ".mhd");
			StoreImage(result->svmProbImages[i], probFileName, true);
		}
		NormalizeProbabilities(result->erwProbImages);
		for (int i = 0; storeERWProbabilities && i<result->erwProbImages.size(); ++i)
		{
			QString probFileName(outputDir + "\\prob" + QString::number(i) + ".mhd");
			StoreImage<ProbabilityImageType>(result->erwProbImages[i], probFileName);
		}
	}
	catch (std::exception & e)
	{
		DEBUG_LOG(QString("Segmentation failed: %1").arg(e.what()));
		return 1;
	}
}

