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
#include "SegmentationRunner.h"

#include "iAConsole.h"
#include "iAImageCoordinate.h"
#include "iAMMSegParameter.h"
#include "iAModality.h"
#include "iANormalizerImpl.h"
#include "iAPCA.h"
#include "iASpectraDistance.h"
#include "SVMImageFilter.h"

#include "iAConnector.h"
#include <itkGPUGradientAnisotropicDiffusionImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>

#include <vtkImageAccumulate.h>
#include <vtkImageData.h>

SegmentationRunner::SegmentationRunner(
	QSharedPointer<iAModalityList const> modalities,
	SeedsPointer seeds,
	iAMMSegParameter const & parameter,
	RWMode mode
) :
	m_modalities(modalities),
	m_seeds(seeds),
	m_param(parameter),
	m_aborted(false),
	m_mode(mode)
{
}

void SegmentationRunner::run()
{
	try
	{
		if (m_param.doGAD())
		{
			typedef itk::Image< double, 3 > InputImageType;

			for (int iM = 0; iM < m_modalities->size(); iM++) {
				iAConnector* m_conn = new iAConnector;
				m_conn->SetImage(m_modalities->Get(iM)->GetImage());

				typedef itk::GradientAnisotropicDiffusionImageFilter< InputImageType, InputImageType > GGADIFType;
				typename GGADIFType::Pointer filter = GGADIFType::New();

				filter->SetNumberOfIterations(m_param.gad_iter());
				filter->SetTimeStep(m_param.gad_step());
				filter->SetConductanceParameter(m_param.gad_cond());
				filter->SetInput(dynamic_cast<InputImageType *>(m_conn->GetITKImage()));

				filter->Update();
				m_conn->SetImage(filter->GetOutput());
				m_modalities->Get(iM)->SetImage(m_conn->GetVTKImage());
				m_modalities->Get(iM)->GetImage()->Modified();
				filter->ReleaseDataFlagOn();
			}
		}


		if (m_mode == Extended)
		{
			iAPerformanceTimer svmPerfTimer;
			svmPerfTimer.start();

			DEBUG_LOG(QString("SVM (c=%1, gamma=%2, channels=%3) START")
				.arg(m_param.svm_c())
				.arg(m_param.svm_gamma())
				.arg(m_param.svm_channels()));
			SVMImageFilter svm(
				m_param.svm_c(),
				m_param.svm_gamma(),
				m_modalities,
				m_seeds,
				m_param.svm_channels());
			svm.Run();

			double svmTime = svmPerfTimer.elapsed();

			DEBUG_LOG(QString("SVM Finished in %1 seconds.").arg(svmTime));

			DEBUG_LOG(QString("ERW (beta=%1, gamma=%2, max_iter=%3) START")
				.arg(m_param.erw_beta()).arg(m_param.erw_gamma()).arg(m_param.erw_maxIter()));
			for (int i = 0; i < m_param.modalityCount(); ++i)
			{
				DEBUG_LOG(QString("   modality %1 (weight=%2, metric=%3, pcaDim=%4)")
					.arg(i + 1)
					.arg(m_param.modalityWeight(i))
					.arg(m_param.distanceFuncIdx(i))
					.arg(m_param.pcaDim(i)));
			}

			QSharedPointer<iAGaussianNormalizer> gauss(new iAGaussianNormalizer);
			gauss->SetBeta(m_param.erw_beta());

			QSharedPointer<QVector<iARWInputChannel> > inputChannels(new QVector<iARWInputChannel>());
			for (int i = 0; i < m_param.modalityCount(); ++i)
			{
				iARWInputChannel input;
				QSharedPointer<iAModality const> modI = m_modalities->Get(i);
				if (modI->GetData()->channelCount() != m_param.pcaDim(i))
				{
					iAPCA pca(modI->GetData());
					input.image = pca.GetReduced(modI->GetConverter(), m_param.pcaDim(i));
				}
				else
				{
					input.image = modI->GetData();
				}
				input.distanceFunc = m_param.distanceFunc(i);
				input.normalizeFunc = gauss;
				input.weight = m_param.modalityWeight(i);
				inputChannels->push_back(input);
			}
			QSharedPointer<iAModality const> mod0 = m_modalities->Get(0);

			iAExtendedRandomWalker* extendedRandomWalker = new iAExtendedRandomWalker(
				mod0->GetWidth(), mod0->GetHeight(), mod0->GetDepth(), mod0->GetSpacing(),
				inputChannels,
				svm.GetResult(),
				m_param.erw_gamma(),
				m_param.erw_maxIter()
			);

			//connect(extendedRandomWalker, SIGNAL(finished()), this, SLOT(erwFinished()) );
			//extendedRandomWalker->start();

			m_erwPerfTimer.start();
			extendedRandomWalker->run();
			double erwTime = m_erwPerfTimer.elapsed();

			DEBUG_LOG(QString("ERW Finished in %1 seconds!").arg(QString::number(erwTime)));

			m_result = extendedRandomWalker->GetResult();

			auto svmResult = svm.GetResult();

			double maxMinMaxDist = std::numeric_limits<double>::max();
			const double EPSILON = 0.01;
			for (int i = 0; i < svmResult->size(); ++i)
			{
				m_result->svmProbImages.push_back(svmResult->at(i));
				vtkSmartPointer<vtkImageAccumulate> acc = vtkSmartPointer<vtkImageAccumulate>::New();
				acc->SetInputData(svmResult->at(i));
				acc->Update();
				//DEBUG_LOG(QString("Modality %1: Min=%2, Max=%3").arg(i).arg(*acc->GetMin()).arg(*acc->GetMax()));
				if (abs(*acc->GetMax() - *acc->GetMin()) < maxMinMaxDist)
				{
					maxMinMaxDist = *acc->GetMax() - *acc->GetMin();
				}

			}
			if (maxMinMaxDist < EPSILON)
			{
				DEBUG_LOG(QString("WARNING: Maximum Min-Max distance small: %1").arg(maxMinMaxDist));
			}

			delete extendedRandomWalker;
		}

		else
		{
			QSharedPointer<iAGaussianNormalizer> gauss(new iAGaussianNormalizer);
			gauss->SetBeta(m_param.erw_beta());

			QSharedPointer<QVector<iARWInputChannel> > inputChannels(new QVector<iARWInputChannel>());
			for (int i = 0; i < m_param.modalityCount(); ++i)
			{
				iARWInputChannel input;
				QSharedPointer<iAModality const> modI = m_modalities->Get(i);
				if (modI->GetData()->channelCount() != m_param.pcaDim(i))
				{
					iAPCA pca(modI->GetData());
					input.image = pca.GetReduced(modI->GetConverter(), m_param.pcaDim(i));
				}
				else
				{
					input.image = modI->GetData();
				}
				input.distanceFunc = m_param.distanceFunc(i);
				input.normalizeFunc = gauss;
				input.weight = m_param.modalityWeight(i);
				inputChannels->push_back(input);
			}

			if (inputChannels->size() == 0) {
				DEBUG_LOG("Need at least one input channel!");
				return;
			}

			QSharedPointer<iAModality const> mod0 = m_modalities->Get(0);

			iARandomWalker* randomWalker = new iARandomWalker(
				mod0->GetWidth(), mod0->GetHeight(), mod0->GetDepth(), mod0->GetSpacing(),
				inputChannels,
				m_seeds,
				m_param.erw_maxIter()
			);

			//connect(extendedRandomWalker, SIGNAL(finished()), this, SLOT(erwFinished()) );
			//extendedRandomWalker->start();

			m_erwPerfTimer.start();
			randomWalker->run();

			double erwTime = m_erwPerfTimer.elapsed();
			/*
			if (randomWalker->GetError() != "")
			{
				DEBUG_LOG(QString("ERW Error: %1")
					.arg(randomWalker->GetError()));
				m_aborted = true;
				return;
			}
			*/
			DEBUG_LOG(QString("ERW Finished in %1 seconds!").arg(QString::number(erwTime)));
			m_result = randomWalker->GetResult();
			delete randomWalker;
		}
	}
	catch (std::exception& e)
	{
		DEBUG_LOG(QString("An exception has occured: %1").arg(e.what()));
		m_aborted = true;
	}
}

QSharedPointer<iARWResult> SegmentationRunner::GetResult()
{
	return m_result;
}
