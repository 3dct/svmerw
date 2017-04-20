/************************************  SVMERW  ************************************** *
* **  Command line tool for segmenting multi-modal 3D data with a Support Vector   ** *
* **  Machine (SVM) and the Extended Random Walker (ERW)                           ** *
* *********************************************************************************** *
* Copyright (C) 2017 Bernhard Fröhler                                                 *
* *********************************************************************************** *
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
* *********************************************************************************** *
* Contact: FH OÖ Forschungs & Entwicklungs GmbH, Campus Wels, CT-Gruppe,              *
*       Stelzhamerstraße 23, 4600 Wels / Austria, Email: bernhard.froehler@fh-wels.at *
* ************************************************************************************/
#include "iARandomWalker.h"

#include "iAImageGraph.h"
#include "iAGraphWeights.h"
#include "iASpectraDistance.h"
#include "iAMathUtility.h"
#include "iAConsole.h"

#include <itkImage.h>

#include <vtkImageData.h>

#ifdef USE_EIGEN

#include <Eigen/Core>
#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> MatrixType;
typedef Eigen::VectorXd VectorType;

#else // using VNL

#include <vnl/vnl_vector.h>
#include <vnl/vnl_sparse_matrix.h>
#include <vnl/vnl_sparse_matrix_linear_system.h>
#include <vnl/algo/vnl_sparse_lu.h>
#include <vnl/algo/vnl_lsqr.h>
#include <vnl/algo/vnl_conjugate_gradient.h>
//<vnl/vnl_conjugate_gradient.h>
typedef vnl_sparse_matrix<double> MatrixType;
typedef vnl_vector<double> VectorType;

#endif

#include <QMap>
#include <QSet>

namespace
{
	
	// TODO: TEST!
	/*
	bool SeedsContainIdx(SeedsPointer seeds,
		int vertexIdx, iAImageCoordConverter const & conv)
	{
		for (int label = 0; label < seeds->size(); ++label)
		{
			for (iAImageCoordinate coord : seeds->at(label))
			{
				if (conv.GetIndexFromCoordinates(coord) == vertexIdx)
				{
					return true;
				}
			}
		}
		return false;
	}
	*/

	typedef QMap<iAVertexIndexType, iAVertexIndexType> IndexMap;

	void CreateLaplacianPart( MatrixType & output,
		IndexMap const & rowIndices,
		IndexMap const & colIndices,
		QSharedPointer<iAImageGraph> imageGraph,
		QSharedPointer<iAGraphWeights const> finalWeight,
		QVector<double> const & vertexWeightSum,
		int vertexCount,
		bool noSameIndices = false)
	{
#ifdef USE_EIGEN
		output.reserve(Eigen::VectorXi::Constant(rowIndices.size(), 7) );
#endif
		// edge weights:
		for (iAEdgeIndexType edgeIdx = 0; edgeIdx < imageGraph->GetEdgeCount(); ++edgeIdx)
		{
			iAEdgeType const & edge = imageGraph->GetEdge(edgeIdx);
			if (rowIndices.contains(edge.first) && colIndices.contains(edge.second))
			{
				iAVertexIndexType newRowIdx = rowIndices[edge.first];
				iAVertexIndexType newColIdx = colIndices[edge.second];
#ifdef USE_EIGEN
				output.insert(newRowIdx, newColIdx) = -finalWeight->GetWeight(edgeIdx);
#else
				output(newRowIdx, newColIdx) = -finalWeight->GetWeight(edgeIdx);
#endif
			}
			if (rowIndices.contains(edge.second) && colIndices.contains(edge.first))
			{
				iAVertexIndexType newRowIdx = rowIndices[edge.second];
				iAVertexIndexType newColIdx = colIndices[edge.first];
#ifdef USE_EIGEN
				output.insert(newRowIdx, newColIdx) = -finalWeight->GetWeight(edgeIdx);
#else
				output(newRowIdx, newColIdx) = -finalWeight->GetWeight(edgeIdx);
#endif
			}
		}
		if (noSameIndices)
		{ // optimization: if rowIndices and colIndices don't share any elements, there is no diagonal entry
			return;
		}
		// sum of incident edge weights at diagonal:
		for (iAVertexIndexType vertexIdx = 0; vertexIdx < vertexCount; ++vertexIdx)
		{
			if (rowIndices.contains(vertexIdx) && colIndices.contains(vertexIdx))
			{
				iAVertexIndexType newRowIdx = rowIndices[vertexIdx];
				iAVertexIndexType newColIdx = colIndices[vertexIdx];
#ifdef USE_EIGEN
				output.insert(newRowIdx, newColIdx) = vertexWeightSum[vertexIdx];
#else
				output(newRowIdx, newColIdx) = vertexWeightSum[vertexIdx];
#endif
			}
		}
	}

	void SetIndexMapValues(ProbabilityImagePointer image,
		VectorType const & values,
		IndexMap const & indexMap,
		iAImageCoordConverter const & conv)
	{
		for (IndexMap::const_iterator it = indexMap.begin(); it != indexMap.end(); ++it)
		{
#ifdef USE_EIGEN
			double imgVal = values.coeff(it.value());
#else
			double imgVal = values[it.value()];
#endif
			iAImageCoordinate coord = conv.GetCoordinatesFromIndex(it.key());
			ProbabilityImageType::IndexType pixelIndex;
			pixelIndex[0] = coord.x;
			pixelIndex[1] = coord.y;
			pixelIndex[2] = coord.z;
			if (imgVal < 0 || imgVal > 1 || std::isinf(imgVal) ||std::isnan(imgVal))
			{
				/*
				DebugOut() << "Invalid pixel value at ("
					<<   "x=" << pixelIndex[0]
					<< ", y=" << pixelIndex[1]
					<< ", z=" << pixelIndex[2]
					<< "):  " << imgVal
					<< std::endl;
				*/
				imgVal = 0;
			}
			image->SetPixel(pixelIndex,	imgVal);
		}
	}

	ProbabilityImagePointer CreateProbabilityImage(
		iAImageCoordConverter const & conv,
		double const spacing[3])
	{
		ProbabilityImageType::IndexType start;
		start[0] = start[1] = start[2] = 0;
 
		ProbabilityImageType::SizeType size;
		size[0] = conv.GetWidth();
		size[1] = conv.GetHeight();
		size[2] = conv.GetDepth();
		
		ProbabilityImageType::RegionType region;
		region.SetSize(size);
		region.SetIndex(start);
 
		ProbabilityImagePointer image = ProbabilityImageType::New();
		image->SetRegions(region);
		image->Allocate();

		image->SetSpacing(spacing);
		return image;
	}

	LabelImagePointer CreateLabelImage(
		iAImageCoordConverter const & conv,
		double const  spacing[3])
	{
		LabelImageType::IndexType start;
		start[0] = start[1] = start[2] = 0;
 
		LabelImageType::SizeType size;
		size[0] = conv.GetWidth();
		size[1] = conv.GetHeight();
		size[2] = conv.GetDepth();
		
		LabelImageType::RegionType region;
		region.SetSize(size);
		region.SetIndex(start);
 
		LabelImagePointer image = LabelImageType::New();
		image->SetRegions(region);
		image->Allocate();
		image->SetSpacing(spacing);
		return image;
	}
	
	LabelImagePointer CreateLabelImage(
		iAImageCoordConverter const & conv,
		double const spacing[3],
		QVector<ProbabilityImagePointer> const & probabilityImages,
		int labelCount)
	{
		// create labelled image (as value at k = arg l max(p_l^k) for each pixel k)
		LabelImagePointer pImg = CreateLabelImage(conv, spacing);
		for (int x=0; x<conv.GetWidth(); ++x)
		{
			for (int y=0; y<conv.GetHeight(); ++y)
			{
				for (int z=0; z<conv.GetDepth(); ++z)
				{
					double maxProb =0;
					int maxProbLabel = -1;
					ProbabilityImageType::IndexType idx;
					idx[0] = x;
					idx[1] = y;
					idx[2] = z;
					for (int l=0; l<labelCount; ++l)
					{
						double prob = probabilityImages[l]->GetPixel(idx);
						if (prob >= maxProb)
						{
							maxProb = prob;
							maxProbLabel = l;
						}
					}
					pImg->SetPixel(idx, maxProbLabel );
				}
			}
		}
		return pImg;
	}
}

iARandomWalker::iARandomWalker(
		iAVoxelIndexType width,
		iAVoxelIndexType height,
		iAVoxelIndexType depth,
		double const spacing[3],
		QSharedPointer<QVector<iARWInputChannel> > inputChannels,
		SeedsPointer seeds,
		int maxIter
):
	m_imageGraph(new iAImageGraph(width, height, depth, iAImageCoordinate::ColRowDepMajor)),
	m_inputChannels(inputChannels),
	m_vertexCount(width*height*depth),
	m_minLabel(std::numeric_limits<iALabelType>::max()),
	m_maxLabel(std::numeric_limits<iALabelType>::lowest()),
	m_seeds(seeds),
	m_maxIterations(maxIter)
{
	for(int i=0; i<3; ++i) m_spacing[i] = spacing[i];
	assert((*m_inputChannels)[0].image->size() == width*height*depth);
	m_minLabel = 0;
	m_maxLabel = seeds->size() - 1;
}

void iARandomWalker::run()
{
	if (m_inputChannels->size() == 0)
	{
		DEBUG_LOG("Input Channels must not be empty!");
		return;
	}
	if (m_seeds->size() == 0)
	{
		DEBUG_LOG("Seeds must not be empty!");
		return;
	}

	if (m_minLabel != 0)
	{
		DEBUG_LOG("Labels must start at 0");
		return;
	}

	// iATimeGuard perf("RW Start");

	IndexMap seedMap;
	QSet<int> labelSet;
	int seedIdx = 0;
	for (int label=0; label < m_seeds->size(); ++label)
	{
		for(Seed seed : m_seeds->at(label))
		{
			seedMap.insert(m_imageGraph->GetConverter().GetIndexFromCoordinates(seed.coord()), seedIdx);
			labelSet.insert(label);
			seedIdx++;
		}
	}
	int labelCount = labelSet.size();
	if (m_maxLabel != labelCount-1)
	{
		DEBUG_LOG("Labels must be consecutive from 0 .. maxLabel!");
		return;
	}
	//perf.time("RW: Labels set up");

	QVector<QSharedPointer<iAGraphWeights> > graphWeights(m_inputChannels->size());
	QVector<double> weightsForChannels(m_inputChannels->size());
	for (int i=0; i<m_inputChannels->size(); ++i)
	{
		QSharedPointer<iAGraphWeights> currentMeasureWeight =
			CalculateGraphWeights(*m_imageGraph, *(*m_inputChannels)[i].image, *(*m_inputChannels)[i].distanceFunc);
		graphWeights[i] = currentMeasureWeight;
		weightsForChannels[i] = (*m_inputChannels)[i].weight;
	}
	//perf.time("RW: weights calculated");

	for (int i=0; i<m_inputChannels->size(); ++i)
	{
		graphWeights[i]->Normalize((*m_inputChannels)[i].normalizeFunc);
	}
	// perf.time("RW: weights normalized");

	QSharedPointer<const iAGraphWeights > finalWeight =
		CombineGraphWeights(graphWeights, weightsForChannels);
	// perf.time("RW: weights combined");

	QVector<double> vertexWeightSum(m_vertexCount);
	for (iAEdgeIndexType edgeIdx = 0; edgeIdx < m_imageGraph->GetEdgeCount(); ++edgeIdx)
	{
		iAEdgeType const & edge = m_imageGraph->GetEdge(edgeIdx);
		vertexWeightSum[edge.first]  += finalWeight->GetWeight(edgeIdx);
		vertexWeightSum[edge.second] += finalWeight->GetWeight(edgeIdx);
	}
	// perf.time("RW: vertex weight sum");

	IndexMap unlabeledMap;
	for(iAVertexIndexType vertexIdx=0, newIdx=0;
		vertexIdx < m_vertexCount; ++vertexIdx)
	{
		if(!seedMap.contains(vertexIdx)) {
			unlabeledMap.insert(vertexIdx, newIdx);
			++newIdx;
		}
	}
	int seedCount  = seedMap.size();
	
	MatrixType A(m_vertexCount-seedCount,  m_vertexCount-seedCount);
	CreateLaplacianPart(A, unlabeledMap, unlabeledMap, m_imageGraph, finalWeight, vertexWeightSum, m_vertexCount);
#ifdef USE_EIGEN
	A.makeCompressed();
#endif
	// perf.time("RW: Laplacian generated");
	
	MatrixType BT(m_vertexCount-seedCount, seedCount);
	CreateLaplacianPart(BT, unlabeledMap, seedMap, m_imageGraph, finalWeight, vertexWeightSum, m_vertexCount, true);
	BT = -BT;
#ifdef USE_EIGEN
	BT.makeCompressed();
#endif
	// perf.time("RW: Other Laplacian generated");

	// Conjugate Gradient: very fast and small memory usage!
	// Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::ColMajor> > solver;

#ifdef USE_EIGEN
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > solver;
	solver.analyzePattern(A);
	//solver.setMaxIterations(m_maxIterations);
	std::string error = solver.lastErrorMessage();
	if (error != "")
	{
		DEBUG_LOG(QString(error.c_str()));
		return;
	}
	solver.factorize(A);
	error = solver.lastErrorMessage();
	if (error != "")
	{
		DEBUG_LOG(QString(error.c_str()));
		return;
	}
#else
	vnl_sparse_lu linear_solver(A, vnl_sparse_lu::quiet);
#endif

	// BiCGSTAB: uses a bit more memory, but is a bit faster!
	// Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::ColMajor> > solver;
	// solver.compute(A);
	// perf.time("RW: Solver compute executed");

	m_result = QSharedPointer<iARWResult>(new iARWResult);
	for (int i=0; i<labelCount; ++i)
	{
		VectorType boundary(seedCount);
		int seedIdx = 0;
		for (int label = 0; label < m_seeds->size(); ++label)
		{
			for (Seed coord : m_seeds->at(label))
			{
				boundary[seedIdx] = label == i;
				seedIdx++;
			}
		}
		VectorType b(m_vertexCount-seedCount);
#ifdef USE_EIGEN
		b = BT * boundary;
#else
		BT.mult(boundary, b);
#endif
		VectorType x(m_vertexCount-seedCount);
#ifdef USE_EIGEN

		x = solver.solve(b);
#else
		linear_solver.solve(b,&x);
#endif
		// perf.time("RW: part solved");

		// put values into probability image
		ProbabilityImagePointer pImg = CreateProbabilityImage(m_imageGraph->GetConverter(), m_spacing);
		SetIndexMapValues(pImg, x, unlabeledMap, m_imageGraph->GetConverter());
		SetIndexMapValues(pImg, boundary, seedMap, m_imageGraph->GetConverter());

		m_result->erwProbImages.push_back(pImg);
	}
	m_result->labelledImage = CreateLabelImage(m_imageGraph->GetConverter(), m_spacing, m_result->erwProbImages, labelCount);
}


iAExtendedRandomWalker::iAExtendedRandomWalker(
	iAVoxelIndexType width,
	iAVoxelIndexType height,
	iAVoxelIndexType depth,
	double const spacing[3],
	QSharedPointer<QVector<iARWInputChannel> > inputChannels,
	QSharedPointer<QVector<PriorModelImagePointer> > priorModel,
	double priorModelWeight,
	int maxIterations
):
	m_imageGraph(new iAImageGraph(width, height, depth, iAImageCoordinate::ColRowDepMajor)),
	m_inputChannels(inputChannels),
	m_vertexCount(width*height*depth),
	m_labelCount(0),
	m_priorModel(priorModel),
	m_priorModelWeight(priorModelWeight),
	m_maxIterations(maxIterations)
{
	for(int i=0; i<3; ++i) m_spacing[i] = spacing[i];
	int inputChannelSize = (*m_inputChannels)[0].image->size();
	assert(inputChannelSize == width*height*depth);
	m_labelCount = m_priorModel->size();
	assert(m_labelCount > 0);
	assert(m_priorModelWeight != 0);
	/*
	PriorModelImageType::RegionType reg = (*m_priorModel)[0]->GetLargestPossibleRegion();
	int priorWidth = reg.GetSize()[0];
	int priorHeight = reg.GetSize()[1];
	int priorDepth = reg.GetSize()[2];
	*/
	int * extent = (*m_priorModel)[0]->GetExtent();
	int priorWidth = extent[1] - extent[0] + 1;
	int priorHeight = extent[3] - extent[2] + 1;
	int priorDepth = extent[5] - extent[4] + 1;
	assert(width == priorWidth && height == priorHeight && depth == priorDepth);
}


iAExtendedRandomWalker::iAExtendedRandomWalker(
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
):
	m_imageGraph(new iAImageGraph(width, height, depth, iAImageCoordinate::ColRowDepMajor)),
	m_inputChannels(inputChannels),
	m_vertexCount(width*height*depth),
	m_labelCount(0),
	m_priorModel(priorModel),
	m_priorModelWeight(priorModelWeight),
	m_maxIterations(maxIterations)
{
	for(int i=0; i<3; ++i) m_spacing[i] = spacing[i];
	assert((*m_inputChannels)[0].image->size() == width*height*depth);
	m_labelCount = m_priorModel->size();
	assert(m_labelCount > 0);
	assert(m_priorModelWeight != 0);
	/*
	PriorModelImageType::RegionType reg = (*m_priorModel)[0]->GetLargestPossibleRegion();
	int priorWidth = reg.GetSize()[0];
	int priorHeight = reg.GetSize()[1];
	int priorDepth = reg.GetSize()[2];
	*/
	int * extent = (*m_priorModel)[0]->GetExtent();
	int priorWidth = extent[1] - extent[0] + 1;
	int priorHeight = extent[3] - extent[2] + 1;
	int priorDepth = extent[5] - extent[4] + 1;
	assert(width == priorWidth && height == priorHeight && depth == priorDepth);
}

const double EPSILON = 1e-6;
		
void iAExtendedRandomWalker::run()
{
	QVector<QSharedPointer<iAGraphWeights> > graphWeights(m_inputChannels->size());
	QVector<double> weightsForChannels(m_inputChannels->size());
	for (int i=0; i<m_inputChannels->size(); ++i)
	{
		QSharedPointer<iAGraphWeights> currentMeasureWeight =
			CalculateGraphWeights(*m_imageGraph, *(*m_inputChannels)[i].image, *(*m_inputChannels)[i].distanceFunc);
		graphWeights[i] = currentMeasureWeight;
		weightsForChannels[i] = (*m_inputChannels)[i].weight;
	}
	// perf.time("ERW: weights calculated");

	for (int i=0; i<m_inputChannels->size(); ++i)
	{
		graphWeights[i]->Normalize((*m_inputChannels)[i].normalizeFunc);
	}
	// perf.time("ERW: weights normalized");

	QSharedPointer<iAGraphWeights const> finalWeight =
		CombineGraphWeights(graphWeights, weightsForChannels);
	// perf.time("ERW: weights combined");

	QVector<double> vertexWeightSum(m_vertexCount);
	for (iAEdgeIndexType edgeIdx = 0; edgeIdx < m_imageGraph->GetEdgeCount(); ++edgeIdx)
	{
		iAEdgeType const & edge = m_imageGraph->GetEdge(edgeIdx);
		vertexWeightSum[edge.first]  += finalWeight->GetWeight(edgeIdx);
		vertexWeightSum[edge.second] += finalWeight->GetWeight(edgeIdx);
	}
	// perf.time("ERW: vertex weight sums");
	
	//bool priorNormalized = true;
	// add priors into vertexWeightSum:
	// if my thinking is correct it should be enough to add the weight factor to each entry,
	// since for one voxel, the probabilities for all labels should add up to 1!
	for (iAVoxelIndexType voxelIdx = 0; voxelIdx < m_vertexCount; ++ voxelIdx)
	{
		double sum = 0;
		
		//PriorModelImageType::IndexType idx;
		iAImageCoordinate coord = m_imageGraph->GetConverter().GetCoordinatesFromIndex(voxelIdx);
		/*
		idx[0] = coord.x;
		idx[1] = coord.y;
		idx[2] = coord.z;
		*/
		for (int labelIdx = 0; labelIdx < m_priorModel->size(); ++labelIdx)
		{
			//sum += (*m_priorModel)[labelIdx]->GetPixel(idx);
			sum += (*m_priorModel)[labelIdx]->GetScalarComponentAsDouble(coord.x, coord.y, coord.z, 0);
		}
		assert (std::abs(sum-1.0) < EPSILON);
		if (std::abs(sum-1.0) >= EPSILON)
		{
			//priorNormalized = false;
			DEBUG_LOG(QString("Prior Model not normalized at (x=%1, y=%2, z=%3): %4").
				arg(coord.x).arg(coord.y).arg(coord.z).arg(sum));
		}
		vertexWeightSum[voxelIdx] += (m_priorModelWeight * sum);
	}
	//if (!priorNormalized)
	//{
	//	DebugOut() << "Prior Model not normalized." << std::endl;
	//}
	// perf.time("ERW: priors added to vertex weight");

	IndexMap fullMap;
	for(iAVertexIndexType vertexIdx=0, newIdx=0;
		vertexIdx < m_vertexCount; ++vertexIdx)
	{
		fullMap.insert(vertexIdx, vertexIdx);
	}

	MatrixType A(m_vertexCount,  m_vertexCount);
	CreateLaplacianPart(A, fullMap, fullMap, m_imageGraph, finalWeight, vertexWeightSum, m_vertexCount);
#ifdef USE_EIGEN
	A.makeCompressed();

	// perf.time("ERW: Laplacian generated");

	// Sparse LU: very slow!
	/*
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > solver;
	solver.analyzePattern(A);
	solver.factorize(A);
	*/

	// Conjugate Gradient: very fast and small memory usage!
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::ColMajor> > solver;
	// TODO:
	//    - check if matrix is positive-definite
	//    - find a good maximum number of iterations!
	solver.setMaxIterations(m_maxIterations);
	solver.compute(A);

	// BiCGSTAB: seems to be a bit faster than Conjugate Gradient, but using more memory
	/*
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::ColMajor> > solver;
	solver.setMaxIterations(m_maxIterations);
	solver.compute(A);
	*/

	// std::string error = solver.lastErrorMessage();
	// TODO: use this somewhere?

//#else
	// in case sparse LU should be used, do this:
//	vnl_sparse_lu linear_solver(A, vnl_sparse_lu::quiet);
#endif

	// perf.time("ERW: solver compute done");

	m_result = QSharedPointer<iARWResult>(new iARWResult);

	for (int i=0; i<m_labelCount; ++i)
	{
		VectorType priorForLabel(m_vertexCount);
		// fill from image
		for (iAVoxelIndexType voxelIdx = 0; voxelIdx < m_vertexCount; ++ voxelIdx)
		{
			//PriorModelImageType::IndexType idx;
			iAImageCoordinate coord = m_imageGraph->GetConverter().GetCoordinatesFromIndex(voxelIdx);
			/*
			idx[0] = coord.x;
			idx[1] = coord.y;
			idx[2] = coord.z;
			*/
			priorForLabel[voxelIdx] = (*m_priorModel)[i]->GetScalarComponentAsDouble(coord.x, coord.y, coord.z, 0);
		}

		VectorType x(m_vertexCount);
#ifdef USE_EIGEN
		x = solver.solve(priorForLabel);
#else
		vnl_sparse_matrix_linear_system<double> problem(A, priorForLabel);
		vnl_lsqr solver(problem);
		int returnCode = solver.minimize(x);
		// in case sparse LU should be used, do this instead:
		//linear_solver.solve(priorForLabel, &x);
#endif
		// perf.time("ERW: part solved");
		// put values into probability image
		ProbabilityImagePointer pImg = CreateProbabilityImage(m_imageGraph->GetConverter(), m_spacing);
		SetIndexMapValues(pImg, x, fullMap, m_imageGraph->GetConverter());
		m_result->erwProbImages.push_back(pImg);
	}
	// create labelled image (as value at k = arg l max(p_l^k) for each pixel k)
	m_result->labelledImage = CreateLabelImage(m_imageGraph->GetConverter(), m_spacing, m_result->erwProbImages, m_labelCount);
}

QSharedPointer<iARWResult> iAExtendedRandomWalker::GetResult()
{
	return m_result;
}

QSharedPointer<iARWResult> iARandomWalker::GetResult()
{
	return m_result;
}
