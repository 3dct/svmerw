/************************************  SVMERW  ****************************************
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
**************************************************************************************/
#include "iAToolsITK.h"

#include "iAITKIO.h"
#include "iATypedCallHelper.h"


itk::ImageIOBase::IOComponentType GetITKScalarPixelType(ImagePointer image)
{
	itk::ImageIOBase::IOComponentType result = itk::ImageIOBase::UNKNOWNCOMPONENTTYPE;
	ImageBaseType * imagePtr = image.GetPointer();

	if (dynamic_cast<itk::Image< unsigned char, m_DIM> *>(imagePtr))
		result = itk::ImageIOBase::UCHAR;
	else if (dynamic_cast<itk::Image< char, m_DIM> *>(imagePtr))
		result = itk::ImageIOBase::CHAR;
	else if (dynamic_cast<itk::Image< short, m_DIM> *>(imagePtr))
		result = itk::ImageIOBase::SHORT;
	else if (dynamic_cast<itk::Image< unsigned short, m_DIM> *>(imagePtr))
		result = itk::ImageIOBase::USHORT;
	else if (dynamic_cast<itk::Image< int, m_DIM> *>(imagePtr))
		result = itk::ImageIOBase::INT;
	else if (dynamic_cast<itk::Image< unsigned int, m_DIM> *>(imagePtr))
		result = itk::ImageIOBase::UINT;
	else if (dynamic_cast<itk::Image< long, m_DIM> *>(imagePtr))
		result = itk::ImageIOBase::LONG;
	else if (dynamic_cast<itk::Image< unsigned long, m_DIM> *>(imagePtr))
		result = itk::ImageIOBase::ULONG;
	else if (dynamic_cast<itk::Image< float, m_DIM> *>(imagePtr))
		result = itk::ImageIOBase::FLOAT;
	else if (dynamic_cast<itk::Image< double, m_DIM> *>(imagePtr))
		result = itk::ImageIOBase::DOUBLE;

	return result;
}


template <class T>
void alloc_image_tmpl(ImagePointer otherImg, ImagePointer & result)
{
	typedef itk::Image<T, m_DIM > ImageType;
	typedef typename ImageType::Pointer ImagePointer;

	ImagePointer image = ImageType::New();
	typename ImageType::RegionType reg(
		otherImg->GetLargestPossibleRegion().GetIndex(),
		otherImg->GetLargestPossibleRegion().GetSize()
	);
	image->SetRegions(reg);
	image->Allocate();
	image->FillBuffer(0);
	image->SetSpacing(otherImg->GetSpacing());
	result = image;
}


ImagePointer AllocateImage(ImagePointer img)
{
	ImagePointer result;
	ITK_TYPED_CALL(alloc_image_tmpl, GetITKScalarPixelType(img), img, result);
	return result;
}


void StoreImage(ImagePointer image, QString const & filename, bool useCompression)
{
	iAITKIO::writeFile(filename, image, GetITKScalarPixelType(image), useCompression);
}

template <class TImage>
void GetITKPixel2(double & result, TImage* image, typename TImage::IndexType idx)
{
	result = image->GetPixel(idx);
}

template <class T>
void GetITKPixel(double & result, ImagePointer img, ImageBaseType::IndexType idx)
{
	typedef itk::Image<T, m_DIM > ImageType;
	GetITKPixel2(result, dynamic_cast<ImageType*>(img.GetPointer()), idx);
}


double GetITKPixel(ImagePointer img, ImageBaseType::IndexType idx)
{
	double result;
	ITK_TYPED_CALL(GetITKPixel, GetITKScalarPixelType(img), result, img,  idx);
	return result;
}

template <class TImage>
void SetITKPixel2(double value, TImage* image, typename TImage::IndexType idx)
{
	image->SetPixel(idx, value);
}

template <class T>
void SetITKPixel(double value, ImagePointer img, ImageBaseType::IndexType idx)
{
	typedef itk::Image<T, m_DIM > ImageType;
	SetITKPixel2(value, dynamic_cast<ImageType*>(img.GetPointer()), idx);
}


void SetITKPixel(ImagePointer img, ImageBaseType::IndexType idx, double value)
{
	ITK_TYPED_CALL(SetITKPixel, GetITKScalarPixelType(img), value, img, idx);
}
