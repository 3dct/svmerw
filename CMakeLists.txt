#*************************************  SVMERW  ****************************************
#* **  Command line tool for segmenting multi-modal 3D data with a Support Vector   ** *
#* **  Machine (SVM) and the Extended Random Walker (ERW)                           ** *
#***************************************************************************************
#* Copyright (C) 2017 Bernhard Fröhler                                                 *
##***************************************************************************************
#* This program is free software: you can redistribute it and/or modify it under the   *
#* terms of the GNU General Public License as published by the Free Software           *
#* Foundation, either version 3 of the License, or (at your option) any later version. *
#*                                                                                     *
#* This program is distributed in the hope that it will be useful, but WITHOUT ANY     *
#* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A     *
#* PARTICULAR PURPOSE.  See the GNU General Public License for more details.           *
#*                                                                                     *
#* You should have received a copy of the GNU General Public License along with this   *
#* program.  If not, see http://www.gnu.org/licenses/                                  *
#***************************************************************************************
#* Contact: FH OÖ Forschungs & Entwicklungs GmbH, Campus Wels, CT-Gruppe,              *
#*       Stelzhamerstraße 23, 4600 Wels / Austria, Email: bernhard.froehler@fh-wels.at *
#***************************************************************************************
cmake_minimum_required(VERSION 3.5.0)

SET (PROJECT_NAME "svmerw")

PROJECT(${PROJECT_NAME})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})
ADD_DEFINITIONS (-DITK4)
SET( ITK_LIBRARIES
	ITKCommon
	ITKIOImageBase
	ITKIOBioRad		ITKIOBMP		ITKIOGDCM		ITKIOGE
	ITKIOGIPL		ITKIOHDF5		ITKIOIPL		ITKIOJPEG
	ITKIOLSM		ITKIOMeta		ITKIOMRC		ITKIONIFTI
	ITKIONRRD		ITKIOPNG		ITKIOSiemens	ITKIOSpatialObjects
	ITKIOStimulate	ITKIOTIFF		ITKIOVTK
	ITKVTK
)
IF (SCIFIO_LOADED)
	SET(ITK_LIBRARIES ${ITK_LIBRARIES} SCIFIO)
ENDIF(SCIFIO_LOADED)

SET(CMAKE_AUTOMOC ON)
FIND_PACKAGE(Qt5 COMPONENTS Core REQUIRED)

FIND_PACKAGE(VTK REQUIRED)
#SET (VTK_LIBRARIES
	#vtkImagingStatistics
#)
#include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

SET (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
FIND_PACKAGE(Eigen3)
IF(EIGEN3_FOUND)
	ADD_DEFINITIONS( -DUSE_EIGEN )
	INCLUDE_DIRECTORIES( ${EIGEN3_INCLUDE_DIR} )
ENDIF(EIGEN3_FOUND)


SET(QT_LIBRARIES ${Qt5Core_LIBRARIES})

FILE( GLOB MySRCS "*.cpp" "*.h" )
ADD_EXECUTABLE(${PROJECT_NAME} ${MySRCS} )

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${VTK_LIBRARIES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${QT_LIBRARIES})