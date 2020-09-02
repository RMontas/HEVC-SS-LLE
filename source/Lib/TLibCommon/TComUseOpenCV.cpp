/**
 \file     TComUseOpenCV.cpp
 \brief    functions for show images with opencv
 */

#include <stdint.h>
#include <cassert>
#include <vector>

#include "TComUseOpenCV.h"

#if( HAVE_OPENCV == 1 )
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/contrib/contrib.hpp>
#include <opencv2/imgproc/types_c.h>
using namespace cv;
#endif




TUseOpenCV::TUseOpenCV()
{
#if( HAVE_OPENCV == 1 )
  m_bHaveOpenCV = 1;
#else
  m_bHaveOpenCV = 0;
#endif
}

#if( HAVE_OPENCV == 1 )

static void copyData( Pel* src, Int iSrcStride, UChar* dst, Int iDstStride, Int iWidth, Int iHeight,
    Int bIs16bit )
{
  if ( !bIs16bit )
  {
    for ( Int y = 0; y < iHeight; y++ )
    {
      for ( Int x = 0; x < iWidth; x++ )
      {
        dst[x] = (UChar) src[x];
      }
      src += iSrcStride;
      dst += iDstStride;
    }
  }
//  else
//  {
//    for ( Int y = 0; y < pPicYuv->getHeight(); y++ )
//    {
//      for ( Int x = 0; x < pPicYuv->getWidth(); x++ )
//      {
//        dst[2 * x] = src[x] & 0xff;
//        dst[2 * x + 1] = (src[x] >> 8) & 0xff;
//      }
//      src += pPicYuv->getStride();
//    }
//  }
}

static Mat convertTComPicYuv2CVMat( TComPicYuv* pPicYuv, Bool bHaveChroma )
{
  Mat acMatYCbCr[bHaveChroma ? 3 : 1], cMatRGB;

  Int Width = pPicYuv->getWidth();
  Int Height = pPicYuv->getHeight();

  acMatYCbCr[0].create( Height, Width, CV_8UC1);

  copyData( pPicYuv->getLumaAddr(), pPicYuv->getStride(), acMatYCbCr[0].data, acMatYCbCr[0].cols,
      Width, Height, false );

  if ( bHaveChroma )
  {
    Mat acMatCbCrHalf[2], cMatYUV;

    acMatCbCrHalf[0].create( Height / 2, Width / 2, CV_8UC1);
    acMatCbCrHalf[1].create( Height / 2, Width / 2, CV_8UC1);
    copyData( pPicYuv->getCbAddr(), pPicYuv->getCStride(), acMatCbCrHalf[0].data,
        acMatCbCrHalf[0].cols, Width >> 1, Height >> 1, false );

    copyData( pPicYuv->getCrAddr(), pPicYuv->getCStride(), acMatCbCrHalf[1].data,
        acMatCbCrHalf[1].cols, Width >> 1, Height >> 1, false );

    resize( acMatCbCrHalf[0], acMatYCbCr[1], Size( Width, Height ), 0, 0, INTER_CUBIC );
    resize( acMatCbCrHalf[1], acMatYCbCr[2], Size( Width, Height ), 0, 0, INTER_CUBIC );

    merge( acMatYCbCr, 3, cMatYUV );

    cvtColor( cMatYUV, cMatRGB, CV_YCrCb2RGB, 0 );
  }
  else
  {
    cMatRGB = acMatYCbCr[0].clone();
  }
  return cMatRGB;
}

#endif


void TUseOpenCV::showTComPicYuv( TComPicYuv* pPicYuv, char *WindowName )
{
#if( HAVE_OPENCV == 1 )
  Mat cMatcvPic = convertTComPicYuv2CVMat( pPicYuv, true );
  //namedWindow( WindowName, WINDOW_AUTOSIZE);
  namedWindow( WindowName, WINDOW_NORMAL);
  imshow( WindowName, cMatcvPic );
  waitKey( 10 );
//  destroyWindow( WindowName );
#endif
  return;
}

void TUseOpenCV::saveTComPicYuv( TComPicYuv* pPicYuv, char *FileName )
{
#if( HAVE_OPENCV == 1 )
  Mat cMatcvPic = convertTComPicYuv2CVMat( pPicYuv, false );
  imwrite( FileName, cMatcvPic);
#endif
  return;
}

int TUseOpenCV::getAverageDisparity( TComPicYuv* pPicYuvLeft, TComPicYuv* pPicYuvRight )
{
  Int iMeanDisp = 0;
#if( HAVE_OPENCV == 1 )
  StereoVar var;
  Mat cMatcvPicLeft = convertTComPicYuv2CVMat( pPicYuvLeft, false );
  Mat cMatcvPicRight = convertTComPicYuv2CVMat( pPicYuvRight, false );
  Int numberOfDisparities = 128;

  var.levels = 3;                                 // ignored with USE_AUTO_PARAMS
  var.pyrScale = 0.5;                             // ignored with USE_AUTO_PARAMS
  var.nIt = 25;
  var.minDisp = -numberOfDisparities;
  var.maxDisp = 0;
  var.poly_n = 3;
  var.poly_sigma = 0.0;
  var.fi = 15.0f;
  var.lambda = 0.03f;
  var.penalization = var.PENALIZATION_TICHONOV;   // ignored with USE_AUTO_PARAMS
  var.cycle = var.CYCLE_V;                        // ignored with USE_AUTO_PARAMS
  var.flags = var.USE_SMART_ID | var.USE_AUTO_PARAMS | var.USE_INITIAL_DISPARITY | var.USE_MEDIAN_FILTERING ;

  Mat disp;
  var( cMatcvPicLeft, cMatcvPicRight, disp);

  iMeanDisp = mean( disp ).val[0] / (256/128);

//  namedWindow( "Disparity", CV_WINDOW_NORMAL);
//  imshow( "Disparity", disp );
//  waitKey( 0 );
//  destroyWindow( "Disparity" );

#endif
  return iMeanDisp;
}

