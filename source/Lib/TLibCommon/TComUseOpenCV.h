/**
 \file     TComUseOpenCV.h
 \brief    functions for show images with opencv
 */

#ifndef __TCOMUSEOPENCV_H__
#define __TCOMUSEOPENCV_H__

#include "TComPicYuv.h"


#ifndef HAVE_OPENCV
#define HAVE_OPENCV 0
#endif

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// Class for use OpenCV functions
class TUseOpenCV
{
private:
  Bool m_bHaveOpenCV;
public:
  TUseOpenCV();
  void showTComPicYuv( TComPicYuv* pPicYuv, char *WindowName );
  void saveTComPicYuv( TComPicYuv* pPicYuv, char *FileName );
  int getAverageDisparity( TComPicYuv* pPicYuvLeft, TComPicYuv* pPicYuvRight );
protected:
};

#endif //  __TCOMUSEOPENCV_H__
