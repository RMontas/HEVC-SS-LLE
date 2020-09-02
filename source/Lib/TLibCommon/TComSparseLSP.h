/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.  
 *
 * Copyright (c) 2010-2013, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     TComSparseLSP.h
    \brief    sparse-lsp class (header)
*/

#ifndef __TCOMSPARSELSP__
#define __TCOMSPARSELSP__


// Include files
#include "TComPic.h"

/*! \def Motion_vector
    \brief struct Motion_vector is used to store a motion vector */
typedef struct MVector{
        short mv_x;                     // Motion vector X component
        short mv_y;                     // Motion vector Y component
} Motion_vector;

#define TINY 1.0e-20

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// sparse-lsp class
class TComSparseLSP
{
private:
  //TComYuv   m_acYuvPred[2];
  //TComYuv   m_cYuvPredTemp;
  int m_imXsize;
  int m_imYsize;
  Pel **CausalImage;
  Pel **DstBlock;
  
  int m_availparts[4];
  int m_T;
  Motion_vector *m_MotionVector;
  float *m_TM_errors=NULL;
  int m_maxNumTemplates;
  int tm_search_success;


  int template_matching_search(Pel **curImg, Pel **refImg, int imgWidth, int imgHeight, int width, int height, int posX, int posY, int dtop, int dleft, int dright, Motion_vector **vectorsTM, float **errorsTM, int *numtemplates, int T, int availparts[4]);
  int tm_search_for_lsp(Pel **curImg, Pel **refImg, int imgWidth, int imgHeight, int width, int height, int posX, int posY, int dtop, int dleft, int dright, Motion_vector **vectorsTM, float **errorsTM, int *numtemplates, int T, int availparts[4]);
  float local_linear_embedding_function(Pel **data_recon, int imgWidth, int imgHeight, Motion_vector *vectorsTM, double *lle_coefs, int posX, int posY, int width, int height, int k_NN, int T, int availparts[4]);
  void solve_lsp_problem(double **C,double *y,double *a, int m, int n);
  void solve_lle_problem(double **Z,double *Xi,double *a, int m, int n);
  int ludcmp(double **a,int n, int *indx,double *d);
  void lubksb(double **a,int n,int *indx,double b[]);
  void getTMpels_inColumn(Pel **img, Pel *buffer, int xpel, int ypel, int xblksize, int yblksize, int T, int avail[4]);
  void getTemplatePels_pel(Pel **img, Pel **buffer, int xpel, int ypel, int xblksize, int yblksize, int T, int avail[4]);
  void getTemplatePels_uchar(unsigned char **img, Pel **buffer, int xpel, int ypel, int xblksize, int yblksize, int T, int avail[4]);
  float tmQuadError(Pel *buffer1, Pel *buffer2, int tam, float min_err);
  inline int TMQuadError_all(Pel **cur_buffer, Pel **ref_buffer, int xcurpel, int ycurpel, int xrefpel, int yrefpel, int xblksize, int yblksize, int T, int avail[4]);
  void initTMpels(Pel *buffer, int tam, int value);
  Pel **allocTMpels(int *fraction_sizes, int numfractions);
  void freeTMpels(Pel **buffer, int numfractions);
  Pel **pelmatrix(int nr, int nc);
  double **doublematrix(long nr, long nc);
  void free_pelmatrix(Pel **matrix, int nr, int nc);
  void free_doublematrix(double **matrix, int nr, int nc);
  int xpredSparseLSP ( Pel **Image_out, Pel **pred_block, int imgWidth, int imgHeight, int blkWidth, int blkHeight, int xblock, int yblock );
  int xpredLLE ( Pel **Image_out, Pel **pred_block, int imgWidth, int imgHeight, int blkWidth, int blkHeight, int xblock, int yblock, int uiMode );

public:
  TComSparseLSP();
  virtual ~TComSparseLSP();

  void initCausalImage(Pel* piSrc, Int imgWidth, Int imgHeight, Int iSrcStride);
  void initDstBlock( Pel *piPred, Int iWidth, Int iHeight, Int iPredStride );
  void destroyBuffers();
  void executeTemplateSearch(int blkWidth, int blkHeight, int xblock, int yblock);
  int getTMstatus();

  int predSparseLSP( int blkWidth, int blkHeight, int xblock, int yblock, int uiMode );
  int predLLE( int blkWidth, int blkHeight, int xblock, int yblock, int uiMode );


};

//! \}

#endif // __TCOMSPARSELSP__
