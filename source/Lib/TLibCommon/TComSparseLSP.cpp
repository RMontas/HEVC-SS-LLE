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

/** \file     TComSparseLSP.cpp
    \brief    sparse-lsp class
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <malloc.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <unistd.h>              /*to use write(), read() and close() with HP*/
#include <limits.h>
#include <float.h>

#include "TComSparseLSP.h"

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Constructor / destructor / initialize
// ====================================================================================================================

TComSparseLSP::TComSparseLSP()
{
  m_imXsize=0;
  m_imYsize=0;
  CausalImage=NULL;
  DstBlock=NULL;
  m_MotionVector=NULL;
  m_TM_errors=NULL;
  tm_search_success=-1;
}

TComSparseLSP::~TComSparseLSP()
{
  
}

void TComSparseLSP::initCausalImage(Pel* piSrc, Int imgWidth, Int imgHeight, Int iSrcStride)
{
  m_imXsize = imgWidth;
  m_imYsize = imgHeight;

  CausalImage=(Pel **)xMalloc(Pel *,imgHeight);

  for(Int i=0; i<imgHeight; i++)
  {
    CausalImage[i]=(Pel *)piSrc;
    piSrc+=iSrcStride;
  }
}

void TComSparseLSP::initDstBlock(Pel *piPred, Int iWidth, Int iHeight, Int iPredStride )
{
  DstBlock=(Pel **)xMalloc(Pel *,iHeight);

  for(Int i=0; i<iHeight; i++)
  {
    DstBlock[i]=(Pel *)piPred;
    piPred+=iPredStride;
  }
}

void TComSparseLSP::executeTemplateSearch(int blkWidth, int blkHeight, int xblock, int yblock)
{
  int imgWidth  = m_imXsize;
  int imgHeight = m_imYsize;
  Pel **Image_out = CausalImage;

  m_availparts[0]=1;
  m_availparts[1]=1;
  m_availparts[2]=1;
  m_availparts[3]=0;
  m_T=4;
  m_maxNumTemplates=8;
  m_MotionVector=NULL;
  m_TM_errors=NULL;

  tm_search_success=template_matching_search(Image_out, Image_out, imgWidth, imgHeight, blkWidth, blkHeight, xblock, yblock, 128, 128, 64, &m_MotionVector, &m_TM_errors, &m_maxNumTemplates, m_T, m_availparts);
}


void TComSparseLSP::destroyBuffers()
{
  if(CausalImage)
    xFree(CausalImage);
  CausalImage=NULL;

  if(DstBlock)
    xFree(DstBlock);
  DstBlock=NULL;

  if(m_MotionVector)
    xFree(m_MotionVector);
  m_MotionVector=NULL;

  if(m_TM_errors)
    xFree(m_TM_errors);
  m_TM_errors=NULL;
}

int TComSparseLSP::predLLE( int blkWidth, int blkHeight, int xblock, int yblock, int uiMode )
{
  int imgWidth  = m_imXsize;
  int imgHeight = m_imYsize;
  Pel **Image_out = CausalImage;
  Pel **dstBlock = DstBlock;

  assert(Image_out);
  assert(dstBlock);

  if( (!Image_out) || (!dstBlock) )
    return 0;

  return xpredLLE ( Image_out, dstBlock, imgWidth, imgHeight, blkWidth, blkHeight, xblock, yblock, uiMode );
}

int TComSparseLSP::getTMstatus()
{
  return tm_search_success;
}

int TComSparseLSP::xpredLLE ( Pel **Image_out, Pel **pred_block, int imgWidth, int imgHeight, int blkWidth, int blkHeight, int xblock, int yblock, int uiMode )
{
  int i, j, k;
  int availparts[4], T;
  Motion_vector *MvVector=m_MotionVector;
  float *lle_errors=NULL, min_lle_error;
  float *errorsTM=m_TM_errors;
  int numtemplates, maxK_NN, best_lle_k=-1;
  double *lle_coefs=NULL, *best_lle_coefs=NULL, *aux_double_ptr=NULL, prediction;

  T = m_T;

  availparts[0]=m_availparts[0];
  availparts[1]=m_availparts[1];
  availparts[2]=m_availparts[2];
  availparts[3]=m_availparts[3];

  maxK_NN=(uiMode+1)>>2;

  numtemplates=maxK_NN;

  if(tm_search_success==1)
  {

    for(i=0; i<numtemplates; i++)
    {
      if(errorsTM[i]==FLT_MAX)
        break;
    }
    if(maxK_NN>i)
      maxK_NN=i;

    lle_coefs=(double*)malloc(sizeof(double)*maxK_NN);
    best_lle_coefs=(double*)malloc(sizeof(double)*maxK_NN);
    lle_errors=(float*)malloc(sizeof(float)*maxK_NN);

    min_lle_error=FLT_MAX;
  //  for(k=1; k<=maxK_NN; k++)
  //  {
      k=maxK_NN;
      lle_errors[k-1] = local_linear_embedding_function(Image_out,imgWidth,imgHeight,MvVector,lle_coefs,xblock,yblock,blkWidth,blkHeight,k,T,availparts);

      if(lle_errors[k-1]<min_lle_error)
      {
        min_lle_error=lle_errors[k-1];
        best_lle_k=k;
        aux_double_ptr=best_lle_coefs;
        best_lle_coefs=lle_coefs;
        lle_coefs=aux_double_ptr;
      }
//    }

    for (i=0; i<blkHeight; i++)
    {
      for (j=0; j<blkWidth; j++)
      {
        prediction=0.0;
        for(k=0;k<best_lle_k;k++)
          prediction+=best_lle_coefs[k]*(double)(Image_out[yblock+i+MvVector[k].mv_y][xblock+j+MvVector[k].mv_x]);

        pred_block[i][j] = Clip3(0,255,(int)rint(prediction));
        Image_out[yblock+i][xblock+j]=pred_block[i][j];
      }
    }

    for(i=0;i<blkHeight;i++)
      for(j=0;j<blkWidth;j++)
        Image_out[yblock+i][xblock+j] = -1;

    free(lle_coefs);
    free(lle_errors);
    free(best_lle_coefs);

    return 1;
  }
  return 0;
}


int TComSparseLSP::tm_search_for_lsp(Pel **curImg, Pel **refImg, int imgWidth, int imgHeight, int width, int height, int posX, int posY, int dtop, int dleft, int dright, Motion_vector **vectorsTM, float **errorsTM, int *numtemplates, int T, int availparts[4])
{
  int i, j, err_tmp;
  float err_total=FLT_MAX;
//  Pel **orig_pels, **tm_pels;
  int tam[4], tam_total;
  int avail[4]={1,1,1,1};
  int blockX=width;
  int blockY=height;
  int imgsizeX=imgWidth;
  int imgsizeY=imgHeight;
  float *tm_errors=NULL, minError;
  Motion_vector *tm_vectors=NULL, bestMv;
  int ind, total_ind, best_index;

  /*  2  3  4
   *  1  B
   */

  tam[0]=T*blockY;
  tam[1]=T*T;
  tam[2]=T*blockX;
  tam[3]=T*blockX;

  //clean-up any garbish that may be in the buffer
  for(j=0;j<blockY;j++)
     for(i=0;i<blockX;i++)
        refImg[posY+j][posX+i] = -1;


  if( (posY-dtop-T<0) || (posX-dleft-T<0) )
  {
     while( (posY-dtop-T) < 0 )
        dtop--;

     while( (posX-dleft-T) < 0 )
        dleft--;
  }

  if( (posX+blockX+dright>imgsizeX) )
  {
     while( (posX+blockX+dright) > imgsizeX )
       dright--;
  }

  // dright should be greater or equal to zero
  assert(dright>=0);

  if(dtop<=0 && dleft<=0)
  {
          dtop=T;
          dleft=T;
          avail[0]=0;
          avail[1]=0;
          avail[2]=0;
          avail[3]=0;
          return 0;
  }
  if(dtop<0)
  {
          dtop=T;
          avail[0]=1;
          avail[1]=0;
          avail[2]=0;
          avail[3]=0;
          return 0;
  }
  if(dleft<0)
  {
          dleft=T;
          avail[0]=0;
          avail[1]=0;
          avail[2]=1;
          avail[3]=1;
          return 0;
  }

  avail[0]&=availparts[0];
  avail[1]&=availparts[1];
  avail[2]&=availparts[2];
  avail[3]&=availparts[3];
  if(!avail[0] && !avail[1] && !avail[2] && !avail[3])
  {
      return 0;
  }

  if(availparts[3]==1)
  {
    if( (posX+blockX*2+dright>imgsizeX) || ( refImg[posY-1][posX+blockX*2-1] == -1) )
    {
          return 0;
    }
  }

  tam_total=0;
  if(avail[0])
    tam_total+=tam[0];
  if(avail[1])
      tam_total+=tam[1];
  if(avail[2])
      tam_total+=tam[2];
  if(avail[3])
      tam_total+=tam[3];


  total_ind=(dleft+dright+1)*(dtop+1);
  tm_errors=(float*)malloc(sizeof(float)*total_ind);
  tm_vectors=(Motion_vector*)malloc(sizeof(Motion_vector)*total_ind);

  for(i=0; i<total_ind; i++){
    tm_errors[i]=FLT_MAX;
    tm_vectors[i].mv_x=0;
    tm_vectors[i].mv_y=0;
  }

  err_tmp=0;
  err_total=FLT_MAX;

  for(i=0; i>=-dtop; i--){
    for(j=-dleft; j<=dright; j++)
    {
        if(i==0 && j==0)
        {
          continue;
        }

        //Check bottom-right pixel of the block associated to the template
        if( (refImg[posY+i+blockY-1][posX+j+blockX-1] == -1) && (j>0))
        {
          continue;
        }

        err_tmp = TMQuadError_all(curImg, refImg, posX, posY, posX+j, posY+i, blockX, blockY, T, avail);

        err_total=(float)err_tmp;///(float)tam_total;

        ind=(i*(-1))*(dright+dleft+1)+(j+dleft);

        tm_errors[ind]=err_total;
        tm_vectors[ind].mv_x=j;
        tm_vectors[ind].mv_y=i;

    }
  }

  if(err_total==FLT_MAX)
  {
    free(tm_errors);
    free(tm_vectors);
    return 0;
  }

  if(total_ind<(*numtemplates))
    (*numtemplates)=total_ind;
  *errorsTM=(float*)malloc(sizeof(float)*(*numtemplates));
  *vectorsTM=(Motion_vector*)malloc(sizeof(Motion_vector)*(*numtemplates));

  for(i=0; i<(*numtemplates); i++)
  {
    minError=FLT_MAX;
    bestMv.mv_x=0;
    bestMv.mv_y=0;
    best_index=0;

    for(j=0; j<total_ind; j++)
    {
      if(tm_errors[j]<minError)
      {
        minError=tm_errors[j];
        bestMv=tm_vectors[j];
        best_index=j;
      }
    }
    (*errorsTM)[i]=minError;
    (*vectorsTM)[i]=bestMv;
    tm_errors[best_index]=FLT_MAX;
  }

  free(tm_errors);
  free(tm_vectors);

  return 1;
}



int TComSparseLSP::template_matching_search(Pel **curImg, Pel **refImg, int imgWidth, int imgHeight, int width, int height, int posX, int posY, int dtop, int dleft, int dright, Motion_vector **vectorsTM, float **errorsTM, int *numtemplates, int T, int availparts[4])
{
  int i, j, err_tmp;
  float err_total=FLT_MAX;
//  Pel **orig_pels, **tm_pels;
  int tam[4], tam_total;
  int avail[4]={1,1,1,1};
  int blockX=width;
  int blockY=height;
  int imgsizeX=imgWidth;
  int imgsizeY=imgHeight;
  float *tm_errors=NULL, minError;
  Motion_vector *tm_vectors=NULL, bestMv;
  int ind, total_ind, best_index;

  /*  2  3  4
   *  1  B
   */

  tam[0]=T*blockY;
  tam[1]=T*T;
  tam[2]=T*blockX;
  tam[3]=T*blockX;

  //clean-up any garbish that may be in the buffer
  for(j=0;j<blockY;j++)
     for(i=0;i<blockX;i++)
        refImg[posY+j][posX+i] = -1;


  if( (posY-dtop-T<0) || (posX-dleft-T<0) )
  {
     while( (posY-dtop-T) < 0 )
        dtop--;

     while( (posX-dleft-T) < 0 )
        dleft--;
  }

  if( (posX+blockX+dright>imgsizeX) )
  {
     while( (posX+blockX+dright) > imgsizeX )
       dright--;
  }

  // dright should be greater or equal to zero
  assert(dright>=0);

  if(dtop<=0 && dleft<=0)
  {
          dtop=T;
          dleft=T;
          avail[0]=0;
          avail[1]=0;
          avail[2]=0;
          avail[3]=0;
          return 0;
  }
  if(dtop<0)
  {
          dtop=T;
          avail[0]=1;
          avail[1]=0;
          avail[2]=0;
          avail[3]=0;
          return 0;
  }
  if(dleft<0)
  {
          dleft=T;
          avail[0]=0;
          avail[1]=0;
          avail[2]=1;
          avail[3]=1;
          return 0;
  }

  avail[0]&=availparts[0];
  avail[1]&=availparts[1];
  avail[2]&=availparts[2];
  avail[3]&=availparts[3];
  if(!avail[0] && !avail[1] && !avail[2] && !avail[3])
  {
      return 0;
  }

  if(availparts[3]==1)
  {
    if( (posX+blockX*2+dright>imgsizeX) || ( refImg[posY-1][posX+blockX*2-1] == -1) )
    {
          return 0;
    }
  }

  tam_total=0;
  if(avail[0])
    tam_total+=tam[0];
  if(avail[1])
      tam_total+=tam[1];
  if(avail[2])
      tam_total+=tam[2];
  if(avail[3])
      tam_total+=tam[3];


  total_ind=(dleft+dright+1)*(dtop+1);
  tm_errors=(float*)malloc(sizeof(float)*total_ind);
  tm_vectors=(Motion_vector*)malloc(sizeof(Motion_vector)*total_ind);

  for(i=0; i<total_ind; i++){
    tm_errors[i]=FLT_MAX;
    tm_vectors[i].mv_x=0;
    tm_vectors[i].mv_y=0;
  }


  err_tmp=0;
  err_total=FLT_MAX;

  for(i=0; i>=-dtop; i--){
    for(j=-dleft; j<=dright; j++)
    {
        if(i==0 && j==0)
        {
          continue;
        }

        //Check bottom-right pixel of the block associated to the template
        if( refImg[posY+i+blockY-1][posX+j+blockX-1] == -1)
        {
          continue;
        }

        err_tmp = TMQuadError_all(curImg, refImg, posX, posY, posX+j, posY+i, blockX, blockY, T, avail);

        err_total=(float)err_tmp;///(float)tam_total;

        ind=(i*(-1))*(dright+dleft+1)+(j+dleft);

        tm_errors[ind]=err_total;
        tm_vectors[ind].mv_x=j;
        tm_vectors[ind].mv_y=i;

    }
  }

  if(err_total==FLT_MAX)
  {
    free(tm_errors);
    free(tm_vectors);
    return 0;
  }

  if(total_ind<(*numtemplates))
    (*numtemplates)=total_ind;
  *errorsTM=(float*)malloc(sizeof(float)*(*numtemplates));
  *vectorsTM=(Motion_vector*)malloc(sizeof(Motion_vector)*(*numtemplates));

  for(i=0; i<(*numtemplates); i++)
  {
    minError=FLT_MAX;
    bestMv.mv_x=0;
    bestMv.mv_y=0;
    best_index=0;

    for(j=0; j<total_ind; j++)
    {
      if(tm_errors[j]<minError)
      {
        minError=tm_errors[j];
        bestMv=tm_vectors[j];
        best_index=j;
      }
    }
    (*errorsTM)[i]=minError;
    (*vectorsTM)[i]=bestMv;
    tm_errors[best_index]=FLT_MAX;
  }

  free(tm_errors);
  free(tm_vectors);

  return 1;
}


float TComSparseLSP::local_linear_embedding_function(Pel **data_recon, int imgWidth, int imgHeight, Motion_vector *vectorsTM, double *lle_coefs, int posX, int posY, int width, int height, int k_NN, int T, int availparts[4])
{
  int i, j, k, l;
  Pel *ref_tm, **nearest_patches;
  int tam[4], tam_total;
  int blockX=width;
  int blockY=height;
  int imgsizeX=imgWidth;
  int imgsizeY=imgHeight;
  int offsetX, offsetY;
  Motion_vector mv;
  double **Z, *Xi, *coefs, coefmax, prediction;
  Pel *pred;
  int lle_error;

        /*  2  3  4
         *  1  B
         */
  tam[0]=T*blockY;
  tam[1]=T*T;
  tam[2]=T*blockX;
  tam[3]=T*blockX;

  tam_total=0;
  if(availparts[0])
    tam_total+=tam[0];
  if(availparts[1])
    tam_total+=tam[1];
  if(availparts[2])
    tam_total+=tam[2];
  if(availparts[3])
    tam_total+=tam[3];

  //clean-up any garbish that may be in the buffer
  for(j=0;j<blockY;j++)
    for(i=0;i<blockX;i++)
      data_recon[posY+j][posX+i] = -1;

  coefs=(double*)malloc(sizeof(double)*k_NN);

  ref_tm=(Pel*)malloc(sizeof(Pel)*tam_total);
  Xi=(double*)malloc(sizeof(double)*tam_total);
  nearest_patches=pelmatrix(k_NN,tam_total);
  Z=doublematrix(k_NN,tam_total);

  getTMpels_inColumn(data_recon, ref_tm, posX, posY, blockX, blockY, T, availparts);
  for(i=0; i<tam_total ;i++)
    Xi[i]=(double)ref_tm[i];

  for(k=0; k<k_NN; k++)
  {
    offsetX=vectorsTM[k].mv_x;
    offsetY=vectorsTM[k].mv_y;
    getTMpels_inColumn(data_recon, nearest_patches[k], posX+offsetX, posY+offsetY, blockX, blockY, T, availparts);
  }

  /*
#ifdef DEBUG_LLE
  printf("\n\n ------------> Use LLE:\n");
#endif
  for(i=0; i<k_NN ;i++)
    for(j=0; j<tam_total ;j++)
    {
      Z[i][j]=(double)nearest_patches[i][j];
    }
  solve_lle_problem(Z,Xi,coefs,k_NN,tam_total);
*/
#ifdef DEBUG_LLE
  printf("\n\n ------------> Use LSP:\n");
#endif
  for(i=0; i<k_NN ;i++)
    for(j=0; j<tam_total ;j++)
    {
      Z[i][j]=(double)nearest_patches[i][j];
    }
  solve_lle_problem(Z,Xi,coefs,k_NN,tam_total);

  coefmax=0.0;
  for(l=0;l<k_NN;l++)
    if(fabs(coefs[l])>coefmax)
      coefmax=fabs(coefs[l]);
  if (coefmax>2.0)
  {
    for(l=0;l<k_NN;l++)
      coefs[l]=1.0/k_NN;
  }

  pred=(Pel*)malloc(sizeof(Pel)*tam_total);

  for (i=0; i<tam_total; i++)
  {
    prediction=0;
    for(k=0;k<k_NN;k++)
      prediction+=coefs[k]*(double)(nearest_patches[k][i]);

    pred[i] = Clip3((Pel)0,(Pel)255,(Pel)rint(prediction));
  }

  lle_error=0;
  for (i=0; i<tam_total; i++)
  {
    lle_error+=(Xi[i]-pred[i])*(Xi[i]-pred[i]);
  }

  // Put coefficients in input pointer variable
  for(i=0; i<k_NN ;i++){
    lle_coefs[i]=coefs[i];
  }

  free_pelmatrix(nearest_patches,k_NN,tam_total);
  free_doublematrix(Z,k_NN,tam_total);
  free(coefs);
  free(Xi);
  free(ref_tm);
  free(pred);

  return (float)(lle_error);
}


void TComSparseLSP::solve_lsp_problem(double **C,double *y,double *a, int m, int n)
{
  int i,j,k,*index,flag;
  double **cov,d,*col, sum_a;

#ifdef DEBUG_LLE
  printf("\nC:\n");
  for(i=0; i<m ;i++){
    for(j=0; j<n ;j++)
      printf(" %f", C[i][j]);
    printf("\n");
  }
  printf("\ny:\n");
  for(j=0; j<n ;j++)
    printf(" %f", y[j]);
  printf("\n");
#endif

  index=(int *)malloc((m+1)*sizeof(int));
  col=(double *)malloc((m+1)*sizeof(double));
  cov=(double **)malloc((m+1)*sizeof(double *));
  for (i=0;i<m+1;i++)
    cov[i]=(double *)malloc((m+1)*sizeof(double));
  for(i=0;i<m;i++)
    for(j=i;j<m;j++)
    {
      cov[i+1][j+1]=0;
      for(k=0;k<n;k++)
        cov[i+1][j+1]+=C[i][k]*C[j][k];
    }
  for(i=0;i<m;i++)
    for(j=0;j<=i;j++)
      cov[i+1][j+1]=cov[j+1][i+1];
  for(i=0;i<m;i++)
  {
    col[i+1]=0;
    for(j=0;j<n;j++)
      col[i+1]+=C[i][j]*y[j];
  }
#ifdef DEBUG_LLE
  printf("\ncov:\n");
  for(i=1; i<=m ;i++){
    for(j=1; j<=m ;j++)
      printf(" %f", cov[i][j]);
    printf("\n");
  }
  printf("\ncol:\n");
  for(i=1; i<=m ;i++)
    printf(" %f", col[i]);
  printf("\n");

#endif

  flag=ludcmp(cov,m,index,&d);
#ifdef DEBUG_LLE
  printf("\nflag:%d", flag);
  printf("\ncov2:\n");
  for(i=1; i<=m ;i++){
    for(j=1; j<=m ;j++)
      printf(" %f", cov[i][j]);
    printf("\n");
  }
#endif
  if(flag==0)
  {
    lubksb(cov,m,index,col);
    sum_a=0;
    for(i=1;i<=m;i++){
      a[i-1]=col[i];
      sum_a+=col[i];
    }
    for(i=0;i<m;i++)
      a[i]=a[i]/sum_a;
  }
  else
    for(i=1;i<=m;i++)
      a[i-1]=1.0/m;

#ifdef DEBUG_LLE
  printf("\na:\n");
  for(j=0; j<m ;j++)
    printf(" %f", a[j]);
  printf("\n");
#endif
  free(index);free(col);
  for(i=0;i<m+1;i++)
    free(cov[i]);
  free(cov);
}

/*!
 ********************************************************************************************
 *  \brief
 *      Solve LLE optimization problem based on Z' matrix and Xi vector.
 *      See: http://www.cs.nyu.edu/~roweis/lle/algorithm.html
 *
 *      Xi is original template
 *      m - number of rows of Z matrix
 *      n - number of columns of Z matrix
 *
 ********************************************************************************************
*/
void TComSparseLSP::solve_lle_problem(double **Z,double *Xi,double *a, int m, int n)
{
  int i,j,k,*index,flag;
  double **cov,d,*col, sum_a;

#ifdef DEBUG_LLE
  printf("\nZ:\n");
  for(i=0; i<m ;i++){
    for(j=0; j<n ;j++)
      printf(" %f", Z[i][j]);
    printf("\n");
  }
  printf("\nXi:\n");
  for(j=0; j<n ;j++)
    printf(" %f", Xi[j]);
  printf("\n");
#endif

  index=(int *)malloc((m+1)*sizeof(int));
  col=(double *)malloc((m+1)*sizeof(double));
  cov=(double **)malloc((m+1)*sizeof(double *));
  for (i=0;i<m+1;i++)
    cov[i]=(double *)malloc((m+1)*sizeof(double));
  for(i=0; i<m ;i++)
    for(j=0; j<n ;j++)
    {
      Z[i][j]=Z[i][j]-Xi[j];
    }
  for(i=0;i<m;i++)
    for(j=i;j<m;j++)
    {
      cov[i+1][j+1]=0;
      for(k=0;k<n;k++)
        cov[i+1][j+1]+=Z[i][k]*Z[j][k];
    }
  for(i=0;i<m;i++)
    for(j=0;j<=i;j++)
      cov[i+1][j+1]=cov[j+1][i+1];

#ifdef DEBUG_LLE
  printf("\nZ2:\n");
  for(i=0; i<m ;i++){
    for(j=0; j<n ;j++)
      printf(" %f", Z[i][j]);
    printf("\n");
  }
  printf("\ncov:\n");
  for(i=1; i<=m ;i++){
    for(j=1; j<=m ;j++)
      printf(" %f", cov[i][j]);
    printf("\n");
  }
#endif

  for(i=0;i<m;i++)
  {
    col[i+1]=1;
  }
  flag=ludcmp(cov,m,index,&d);
/*  for(i=1; i<=m; i++)
    if(cov[i][i]==0.0)
    {
#ifdef DEBUG_LLE
      printf("\nIt wasn't full rank.\nflag:%d", flag);
      printf("\ncov2b:\n");
      for(i=1; i<=m ;i++){
        for(j=1; j<=m ;j++)
          printf(" %f", cov[i][j]);
        printf("\n");
      }
#endif
      for(j=1; j<=m; j++)
        cov[j][j]+=TINY;
      flag=ludcmp(cov,m,index,&d);
      break;
    }
*/
#ifdef DEBUG_LLE
  printf("\nflag:%d", flag);
  printf("\ncov2:\n");
  for(i=1; i<=m ;i++){
    for(j=1; j<=m ;j++)
      printf(" %f", cov[i][j]);
    printf("\n");
  }
#endif
  if(flag==0)
  {
    lubksb(cov,m,index,col);
    sum_a=0;
    for(i=1;i<=m;i++){
      a[i-1]=col[i];
      sum_a+=col[i];
    }
    for(i=0;i<m;i++)
      a[i]=a[i]/sum_a;
  }
  else
    for(i=1;i<=m;i++)
      a[i-1]=1.0/m;

#ifdef DEBUG_LLE
  printf("\na:\n");
  for(j=0; j<m ;j++)
    printf(" %f", a[j]);
  printf("\n");
#endif
  free(index);free(col);
  for(i=0;i<m+1;i++)
    free(cov[i]);
  free(cov);
}


/*!
 ********************************************************************************************
 * \brief
 *    LU decomposition : write a matrix as a product of two matrices, L (lower) and U (upper)
 *    Crout's algorithm is the decomposition in place. The follwoing algorithm doesn't actually
 *    decompose the matrix A in LU form; it rather decompose a rowwise permutation of the input
 *    matrix. The permutation is recorded in the index output vector
 *
 *  \param input
 *  double matrix a
 *  size of matrix n (A[nxn])
 *  \param output
 *  double matrix a (decomposed matrix)
 *  int vector index (permutation is recorded here)
 *  d = +1 (even rows interchanged), -1 (odd rows interchanged)
 *
 ********************************************************************************************
*/
int TComSparseLSP::ludcmp(double **a,int n, int *indx,double *d)
{
int i,j,k;
int imax = 0;
double big,dum,sum,temp;
double *vv;

  vv=(double *)malloc((n+1)*sizeof(double));
  *d=1.0;
  for(i=1;i<=n;i++)
  {
    big=0.0;
    for(j=1;j<=n;j++)
      if((temp=fabs(a[i][j]))>big)
        big=temp;
    if(big==0.0)
    {
      free(vv);
      return(1);
    }
    vv[i]=1.0/big;
  }
  for(j=1;j<=n;j++)
  {
    for(i=1;i<j;i++)
    {
      sum=a[i][j];
      for(k=1;k<i;k++)
        sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for(i=j;i<=n;i++)
    {
      sum=a[i][j];
      for(k=1;k<j;k++)
        sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
      if((dum=vv[i]*fabs(sum))>=big)
      {
        big=dum;
        imax=i;
      }
    }
    if(j!=imax)
    {
      for(k=1;k<=n;k++)
      {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d=-(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if(a[j][j]==0.0)
      a[j][j]=TINY;
    if(j!=n)
    {
      dum=1.0/(a[j][j]);
      for(i=j+1;i<=n;i++)
        a[i][j]*=dum;
    }
  }
  free(vv);
  return(0);
}

/*!
 ********************************************************************************************
 *  \brief
 *      The lubksb algorithm (for forward substitution and backward substitution) solves the
 *      set of N linear equations Ax=b, where A is the output of the ludcmp algorithm.
 *
 * \param input:
 *  double matrix a (output from ludcmp, it is not altered by the algorithm)
 *  size n of matrix (A[nxn])
 *  int vector indx (permutations)
 *  double vector b
 *   \param output:
 *  double vector b (includes the solution of the matrix
 ********************************************************************************************
*/
void TComSparseLSP::lubksb(double **a,int n,int *indx,double b[])
{
int i,ii=0,ip,j;
double sum;

  for(i=1;i<=n;i++)
  {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if(ii)
      for(j=ii;j<=i-1;j++)
        sum-=a[i][j]*b[j];
    else if
      (sum) ii=i;
    b[i]=sum;
  }
  for(i=n;i>=1;i--)
  {
    sum=b[i];
    for(j=i+1;j<=n;j++)
      sum-=a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


void TComSparseLSP::getTMpels_inColumn(Pel **img, Pel *buffer, int xpel, int ypel, int xblksize, int yblksize, int T, int avail[4])
{
  int i, j, ind;

  ind=0;

  if(avail[0])
  {
    for(i=ypel; i<ypel+yblksize; i++){
      for(j=xpel-T; j<xpel; j++)
      {
        buffer[ind]=img[i][j];
        ind++;
      }
    }
  }
  if(avail[1])
  {
    for(i=ypel-T; i<ypel; i++){
      for(j=xpel-T; j<xpel; j++)
      {
        buffer[ind]=img[i][j];
        ind++;
      }
    }
  }
  if(avail[2])
  {
    for(i=ypel-T; i<ypel; i++){
      for(j=xpel; j<xpel+xblksize; j++)
      {
        buffer[ind]=img[i][j];
        ind++;
      }
    }
  }
  if(avail[3])
  {
    for(i=ypel-T; i<ypel; i++){
      for(j=xpel+xblksize; j<xpel+xblksize*2; j++)
      {
        buffer[ind]=img[i][j];
        ind++;
      }
    }
  }
}


void TComSparseLSP::getTemplatePels_pel(Pel **img, Pel **buffer, int xpel, int ypel, int xblksize, int yblksize, int T, int avail[4])
{
  int i, j, ind;

  if(avail[0])
  {
          ind=0;
          for(i=ypel; i<ypel+yblksize; i++){
                for(j=xpel-T; j<xpel; j++)
                {
                  buffer[0][ind]=img[i][j];
        //      printf(" %d", Img[i][j]);
                  ind++;
                }
          }
  }

  if(avail[1])
  {
          ind=0;
          for(i=ypel-T; i<ypel; i++){
                for(j=xpel-T; j<xpel; j++)
                {
                        buffer[1][ind]=img[i][j];
        //      printf(" %d", Img[i][j]);
                        ind++;
                }
          }
  }

  if(avail[2])
  {
          ind=0;
          for(i=ypel-T; i<ypel; i++){
                for(j=xpel; j<xpel+xblksize; j++)
                {
                        buffer[2][ind]=img[i][j];
        //      printf(" %d", Img[i][j]);
                        ind++;
                }
          }
  }

  if(avail[3])
  {
          ind=0;
          for(i=ypel-T; i<ypel; i++){
                for(j=xpel+xblksize; j<xpel+xblksize*2; j++)
                {
                        buffer[3][ind]=img[i][j];
          //      printf(" %d", Img[i][j]);
                        ind++;
                }
          }
  }
//  printf("\n");
}


inline int TComSparseLSP::TMQuadError_all(Pel **cur_buffer, Pel **ref_buffer, int xcurpel, int ycurpel, int xrefpel, int yrefpel, int xblksize, int yblksize, int T, int avail[4])
{
  int i, j;
  int error=0;

  if(avail[0])
  {
      for(i=0; i<yblksize; i++){
            for(j=-T; j<0; j++)
            {
                error+=(cur_buffer[ycurpel+i][xcurpel+j]-ref_buffer[yrefpel+i][xrefpel+j])*(cur_buffer[ycurpel+i][xcurpel+j]-ref_buffer[yrefpel+i][xrefpel+j]);
            }
      }
  }

  if(avail[1])
  {
      for(i=-T; i<0; i++){
            for(j=-T; j<0; j++)
            {
                error+=(cur_buffer[ycurpel+i][xcurpel+j]-ref_buffer[yrefpel+i][xrefpel+j])*(cur_buffer[ycurpel+i][xcurpel+j]-ref_buffer[yrefpel+i][xrefpel+j]);
            }
      }
  }

  if(avail[2])
  {
      for(i=-T; i<0; i++){
            for(j=0; j<xblksize; j++)
            {
                error+=(cur_buffer[ycurpel+i][xcurpel+j]-ref_buffer[yrefpel+i][xrefpel+j])*(cur_buffer[ycurpel+i][xcurpel+j]-ref_buffer[yrefpel+i][xrefpel+j]);
            }
      }
  }

  if(avail[3])
  {
      for(i=-T; i<0; i++){
            for(j=xblksize; j<xblksize*2; j++)
            {
                error+=(cur_buffer[ycurpel+i][xcurpel+j]-ref_buffer[yrefpel+i][xrefpel+j])*(cur_buffer[ycurpel+i][xcurpel+j]-ref_buffer[yrefpel+i][xrefpel+j]);
            }
      }
  }

  return error;
}


void TComSparseLSP::getTemplatePels_uchar(unsigned char **img, Pel **buffer, int xpel, int ypel, int xblksize, int yblksize, int T, int avail[4])
{
  int i, j, ind;

  if(avail[0])
  {
          ind=0;
          for(i=ypel; i<ypel+yblksize; i++){
                for(j=xpel-T; j<xpel; j++)
                {
                  buffer[0][ind]=(Pel)img[i][j];
        //      printf(" %d", Img[i][j]);
                  ind++;
                }
          }
  }

  if(avail[1])
  {
          ind=0;
          for(i=ypel-T; i<ypel; i++){
                for(j=xpel-T; j<xpel; j++)
                {
                        buffer[1][ind]=(Pel)img[i][j];
        //      printf(" %d", Img[i][j]);
                        ind++;
                }
          }
  }

  if(avail[2])
  {
          ind=0;
          for(i=ypel-T; i<ypel; i++){
                for(j=xpel; j<xpel+xblksize; j++)
                {
                        buffer[2][ind]=(Pel)img[i][j];
        //      printf(" %d", Img[i][j]);
                        ind++;
                }
          }
  }

  if(avail[3])
  {
          ind=0;
          for(i=ypel-T; i<ypel; i++){
                for(j=xpel+xblksize; j<xpel+xblksize*2; j++)
                {
                        buffer[3][ind]=(Pel)img[i][j];
          //      printf(" %d", Img[i][j]);
                        ind++;
                }
          }
  }
//  printf("\n");
}



float TComSparseLSP::tmQuadError(Pel *buffer1, Pel *buffer2, int tam, float min_err)
{
  int i;
  float err=0;

  for(i=0; i<tam && err<=min_err; i++){
      err+=(float)((buffer1[i]-buffer2[i])*(buffer1[i]-buffer2[i]));
//      assert(buffer1[i]>=0);
//      assert(buffer2[i]>=0);
  }

  return err;
}

void TComSparseLSP::initTMpels(Pel *buffer, int tam, int value)
{
        int i;
        for(i=0; i<tam; i++)
        {
                buffer[i]=value;
        }
}

Pel **TComSparseLSP::allocTMpels(int *fraction_sizes, int numfractions)
{
        int i;
        Pel **buffer;

        buffer=(Pel **)calloc( numfractions, sizeof(Pel*) );

        for(i=0; i<numfractions; i++)
                buffer[i]=(Pel*)calloc( fraction_sizes[i], sizeof(Pel) );

        return buffer;
}

void TComSparseLSP::freeTMpels(Pel **buffer, int numfractions)
{
        int i;
        for(i=0; i<numfractions; i++)
          free(buffer[i]);

        free(buffer);
}

/*!
 ********************************************************************************************
 * \brief
 *    Allocates memory for a matrix of variables of type int
 *
 * \param nr
 *    number of rows
 * \param nc
 *     number of columns
 *
 * \return
 *    a poiter to a int matrix (int **)
 *
 ********************************************************************************************
*/
Pel **TComSparseLSP::pelmatrix(int nr, int nc)
{
  int i;
  Pel **m;

  m=(Pel **)malloc((unsigned) (nr)*sizeof(Pel *));
  if (!m){
  printf("intmatrix() - allocation failure 1 \n");
  exit(1);
  }
  for(i=0;i<nr;i++)
  {
  m[i]=(Pel *)calloc((unsigned)(nc),sizeof(Pel));
  if (!m[i]){
    printf("intmatrix() - allocation failure 2 \n");
    exit(1);
      }
  }
  return m;
}

/*!
 ********************************************************************************************
 * \brief
 *    Allocates memory for a matrix of variables of type double
 *
 * \param nr
 *    number of rows
 * \param nc
 *     number of columns
 *
 * \return
 *    a poiter to a float matrix (double **)
 *
 ********************************************************************************************
*/
double **TComSparseLSP::doublematrix(long nr, long nc)
{
  long i;
  double **m;

  m=(double **)malloc(nr*nc*sizeof(double *));
  if (!m)
  {
  printf("doublematrix() - allocation failure 1 \n");
  exit(1);
  }
  for(i=0;i<nr;i++)
  {
  m[i]=(double *)calloc(nc, sizeof(double));
  if (!m[i])
  {
    printf("doublematrix() - allocation failure 2 \n");
    exit(1);
  }
  }
  return m;
}

/*!
 ********************************************************************************************
 * \brief
 *    Deallocates memory for a matrix of variables of type int
 *
 * \param matrix
 *   matrix to be deleted
 * \param nr
 *    number of rows
 * \param nc
 *    number of columns
 *
 * \return
 *    none.
 *
 ********************************************************************************************
*/
void TComSparseLSP::free_pelmatrix(Pel **matrix, int nr, int nc)
{
  int i;

  for (i=0; i<nr; i++)
  free(matrix[i]);
  free(matrix);
}

/*!
 ********************************************************************************************
 * \brief
 *    Deallocates memory for a matrix of variables of type double
 *
 * \param matrix
 *   matrix to be deleted
 * \param nr
 *    number of rows
 * \param nc
 *    number of columns
 *
 * \return
 *    none.
 *
 ********************************************************************************************
*/
void TComSparseLSP::free_doublematrix(double **matrix, int nr, int nc)
{
  int i;

  for (i=0; i<nr; i++)
  free(matrix[i]);
  free(matrix);
}

//! \}
