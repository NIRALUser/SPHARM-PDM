// -*- C++ -*-
/****************************************************************************
 *
 *  class TVoxelVolume
 *
 ****************************************************************************
 *
 *  File         :  TVoxelVolume.cc
 *  Type         :
 *  Purpose      :
 *
 ****************************************************************************
 *
 *  started      :  4  Jan 95   Robert Schweikert
 *  last change  :  1  Aug 2005   Martin Styner
 *
 ****************************************************************************/

#include <iostream>
#include "TVoxelVolume.h"
#include <string.h>

TVoxelVolume::TVoxelVolume()
{
  dim_[0] = 0;
  dim_[1] = 0;
  dim_[2] = 0;
  data_   = NULL;
  size_   = 0;
}

TVoxelVolume::TVoxelVolume(int xs, int ys, int zs, const unsigned char* data)
{
  dim_[0] = xs;
  dim_[1] = ys;
  dim_[2] = zs;
  size_   = dim_[0] * dim_[1] * dim_[2];

  data_ = new unsigned char[dim_[0] * dim_[1] * dim_[2]];
  if( data == 0 )
    {
    this->Fill(0);
    }
  else
    {
    memcpy(data_, data, sizeof(unsigned char) * dim_[0] * dim_[1] * dim_[2]);
    }
}

TVoxelVolume::TVoxelVolume(const int dim[3], const unsigned char* data)
{
  dim_[0] = dim[0];
  dim_[1] = dim[1];
  dim_[2] = dim[2];
  size_   = dim_[0] * dim_[1] * dim_[2];

  data_ = new unsigned char[dim_[0] * dim_[1] * dim_[2]];
  if( data == 0 )
    {
    this->Fill(0);
    }
  else
    {
    memcpy(data_, data, sizeof(unsigned char) * dim_[0] * dim_[1] * dim_[2]);
    }
}

TVoxelVolume::TVoxelVolume(const TVoxelVolume& vox)
{
  dim_[0] = vox.dim_[0];
  dim_[1] = vox.dim_[1];
  dim_[2] = vox.dim_[2];
  size_   = dim_[0] * dim_[1] * dim_[2];

  data_ = new unsigned char[dim_[0] * dim_[1] * dim_[2]];
  memcpy(this->data_, vox.data_, sizeof(unsigned char) * dim_[0] * dim_[1] * dim_[2]);
}

TVoxelVolume::~TVoxelVolume()
{
  delete [] data_;
}

TVoxelVolume & TVoxelVolume::operator=(const TVoxelVolume& src)
{
  for( long i = 0; i < size_; i++ )
    {
    data_[i] = src.data_[i];
    }
  return *this;
}

TVoxelVolume TVoxelVolume::operator+(const TVoxelVolume& op2)
{
  TVoxelVolume res(dim_, 0);

  for( long i = 0; i < size_; i++ )
    {
    res.data_[i] = data_[i] + op2.data_[i];
    }
  return res;
}

TVoxelVolume & TVoxelVolume::operator+=(const TVoxelVolume& src)
{
  for( long i = 0; i < size_; i++ )
    {
    data_[i] += src.data_[i];
    }
  return *this;
}

TVoxelVolume TVoxelVolume::operator-(const TVoxelVolume& op2)
{
  TVoxelVolume res(dim_, 0);

  for( long i = 0; i < size_; i++ )
    {
    res.data_[i] = data_[i] - op2.data_[i];
    }
  return res;

}

TVoxelVolume & TVoxelVolume::operator-=(const TVoxelVolume& src)
{
  for( long i = 0; i < size_; i++ )
    {
    data_[i] -= src.data_[i];
    }
  return *this;
}

TVoxelVolume TVoxelVolume::operator*(const TVoxelVolume& op2)
{
  TVoxelVolume res(dim_, 0);

  for( long i = 0; i < size_; i++ )
    {
    res.data_[i] = data_[i] * op2.data_[i];
    }
  return res;
}

TVoxelVolume TVoxelVolume::operator*(const unsigned char& val)
{
  TVoxelVolume res(dim_, 0);

  for( long i = 0; i < size_; i++ )
    {
    res.data_[i] = data_[i] * val;
    }
  return res;
}

TVoxelVolume & TVoxelVolume::operator*=(const TVoxelVolume& src)
{
  for( long i = 0; i < size_; i++ )
    {
    data_[i] *= src.data_[i];
    }
  return *this;
}

TVoxelVolume & TVoxelVolume::operator*=(const unsigned  char& val)
{
  for( long i = 0; i < size_; i++ )
    {
    data_[i] *= val;
    }
  return *this;
}

void TVoxelVolume::Fill(unsigned char val)
{
  for( long i = 0; i < size_; data_[i++] = val )
    {
    ;
    }
}

void TVoxelVolume::SetBoundary(unsigned char val)
{
  long x, y, z;

  for( x = 0; x < dim_[0]; x++ )
    {
    for( y = 0; y < dim_[1]; y++ )
      {
      VOXEL(x, y, 0) = VOXEL(x, y, dim_[2] - 1) = val;
      }
    }
  for( y = 0; y < dim_[1]; y++ )
    {
    for( z = 1; z < dim_[2] - 1; z++ )
      {
      VOXEL(0, y, z) = VOXEL(dim_[0] - 1, y, z) = val;
      }
    }
  for( z = 1; z < dim_[2] - 1; z++ )
    {
    for( x = 1; x < dim_[0] - 1; x++ )
      {
      VOXEL(x, 0, z) = VOXEL(x, dim_[1] - 1, z) = val;
      }
    }
}

void TVoxelVolume::Binarize(unsigned char thresh)
// VOXEL(x,y,z) <= thresh -> 0
// VOXEL(x,y,z) >  thresh -> 255
{
  long x, y, z;

  for( z = 0; z < dim_[2]; z++ )
    {
    for( y = 0; y < dim_[1]; y++ )
      {
      for( x = 0; x < dim_[0]; x++ )
        {
        VOXEL(x, y, z) = -( VOXEL(x, y, z) > thresh );
        }
      }
    }
}

#define VOXEL2(x, y, z) data2[(x) + newxs * ( (y) + newys * (z) )]

void      TVoxelVolume::DimensionsChanged(int newxs, int newys, int newzs)
{
  if( newxs == dim_[0] && newys == dim_[1] && newzs == dim_[2] )
    {
    return;
    }
  else
    {
    unsigned char *data2 = new unsigned char[newxs * newys * newzs];

    long x, y, z;
    for( z = 0; z < newzs; z++ )
      {
      for( y = 0; y < newys; y++ )
        {
        for( x = 0; x < newxs; x++ )
          {
          if( x < dim_[0] && y < dim_[1] && z < dim_[2] )
            {
            VOXEL2(x, y, z) = VOXEL(x, y, z);
            }
          }
        }
      }

    delete data_;
    data_ = data2;
    dim_[0] = newxs; dim_[1] = newys; dim_[2] = newzs;
    size_ = dim_[0] * dim_[1] * dim_[2];
    return;
    }
}
