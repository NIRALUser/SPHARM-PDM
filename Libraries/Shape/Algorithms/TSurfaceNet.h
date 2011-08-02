// -*- C++ -*-
/****************************************************************************
 *
 *  class TSurfaceNet
 *
 ****************************************************************************
 *
 *  File         :  TSurfaceNet.h
 *  Type         :
 *  Purpose      :  stores the surface net data structure
 *
 ****************************************************************************
 *
 *  started      :  1  Jan 95     Robert Schweikert
 *  last change  :  1  Aug 2005   Martin Styner
 *
 ****************************************************************************/

#ifndef TSURFACENET_H
#define TSURFACENET_H

#include <iostream>

#include "TVoxelVolume.h"

class TVertex
{
public:
  TVertex();
  TVertex(int wh[3], int count, int neighb[14]);

  TPoint3D GetCoords3D()
  {
    return TPoint3D(wh_);
  }

  TPointFloat3D GetCoordsFloat3D()
  {
    return TPointFloat3D(wh_);
  }

  friend int operator<(const TVertex& v1, const TVertex& v2)
  {
    return 0;
  }

  friend int operator>(const TVertex& v1, const TVertex& v2)
  {
    return 0;
  }

  friend int operator==(const TVertex& v1, const TVertex& v2);

  friend std::ostream & operator<<(std::ostream& os, const TVertex& v);

  friend float vdist(const TVertex& v1, const TVertex& v2)
  {
    return sqrt( (double) (v1.wh_[0] - v2.wh_[0]) * (v1.wh_[0] - v2.wh_[0])
                 + (v1.wh_[1] - v2.wh_[1]) * (v1.wh_[1] - v2.wh_[1])
                 + (v1.wh_[2] - v2.wh_[2]) * (v1.wh_[2] - v2.wh_[2]) );
  }

  int wh_[3];
  int count_;
  int neighb_[14];
};

class TVertexArray
{
public:
  TVertexArray(int sz = 10)
  {
    init(0, sz);
  }

  TVertexArray(TVertex *ar, int sz)
  {
    init(ar, sz);
  }

  TVertexArray(const TVertexArray& iA)
  {
    init(iA.ia_, iA.size_);
  }

  ~TVertexArray()
  {
    delete [] ia_;
  }

  TVertexArray & operator=(const TVertexArray &);

  int size()
  {
    return size_;
  }
  int length()
  {
    return size_;
  }

  TVertex & operator[](int ix)
  {
    return ia_[ix];
  }

  const TVertex & operator[](int ix) const
  {
    return ia_[ix];
  }

  int find(TVertex) const;

  void print(std::ostream & = std::cout) const;

  friend std::ostream & operator<<(std::ostream& os, const TVertexArray& ar);

protected:
  void init(const TVertex *, int);

  int      size_;
  TVertex* ia_;
};
class TVertexItem
{
public:
  TVertexItem(TVertex val, TVertexItem* next = 0) : val_(val), next_(next)
  {
  }

  TVertex      val_;
  TVertexItem* next_;
};

class TVertexList
{
public:
  TVertexList() : last_(0)
  {
  }

  ~TVertexList();

  TVertexItem * succ(TVertexItem* it) const
  {
    return it == last_ ? 0 : it->next_;
  }

  TVertexItem * cyclic_succ(TVertexItem* it) const
  {
    return it->next_;
  }

  TVertexItem * first() const
  {
    return last_ == 0 ? 0 : last_->next_;
  }

  TVertexItem * last() const
  {
    return last_;
  }

  TVertexItem * search(TVertex val) const;

  int size()   const;

  int length() const
  {
    return size();
  }

  int empty()  const
  {
    return last_ == 0 ? 1 : 0;
  }

  int contains(TVertex val) const;

  void insert(TVertex val);

  void append(TVertex val);

  TVertex pop();

  void push(TVertex val)
  {
    insert(val);
  };
  void insert_after(TVertex val, TVertexItem* it);

  void remove(TVertexItem* it);

  void print(std::ostream & = std::cout) const;

  friend std::ostream & operator<<(std::ostream& os, const TVertexList& li);

private:
  TVertexItem* last_;
};

#define forall_vertexlistitems(a, l) for( a = (l).first(); a; a = (l).succ(a) )

class TNetExtractor
{
public:
  TNetExtractor(const TVoxelVolume* vol) : vol_(vol)
  {
  }

  TVertexList * extractnet();

  struct t1entry
    {
    int nbc;
    int count;
    int nbcs[4];
    };

  struct t2entry
    {
    int nbc;
    int count;
    int dir[7];
    int ne[7];
    };

  struct numentry
    {
    int nums[8];
    };

  enum { W, E, S, N, B, T };
private:
  int direct_translate(int dir, int ne);

  int diag_translate(int dir1, int dir2, int ne);

  void encode(int nx, int ny, int nz, int& number, unsigned char* nbc, numentry* num, int z);

  const TVoxelVolume* vol_;                   // pointer to object volume
};

class TSurfaceNet
{
public:
  TSurfaceNet(TVoxelVolume& vol, int thresh);
  ~TSurfaceNet()
  {
    delete va_;
  }

  int length() const
  {
    return va_->length();
  }

  int size() const
  {
    return va_->length();
  }

  int nface() const;

  // gets first vertex that fits
  int GetIndex(const TPoint3D& pt) const;

  // given index, get corresponding vertex
  TVertex GetVertex(int index) const
  {
    return (*va_)[index];
  }

  // distance between two vertices
  float dist(int inx1, int inx2)
  {
    return vdist( (*va_)[inx1], (*va_)[inx2] );
  }

  // output
  void read(std::istream& is);

  void print(std::ostream& os) const;

  friend std::ostream & operator<<(std::ostream& os, const TSurfaceNet& sn);

private:
  // this creates a TNetExtractor object to do the dirty work
  void vol2net(TVoxelVolume& vol);

  TVertexArray* va_;                      // vertex data
  unsigned char thresh_;                  // threshold value
};

#endif
