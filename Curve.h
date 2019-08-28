#pragma once

#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <random>

#include "LinearAlgebra.h"

typedef int SQuadTreeDataType;

struct SRect
{
	double x1, y1, x2, y2; // x1 < x2, y1 < y2
	SRect() {}
	SRect(double x1_, double y1_, double x2_, double y2_) : x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}

	bool IsInside(const rvector& pt) const;
	bool IsOverlapping(const SRect& rc) const;
};

// quad tree to store grid points
struct SQuadTreeNode
{
	const int maxpoints = 20;
	SRect m_rc; // rect
	SQuadTreeNode *n[4]; // children
	std::vector<std::pair<rvector, SQuadTreeDataType>> pts; // points

	SQuadTreeNode();

	bool AddPoint(rvector &p, SQuadTreeDataType data);
	void Query(const SRect& rc, std::vector<SQuadTreeDataType>& out) const;
	void Clear();
};

struct SGrid
{
	std::uniform_real_distribution<> m_udist;
	std::mt19937 m_rng;

	SQuadTreeNode m_quadtree;
	rrvector p; // points
	
	void Create(int N, const SRect& rc, bool bSquareGrid);
	void Create(std::string fname);
	void QueryRect(const SRect& rc, std::vector<int>& out) const; // return points in the rect
	void QueryTriangle(int i1, int i2, int i3, std::vector<int>& out) const; // return points in the triangle, not including the vertices
	void BarycentricCoords(int i1, int i2, int i3, int i, double &l1, double &l2, double &l3) const;
	void Save(std::string fname);

	double drand(double x, double X);

	std::vector<int> ReduceTriangle(int p1, int p2, int p3);
	SGrid();
	~SGrid();
};


struct SCurve
{
	SGrid m_grid;
	std::vector<int> m_pts;
	std::vector<int> m_nextpts;
	rrvector  m_rpts;
	rrvector  m_nextrpts;
	
	bool m_nextStepReady;
	
	SCurve() : m_nextStepReady(false) {}

	void Step();
	void CalcNextStep();
	void Load(const rrvector& pts, std::string gridfname, int nGridPoints, bool bSquareGrid);
	void SaveGrid(std::string fname);
};
