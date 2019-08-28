#pragma once

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

typedef std::vector<double> rvector;
typedef std::vector<rvector> rrvector;
typedef std::vector<rrvector> rrrvector;
typedef std::vector<rrrvector> rrrrvector;

double DotProduct(const rvector& V1, const rvector& V2);
double VectorProduct(rvector& V1, rvector& V2);
void CrossProduct(rvector& V1, rvector& V2, rvector& out);
void Normalise(rvector& V);
double Norm(const rvector& V);
double SqNorm(const rvector& V);
double Norm(rrvector& V);
double SqNorm(std::vector<rvector >& V);
void Transpose(rrvector & V);
double DistPtSegment(const rvector &pt, const rvector& p1, const rvector &p2);
bool IntersectLines(const rvector& p1, const rvector& p2, const rvector& p3, const rvector& p4, rvector& out);
void CircumCircle(const rvector& A, const rvector& B, const rvector& C, rvector& center, double& r);
bool IntCircles(const rvector& c1, const rvector& c2, double r1, double r2, rvector& out1, rvector& out2);
double Area(const rrvector& pts);
double Length(const rrvector& pts);

inline void operator *=(rvector& v, double d)
{
	for (size_t i = 0; i < v.size(); ++i)
		v[i] *= d;
}

inline void operator /=(rvector& v, double d)
{
	for (size_t i = 0; i < v.size(); ++i)
		v[i] /= d;
}

inline rvector operator *(const rvector& v, double d)
{
	rvector r(v.size());
	for (size_t i = 0; i < r.size(); ++i)
		r[i] = v[i] * d;
	return r;
}

inline rvector operator /(rvector& v, double d)
{
	rvector r(v.size());
	for (size_t i = 0; i < r.size(); ++i)
		r[i] = v[i] / d;
	return r;
}

inline rvector operator +(const rvector& v1, const rvector& v2)
{
	rvector r(v1.size());
	for (size_t i = 0; i < r.size(); ++i)
		r[i] = v1[i] + v2[i];
	return r;
}

inline rvector operator -(const rvector& v1, const rvector& v2)
{
	rvector r(v1.size());
	for (size_t i = 0; i < r.size(); ++i)
		r[i] = v1[i] - v2[i];
	return r;
}

inline std::vector<rvector > operator*(rrvector & M1, rrvector& M2)
{
	rrvector r(M1.size(), rvector(M2[0].size(), 0));
	for (size_t i = 0; i < M1.size(); ++i)
		for (size_t j = 0; j < M2[0].size(); ++j)
			for (size_t k = 0; k < M1[0].size(); ++k)
				r[i][j] += M1[i][k] * M2[k][j];

	return r;
}

inline std::vector<rvector > operator+(const rrvector& M1, const rrvector& M2)
{
	std::vector<rvector > r(M1.size(), rvector(M2[0].size(), 0));
	for (size_t i = 0; i < r.size(); ++i)
		for (size_t j = 0; j < r[0].size(); ++j)
			r[i][j] = M1[i][j] + M2[i][j];

	return r;
}

inline rrvector operator-(const rrvector& M1, const rrvector& M2)
{
	std::vector<rvector > r(M1.size(), rvector(M2[0].size(), 0));
	for (size_t i = 0; i < r.size(); ++i)
		for (size_t j = 0; j < r[0].size(); ++j)
			r[i][j] = M1[i][j] - M2[i][j];

	return r;
}

inline std::vector<rvector > operator*(rrvector& M, double d)
{
	std::vector<rvector > r(M);
	for (size_t i = 0; i < r.size(); ++i)
		r[i] = r[i] * d;

	return r;
}