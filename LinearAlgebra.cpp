#include "LinearAlgebra.h"

void Normalise(rvector& V)
{
	double N = Norm(V);

	for (size_t i = 0; i < V.size(); ++i)
		V[i] /= N;
}

double SqNorm(const rvector& V)
{
	double N = 0;
	for (size_t i = 0; i < V.size(); ++i)
		N += V[i] * V[i];

	return N;
}

double Norm(const rvector& V)
{
	double N = 0;
	for (size_t i = 0; i < V.size(); ++i)
		N += V[i] * V[i];
	N = sqrt(N);

	return N;
}

double Norm(rrvector& V)
{
	double N = 0;
	for (size_t i = 0; i < V.size(); ++i)
	{
		for (size_t j = 0; j < V[i].size(); ++j)
			N += V[i][j] * V[i][j];
	}
	N = sqrt(N);

	return N;
}

double SqNorm(std::vector<rvector >& V)
{
	double N = 0;
	for (size_t i = 0; i < V.size(); ++i)
	{
		for (size_t j = 0; j < V[i].size(); ++j)
			N += V[i][j] * V[i][j];
	}

	return N;
}

double DotProduct(const rvector& V1, const rvector& V2)
{
	double r = 0;
	for (size_t i = 0; i < V1.size(); ++i)
		r += V1[i] * V2[i];

	return r;
}

double VectorProduct(rvector& V1, rvector& V2)
{
	return V1[0] * V2[1] - V1[1] * V2[0];
}

void CrossProduct(rvector& V1, rvector& V2, rvector& out)
{
	out[0] = V1[1] * V2[2] - V2[1] * V1[2];
	out[1] = V1[2] * V2[0] - V2[2] * V1[0];
	out[2] = V1[0] * V2[1] - V2[0] * V1[1];
}

void Transpose(rrvector & V)
{
	rrvector out(V[0].size(), rvector(V.size()));

	for (size_t i = 0; i < out.size(); ++i)
		for (size_t j = 0; j < out[0].size(); ++j)
			out[i][j] = V[j][i];

	V = out;
}

double DistPtSegment(const rvector &pt, const rvector& p1, const rvector &p2)
{
	rvector vLine = p2 - p1;
	rvector vToPoint = pt - p1;

	double sc = DotProduct(vLine, vToPoint);
	rvector vToClosest = vLine * sc / SqNorm(vLine);

	if( (Norm(vToClosest) < Norm(vLine)) && (sc >= 0) )
		return Norm(pt - (p1 + vToClosest));
	
	double d1 = Norm(pt - p1);
	double d2 = Norm(pt - p2);

	return (d1 < d2) ? d1 : d2;
}


bool IntersectLines(const rvector& p1, const rvector& p2, const rvector& p3, const rvector& p4, rvector& out)
{
	rvector v1 = p2 - p1;
	rvector v3 = p4 - p3;

	if (fabs((v1[1] * v3[0] - v3[1] * v1[0])) < 1e-5)
		return false;

	double x = (p3[1] * v3[0] + v3[1] * p1[0] - v3[1] * p3[0] - p1[1] * v3[0]) / (v1[1] * v3[0] - v3[1] * v1[0]);

	out = p1 + v1*x;

	return true;
}

void CircumCircle(const rvector& A, const rvector& B, const rvector& C, rvector& center, double& r)
{
	double a = Norm(B - C);
	double b = Norm(A - C);
	double c = Norm(A - B);
	double x = a*a*(b*b + c*c - a*a);
	double y = b*b*(a*a + c*c - b*b);
	double z = c*c*(b*b + a*a - c*c);

	center = (A * x + B * y + C * z) / (x + y + z);
	r = Norm(center - A);
}

bool IntCircles(const rvector& c1, const rvector& c2, double r1, double r2, rvector& out1, rvector& out2)
{
	double x0 = c1[0];
	double y0 = c1[1];
	double x1 = c2[0];
	double y1 = c2[1];

	double d = sqrt( (x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) );
	double a = (r1*r1 - r2*r2 + d*d) / (2 * d);

	double h = (r1*r1 - a*a);
	if (h < 0)
		return false;
	h = sqrt(h);
	
	double x2 = x0 + a*(x1 - x0) / d;
	double y2 = y0 + a*(y1 - y0) / d;

	out1[0] = x2 + h*(y1 - y0) / d;
	out1[1] = y2 - h*(x1 - x0) / d;

	out2[0] = x2 - h*(y1 - y0) / d;
	out2[1] = y2 + h*(x1 - x0) / d;

	return true;
}

double Area(const rrvector& pts)
{
	double A = 0;
	for (size_t i = 0; i < pts.size(); ++i)
	{
		size_t i1 = (i + 1) % pts.size();
		A += (pts[i1][0] - pts[i][0])*(pts[i1][1] + pts[i][1]);
	}
	return A / 2.0;
}

double Length(const rrvector& pts)
{
	double L = 0;
	for (size_t i = 0; i < pts.size(); ++i)
	{
		size_t i1 = (i + 1) % pts.size();
		L += Norm(pts[i] - pts[i1]);
	}
	return L;
}