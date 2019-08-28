#include "Curve.h"



bool SRect::IsInside(const rvector& pt) const
{
	if (pt[0] < x1)
		return false;
	if (pt[0] > x2)
		return false;
	if (pt[1] < y1)
		return false;
	if (pt[1] > y2)
		return false;

	return true;
}

bool SRect::IsOverlapping(const SRect& rc) const
{
	return (x1 < rc.x2) && (x2 > rc.x1) && (y1 < rc.y2) && (y2 > rc.y1);
}

SQuadTreeNode::SQuadTreeNode()
{
	for (int i = 0; i < 4; ++i)
		n[i] = 0;
}

bool SQuadTreeNode::AddPoint(rvector &p, SQuadTreeDataType data)
{
	if ((p[0] < m_rc.x1) || (p[0] > m_rc.x2) || (p[1] < m_rc.y1) || (p[1] > m_rc.y2))
		return false;
	
	if (n[0] != 0)
	{
		for (int i = 0; i < 4; ++i)
			if (n[i]->AddPoint(p, data))
				return true;
	}
	
	pts.push_back(std::make_pair(p, data));

	if (pts.size() > maxpoints)
	{
		for(int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
			{
				double tol = 1e-6;
				n[i * 2 + j] = new SQuadTreeNode;
				n[i * 2 + j]->m_rc = SRect(
					m_rc.x1 + i*(m_rc.x2 - m_rc.x1) / 2.0 - tol,
					m_rc.y1 + j*(m_rc.y2 - m_rc.y1) / 2.0 - tol,
					m_rc.x1 + (i + 1)*(m_rc.x2 - m_rc.x1) / 2.0 + tol,
					m_rc.y1 + (j + 1)*(m_rc.y2 - m_rc.y1) / 2.0 + tol);
			}

		while (!pts.empty())
		{
			for (int i = 0; i < 4; ++i)
				if (n[i]->AddPoint(pts.back().first, pts.back().second))
					break;
			pts.pop_back();
		}
		pts.shrink_to_fit();
	}

	return true;
}

void SQuadTreeNode::Query(const SRect& rc, std::vector<SQuadTreeDataType>& out) const
{
	if (n[0] == 0)
	{
		for (auto &p : pts)
		{
			if(rc.IsInside(p.first))
				out.push_back(p.second);
		}
		return;
	}

	for (int i = 0; i < 4; ++i)
	{
		if(m_rc.IsOverlapping(rc))
			n[i]->Query(rc, out);
	}
}

void SQuadTreeNode::Clear()
{
	if (n[0] != 0)
	{
		for (int i = 0; i < 4; ++i)
		{
			n[i]->Clear();
			delete n[i];
			n[i] = 0;
		}
	}
}

double SGrid::drand(double x, double X)
{
	double d = m_udist(m_rng);

	return d*(X - x) + x;
}

SGrid::SGrid() : m_udist(0.0 + 1e-8, 1.0 - 1e-8)
{
	double seedval = time(NULL);
	m_rng.seed(seedval);
}

SGrid::~SGrid()
{
	m_quadtree.Clear();
}

void SGrid::Create(int N, const SRect& rc, bool bSquareGrid)
{
	
	
	p.clear();
	m_quadtree.Clear();

	m_quadtree.m_rc = rc;

	if (!bSquareGrid)
	{
		for (int i = 0; i < N; ++i)
		{
			p.push_back({ drand(rc.x1, rc.x2), drand(rc.y1, rc.y2) });
			m_quadtree.AddPoint(p[i], i);
		}
	}
	else
	{
		double d = (rc.x2 - rc.x1) / (rc.y2 - rc.y1);
		double nx = sqrt(N * d);
		double s = (rc.x2 - rc.x1) / nx;

		for(int x = 0; x < nx; ++x)
			for (int y = 0; y < nx / d; ++y)
			{
				p.push_back({ 
					rc.x1 + x*s + drand(-s*0.02, s*0.02),
					rc.y1 + y*s + drand(-s*0.02, s*0.02)
				});
				m_quadtree.AddPoint(p.back(), p.size()-1);
			}
	}
}

void SGrid::Create(std::string fname)
{
	p.clear();
	m_quadtree.Clear();

	std::ifstream f;
	f.setf(std::ios::scientific);
	f.open(fname);

	while (f.good())
	{
		double x, y;
		f >> x;
		f >> y;
		p.push_back({ x, y });
	}

	f.close();

	double x1(1e10), x2(-1e10), y1(1e10), y2(-1e10);
	for (auto &pt : p)
	{
		x1 = std::min<double>(x1, pt[0]);
		x2 = std::max<double>(x2, pt[0]);
		y1 = std::min<double>(y1, pt[1]);
		y2 = std::max<double>(y2, pt[1]);
	}

	m_quadtree.m_rc = SRect(x1, y1, x2, y2);

	for(size_t i = 0; i < p.size(); ++i)
		m_quadtree.AddPoint(p[i], i);
}

void SGrid::QueryRect(const SRect& rc, std::vector<int>& out) const
{
	m_quadtree.Query(rc, out);
}

// barycentric coordinates of i in the triangle i1, i2, i3 are stored in l1, l2, l3 
void SGrid::BarycentricCoords(int i1, int i2, int i3, int i, double &l1, double &l2, double &l3) const
{
	double A = (p[i2][1] - p[i3][1])*(p[i][0] - p[i3][0]) + (p[i3][0] - p[i2][0])*(p[i][1] - p[i3][1]);
	double B = (p[i2][1] - p[i3][1])*(p[i1][0] - p[i3][0]) + (p[i3][0] - p[i2][0])*(p[i1][1] - p[i3][1]);
	double C = (p[i3][1] - p[i1][1])*(p[i][0] - p[i3][0]) + (p[i1][0] - p[i3][0])*(p[i][1] - p[i3][1]);
	double D = (p[i2][1] - p[i3][1])*(p[i1][0] - p[i3][0]) + (p[i3][0] - p[i2][0])*(p[i1][1] - p[i3][1]);

	l1 = A / B;
	l2 = C / D;
	l3 = 1.0 - l1 - l2;
}

// returns the grid points in the triangle i1, i2, i3 except for vertices
void SGrid::QueryTriangle(int i1, int i2, int i3, std::vector<int>& out) const
{
	double x1 = std::min<double>(std::min<double>(p[i1][0], p[i2][0]), p[i3][0]);
	double x2 = std::max<double>(std::max<double>(p[i1][0], p[i2][0]), p[i3][0]);
	double y1 = std::min<double>(std::min<double>(p[i1][1], p[i2][1]), p[i3][1]);
	double y2 = std::max<double>(std::max<double>(p[i1][1], p[i2][1]), p[i3][1]);

	std::vector<int> v;
	QueryRect(SRect(x1, y1, x2, y2), v);
	double l1, l2, l3;
	double tol = 1e-6;
	for (auto &i : v)
	{
		BarycentricCoords(i1, i2, i3, i, l1, l2, l3);
		if( (l1 > tol) && (l2 > tol) && (l3 > tol) )
			out.push_back(i);
	}
}

// if p1, A1, ..., Ak, p3 is the convex hull of all the grid points in the triangle p1, p2, p3 except for p2
// then ReduceTriangle returns A1, ..., Ak
std::vector<int> SGrid::ReduceTriangle(int p1, int p2, int p3)
{
	std::vector<int> out;
		
	std::map<double, int > ptsm; 
		
	// find all points in the triangle
	std::vector<int> trpts;
	QueryTriangle(p1, p2, p3, trpts);
	double l1, l2, l3;
	for (auto &i : trpts)
	{
		BarycentricCoords(p1, p2, p3, i, l1, l2, l3);
		ptsm.insert(std::make_pair(l3/(l1+l3), i));
	}

	if (ptsm.empty())
		return out;

	std::list<std::pair<int, bool>> pts; // <pt, toCheck>
	std::vector<std::list<std::pair<int, bool>>::iterator> toCheck;
	for (auto &p : ptsm)
	{
		pts.push_back(std::make_pair(p.second, true));
		toCheck.push_back(--pts.end());
	}

	std::list<std::pair<int, bool>>::iterator it, it1, it2;
	while (!toCheck.empty())
	{
		it = toCheck.back();
		toCheck.pop_back();
		it->second = false;

		int i1, i2;

		it2 = it;
		++it2;
		if (it2 != pts.end())
			i2 = it2->first;
		else
			i2 = p3;

		if (it != pts.begin())
		{
			it1 = it;
			--it1;
			i1 = it1->first;
		}
		else
		{
			i1 = p1;
			it1 = pts.end();
		}

		BarycentricCoords(i1, p2, i2, it->first, l1, l2, l3);
		if ((l1 < 0) || (l2 < 0) || (l3 < 0))
		{
			pts.erase(it);
			if (it1 != pts.end())
			{
				if (it1->second == false)
					toCheck.push_back(it1);
				it1->second = true;
			}
			if (it2 != pts.end())
			{
				if (it2->second == false)
					toCheck.push_back(it2);
				it2->second = true;
			}
		}
	}

	for (auto &p : pts)
		out.push_back(p.first);

	return out;
}

void SGrid::Save(std::string fname)
{
	std::fstream f;
	f.open(fname, std::fstream::out | std::fstream::trunc);

	for (auto &pt : p)
	{
		f << pt[0] << " " << pt[1] << std::endl;
	}

	f.close();
}

void SCurve::Step()
{
	CalcNextStep();
	m_pts = m_nextpts;
	m_nextStepReady = false;
	m_rpts = m_nextrpts;
}

void SCurve::CalcNextStep()
{
	static int counter = 0;
	counter++;

	if (m_nextStepReady)
		return;

	m_nextStepReady = true;

	if (m_pts.size() < 3)
	{
		m_nextpts.clear();
		return;
	}

	struct SNail
	{
		int p; // point
		double a; // angle
		bool toCheck;
		SNail(int p_, double a_, bool toCheck_) : p(p_), a(a_), toCheck(toCheck_) {}
	};

	std::list<SNail> N; // list of nails
	std::vector<std::list<SNail>::iterator> toCheck;
	for (auto &i : m_pts)
	{
		N.push_back(SNail(i, 0, true));
		toCheck.push_back(--N.end());
	}

	double Pi = 3.1415926535;
	while (true)
	{
		if (toCheck.empty())
			break;

		std::list<SNail>::iterator it2 = toCheck.back(); // nail
		toCheck.pop_back();
		it2->toCheck = false;

		if (fabs(it2->a) > Pi)
			continue;
		
		std::list<SNail>::iterator it1 = it2; // prev nail
		if (it1 == N.begin())
			it1 = --N.end();
		else
			--it1;
		std::list<SNail>::iterator it3 = it2; // next nail
		++it3;
		if (it3 == N.end())
			it3 = N.begin();

		if ((it2->p == it1->p) || (it2->p == it3->p)) //zero length edge
		{
			N.clear();
			break;
		}

		std::vector<int> v = m_grid.ReduceTriangle(it1->p, it2->p, it3->p);
		
		// update angles
		{
			rvector vnew;
			if (!v.empty())
				vnew = m_grid.p[v[0]] - m_grid.p[it1->p];
			else
				vnew = m_grid.p[it3->p] - m_grid.p[it1->p];
			rvector vold = m_grid.p[it2->p] - m_grid.p[it1->p];
			double s = VectorProduct(vold, vnew);
			double a = acos(DotProduct(vnew, vold) / (Norm(vnew)*Norm(vold)));
			if (s > 0)
				it1->a += a;
			else
				it1->a -= a;
		}
		{
			rvector vnew;
			if(!v.empty())
				vnew = m_grid.p[v.back()] - m_grid.p[it3->p];
			else
				vnew = m_grid.p[it1->p] - m_grid.p[it3->p];
			rvector vold = m_grid.p[it2->p] - m_grid.p[it3->p];
			double s = VectorProduct(vnew, vold);
			double a = acos(DotProduct(vnew, vold) / (Norm(vnew)*Norm(vold)));
			if (s > 0)
				it3->a += a;
			else
				it3->a -= a;
		}
		if (!it1->toCheck)
		{
			it1->toCheck = true;
			toCheck.push_back(it1);
		}
		if (!it3->toCheck)
		{
			it3->toCheck = true;
			toCheck.push_back(it3);
		}

		// erase it2
		N.erase(it2);
		
		// insert
		std::list<SNail>::iterator itInsertBefore(it3);
		for (int i = v.size() - 1; i >= 0; --i)
		{
			SNail n(v[i], 0, false);
			rvector v1, v2;
			if (i > 0)
				v1 = m_grid.p[v[i - 1]] - m_grid.p[v[i]];
			else
				v1 = m_grid.p[it1->p] - m_grid.p[v[i]];
			if (i + 1 < v.size())
				v2 = m_grid.p[v[i + 1]] - m_grid.p[v[i]];
			else
				v2 = m_grid.p[it3->p] - m_grid.p[v[i]];

			double s = VectorProduct(v2, v1);
			double a = acos(DotProduct(v2, v1) / (Norm(v2)*Norm(v1)));
			a = 2 * Pi - a;
			if (s > 0)
				n.a = a;
			else
				n.a = -a;

			itInsertBefore = N.insert(itInsertBefore, n);
		}

		// if we inserted nothing
		if (v.empty())
		{
			if (it1->p == it3->p)
			{
				it1->a += it3->a;
				if (it3->toCheck)
				{
					for (auto &c : toCheck)
					{
						if (c == it3)
						{
							std::swap(toCheck.back(), c);
							toCheck.pop_back();
							break;
						}
					}
				}
				N.erase(it3);
			}
		}
	}

	m_nextpts.clear();
	m_nextrpts.clear();
	for (auto &n : N)
	{
		m_nextpts.push_back(n.p);
		m_nextrpts.push_back(m_grid.p[n.p]);
	}
}

void SCurve::SaveGrid(std::string fname)
{
	m_grid.Save(fname);
}

// loads from pts
void SCurve::Load(const rrvector& pts, std::string gridfname, int nGridPoints, bool bSquareGrid)
{
	double x1(1e10), x2(-1e10), y1(1e10), y2(-1e10);
	for (auto &p : pts)
	{
		x1 = std::min<double>(x1, p[0]);
		x2 = std::max<double>(x2, p[0]);
		y1 = std::min<double>(y1, p[1]);
		y2 = std::max<double>(y2, p[1]);
	}

	x1 -= (x2 - x1)*0.2;
	x2 += (x2 - x1)*0.2;
	y1 -= (y2 - y1)*0.2;
	y2 += (y2 - y1)*0.2;

	if(gridfname == "")
		m_grid.Create(nGridPoints, SRect(x1, y1, x2, y2), bSquareGrid);
	else
		m_grid.Create(gridfname);


	std::vector<int> tmp;
	for (auto &p : pts)
	{
		std::vector<int> v;
		double rx1 = (x2 - x1) / sqrt(nGridPoints);
		double ry1 = (y2 - y1) / sqrt(nGridPoints);
		while (v.empty())
		{
			m_grid.QueryRect(SRect(p[0] - rx1, p[1] - ry1, p[0] + rx1, p[1] + ry1), v);
			rx1 *= 2;
			ry1 *= 2;
		}

		int besti = v[0];
		for (auto &i : v)
		{
			double d1 = Norm(p - m_grid.p[besti]);
			double d2 = Norm(p - m_grid.p[i]);
			if(d1 > d2)
				besti = i;
		}
		tmp.push_back(besti);
	}

	// clear duplicates
	m_pts.clear();
	while (tmp.front() == tmp.back())
		tmp.pop_back();
	for (auto &i : tmp)
	{
		if (!m_pts.empty())
			if (m_pts.back() == i)
				continue;
		m_pts.push_back(i);
	}

	m_rpts.clear();
	for (auto &p : m_pts)
		m_rpts.push_back(m_grid.p[p]);

	m_nextStepReady = false;
}
