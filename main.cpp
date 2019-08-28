#include <time.h>

#include "Curve.h"


int main()
{
	// starting length
	double StaringLength = 0;

	// default
	std::string gridfname = "";
	std::string savegridfname = "obstacles.txt";
	std::string curvefname = "curve.txt";
	int nGridPoints = 10000;
	double dLengthLimit = 0.7;
		
	// load settings
	{
		
		std::ifstream fsettings;
		fsettings.open("settings.txt");
		char tmp[1000];
		std::string str;

		// obstacles file name, leave "" for random obstacles
		fsettings.getline(tmp, 1000);
		str = std::string(tmp);
		gridfname = str.substr(str.find('=')+2, str.length());
		// number of random obstacles (ignored, if the obstacles are loaded from a file)
		fsettings.getline(tmp, 1000);
		str = std::string(tmp);
		nGridPoints = atoi(str.substr(str.find('=') + 2, str.length()).c_str());
		// save obstacles to
		fsettings.getline(tmp, 1000);
		str = std::string(tmp);
		savegridfname = str.substr(str.find('=') + 2, str.length());
		// curve file name
		fsettings.getline(tmp, 1000);
		str = std::string(tmp);
		curvefname = str.substr(str.find('=') + 2, str.length());
		// go until length
		fsettings.getline(tmp, 1000);
		str = std::string(tmp);
		dLengthLimit = atof(str.substr(str.find('=') + 2, str.length()).c_str());
	
		fsettings.close();
	}

	SCurve curve;

	time_t timeStart;
	time(&timeStart);
	// load curve
	{
		std::ifstream fcurve;
		fcurve.setf(std::ios::scientific);
		fcurve.open(curvefname);

		rrvector pts;

		while (fcurve.good())
		{
			double x, y;
			fcurve >> x;
			fcurve >> y;
			pts.push_back({ x, y });
		}

		StaringLength = Length(pts);
		curve.Load(pts, gridfname, nGridPoints, false);

		fcurve.close();
	}

	// save grid
	if(!savegridfname.empty())
	{
		curve.SaveGrid(savegridfname);
	}

	// steps
	for (int n = 1; true; ++n)
	{
		curve.Step();

		std::cout << "Step: " << n << std::endl;

		double L = Length(curve.m_rpts);

		// steps output
		if (L < dLengthLimit*StaringLength)
		{
			{
				std::fstream fout;
				std::stringstream ss;
				ss << "result.txt";
				fout.open(ss.str(), std::fstream::out | std::fstream::trunc);

				for (auto &p : curve.m_rpts)
					fout << p[0] << " " << p[1] << std::endl;

				fout.close();
			}
			break;
		}
	}
	
	return 0;
}