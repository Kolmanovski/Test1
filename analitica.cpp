#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

int main()
{
	setlocale(LC_ALL, "ru");
	double L = 1.0;
	size_t NX = 10001;
	double dx = L / (NX - 1);
	double x;
	double x_center = 0.5;
	double p, u, rho, E;

	const double pL = 1, rhoL = 1.0, uL = 0.0;
	const double pR = 0.1, rhoR = 0.125, uR = 0.0;

	double pStar = 0.30313;
	double uStar = 0.92745;
	double rhoStarL = 0.42632;
	double rhoStarR = 0.26557;

	const double gamma = 1.4;
	const double t = 0.25;

	std::ofstream fout("Test_1.plt");
	fout << "Variables = \"Position\", \"Pressure\", \"Density\", \"Velocity\", \"Energy\"" << std::endl;

	double Q = sqrt(0.5 * rhoR * ((gamma + 1) * pStar + (gamma - 1) * pR));
	double S_wave = Q / rhoR + uR;

	double x_wave = x_center + S_wave * t;
	double x_contact = x_center + uStar * t;
	
	const double aL = sqrt(gamma*pL/rhoL);
	const double aStar = aL * std::pow(pStar / pL, (gamma - 1) / 2 / gamma);
	
	double S_HL = uL - aL;
	double S_TL = uStar - aStar;

	double x_head = x_center + S_HL * t;
	double x_tail = x_center + S_TL * t;
	
	std::cout << "V_Shock = " << S_wave << std::endl;
	std::cout << "V_Head = " << S_HL << std::endl;
	std::cout << "V_Tail = " << S_TL << std::endl;
	std::cout << "aL = " << aL << std::endl;
	std::cout << "aL* = " << aStar << std::endl;
	std::cout << "X_Contact = " << x_contact << std::endl;
	std::cout << "X_Head = " << x_head << std::endl;
	std::cout << "X_Tail = " << x_tail << std::endl;

	for (size_t i = 0; i < NX; i++)
	{
		x = dx * i;
		if (x >= x_wave)
		{
			u = uR;
			p = pR;
			rho = rhoR;
		}
		if (x >= x_contact && x < x_wave)
		{
			u = uStar;
			p = pStar;
			rho = rhoStarR;
		}
		if (x >= x_tail && x < x_contact)
		{
			u = uStar;
			p = pStar;
			rho = rhoStarL;
		}
		if (x < x_tail && x>= x_head )
		{
			double _x = x - 0.5;
			double temp = 2 / (gamma + 1) + (gamma - 1) / (gamma + 1) / aL * (uL - _x / t);
			u = 2.0 / (gamma + 1) * (aL + uL * (gamma - 1) / 2.0 + _x/t);
			p = pL * std::pow(temp,2*gamma/(gamma - 1));
			rho = rhoL * std::pow(temp, 2/(gamma - 1));
		}
		if (x < x_head)
		{
			u = uL;
			p = pL;
			rho = rhoL;
		}
		
		E = p / (gamma - 1) / rho;

		fout << x << "\t" << p << "\t" << rho << "\t" << u << "\t" << E << std::endl;
	}

	return 0;
}