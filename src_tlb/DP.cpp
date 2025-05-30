#include "DP.h"
#include <limits>
#include <cmath>
#include <iostream>
#include <cstdio>
 
//needs to receive the vector of capacities h_i
double DP::GetCostDP(std::vector<int>& h, int Qtot, std::vector<int>& y) 
{
    int n = h.size();
	std::vector<int> hh = h;
	std::vector< std::vector<double> > dp(n + 1, std::vector<double>(Qtot + 1, 99999));
	std::vector< std::vector<int> > x(n + 1, std::vector<int>(Qtot + 1));

	for (int q = 0; q <= Qtot; q++)
	{
		dp[0][q] = std::pow(q - (int)std::round(0.5*hh[ 0 ]) , 2);
		x[0][q] = q;
	}
	// Dynamic programming
	for (int i = 1; i < n; i++)
		for (int q = 0; q <= Qtot; q++)
			for (int x_val = 0; x_val <= std::min(q, (int)std::round(0.5*hh[i])); x_val++)
			{
				double cost = std::pow(x_val - (int)std::round(0.5*hh[i]), 2) + dp[i - 1][q - x_val];
				if (cost < dp[i][q])
				{
					dp[i][q] = cost;
					x[i][q] = x_val;
				}
			}
	//for (int i = 0; i < n; i++)
	//{
		//printf("i:%2d dp:", i);
		//for (int q = 0; q <= Qtot; q++)
			//printf("%3d ", (int) (dp[i][q]+0.0001) );
		//printf("\n");
	//}
	//for (int i = 0; i < n; i++)
	//{
		//printf("i:%2d x:", i);
		//for (int q = 0; q <= Qtot; q++)
			//printf("%3d ", x[i][q]);
		//printf("\n");
	//}

	int qmin = Qtot;
	double dpmin = dp[n-1][Qtot];
	int xmin = x[n-1][Qtot];
	for (int q = 0; q <= Qtot; q++)
		if(dp[n-1][q] < dpmin)
		{
			dpmin = dp[n-1][q];
			xmin = x[n-1][q];
			qmin = q;
		}

	// Reconstruction of y and x
	y.resize(n);
	y[n-1] = std::max(0, std::min(xmin, h[n-1]));
	std::min(xmin,h[n-1]);
	int q = qmin - xmin;
	for(int i=n-2;i>=0;i--)
	{
		y[i] = std::max(0, std::min(x[i][q],h[i]));
		q -= x[i][q];
	}

	printf("Best:%d Qtot:%d qmin:%d\n", (int)(dpmin + 0.0001), Qtot, qmin);
	for (int i = 0; i < n; i++)
		printf("i:%d h:%d hl:%d y:%d\n",i, h[i], hh[i], y[i]);
	
	return dpmin;
}
