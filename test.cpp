#include"peg.h"
#include<Windows.h>
void main()
{
	double r0=6371000;
	double h=r0+180000;
	Vector3 r(0,0,r0);
	//Vector3 v(-1000,0,800);
	Vector3 v(-200, 0, 0);
	double max_f = 0;
	double min_acc = 0;
	const double weight1 = 120;
	const double weight2 = 15;
	const double thrust1 = 180 * 10;
	const double thrust2 = 18*10;
	for (int j = 0;j < 500;j++)
	{
		double dw = (double)j*0.1;
		pegas peg;
		peg.add_stage(weight1, weight2 + (weight1-weight2-dw)*0.04 + dw, thrust1, 310);
		peg.add_stage(weight2 + dw, 5, thrust2, 465);
		peg.set_target_orbit(0, 0, h, 7900+520, 0);
		peg.init_peg(r, v, weight1, 0, 398600446148608.0);
		for (int i = 0;i < 100;i++)
		{
			peg.update(0, weight1, r, v);
		}
		double acc = thrust2 / (weight2 + dw) / g0;
		double final_weight = peg.final_weight()-(weight2 + dw - peg.final_weight())*0.08-0.24;
		if (final_weight > max_f)
		{
			max_f = final_weight;
			min_acc = acc;
			printf("i:%lf finalWeight:%lf,tgo:%lf dv:%lf acc:%lf k:%lf\n", 
				(double)(dw + weight2), final_weight, peg.time_to_go(), peg.dv_to_go(), acc,final_weight/weight1);
		}
		else
		{
			break;
		}
	}
	system("pause");
}