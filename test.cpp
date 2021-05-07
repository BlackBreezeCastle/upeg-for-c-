#include"peg.h"
#include<Windows.h>
void main()
{
	double r0=6371000;
	double h=r0+180000;
	Vector3 r(0,0,r0);
	Vector3 v(-1000,0,500);
	pegas peg;
	peg.add_stage(1000,100,7000,465);
	peg.set_target_orbit(0,0,h,10900,0);
	peg.init_peg(r,v,1000,0,398600446148608.0);
	for(int i=0;i<100;i++)
	{
		peg.update(0,1000,r,v);
		//printf("tgo:%lf\n",peg.time_to_go());
	}
	peg.update(0,1000,r,v);
	system("pause");
}