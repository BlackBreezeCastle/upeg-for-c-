#include"peg.h"

void main()
{
	double h=6371000+180000;
	Vector3 r(0,0,h);
	Vector3 v(-7900,0,0);
	pegas peg;
	peg.add_stage(1000,100,2000,400);
	peg.set_target_orbit(0,0,h,8000,0);
	peg.init_peg(r,v,10000,0,398600446148608.0);
	for(int i=0;i<100;i++)
	{
		peg.update(0,1000,r,v);
	}
	peg.update(0,1000,r,v);
}
