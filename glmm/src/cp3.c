#include "myheader.h"
void cp3(double *etain, int *neta, int *typein, double *cpout)
{
    int leneta=neta[0];
    int type = typein[0];
    int i;

    for(i=0;i<leneta;++i)
    {
	double eta = etain[i];
	switch (type) {
	case 1:
	    cpout[i]= 1/(1+exp(-eta));
	    break;
        case 2:
            cpout[i]= exp(eta);
	    break;
        default:
            error("unrecognized type");
	    break;
        }
    }
}
