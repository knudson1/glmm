
#include <R.h>
#include "myheader.h"


void cpp3(double *etain, int *neta, int *typein, double *cppout)
{
    int leneta=neta[0];
    int type = typein[0];
    int i;


    for(i=0;i<leneta;++i)
    {
	double eta = etain[i];
	switch (type) {
	case 1:
	    cppout[i]= 1/((1+exp(-eta))*(1+exp(eta)));
	    break;
        case 2:
            cppout[i]= exp(eta);
	    break;
        default:
            error("unrecognized type");
	    break;
        }
    }
}
   
