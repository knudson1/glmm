#include "myheader.h"
void cp3(double *etain, int *neta, int *typein, int *ntrials, double *cpout)
{
    int leneta=neta[0];
    int type = typein[0];
    int i;

    for(i=0;i<leneta;++i)
    {
		double eta = etain[i];
		switch (type) {
			case BERNOULLI:
			    cpout[i]= 1/(1+exp(-eta));
			    break;

		    case POISSON:
		        cpout[i]= exp(eta);
			    break;

			case BINOMIAL:
			    cpout[i]= ntrials[i]/(1+exp(-eta));
			    break;

		    default:
		        error("unrecognized type");
                            // clang complains about unreachable
                            // break;
        }
    }
}
