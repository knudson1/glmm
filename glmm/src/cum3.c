
#include <R.h>
#include "myheader.h"


void cum3(double *etain, int *neta, int *typein, int *ntrials, double *wts, double *cumout)
{
    int leneta=neta[0];
    int type = typein[0];
    int i;
    
 
    for(i=0;i<leneta;++i)
    	{
		double eta = etain[i];
		int ntr = ntrials[i];
            double w = wts[i];
		switch (type) {
		case BERNOULLI:
	 	   if(eta>0)
				{
				cumout[0]+= w * (eta+log1p(exp(-eta)));
				break;
				}
   	    	else
				{
				cumout[0]+= w * (log1p(exp(eta)));
				break;
				}

	    case POISSON:
	        cumout[0]+= w * (exp(eta));
		    break;

		case BINOMIAL:
	 	   if(eta>0)
				{
				cumout[0]+= w * (ntr*(eta+log1p(exp(-eta))));
				break;
				}
   	    	else
				{
				cumout[0]+= w * (ntr*log1p(exp(eta)));
				break;
				}

	    default:
	        error("unrecognized type");
			break;
        }
    }
}

