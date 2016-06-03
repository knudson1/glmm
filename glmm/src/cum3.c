
#include <R.h>
#include "myheader.h"


void cum3(double *etain, int *neta, int *typein, int *ntrials, double *cumout)
{
    int leneta=neta[0];
    int type = typein[0];
    int i;
    
 
    for(i=0;i<leneta;++i)
    	{
		double eta = etain[i];
		switch (type) {
		case BERNOULLI:
	 	   if(eta>0)
				{
				cumout[0]+= eta+log1p(exp(-eta));
				break;
				}
   	    	else
				{
				cumout[0]+= log1p(exp(eta));
				break;
				}

	    case POISSON:
	        cumout[0]+= exp(eta);
		    break;

		case BINOMIAL:
	 	   if(eta>0)
				{
				cumout[0]+= ntrials[i] * (eta+log1p(exp(-eta)));
				break;
				}
   	    	else
				{
				cumout[0]+= ntrials[i] * log1p(exp(eta));
				break;
				}

	    default:
	        error("unrecognized type");
			break;
        }
    }
}

