#include <R.h>
#include <math.h>
#include <Rmath.h>

/* Convert elo rating from the logit scale to the 0,1 scale
it takes a subset of the vector t, and returns the expit transformed value*/
double expit(double t) {
    return 1.0 / (1.0 + exp(-t));
}

void elo_double_change(
    double*K_s,
    double*K_i,
    int*reps,
    int*games,
    int*N,
    int*M,
    double*theta,
    double*delta,
    double*t1,
    double*d1,
    double*t2,
    double*d2,
    double*mT,
    double*A,
    double*B,
    double*cumsum,
    double*bias_v,
    double*var_v
    ){
/*Arguments
----------------------------------------------------------------------------------------
K_s: student learning rate
K_i: item learning rate
reps: number of replications
games: number of games
N: number of persons
M: number of items
theta: true ability of persons
delta: true difficulty of items
t: person ratings {vector[length=N*reps]}
t1: first elo rating of the students {vector[length=N*reps]}
d1: first elo rating of the items {vector[length=M*reps]}
t2: second elo rating of the students {vector[length=N*reps]}
d2: second elo rating of the items {vector[length=M*reps]}
lqT_baseline: lower quantile of the baseline distribution {vector[length=N]}
uqT_baseline: upper quantile of the baseline distribution {vector[length=N]}
A: parameter for adaptive item selection
B: parameter for adaptive item selection
cumsum: cumulative sum of item selection probabilities {vector[length=M+1]}
*/

/*aux variables
----------------------------------------------------------------------------
*/
double Mp = 0; //normalising constant for adaptive item selection
double p = 0; //random number for sampling
int j =0; //item index
double L = 0; //probability of correct response
double dif = 0; //differences in ratings for NKM
int x = 0; //observed response
double e = 0; //expected accuracy
double md1 = 0.0;
double md2 = 0.0; // mean of the item ratings in a given game and interation

GetRNGstate();
/*loops*/
for(int g=0;g<games[0];g++){
for(int r=0;r<reps[0];r++){
    for(int i=0;i<N[0];i++){

            /*specify item selection probabilities*/
            Mp = 0;
            cumsum[0] = 0.0; //do we need this?
            if(g % 2==1){
                for(int s=0;s<M[0];s++){
                    dif = t2[i+r*N[0]] - d2[s+r*M[0]];
                    L = exp(dif*dif*A[0]+dif*B[0]);
                    Mp += L;
                    cumsum[s+1] = cumsum[s] + L;
                }
            } 
            
            if(g % 2 == 0) {
                for(int s=0;s<M[0];s++){
                    double dif = t1[i+r*N[0]] - d1[s+r*M[0]];
                    double L = exp(dif*dif*A[0]+dif*B[0]);
                    Mp += L;
                    cumsum[s+1] = cumsum[s] + L;
                }
            }
            
            /* sample an item given the selection probabilities*/
            p = runif(0,Mp);
            j=0;
            for(int s=1;s<M[0];s++){
                if(p>cumsum[s]){
                    j=j+1;
                }
            }

            /*compute true probability*/
            L = 1/(1+exp((delta[j]-theta[i+g*N[0]])));
            /*generate the observed response*/
            p = runif(0,1.00);
            x = 1*(L>p);
            /*compute the expected accuracy*/
            if(g % 2 == 1){
                e = 1/(1+exp(d1[j+r*M[0]]-t1[i+r*N[0]]));
                /* update the ratings*/
                t1[i+r*N[0]] = t1[i+r*N[0]] + K_s[0]*(x-e);
                d1[j+r*M[0]] = d1[j+r*M[0]] - K_i[0]*(x-e);     
            } else {
                e = 1/(1+exp(d2[j+r*M[0]]-t2[i+r*N[0]]));
                /* update the ratings*/
                t2[i+r*N[0]] = t2[i+r*N[0]] + K_s[0]*(x-e);
                d2[j+r*M[0]] = d2[j+r*M[0]] - K_i[0]*(x-e);
            }
    }
    if(g % 2 == 1){
        md1 = 0.0;
        for(int m = 0; m < M[0]; m++){
            md1 += d1[m+r*M[0]] / M[0];
        }
        /*after each game in each replication rescale values to ensure identifiability*/
        //NOTE: this only works with mean fixed at 0
        for(int i = 0; i < N[0]; i++){
            t1[i+r*N[0]] -=  md1;
            mT[i+g*N[0]] += expit(t1[i+r*N[0]])/reps[0];
        }
        for(int m = 0; m < M[0]; m++){
            d1[m+r*M[0]] -= md1;
        }
    } else {
        md2 = 0.0;
        for(int m = 0; m < M[0]; m++){
            md2 += d2[m+r*M[0]] / M[0];
        }
        /*after each game in each replication rescale values to ensure identifiability*/
        //NOTE: this only works with mean fixed at 0
        for(int i = 0; i < N[0]; i++){
            t2[i+r*N[0]] -=  md2;
            mT[i+g*N[0]] += expit(t2[i+r*N[0]])/reps[0];
        }
        for(int m = 0; m < M[0]; m++){
            d2[m+r*M[0]] -= md2;
        }
    }
}
    /*calculate bias and variance in each game over replications*/
    for(int i = 0; i < N[0]; i++){
        bias_v[i+g*N[0]] = mT[i+g*N[0]] - expit(theta[i+g*N[0]]);
        for(int r = 0; r<reps[0]; r++){
            if(g % 2 == 1){
                var_v[i+g*N[0]] += pow(expit(t1[i+r*N[0]]) - mT[i+g*N[0]],2) / (reps[0]-1);
            } else {
            var_v[i+g*N[0]] += pow(expit(t2[i+r*N[0]]) - mT[i+g*N[0]],2) / (reps[0]-1);
            }
        }
    }
}    
PutRNGstate();

}