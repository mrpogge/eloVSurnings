#include <R.h>
#include <math.h>
#include <Rmath.h>

/* Convert elo rating from the logit scale to the 0,1 scale
it takes a subset of the vector t, and returns the expit transformed value*/
double expit(double t) {
    return 1.0 / (1.0 + exp(-t));
}

/* Compare function for qsort*/
int compare(const void *a, const void *b) {
    double arg1 = *(const double *)a;
    double arg2 = *(const double *)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}


void elo_HT(double*K_s,
            double*K_i,
            int*reps,
            int*games,
            int*N,
            int*M,
            double*theta,
            double*delta,
            double*t,
            double*d,
            double*lqT_baseline,
            double*uqT_baseline,
            int*adaptive,
            double*A,
            double*B,
            double*cumsum,
            int*htT,
            int*is_HT,
            double*mean_rating
            ){
/*
----------------------------------------------------------------------------
K_s: student learning rate
K_i: item learning rate
reps: number of replications
games: number of games
N: number of persons
M: number of items
theta: true ability of persons
delta: true difficulty of items
t: person ratings {vector[length=N*reps]}
d: item ratings {vector[length=M*reps]}
lqT_baseline: lower quantile of the baseline distribution {vector[length=N]}
uqT_baseline: upper quantile of the baseline distribution {vector[length=N]}
adaptive: whether to use adaptive item selection
A: parameter for adaptive item selection
B: parameter for adaptive item selection
cumsum: cumulative sum of item selection probabilities {vector[length=M+1]}
htT: hitting time for each student in the order of the students {vector[length=N]}
is_HTv: vector showing whether each student in each replication reached the hitting time {vector[length=N]}
is_HT: number of students who reached the hitting time
----------------------------------------------------------------------------
*/
     int x=0;
     double e=0;
     double L=0;
     double p=0;
     int j=0;
     double Mp=0;
     double dif=0;
     int sumHT = 0;   

     GetRNGstate();
     for(int g=0;g<games[0];g++){
        for(int r=0;r<reps[0];r++){
            for(int i=0;i<N[0];i++){

                /*specify probabilities for item selection*/
                if(adaptive[0]==1){
                    Mp=0;		
                    for(int s=0;s<M[0];s++){
                        dif=t[i+r*N[0]]-d[s+r*M[0]];
                        L=exp(dif*dif*A[0]+dif*B[0]);
                        Mp=Mp+L;
                        cumsum[s+1]=cumsum[s]+L;
                    }
                }
                if(adaptive[0]==0){
                    Mp=M[0];
                }
                /* sample an item given the selection probabilities*/
                p=runif(0,Mp);
                j=0;
                for(int s=1;s<M[0];s++){
                    if(p>cumsum[s]){
                        j=j+1;
                    }
                }
                /* compute the true probability correct*/
                L=1/(1+exp((delta[j]-theta[i]))); 

                /*generate the observed response*/
                p=runif(0,1.00);
                x=1*(L>p);

                /*compute the expected accuracy based on the current ratings*/
                e=1/(1+exp(d[j+r*M[0]]-t[i+r*N[0]]));
                /* update the ratings*/
                t[i+r*N[0]]=t[i+r*N[0]]+K_s[0]*(x-e);
                d[j+r*M[0]]=d[j+r*M[0]]-K_i[0]*(x-e);

                /*calculate the mean of the ratings*/
                mean_rating[i+g*N[0]] += t[i+r*N[0]] * 1.00 / reps[0];                
            }
        }
        /*check hitting time (HT) after each game*/
        sumHT = 0;
        for(int i=0;i<N[0];i++){
            if(is_HT[i] != 1){
                if(mean_rating[i+g*N[0]] > lqT_baseline[i] && mean_rating[i+g*N[0]] < uqT_baseline[i]){
                    htT[i] = g;
                    is_HT[i]= 1;
                }
            }
            sumHT += is_HT[i];
        }
        if(sumHT == N[0]){
            break;
        }
    }
     PutRNGstate();
}

void elo_double_HT(
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
            double*lqT_baseline,
            double*uqT_baseline,
            double*A,
            double*B,
            double*cumsum,
            int*htT,
            int*is_HT,
            double*mean_rating
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
htT: hitting time for each student in the order of the students {vector[length=N]}
is_HTv: vector showing whether each student reached the hitting time {vector[length=N]}
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
int sumHT = 0;
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
                    if(g % 2 == 1){
                        for(int s=0;s<M[0];s++){
                            dif = t2[i+r*N[0]] - d2[s+r*M[0]];
                            L = exp(dif*dif*A[0]+dif*B[0]);
                            Mp += L;
                            cumsum[s+1] = cumsum[s] + L;
                        }
                    } else {
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
                    L = 1/(1+exp((delta[j]-theta[i])));
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
                    mean_rating[i+g*N[0]] += expit(t1[i+r*N[0]])/reps[0];
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
                    mean_rating[i+g*N[0]] += expit(t2[i+r*N[0]])/reps[0];
                }
                for(int m = 0; m < M[0]; m++){
                    d2[m+r*M[0]] -= md2;
                }
            }
              
        }
        
        /*check hitting time (HT) after each game*/
        sumHT = 0;
        for(int i=0;i<N[0];i++){
            if(is_HT[i] != 1){
                if(mean_rating[i+g*N[0]] > lqT_baseline[i] && mean_rating[i+g*N[0]] < uqT_baseline[i]){
                    htT[i] = g;
                    is_HT[i]= 1;
                }
            }
            sumHT += is_HT[i];
        }
        if(sumHT == N[0]){
            break;
        }

    }    
    PutRNGstate();

}

void elo_HT(double*K_s,
    double*K_i,
    int*reps,
    int*games,
    int*N,
    int*M,
    double*theta,
    double*delta,
    double*t,
    double*d,
    double*lqT_baseline,
    double*uqT_baseline,
    int*adaptive,
    double*A,
    double*B,
    double*cumsum,
    int*htT,
    int*is_HT,
    double*mean_rating
    ){
/*
----------------------------------------------------------------------------
K_s: student learning rate
K_i: item learning rate
reps: number of replications
games: number of games
N: number of persons
M: number of items
theta: true ability of persons
delta: true difficulty of items
t: person ratings {vector[length=N*reps]}
d: item ratings {vector[length=M*reps]}
lqT_baseline: lower quantile of the baseline distribution {vector[length=N]}
uqT_baseline: upper quantile of the baseline distribution {vector[length=N]}
adaptive: whether to use adaptive item selection
A: parameter for adaptive item selection
B: parameter for adaptive item selection
cumsum: cumulative sum of item selection probabilities {vector[length=M+1]}
htT: hitting time for each student in the order of the students {vector[length=N]}
is_HTv: vector showing whether each student in each replication reached the hitting time {vector[length=N]}
is_HT: number of students who reached the hitting time
----------------------------------------------------------------------------
*/
int x=0;
double e=0;
double L=0;
double p=0;
int j=0;
double Mp=0;
double dif=0;
int sumHT = 0;   

GetRNGstate();
for(int g=0;g<games[0];g++){
for(int r=0;r<reps[0];r++){
    for(int i=0;i<N[0];i++){

        /*specify probabilities for item selection*/
        if(adaptive[0]==1){
            Mp=0;		
            for(int s=0;s<M[0];s++){
                dif=t[i+r*N[0]]-d[s+r*M[0]];
                L=exp(dif*dif*A[0]+dif*B[0]);
                Mp=Mp+L;
                cumsum[s+1]=cumsum[s]+L;
            }
        }
        if(adaptive[0]==0){
            Mp=M[0];
        }
        /* sample an item given the selection probabilities*/
        p=runif(0,Mp);
        j=0;
        for(int s=1;s<M[0];s++){
            if(p>cumsum[s]){
                j=j+1;
            }
        }
        /* compute the true probability correct*/
        L=1/(1+exp((delta[j]-theta[i]))); 

        /*generate the observed response*/
        p=runif(0,1.00);
        x=1*(L>p);

        /*compute the expected accuracy based on the current ratings*/
        e=1/(1+exp(d[j+r*M[0]]-t[i+r*N[0]]));
        /* update the ratings*/
        t[i+r*N[0]]=t[i+r*N[0]]+K_s[0]*(x-e);
        d[j+r*M[0]]=d[j+r*M[0]]-K_i[0]*(x-e);

        /*calculate the mean of the ratings*/
        mean_rating[i+g*N[0]] += t[i+r*N[0]] * 1.00 / reps[0];                
    }
}
/*check hitting time (HT) after each game*/
sumHT = 0;
for(int i=0;i<N[0];i++){
    if(is_HT[i] != 1){
        if(mean_rating[i+g*N[0]] > lqT_baseline[i] && mean_rating[i+g*N[0]] < uqT_baseline[i]){
            htT[i] = g;
            is_HT[i]= 1;
        }
    }
    sumHT += is_HT[i];
}
if(sumHT == N[0]){
    break;
}
}
PutRNGstate();
}

void elo_simple_HT(
    double*K_s,
    int*reps,
    int*games,
    int*N,
    int*M,
    double*theta,
    double*delta,
    double*t,
    double*d,
    double*lqT_baseline,
    double*uqT_baseline,
    double*A,
    double*B,
    double*cumsum,
    int*htT,
    int*is_HT,
    double*mean_rating
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
htT: hitting time for each student in the order of the students {vector[length=N]}
is_HTv: vector showing whether each student reached the hitting time {vector[length=N]}
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
int sumHT = 0;

GetRNGstate();
/*loops*/
for(int g=0;g<games[0];g++){
for(int r=0;r<reps[0];r++){
    for(int i=0;i<N[0];i++){
            
            /*specify item selection probabilities*/
            Mp = 0;
            cumsum[0] = 0.0; //do we need this?
            for(int s=0;s<M[0];s++){
                dif = t[i+r*N[0]] - delta[s];
                L = exp(dif*dif*A[0]+dif*B[0]);
                Mp += L;
                cumsum[s+1] = cumsum[s] + L;
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
            L = 1/(1+exp((delta[j]-theta[i])));
            /*generate the observed response*/
            p = runif(0,1.00);
            x = 1*(L>p);
            /*compute the expected accuracy*/
             e = 1/(1+exp(delta[j]-t1[i+r*N[0]]));
            /* update the ratings*/
            t[i+r*N[0]] = t[i+r*N[0]] + K_s[0]*(x-e);
            mean_rating[i+g*N[0]] += expit(t[i+r*N[0]])/reps[0];    
    }       
}

/*check hitting time (HT) after each game*/
sumHT = 0;
for(int i=0;i<N[0];i++){
    if(is_HT[i] != 1){
        if(mean_rating[i+g*N[0]] > lqT_baseline[i] && mean_rating[i+g*N[0]] < uqT_baseline[i]){
            htT[i] = g;
            is_HT[i]= 1;
        }
    }
    sumHT += is_HT[i];
}
if(sumHT == N[0]){
    break;
}

}    
PutRNGstate();

}