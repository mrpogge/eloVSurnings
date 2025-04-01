#include <R.h>
#include <math.h>
#include <Rmath.h>

/* Convert elo rating from the logit scale to the 0,1 scale
it takes a subset of the vector t, and returns the expit transformed value*/
double expit(double t) {
    return 1.0 / (1.0 + exp(-t));
}


void urnings_change(int* K_s,
                int* K_i,
                int* reps,
                int* games,
                int* N,
                int* M,
                double* theta,
                double* delta,
                int* t,
                int* d,
                int* adaptive,
                double* P,
                double* cumsum,
                int* paired,
                int* queue,
                int* P_queue_norm,
                double* mean_rating,
                double* bias_v,
                double* var_v
                )
{
/*
----------------------------------------------------------------------------
K_s: person urn size
K_i: item urn size
reps: number of replications
games: number of games
N: number of persons
M: number of items
theta: person ability {vector[N]}
delta: item difficulty {vector[M]}
t: player urnings {vector[N*reps]}
d: item urnings {vector[M*reps]}
seT: squared error of person ratings {vector[N*reps]}
seT_sys: squared error of the system {vector[reps]}
baseline: baseline MSE for each student in the order of the students {vector[N]}
baseline_sys: system level baseline: {vector[1]}
adaptive: whether to use adaptive item selection
P: item selection probabilities {vector[(K_s+1)*(K_i+1)]}
cumsum: cumulative sum of item selection probabilities {vector[M+1]}
paired: whether to use paired updates
htT: hitting time for each student in the order of the students {vector[N]}
htT_sys: hitting time for the system {vector[1]}
is_HTv: vector showing whether each student reached the hitting time {vector[N]}
is_HT: number of students who reached the hitting time
----------------------------------------------------------------------------*/

/*aux variables
---------------------------------------------------------------------------*/
    double p=0.00;
    int j=0;
    double Mp=0;
    //double Cp=0;
    int x=0;
    int y=0;
    double Correct=0;
    double Incorrect=0;
    double L=0;

    //aux variables for MH step
    int oldU=0;
    int oldV=0;
    int newU=0;
    int newV=0;
    double Mpp = 0; //normalizing constant for the selection probability with the proposed updates

    //aux variables for paired updates
    int is_c = 0; // indicator whether an item is a candidate for a paired update 
    int j_c = 0; // index of the candidate item

/*loop
--------------------------------------------------------------------------*/
    GetRNGstate();
    /*loop over games*/
    for(int g=0; g<games[0];g++){
        /*loop over replications*/
        for(int r=0;r<reps[0];r++){
            /*loop over persons*/
            for(int i=0;i<N[0];i++){
                
                /*random item selection*/
                if(adaptive[0]==0){
                    p=runif(0,1.00*M[0]);
                    j=0;
                    for(int s=1;s<M[0];s++){
                        if(p>s){
                            j=j+1;
                        }
                    }
                }

                /*adaptive item selection*/
                if(adaptive[0]==1){
                    Mp=0;
                    for(int s=0; s<M[0]; s++){
                        Mp += P[t[i]+d[s]*(K_s[0]+1)];
                        cumsum[s+1] = cumsum[s] + P[t[i]+d[s]*(K_s[0]+1)];
                    }
                    p=runif(0,Mp);
                    j=0;
                    for(int s=1; s<M[0]; s++){
                        if(p>cumsum[s]){
                            j=j+1;
                        }
                    }
                }

                /*save the current person and item values*/
                oldU = t[i+r*N[0]];
                oldV = d[j+r*M[0]];

                /*generate the observed accuracy*/
                L=1/(1+exp((delta[j]-theta[i+g*N[0]])));
                p=runif(0,1.00);
                x=1*(L>p);

                /*augment the urns*/
                t[i+r*N[0]] +=  x;
                d[j+r*M[0]] += (1-x);

                /*compute the probability x = 1*/
                Correct = t[i+r*N[0]]*(K_i[0]+1-d[j+r*M[0]]);
                Incorrect = (K_s[0]+1-t[i+r*N[0]])*d[j+r*M[0]];
                L = Correct/(Correct+Incorrect);

                /*generate the response*/
                p=runif(0,1.00);
                y=1*(L>p);
                /*remove the balls from the urns based on the simulated response*/
                t[i+r*N[0]] -= y;
                d[j+r*M[0]] -= (1-y);

                /*Metropolis step*/ 
                if(adaptive[0] == 1){
                    newU = t[i+r*N[0]];
                    newV = d[j+r*M[0]];
                    t[i+r*N[0]] = oldU;
                    d[j+r*M[0]] = oldV;
                    if(newU != oldU){
                        /*Compute the normalising constant for selection probability with the proposed updates*/
                        Mpp = P[newU+newV*(K_s[0]+1)];
                        for(int s=0;s<K_s[0];s++){
                            if(s!=j){
                                Mpp += P[newU+d[s]*(K_s[0]+1)];
                            }
                        }
                        /*compute the MH acceptance probability*/
                        L = P[newU+newV*(K_s[0]+1)]/Mpp*Mp/P[oldU+oldV*(K_s[0]+1)];
                        /*generate unif(0,1), to decide acceptandce*/
                        p=runif(0,1.00);
                        if(p<L){
                            t[i+r*N[0]] = newU;
                            d[j+r*M[0]] = newV;
                        }
                    }
                }
                /*Paired updated*/
                if(paired[0] ==1){
                    /*set the item j to the old value*/
                    newV = d[j+r*M[0]];
                    d[j+r*M[0]] = oldV;

                    /*if the new and old are different than we enter paired update*/
                    if(newV != oldV){
                        if(newV > oldV){
                            /*Find items that are possible to update*/
                            P_queue_norm[r] = 0;
                            cumsum[0] = 0.0;
                            for(int s=0; s<M[0]; s++){
                                is_c=1*(queue[s+r*M[0]] == -1)*(d[s+r*M[0]] != 0)*(s!=j);
                                P_queue_norm[r] += is_c;
                                cumsum[s+1] = cumsum[s] + is_c;
                            }
                            /*If there is no candidate add the item to the queue*/
                            if(P_queue_norm[r] == 0){
                                queue[j+r*M[0]] = 1;
                            } else { /*if there is a candidate select a random candidate and perform the update*/
                                p=runif(0,1.00*P_queue_norm[r]);
                                j_c =0;
                                for(int s=1; s<M[0]; s++){
                                    if(p>cumsum[s]){
                                        j_c=j_c+1;
                                    }
                                }
                                queue[j_c+r*M[0]] = 0; //remove the candidate from the queue
                                d[j_c+r*M[0]] -= 1; //update the candidate item
                                d[j+r*M[0]] = newV; //update the item
                        }
                    }
                    if(newV < oldV){
                         /*Find items that are possible to update*/
                         P_queue_norm[r] = 0;
                         cumsum[0] = 0.0;
                         for(int s=0; s<M[0]; s++){
                             is_c=1*(queue[s+r*M[0]] == 1)*(d[s+r*M[0]] != K_i[0])*(s!=j);
                             P_queue_norm[r] += is_c;
                             cumsum[s+1] = cumsum[s] + is_c;
                         }
                         /*If there is no candidate add the item to the queue*/
                         if(P_queue_norm[r] == 0){
                             queue[j+r*M[0]] = -1;
                         } else { /*if there is a candidate select a random candidate and perform the update*/
                             p=runif(0,1.00*P_queue_norm[r]);
                             j_c =0;
                             for(int s=1; s<M[0]; s++){
                                 if(p>cumsum[s]){
                                     j_c=j_c+1;
                                 }
                             }
                             queue[j_c+r*M[0]] = 0; //remove the candidate from the queue
                             d[j_c+r*M[0]] += 1; //update the candidate item
                            d[j+r*M[0]] = newV; //update the item
                        }
                    }
                }
            }

                /*calculate the average rating*/
                mean_rating[i+g*N[0]] += t[i+r*N[0]] * 1.00 / reps[0];
            }
        }
       /*calculate bias and variance in each game over replications*/
        for(int i = 0; i < N[0]; i++){
            bias_v[i+g*N[0]] = (mean_rating[i+g*N[0]] * 1.00/K_s[0]) - expit(theta[i+g*N[0]]);
            for(int r = 0; r<reps[0]; r++){
                var_v[i+g*N[0]] += pow(t[i+r*N[0]] * 1.00/K_s[0] - mean_rating[i+g*N[0]] * 1.00/K_s[0],2) / (reps[0]-1);
            }
        }
    }
    PutRNGstate();
}