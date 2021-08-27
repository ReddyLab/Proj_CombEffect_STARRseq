data {
    int <lower=1> NUM_DATA;
    int <lower=1> NUM_FRAG;
    int <lower=1> NUM_SEG;
    int <lower=1> NUM_SEG_PER_FRAG;
    int <lower=1> NUM_REP;
    
    int <lower=0> OUTPUT[NUM_DATA];
    int <lower=1> INPUT[NUM_DATA];
    real <lower=0> PHI;
    
    int L_DNA; // DNA library size
    int L_RNA; // RNA library size
    real shrinkageVar; // variance of shrinkage prior on effect sizes
}

transformed data {
   real l_rna=L_RNA; real l_dna=L_DNA;
   real LIB_RATIO=l_rna/l_dna;
}

parameters {
    real<lower=0.0001> gamma[NUM_SEG]; // effect size of each fragment
}

model {
    // Compute prefix-sum array for log(gamma)
    real PSA[NUM_SEG + 1];
    PSA[1]=0;
    for(i in 2:(NUM_SEG+1)){
        PSA[i]=PSA[i-1]+log(gamma[i-1]);
    }
    
    // shrinkage prior on effect sizes
    for(i in 1:NUM_SEG){
        gamma[i] ~ lognormal(0, shrinkageVar);    
    }
    
    // Likelihood
    for (i in 1:NUM_FRAG) {
        
        // gather effect from corresponding segments
        real gammaProd = exp(PSA[i + NUM_SEG_PER_FRAG] - PSA[i]);
        
        // mean(output) := segment effect x input x library ratio
        for (k in 1:NUM_REP) {
            int idx_data = NUM_FRAG * (i - 1) + k;
            OUTPUT[idx_data] ~ neg_binomial_2(INPUT[idx_data] * LIB_RATIO * gammaProd, PHI);
            
            //print(seg_effect, " ", INPUT[idx_data], " ", INPUT[idx_data] * seg_effect);
        } 
    }
}
