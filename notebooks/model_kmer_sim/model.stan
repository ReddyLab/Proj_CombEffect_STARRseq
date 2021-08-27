data {
    int <lower=1> N_DATA;
    int <lower=1> N_FRAG;
    int <lower=1> N_SEG;
    int <lower=1> N_SEG_PER_FRAG;
    int <lower=1> N_REP;
    
    int <lower=0> OUTPUT[N_DATA];
    int <lower=1> INPUT[N_DATA];
    real <lower=0> PHI;
}

parameters {
    real<lower=0> SEG[N_SEG];
}

model {
    //print("hello");
    //for (j in 1:N_SEG) {
    //    print(SEG[j]);
    //}
    
    // Shrinkage priors on each effect size
    for (j in 1:N_SEG) {
        SEG[j] ~ lognormal(0, 0.3);
    }
    
    // Likelihood
    for (i in 1:N_FRAG) {
    
        // gather effect from corresponding segments
        real seg_effect = 1;
        for (j in i:(i + N_SEG_PER_FRAG - 1)) {
            seg_effect *= SEG[j];
        }
        
        // mean(output) := segment effect x input
        int idx_data;
        for (k in 1:N_REP) {
            idx_data = N_FRAG * (i - 1) + k;
            OUTPUT[idx_data] ~ neg_binomial_2(INPUT[idx_data] * seg_effect, PHI);
            
            //print(seg_effect, " ", INPUT[idx_data], " ", INPUT[idx_data] * seg_effect);
        } 
    }
}
