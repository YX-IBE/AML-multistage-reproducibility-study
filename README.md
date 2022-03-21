# AML Multistage Reproducibility Study

This repository is forked from `gerstung-lab/AML-multistage` for our manuscript *Validating the knowledge bank approach for personalized prediction of survival in acute myeloid leukemia: a reproducibility study*.

### Reproducibility assessment

`/doc/SupplementaryMethodsCode.R` is to present our modifications of the original code by Gerstung *et al.*. We modified the code in order to reproduce their published results in *Gerstung M, Papaemmanuil E, Martincorena I, et al. Precision oncology for acute myeloid leukemia using a knowledge bank approach. Nat Genet. 2017;49(3):332-340. doi:10.1038/ng.3756*.   

For comparison purposes, modifications were only made when coding errors from the original R script prevented the continuation of the process.

### Reanalysis

We tailored the source code and performed a reanalysis for ‘concept testing’ purposes: to clarify the conflation between causal inference and predictive model in the original paper. The results are summarized in **Fig. 3** in our manuscript.

- `/reanalysis/Fig3.R` is to generate **Fig. 3**, which works on the same data in `/data` as did in the original study.

- In addition, we provided a Docker image for our reanalysis: https://hub.docker.com/repository/docker/yujunxu/aml-multistage-reanalysis (built from `/reanalysis/Dockerfile` and `/reanalysis/install_packages.R`, tested on `Debian  11.0 64-bit` and `Ubuntu  18.04 64-bit`)

   - docker pull  
      `docker pull yujunxu/aml-multistage-reanalysis:1.3`
   - docker run  
      `docker run -it yujunxu/aml-multistage-reanalysis:1.3 /bin/bash`
