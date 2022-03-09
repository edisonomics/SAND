# NMR_time_domain_decomposition
Decompose NMR spectra based on fitting time domain data


# introduction

# install
the scripts can be run without installation, though some dependencies are needed.
1. nmrpipe
2. toolbox (follow )
3. To run the decomposition script, the user can communicate with their local systems and dicuss running matlab in mutiple parallel in HPC.

# tutorials
First follow peak_moving_simulaiton workflow `scripts/peak_moving_simulaiton/simu_spec_decomp.m` which works on a single computer though will take a few hours to finish.

A general process of the spectral decomposiiton include preprocess the spectra, decompostion (mostly in HPC), visualize the decomposition results (e.g. `scripts/urine_like_simualtion_phase`). The corresponding running scripts will be named as data_preprocess.m, test_region_separa_deconv_hpc.m, and vis_spec_run.m. To run this process:

1. Run data_preprocess.m in a local folder with necessary library and data in place as in the script.
2. Construct the HPC (sapelo2 for current test) folder:
```
    ./temp,
    ./data (Upload LOCAL_PROJ_FOLDER/res/nmrpipe_dir/archive.zip unzip here),
    ./res,
    runtab.txt (the exploration parameters for each run),
    test_region_separa_deconv_hpc.m (the decomposition script),
    parallel_job_submit.R (submit the decomposition jobs in batch),
    submit.sh (template submission script)
```
3. Load R environment and run parallel_job_submit.R and wait for the submitted job to finish
4. Download the result: whole folder except the content of ./temp and ./data
5. Run vis_spec_run.m

Attention:

1. The user need to ensure Matlab can run in the HPC env in a intensive parallel pattern.
the folder structure is:
      ./temp,
      ./data (unzip archive.zip from nmrpipe_dir here),
      ./res,
      runtab.txt (the exploration parameters for each run),
      test_region_separa_deconv_hpc.m (the decomposition script),
      parallel_job_submit.R (submit the decomposition jobs in batch),
      submit.sh (template submission script)

nmrpipe code please refer to frank.delaglio@nist.gov

# Versions
