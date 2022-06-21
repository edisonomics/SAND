SAND can be run on mutiple environments. We provide examples on local computer, GACRC, NMRBox as examples. We recommend the users to focus on following the example on NMRBox as the most efficient workflow for the general usage.

# Prerequisite

SAND programs can be run without installation, though some dependencies are needed.
1. [NMRPipe](https://www.ibbr.umd.edu/nmrpipe/install.html). Refer to [Frank delaglio](frank.delaglio@nist.gov) for questions on NMRPipe
2. [Metabolomics Toolbox](https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA)
3. Many scripts need HPC and accebility of paralell computing. The user can communicate with their local systems and dicuss running matlab in mutiple parallel in HPC

# Simple first run

The user can test the script here `scripts/peak_moving_simulaiton/simu_spec_decomp.m` which works on a single local computer and will take a few hours to finish.

# Running SAND on GACRC

A general process of the spectral decomposiiton include preprocess the spectra, decompostion (mostly in HPC), visualize the decomposition results (examples in `scripts/urine_like_simualtion_phase`). The corresponding running scripts will be named as data_preprocess.m, test_region_separa_deconv_hpc.m, and vis_spec_run.m. This has been tested on [Sapelo2](https://wiki.gacrc.uga.edu/wiki/Running_Jobs_on_Sapelo2) of GACRC. To run this process:

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

## Attention

1. The user need to ensure Matlab can run in the HPC env in a intensive parallel pattern.

# Running SAND on NMRBox

NMRBox provides virtual environment and efficient support on processing and analyzing NMR data. Data can be efficiently transferred to NMRBox and the analysis can start and end there.
