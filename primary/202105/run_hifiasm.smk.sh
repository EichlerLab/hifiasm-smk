module purge
module load modules modules-init modules-gs/prod modules-eichler/prod 
module load miniconda/4.5.12



#snakemake -s /net/eichler/vol26/7200/software/pipelines/hifiasm/primary/202105/hifiasm.smk --jobname "{rulename}.{jobid}" --drmaa " -l mfree={resources.mem}G -pe serial {threads} -l h_rt={resources.hrs}:00:00 -w n  -V -cwd -e ./log -o ./log  -S /bin/bash" -w 300 --jobs 200 -p

snakemake -s /net/eichler/vol26/7200/software/pipelines/hifiasm/primary/202105/hifiasm.smk --jobname "{rulename}.{jobid}" --drmaa " -l mfree={resources.mem}G -pe serial {threads} -w n  -V -cwd -e ./log -o ./log  -S /bin/bash" -w 300 --jobs 200 -p
