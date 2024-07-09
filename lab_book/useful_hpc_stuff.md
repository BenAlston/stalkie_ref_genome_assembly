
# **Useful HPC Stuff**

Download files directly to HPC: 

~~~
 wget -r https://cgr.liv.ac.uk/pbio/200437_987ddbd846f6fba3/ -P /shared/wright_lab_hpc/Shared/2024_Stalkie_PacBio_HiFi
 # For some reason this sometimes gets an error 404, 
 # fix:
 wget -r https://cgr.liv.ac.uk/pbio/200437_987ddbd846f6fba3/ --header "cgr.liv.ac.uk" -P /shared/wright_lab_hpc/Shared/2024_Stalkie_PacBio_HiFi/
 # couldn't tell you why 
 # wget has also started ignoring the output directory, so had to set the target dir as working directory instead, again the reason escapes me
~~~

Start interactive session in stanage/bessimer: 
~~~
srun --pty bash -i
~~~
Fastdata: 
~~~
/mnt/parscratch/users/$user
~~~
Connect via mac: 
~~~
ssh $username@$cluster.shef.ac.uk
~~~
check module availbity
~~~
module avail |& grep -i $modulename
~~~
Job control: 
~~~
# submit
sbatch $script.sh
squeue -u(user) or -–job(jobnumber)
# Cancel job: 
scancel -u or –-job
# check job resource use
seff $jobid # can use this and trial and error to work out how many resources a job should be using when the developers decide not to give any guidance
~~~
Download from HPC to local machine:
~~~
scp bip23bta@stanage.shef.ac.uk:/mnt/parscratch/users/bip23bta/ref_genomes/whitei/03-QC/purge_dups/1_primary/round_1_man_cutoffs/man_cov_hist_1.png Downloads/
# run command on local machine, not hpc
~~~

Tmux:
~~~
tmux new -s $session_name
tmux attatch -t $session_name
tmux ls
ctrl + b  # to do things
~~~
Disk usage: 
~~~
du -bch $directory # estimates filesize of dir and all subdirs
~~~

Md5sums check:
~~~
d5sum -c /filepath/to/md5checklist.txt
# checklist provided by source, each line looks like: $md5sum  $filepath
# Make sure filepath of data directory matches the provided checklist
~~~

Conda environment
~~~
module load Anaconda3
conda create -n $name bioconda::$package
conda activate $env

# or
conda create -n $name
source activate $name
conda install bioconda::package

# view a list of all installed conda envs:
conda env list

# search available versions
conda search $package

~~~

Find conda env package version
~~~
# activate environment
source activate $env
# list and search 
conda list | grep $env
~~~

search a directory and subdirectory for a file:
~~~
find $dir_name -name '$name e.g., *.fastq.gz'
~~~

some stuff about array jobs
~~~
# good way of easily parralelizing stuff
# hpc documentation on this is poor
${SLURM_ARRAY_TASK_ID} # corresponds to the current array job running, e.g., 1 or 2 for $jobid_1 or $jobid_2
# the problem comes from not being able to define non numeric slurm array ids
# for my purposes, using bash arrays can fix this somewhat:

#!/bin/bash
#SBATCH --job-name=array_test
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --array 0-3%4 # change depending on how many tasks you need doing
#SBATCH -e error.txt

inputs=('in1' 'in2' 'in3' 'in4')

for task in ${inputs[$SLURM_ARRAY_TASK_ID]}
do
echo $task > ${task}_output.txt
done

# this test script essentially matches the slurm array id with an element of an array ive defined ($inputs in this case)
# we can then do a loop for each if needed 
# the UoS hpc documentation on this is shocking

# another (better?) way:
# make file of line sepparated file names or prefixes, e.g.:
cat filenames.txt
A01
A02
A03

# read this as an array (a specific kind of unix object, kinda like vectors in R)

readarray -t array < filenames.txt

# now i can access anything in the array using numbers 0-2, e.g.,
echo "${array[1]}"
A02

# so each sample now correstponds to a different slurm array task id (which go from 0-however many you specify)
FILENAME="${array[$SLURM_ARRAY_TASK_ID]}

# for SLURM_ARRAY_TASK_ID=0:
echo $FILENAME
A01
~~~
Docker/apptainer
~~~
# docker images are like conda environments but less intuative, but deal much better with dependencies
# I tend to use when something isn't available through conda, or the only conda env available is really outdated
# Apptainer is the open source version of docker, and is used on sheffield HPCs

# usage with docker hub (main repository for docker/apptainer containers)
# make the apptainer image (downloads a .sif file to wd):
apptainer pull  docker://ezlabgva/busco:v5.6.1_cv1# equivalent to creating a conda env

# load image. Analagous to <conda activate $env>:
apptainer shell /path/to/$image.sif

# use in a script, according to hpc documentation:
apptainer exec path/to/imgfile.img ls /
# the above doesn't actually work

# but the below method does:
# when running busco, do this:
apptainer exec /users/bip24bta/busco_v5.6.1_cv1.sif busco -i $GENOME \
-l diptera_odb10 \
-o $OUT \
-m genome \
-f \
-c 20

~~~



