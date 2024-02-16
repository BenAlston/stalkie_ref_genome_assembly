
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
Fastdata (because /fastdata/$user is too convenient): 
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
scp bip23bta@stanage.shef.ac.uk:/shared/wright_lab_hpc/Shared/2024_Stalkie_illumina/index.html /home/benalston/
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

Md5sums check: m
~~~
d5sum -c /filepath/to/md5checklist.txt
# checklist provided by source, each line looks like: $md5sum  $filepath
# Make sure filepath of data directory matches the provided checklist
~~~

Conda environment
~~~
module load Anaconda3
conda create -n $name bioconda::$package
Conda activate $env
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
~~~
