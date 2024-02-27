#########################################################################################
#	Script Name: gfa_to_fa
#	Description: converts .gfa files to .fa
#	Author:      Ben Alston
#	Date:        Feb 2024
#########################################################################################

# Usage: bash ${SCRIPT_NAME}.sh [PATH/TO/FILE.gfa] [OUTPUT_FILE_PREFIX]

# Assign variables
filepath=$(dirname $1)
input_file=$1
output_prefix=$2
output_file=$filepath/${2}.fa

# Check if path to file is provided
if [ $# -eq 0 ]; then
    echo "Error: Please provide input file. Usage: bash [script_name].sh [filepath/to/file.gfa] [output_file_prefix]"
    exit 1
fi

# Check if the file exists
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found."
    exit 1
fi

# check output prefix is supplied
if [ -z "$output_prefix" ]; then
    echo "Error: Please provide an output file prefix. Usage: bash [script_name].sh [filepath/to/file.gfa] [output_file_prefix]"
    exit 1
fi

# Process file
echo "Processing file: $input_file"
awk '/^S/{print ">"$2"\n"$3}' $input_file | fold > $output_file
echo "Output: $output_file"
echo "sorted"






