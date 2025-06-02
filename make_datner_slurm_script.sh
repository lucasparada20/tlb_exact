#!/bin/sh

# This script builds TLB scripts to be launched in high-performance computing systems i.e., SLURM
# The model number need to be hard coded

# Define directories
INSTANCE_DIR=~/work/tlb/instances_datner
TARGETS_DIR=~/work/tlb/targets_datner
RESULTS_DIR=~/work/tlb/results
SCRIPTS_DIR=~/work/tlb/scripts
RUN_SCRIPT=run_slurm_scripts.sh
EXECUTABLE=build/exec_tlb

mkdir -p "$SCRIPTS_DIR"

script_count=0

# Loop through each *_nXXX_eXXX.txt file
for instance_file in "$INSTANCE_DIR"/*_n*_e*.txt; do
    [ -f "$instance_file" ] || continue

    filename=$(basename "$instance_file")
    base_name="${filename%.txt}"  # e.g., boston1_n132_e100
    city_prefix=$(echo "$base_name" | cut -d'_' -f1)  # e.g., boston1
    scenario_file="${INSTANCE_DIR}/${city_prefix}_e400.txt"
    re_file="re_${city_prefix}.txt"
    output_file="targets_${city_prefix}.txt"

    # Skip if result already exists
    if [ -f "$RESULTS_DIR/$re_file" ]; then
        echo "Skipping $filename since $re_file exists."
        continue
    fi

    script_name="$SCRIPTS_DIR/${base_name}.sh"
    cat > "$script_name" <<EOF
#!/bin/sh -l
#SBATCH --mem=30G
#SBATCH --time=1:00:00
#SBATCH --account=def-cotej
#SBATCH --cpus-per-task=64
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

$EXECUTABLE mcfp_solver=3 instance_format=trips \\
    instance_file=$INSTANCE_DIR/$filename \\
    scenarios_file_name=$scenario_file \\
    model=2 \\
    re_file=$RESULTS_DIR/$re_file \\
    output_file=$TARGETS_DIR/$output_file \\
	> output/${city_prefix}.out 2> output/${city_prefix}.err
EOF

    chmod +x "$script_name"
    echo "Created script: $script_name"
    ((script_count++))
done

# Create the master script to run them all
RUN_SCRIPT_PATH="./$RUN_SCRIPT"
echo "#!/bin/bash" > "$RUN_SCRIPT_PATH"

for script in "$SCRIPTS_DIR"/*.sh; do
    echo "sbatch $script" >> "$RUN_SCRIPT_PATH"
done

chmod +x "$RUN_SCRIPT_PATH"

echo "Total scripts created: $script_count"
echo "Created SLURM batch submission script: $RUN_SCRIPT_PATH"
