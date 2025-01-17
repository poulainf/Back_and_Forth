#!/bin/bash

# Input: Subname provided as an argument
subname1=$1

# Validate input argument
if [[ -z $subname1 ]]; then
    echo "Error: Please provide a subname as the first argument."
    exit 1
fi

# Define directories
Work_dir="The_work_${subname1}"
Res_dir="The_results_${subname1}"

# Navigate to the working directory
if [[ -d $Work_dir ]]; then
    cd "$Work_dir" || { echo "Error: Cannot enter $Work_dir"; exit 1; }
else
    echo "Error: Working directory $Work_dir does not exist."
    exit 1
fi

mkdir -p "$Res_dir"

echo "Compilation of matrices and files..."

# Function to process matrices
process_matrices() {
    local pattern="$1"
    local output_prefix="$2"
    local matrix_suffix="$3"
    local rate_type="$4"

    echo "Processing $pattern..."
    for file in $(ls -lrt "Final_Combine_*_sub_selected_${pattern}" 2>/dev/null | awk '{print $NF}'); do
        Make_a_combined_matrix_universal.pl "$file" "$matrix_suffix" "$rate_type"
    done

    combine_matrixes.sh . "Combined_${output_prefix}_" "$matrix_suffix"
    cp "Final_Combine_${matrix_suffix}_matrix.txt" "./${Res_dir}/" 2>/dev/null
}

# Process Snp_matrix_192
process_matrices "matrix.txt" "Snp_matrix_192" "Snp_matrix_192" "NA"

# Consolidate results
consolidate_results() {
    local pattern="$1"
    local output_file="$2"

    local first_file
    first_file=$(ls -t "$pattern" 2>/dev/null | head -n 1)

    if [[ -n $first_file ]]; then
        head -n 1 "$first_file" > "./${Res_dir}/${output_file}"
        for file in $pattern; do
            grep -va "#" "$file" >> "./${Res_dir}/${output_file}"
        done
    else
        echo "Warning: No files matching $pattern found for consolidation."
    fi
}

# Consolidate specific result types
consolidate_results "*_compilation_single.txt" "ID_mono_file.txt"
consolidate_results "*_Comparaison_sub_rate.txt" "Comp_sub_rate.txt"
consolidate_results "*_compilation_TxS_mono.txt" "ID_mono_file_TxS.txt"
consolidate_results "*_compilation_sim_scores.txt" "Full_SIM_score.txt"
consolidate_results "*_matrix_POSITION_PING-PONG_full.txt" "Full_Pingage_posage.txt"

# Clean up Full Pingage Posage
if [[ -f "./${Res_dir}/Full_Pingage_posage.txt" ]]; then
    perl -pe "s/\t+/\ /g" "./${Res_dir}/Full_Pingage_posage.txt" > "./${Res_dir}/Full_Pingage_posage_treated.txt"
    rm "./${Res_dir}/Full_Pingage_posage.txt"
fi

consolidate_results "*_compilation_sim_scores2.txt" "Full_SIM_score2.txt"

echo "Finished processing matrices and results."

