# OvFilter Test Dataset

This directory provides a miniature dataset and an automated script to test the OvFilter pipeline.

## Structure
- `sim_data/`: Place your extracted test dataset here (not included in the repository by default).
- `run_test.sh`: An automated bash script that runs the entire OvFilter pipeline (Minimap2 -> Jellyfish -> OvFilter).
- `evaluate_alignments.py`: A Python script that evaluates the output overlaps against the ground-truth MAF file.

## Download Test Dataset
Due to the large file size, the test dataset is hosted externally. Before running the test, please download and extract the dataset into the `test/` directory:

1. Download `sim_data.tar.gz` from Zenodo: `https://doi.org/10.5281/zenodo.20554435`
2. Extract the contents:
   ```bash
   tar -xzvf sim_data.tar.gz
   ```
   This will create a `sim_data/` directory containing subdirectories like `chr1/hifi/`. By default, the `run_test.sh` script is configured to test the `chr1_hifi_0001` dataset.

## How to Run
```bash
cd test
bash run_test.sh
```

## Output
The final verified overlaps will be generated at `data/ovfilter.paf`.
