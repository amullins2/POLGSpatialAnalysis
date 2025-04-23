Here's a suggested structure for your `README.md` file on GitHub that explains how to set up the project, including creating the directory structure, verifying it, and downloading necessary files:

---

# 🧬 Endocrine Spatial Transcriptomics Analysis

This repository contains the pipeline and scripts for analyzing mouse pancreas spatial transcriptomics data generated using the **Visium HD** platform. The goal is to analyze the **endocrine cell population** in the pancreas across several TMAs.

## Project Structure

Below is the required file structure for this project:

```
project/
├── fastq/
│   ├── TMA1_HE_S1_L001_R1_001.fastq.gz
│   └── ...
├── images/
│   ├── TMA1_HE.jpg
│   └── ...
├── jsons/
│   ├── TMA1_HE.json
│   └── ...
├── reference/
│   └── 10x_mm10_reference/
└── outputs/
    ├── spaceranger/
    └── microarray/
```

### Required Directories:
1. **`fastq/`**: Contains the raw FASTQ files for each TMA.
2. **`images/`**: Contains the images of the tissue sections for each TMA.
3. **`jsons/`**: Contains the `.json` files generated from **MicroarrayProcessor**.
4. **`reference/`**: Contains the reference genome files needed for **Space Ranger** (e.g., `10x_mm10_reference/`).
5. **`outputs/`**: Output directory for **Space Ranger** and **MicroarrayProcessor** results.

---

## Setup Instructions

### Step 1: Create the Project Directory Structure

Use the following commands to create the necessary directories for your project.

```bash
mkdir -p project/fastq project/images project/jsons project/reference project/outputs/spaceranger project/outputs/microarray
```

---

### Step 2: Download and Install `spaceranger`

First, you need to download and install **spaceranger**.

```bash
# Download spaceranger tarball
wget -O spaceranger.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-2.1.1.tar.gz"

# Extract the spaceranger tarball
tar -xzvf spaceranger.tar.gz

# Add spaceranger to your PATH
export PATH=$PWD/spaceranger-2.1.1:$PATH
```

If you need to make this permanent, add `export PATH=$PWD/spaceranger-2.1.1:$PATH` to your `.bashrc`.

---

### Step 3: Download 10x Mouse Reference Data

In the `reference/` directory, download the **mm10 reference data**:

```bash
cd project/reference

# Download the 10x mm10 reference
wget -O refdata-gex-mm10-2020-A.tar.gz "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz"

# Extract the reference data
tar -xvzf refdata-gex-mm10-2020-A.tar.gz
```

This will create the necessary `10x_mm10_reference/` directory.

---

## Verifying the File Structure

You can verify the directory structure by running the following Python code:

```python
import os

# Define the directory paths
project_dir = "project"
required_structure = {
    "fastq": ["TMA1_HE_S1_L001_R1_001.fastq.gz"],
    "images": ["TMA1_HE.jpg"],
    "jsons": ["TMA1_HE.json"],
    "reference": ["10x_mm10_reference"],
    "outputs": ["spaceranger", "microarray"]
}

# Function to verify the directory structure
def verify_structure(base_dir, structure):
    for dir_name, files in structure.items():
        dir_path = os.path.join(base_dir, dir_name)
        if os.path.isdir(dir_path):
            print(f"✅ Directory {dir_path} exists.")
            for file in files:
                file_path = os.path.join(dir_path, file)
                if os.path.isfile(file_path):
                    print(f"  ✅ File {file} exists in {dir_name}.")
                else:
                    print(f"  ❌ File {file} is missing in {dir_name}.")
        else:
            print(f"❌ Directory {dir_path} is missing.")

# Run the check
verify_structure(project_dir, required_structure)
```

Alternatively, you can verify the file structure using the following shell command:

```bash
tree project/
```

This will print the directory and file structure for you to manually inspect.

---

## Running the Analysis

After setting up the directory structure, downloading the required files, and verifying everything, you can proceed with running the **Space Ranger** analysis for spatial transcriptomics and other processing steps.

Refer to the provided Jupyter notebook (`endocrine_spatial_pipeline.ipynb`) for step-by-step analysis.

Let me know if you need additional instructions or setup help!

---

### Contributing

Feel free to fork this project, create pull requests, or open issues if you encounter any problems.

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

You can save this content in a file called `README.md` in your repository's root directory. Let me know if you'd like further assistance!
