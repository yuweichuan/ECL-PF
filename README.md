# ECL 3.0 (extension of ECL-PF)
ECL 3.0 is a highly sensitive cross-linking mass spectrometry data analysis tool, developed by Yu's group at HKUST (https://bioinformatics.hkust.edu.hk/). It is the extension of ECL-PF and is capable of analyzing both non-cleavable and cleavable cross-linking data. Additionally, a user-friendly graphical interface version of the Windows system has been designed for convenient usage.
## Usage
### Environment
1. The Windows system is required for ECL 3.0.
2. Fetch all the files in the repository.
3. Install Python (3.6 or above) and add it to the system path. (https://www.python.org/downloads/)
4. Install numpy, scipy, lxml, and pyteomics packages.
```bash
pip install numpy scipy lxml pyteomics
```
5. Test if the packages are successfully installed.
```bash
python
import numpy, scipy, lxml, pyteomics
```
you should observe something like below
```bash
C:\Users\zhouchen>python
Python 3.10.10 (tags/v3.10.10:aad5f6a, Feb  7 2023, 17:20:36) [MSC v.1929 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
>>> import numpy, scipy, lxml, pyteomics
>>>
```
### Quick start
You can start the data analysis by running the UI.exe. The program will output the CSV file for each task.
The input file formats and the parameter settings are listed below.
- MS spectra
  - mzXML format is required. To ensure the correct format, users can use MSConvert (https://proteowizard.sourceforge.io/download.html) to transfer the RAW file into mzXML by choosing the <strong>peak picking</strong> and <strong>zero samples</strong> functions as well as the <strong>32-bit precision</strong> and uncheck <strong>Use zlib compression</strong>. Both MS1 and MS2 spectra need to be retained. Details of manipulation can be found in the Tutorial.pdf under the ECL-PF repo.
- Protein sequence database
  - FASTA format is required. It is a standardized format that you can download from https://www.uniprot.org/.
- Threads
  - Number of threads to run the data.
- Peptide mass
  - The minimum and maximum peptide mass as the restriction during the in silico protein digestion.
- Digestion rule
  - Enzyme to in silico digest the protein.
- Miss cleavage
  - Maximum missed cleaved sites allowed per peptide.
- MS1 tolerance
  - Precursor mass tolerance in MS1 spectra in ppm unit.
- MS2 tolerance
  - Fragment ions tolerance in MS2 spectra in Dalton unit (ppm unit in the cleavable module).
- FDR setting
  - False discovery rate setting under the target-decoy strategy.
- Cross-linker
  - Non-cleavable (cleavable) cross-linker selection for the specific data type. You can add, edit, or delete the cross-linker in the database. After changing the database, you need to click the <strong>update database</strong> button.
- Modifications
  - Fixed modifications and variable modifications for the specific data type. Typically, carbamidomethyl[C] is added for the fixed mods. and oxidation[M] is added for the variable mods. You can add, edit, or delete the modification in the database. After changing the database, you need to click the <strong>update mods.</strong> button.


A video tutorial can be found at https://youtu.be/PpZgbi8V2xI.

And the available testing data can be obtained from https://bioinformatics.hkust.edu.hk/Software/ECL_3.html
## Authors
czhouau@connect.ust.hk Chen Zhou

eeyu@ust.hk Weichuan Yu

## License
[MIT LICENSE]