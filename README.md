# Protein_Derived_Peptides

## How to use:
The 'compare_data_cross_species' function will be the main function to call and generate the image.
The function takes in three string arguments.
1. ref_file_name: the reference file, for example, "Homo sapiens (Human)_Ameloblastin12_length.csv".
2. The file name for 'target', which is a folder containning all of the data you are comparing to. For example, "target".
3. The name for the similarity matrix, for example, "twelveAAHAmat.csv".

1 and 3 needs to be in the directory of the program. 2 can be any folder within the directory.

### An example of the function call:
```
compare_data_cross_species("Homo sapiens (Human)_Ameloblastin12_length.csv", "target", "twelveAAHAmat.csv");
```
**Note: This will only work if you're in match.py.**

### If you're calling this function from another program (.py):
include this line at the beginning 
```
from match import compare_data_cross_species
```
A complete example would be 
```
from match import compare_data_cross_species
compare_data_cross_species("Homo sapiens (Human)_Ameloblastin12_length.csv", "target", "twelveAAHAmat.csv");
```

After processing, the function will generate the images in the running directory.
