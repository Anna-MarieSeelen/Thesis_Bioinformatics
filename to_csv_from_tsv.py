# Python program to convert .tsv file to .csv file
# importing re library
import re
from sys import argv

file=argv[1]
file_output=argv[2]

# reading given tsv file
with open(file, 'r') as myfile:
    with open(file_output, 'w') as csv_file:
        for line in myfile:
            # Replace every tab with comma
            fileContent = re.sub("\t", ",", line)

            # Writing into csv file
            csv_file.write(fileContent)

# output
print("Successfully made csv file")