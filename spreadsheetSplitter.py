import pandas as pd
import ast
import re


def eval_list(values, indice):
    try:
        list_string = values[indice]
        cleaned_list_string = list_string.replace('""', '"')
        evalList = ast.literal_eval(cleaned_list_string)
    except:
        print(evalList)
        evalList = []
    return evalList


# Treat csv as a text file, we get all of the columns that have lists, find max length, add no. rows,
# We do it that way,
input_file = 'tempCSV/Enriched Clusters - SynGO.csv'

output_file = 'SynGO/Ora/transformed_file.csv'

new_rows = []

outputfile = open(output_file, "w")



avoidLines = [['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '\n'],
              ['alg', 'clNo', 'clustsize', 'no.Trubetskoy Broad', 'Enrichment', 'Broad pval', 'Broad padj',
               'no. Trubetskoy Priortised', 'Prior pval', 'Prior padj', 'Molecular Function', 'Molecular padj',
               'Biological Function', 'Biological padj', 'Diseases', 'Diseases padj\n']]

outputfile.write(", ".join(avoidLines[1]))

with open(input_file, 'r') as file:
    # Read the header
    header = file.readline().strip().split(',')

    # Find the indices of the relevant columns
    columns_of_interest = ['Molecular Function', 'Molecular padj', 'Biological Function', 'Biological padj', 'Diseases',
                           'Diseases padj']
    indices_of_interest = [header.index(col) for col in columns_of_interest]
    print(indices_of_interest)
    for line in file:
        # Splits all into different ones, molecular function and biological function (0, 2 in list) will be same
        # length as their respective pvals

        list_pattern = r"\[.*?\]"

        # Find all occurrences of lists in the line
        lists = re.findall(list_pattern, line)

        print(f"no of lists is {len(lists)}")
        # Replace the lists in the line with a placeholder
        for i, sublist in enumerate(lists):
            line = line.replace(sublist, f"__LIST__{i}__")

        # Split the line using commas
        values = line.split(',')

        # Replace the placeholders with the original lists
        new_values = []

        listLen = len(lists)

        for i, value in enumerate(values):
            if any(f"__LIST__{j}__" in value for j in range(listLen)):
                new_values.append(lists.pop(0))
            else:
                new_values.append(value)




        noCols = len(values)
        values = new_values

        if values in avoidLines:
            print("Line avoided")
            outputfile.write(", ".join(values))
            continue

        for indice in indices_of_interest:
            values[indice] = eval_list(values, indice)
            print("1231")
            print(values[indice])



        # WE have the maximum number of rows we need
        maxListLen = max(len(values[indices_of_interest[0]]), len(values[indices_of_interest[2]]))
        maxListLen = max(len(values[indices_of_interest[4]]), maxListLen)

        if maxListLen != 0:
            # Where we will iterate through the lists
            for i in range(0, maxListLen):
                # We will need to iterate through columns,
                line = ""

                # Going through each value for each col
                for colIndice in range(0, noCols):
                    if i != 0 and colIndice not in indices_of_interest:
                        line = line + ", "
                        # We only want add once at top
                        continue

                    elif colIndice not in indices_of_interest:
                        # Must be str as col of interest are the only ones with lists
                        line = line + ", " + str(values[colIndice])
                    else:
                        # Must be col of interest, no matter where in list
                        try:
                            # We try and add through the list there

                            listCol = values[colIndice]

                            appendValue = listCol[i]

                            if type(appendValue) is not float:
                                if ',' in appendValue:
                                    # Important to replace a comma as its going into a csv
                                    appendValue = appendValue.replace(",","<")
                                elif appendValue == "[]" or appendValue == "":
                                    print("1111")
                                    appendValue = ""

                            line = line + ', ' + str(appendValue)
                        except:
                            print("Lists not same length")

                if line.startswith(", "):
                    line = line[2:]
                outputfile.write((line + "\n"))

        else:
            line = ""
            for colIndice in range(0, noCols):
                line = line + ", " + str(values[colIndice])

            if line.startswith(", "):
                line = line[2:]
            outputfile.write((line + "\n"))

outputfile.close()


# RE-GEN PREV SPREADSHEETS DUE TO LACK OF CLEANING with resetting padj etc...