import pandas as pd
import matplotlib.pyplot as plt
from wordcloud import WordCloud
import os
import re

# We want to read in from temp csv the enriched clusters

# columns
# #Index(['alg', 'clNo', 'clustsize', 'no.Trubetskoy Broad', 'Enrichment',
# 'Broad pval', 'Broad padj', 'no. Trubetskoy Priortised', 'Prior pval',
# 'Prior padj', 'Molecular Function', 'Molecular padj',
# 'Biological Function', 'Biological padj', 'Diseases', 'Diseases padj'],
# dtype='object')

def read_csv_files(folder_path):
    combined_df = pd.DataFrame()
    # Set folder_path from string to directoary
    folder_path = os.path.join(os.getcwd(), folder_path)

    files = os.listdir(folder_path)

    for filename in files:
        if filename.endswith(".csv") and filename.startswith("enriched"):
            file_path = os.path.join(folder_path, filename)
            try:
                df1 = pd.read_csv(file_path)
                combined_df = pd.concat([combined_df, df1], ignore_index=True)
            except:
                print(f"Error reading file {file_path}")
    return combined_df


def print_latex(dictonaryWords, dictFrequency, column):
    latex_table = "\\subsubsection{" + column + "}\n"
    latex_table += "\\begin{longtable}[width=\\linewidth]{p{7.24cm} c p{4cm}}\n"
    latex_table += "\\caption{Molecular Function Frequencies and Weights}\n"
    latex_table += "\\label{tab:molecular-function-Full} \\\\\n"
    latex_table += "\\hline\n"
    latex_table += "\\textbf{Molecular Function} & \\textbf{Frequency} & \\textbf{Weights} \\\\ \n"
    latex_table += "\\hline\n"

    # Add both into a dataframe, joined on the keys
    df = pd.DataFrame(list(dictonaryWords.items()), columns=['Word', 'Weight'])
    df2 = pd.DataFrame(list(dictFrequency.items()), columns=['Word', 'Frequency'])
    df = df.merge(df2, on='Word')

    # Sort by weight descending then frequency
    df = df.sort_values(by=['Weight', 'Frequency'], ascending=False)

    for index, row in df.iterrows():
        table_row = f"{row['Word']} & {row['Frequency']} & {row['Weight']} \\\\ \\hline \n".strip('"').strip(
            "'").strip("\'").replace('_', " ")
        latex_table += table_row
    latex_table += "\\hline\n\\caption{Word Frequencies" + network + "}\n\\label{tab:word-frequencies-" + network + column.strip(' ') + "}\n\\end{longtable}"
    print(latex_table)


def createWordCloud(dataframe, column, column_pval):
    # weight the words by the pval and frequency
    # create a dictionary of words
    word_dict = {}
    dict_frequency = {}

    for index, row in dataframe.iterrows():

        try:
            # lots of splits to account for all the different ways, this also removes
            listWords = row[column].replace('"', "'")
            listWords = listWords.strip('[]')

            # Removes first and last char, gets rid of ' at start and end, can't just remove all ' as part of some names
            listWords = listWords[1:-1]
            listWords = listWords.split('\', \'')


            pvalList = row[column_pval].strip('[]').split(', ')
            for i, word in enumerate(listWords):
                pval = float(pvalList[i])


                if word not in word_dict:
                    word_dict[word] = 1 - pval
                    dict_frequency[word] = 1

                else:
                    word_dict[word] = word_dict[word] + 1 - pval
                    dict_frequency[word] = dict_frequency[word] + 1
        except:
            continue


    if column == "Diseases":
        wordcloud = WordCloud(width=500, height=300, colormap="tab10", background_color='white').generate_from_frequencies(
            word_dict)
        plt.figure(figsize=(10, 5))
    else:
        wordcloud = WordCloud(width=700, height=500, colormap="tab10", background_color='white').generate_from_frequencies(
            word_dict)
        plt.figure(figsize=(15, 10))

    plt.tight_layout(pad=0)
    plt.imshow(wordcloud, interpolation="bilinear")
    plt.axis('off')
    # plt.savefig(f'wordcloud{network}{column.strip(" ")}.png')
    wordcloud.to_file(f'wordcloud{network}{column.strip(" ")}.png')

    print_latex(word_dict, dict_frequency, column)


df = read_csv_files("Consensus/algorithmSummary")

# print(df.columns)

# We want to create a word cloud of the enriched clusters, remove those with negative Enrichment
#
#

columns = ['Molecular Function', 'Biological Function', 'Diseases', "SynGO"]
colsPadj = ['Molecular padj', 'Biological padj', 'Diseases padj', "SynGO padj"]

network = "Consensus"

for i in range(0, 3):
    print(f"For Network {network}")
    df = df[df['Enrichment'] == 'Positive']

    print("\\subsection{WordCloud Tables:" + network + "}")
    for j, column in enumerate(columns):
        createWordCloud(df, column, colsPadj[j])
        print("\n-------------------\n")

    print("\n\n\n")

    if i == 0:
        network = "SynGO"
        df = read_csv_files("SynGO/Ora/Enriched/algorithmSummary")

    elif i == 1:
        network = "Full"
        df = read_csv_files("Ora/Enriched/algorithmSummary")
