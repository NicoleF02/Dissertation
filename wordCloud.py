import pandas as pd
import matplotlib.pyplot as plt
from wordcloud import WordCloud
import matplotlib.cm as cm

# We want to read in from temp csv the enriched clusters

# columns
# #Index(['alg', 'clNo', 'clustsize', 'no.Trubetskoy Broad', 'Enrichment',
# 'Broad pval', 'Broad padj', 'no. Trubetskoy Priortised', 'Prior pval',
# 'Prior padj', 'Molecular Function', 'Molecular padj',
# 'Biological Function', 'Biological padj', 'Diseases', 'Diseases padj'],
# dtype='object')

def createWordCloud(dataframe, column, column_pval):
    # weight the words by the pval and frequency
    # create a dictionary of words
    word_dict = {}

    for index, row in dataframe.iterrows():
        try:
            listWords = row[column].strip('[]').split(', ')
            pvalList = row[column_pval].strip('[]').split(', ')
            for i, word in enumerate(listWords):
                pval = float(pvalList[i])
                if word not in word_dict:
                    word_dict[word] = 1 - pval

                else:
                    word_dict[word] = word_dict[word] + 1 - pval
        except:
            continue

    wordcloud = WordCloud(width=800, height=500, colormap="tab10", background_color='white').generate_from_frequencies(
        word_dict)
    plt.figure(figsize=(15, 10))
    plt.imshow(wordcloud, interpolation="bilinear")
    plt.axis('off')
    plt.savefig('wordcloudFullDiseases.png')



df = pd.read_csv('tempCSV/Enriched Clusters - Full.csv')
print(df.columns)

# We want to create a word cloud of the enriched clusters, remove those with negative Enrichment

df = df[df['Enrichment'] == 'Positive']

createWordCloud(df, 'Diseases', 'Diseases padj')
