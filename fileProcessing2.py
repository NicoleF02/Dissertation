import re

def read_gml_file(file_path):
    with open(file_path, 'r') as file:
        data = file.read()
    return data

def parse_data(data):
    gene_entries = re.findall(r'node\s*\[\s*id (\d+).*?wt "(\d+)".*?fc "(\d+)".*?infomap "(\d+)".*?louvain "(\d+)".*?sgG1 "(\d+)".*?sgG2 "(\d+)".*?sgG5 "(\d+)".*?spectral "(\d+)".*?Trubetskoy2022broadcoding "(.*?)".*\s*\]', data, re.DOTALL)

    return gene_entries

def print_table(gene_entries):
    algorithms = ['wt', 'fc', 'infomap', 'louvain', 'sgG1', 'sgG2', 'sgG5', 'spectral']

    # Initialize variables for maximum communities and true broad coding
    max_communities = {alg: 0 for alg in algorithms}
    true_broad_coding = {alg: 0 for alg in algorithms}

    # Process each gene entry
    for entry in gene_entries:
        for i, alg in enumerate(algorithms):
            community_count = int(entry[i + 1])  # Skip gene id
            # Update maximum communities
            max_communities[alg] = max(max_communities[alg], community_count)

        broad_coding = entry[-1]  # Trubetskoy2022broadcoding
        if broad_coding.lower() == 'true':
            for i, alg in enumerate(algorithms):
                true_broad_coding[alg] = max(true_broad_coding[alg], int(entry[i + 1]))  # Skip gene id

    # Print the table
    print(f"{'Algorithm':<10}{'Max Communities':<20}{'True Broad Coding':<20}")
    for alg in algorithms:
        print(f"{alg:<10}{max_communities[alg]:<20}{true_broad_coding[alg]:<20}")


if __name__ == "__main__":
    file_path = "PostsynapticNetwork/FullPSDDBNetwork.gml"  # Replace with the actual file path
    data = read_gml_file(file_path)
    gene_entries = parse_data(data)
    print_table(gene_entries)
