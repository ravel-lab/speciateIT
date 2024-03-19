import argparse
import os

try:
    import pandas as pd
except ImportError:
    print("Required package pandas not available. Please install pandas using 'pip install pandas'.")
    exit()

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Process count table and classification files.")
    parser.add_argument("-c", "--count_table", help="Path to the count table CSV file", required=True)
    parser.add_argument("-s", "--classification", help="Path to the classification TXT file", required=True)
    args = parser.parse_args()

    # Extract file name without extension from the count table path
    count_table_name = os.path.splitext(os.path.basename(args.count_table))[0]

    # Read count table
    try:
        count_table = pd.read_csv(args.count_table)
        count_table = count_table.rename(columns={"Unnamed: 0": "sampleID"})
    except FileNotFoundError:
        print("Count table file not found. Please provide count table using -c.")
        return

    # Read classification file
    try:
        classification = pd.read_csv(args.classification, sep='\t', header=None)
        classification.columns = ['ASV', 'Species', 'Score1', 'Score2']
    except FileNotFoundError:
        print("Classification file not found. Please provide speciateIT classification file using -s.")
        return

    # Create a dictionary with keys = classification['ASV'] and values = classification['Species']
    result_dict = classification[['ASV', 'Species']].set_index('ASV').to_dict()['Species']

    # Rename columns from count_table according to classification file
    df = count_table.rename(columns=result_dict)

    # Create a new column - sum the rows
    df['read_count'] = df.drop('sampleID', axis=1).sum(axis=1)

    # Group columns with the same Species name
    grouped_df = df.groupby(df.columns, axis=1).sum()

    # Calculate column sums
    column_sums = grouped_df.drop(columns=['sampleID', 'read_count']).sum()

    # Sort columns by descending sums
    sorted_columns = column_sums.sort_values(ascending=False).index.tolist()

    # Reorder the columns: 1st is SampleID, 2nd is read_count, others are species sorted by column sums
    new_order = ['sampleID', 'read_count'] + sorted_columns
    grouped_df = grouped_df[new_order]

    # Construct the output file name based on input count table name
    output_file_name = f"{count_table_name}_speciateIT.csv"

    # Save as CSV file
    grouped_df.to_csv(output_file_name, index=False)

    print(f"Output saved to {output_file_name}")

if __name__ == "__main__":
    main()
