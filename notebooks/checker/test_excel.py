import sys
import pandas as pd

input_file = "duptest.xlsx"

output_file = "cleaned_duptest.xlsx"


# Read Excel
df = pd.read_excel(input_file)


print("Original rows:")
print(len(df))


# Remove duplicates
clean_df = df.drop_duplicates()


print("Rows after removing duplicates:")
print(len(clean_df))


# Save new Excel
clean_df.to_excel(
    output_file,
    index=False
)


print("Saved:")
print(output_file)