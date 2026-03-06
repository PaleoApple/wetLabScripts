import pandas as pd
import re

input_file = "SrtA_P94D_25C_c.xpk"
output_file = "SrtA_P94D_25C_c.xlsx"

# Columns to extract (we’ll add “label” ourselves)
columns_to_extract = ["label", "HN.P", "NH.P", "vol"]

# --- Step 1: Read the file ---
with open(input_file, "r", encoding="utf-8") as f:
    lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

# --- Step 2: Find header line ---
header_line = None
for i, line in enumerate(lines):
    if all(col in line for col in ["HN.P", "NH.P", "vol"]):
        header_line = line
        header_index = i
        break

if header_line is None:
    raise ValueError("Could not locate a valid header line in file.")

# Split header fields
headers = re.split(r"\s+", header_line.strip())

# --- Step 3: Extract data lines ---
data_lines = lines[header_index + 1 :]

rows = []
for line in data_lines:
    # Split line by whitespace
    parts = re.split(r"\s+", line.strip())

    # Find the label (the token in curly braces like {A_3})
    label_match = next((p for p in parts if re.match(r"\{[A-Za-z0-9_]+\}", p)), None)
    label_value = label_match.strip("{}") if label_match else None

    # Remove the label token from the list so data aligns with headers
    parts = [p for p in parts if not re.match(r"\{[A-Za-z0-9_]+\}", p)]

    # Align header/data counts (skip malformed lines)
    if len(parts) != len(headers):
        continue

    # Create a dictionary for this row
    row = dict(zip(headers, parts))
    row["label"] = label_value  # add extracted label
    rows.append(row)

# --- Step 4: Build DataFrame ---
df = pd.DataFrame(rows)

# --- Step 5: Convert to numeric where appropriate ---
for col in ["HN.P", "NH.P", "vol"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# --- Step 6: Write to Excel ---
df[columns_to_extract].to_excel(output_file, index=False)

print(f"✅ Extracted {columns_to_extract} → saved to {output_file}")
