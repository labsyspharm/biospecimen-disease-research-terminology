# Tools for the Biospecimen Disease Research Terminology

# Installation


```sh
pip install .
```

## Validate

```sh
brttool -d path-to-files/
```

This imports and validates the serialized terminology files:

- `clinical_properties.csv`, `clinical_vocabularies.csv` 
- also uses the `metadata_properties.csv` file to determine the shape of the other files.

Actions:
1. Validate column structure using fields defined in `metadata_properties.csv`
2. Validate integer and boolean columns
3. Update ordinal column
4. Verify uniqueness of key columns and alternate key columns (using titles instead of keys)
5. Verify vocabulary terms are matched with property terms
6. Outputs a `property_filename`_`summary.xlsx` file that list each discovered
   (resource, namespace) set of properties in separate sheets. Each property is 
   listed as a column header, and each vocabulary is listed in the column values.
