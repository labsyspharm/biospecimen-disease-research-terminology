# Clinical disease terminology for annotation of biospecimen records

A collection of clinical properties and vocabularies, relevant within the 
context of a namespace, typically a disease type being researched, to be used to 
annotate a biospecimen, or a collection of biospecimens represented by a surgery event,
or patient record for biospecimens.

These clinical properties are selected from existing specifications, 
and are intended to be used in annotation of clinical data for research purposes 
in the [LSP Sample Tracking System](https://github.com/labsyspharm/ExperimentTracker).

All data recorded in the Sample Tracker are de-identified and contain no
personally identifiable information or personal health information (PHI). 
Accordingly, no clinical property values described by this terminology should contain such information.

The TSV (tab-separated value) files in this repository detail the clinical
fields (properties) and controlled vocabularies in the collection, as well
as metadata specifying how to annotate the clinical fields and vocabularies using
identifiers and descriptive properties.

## Files

- [clinical_properties.csv](metadata_properties.csv) - the collection of clinical properties (fields) to be used
  to annotate the patients, surgeries, and biospecimens.
- [clinical_vocabularies.csv](clinical_vocabularies.csv) - the collection of controlled vocabularies 
  for the clinical properties.
- [clinical_properties_summary.xlsx](clinical_properties_summary.xlsx) - (READ ONLY) a summary format, 
  provided for information and collaboration with domain experts. 
  - Excel format used to collect clinical property and clinical 
  vocabulary terms from domain experts,
  - Separate sheet for each namespace and resource, with clinical properties listed 
    in the first row and vocabulary terms listed in columns.
  - See workflow section below for more information.
- [metadata_properties.csv](metadata_properties.csv) - the list of metadata properties used to annotate 
  [clinical_properties.csv](metadata_properties.csv) and [clinical_vocabularies.csv](clinical_vocabularies.csv).
- [namespaces.csv](namespaces.csv) - A list of namespaces (disease types)
  used to form the natural key namespace for each term and vocabulary.
- [resources.csv](resources.csv) - A list of resources (patient, surgery, biospecimen)
  used to form the natural key with the namespace for each term and vocabulary.

## Definition and identification of terms

This terminology is organized such that each term is defined by a unique natural key.

A **clinical property term** is identified by a three-part natural key formed by the
combination [`resource`, `namespace`, `key`]

** NOTE: the combination of [`namespace`,`key`] must also be unique (a `namespace`
and `key` combination may not be reused between resources).

- a `resource` is the entity record to be annotated, one of:
  - `patient`
  - `surgery`
  - `biospecimen`
- a `namespace` is the identifier of the disease type being studied for the 
  clinical property or vocabulary, 
  one of the values in [namespaces.csv](namespaces.csv)
- a `key` is an identifier for the term that is unique
  for the `namespace` context. 
- a `key` is typically created by "normalizing" the `title` and follows the rules
  for normalization described below.

Additionally, each term is assigned a data type:
- `string` - a text string
- `integer` - a whole number
- `float` - a decimal number
- `date` - a date in the format YYYY-MM-DD
- `boolean` - a true/false value
- `arraystring` - a comma separated list of text strings
- `arrayint` - a comma separated list of whole numbers

Additionally, each clinical property term will be assigned a unique identifier when
it is registered in the LSP Sample Tracker, and will contain references to external
identifiers.

A **clinical vocabulary term** is identified by a four part 
natural key formed by the combination [`resource`, `namespace`, `field_key`, `key`]
where:
- `field_key` is the key of the **clinical property term** (field).
- `key` is an identifier for the vocabulary term that is 
  unique for the clinical property [`resource`, `namespace`, `field_key`] context. 
- a `key` is typically created by "normalizing" the `title` and follows the rules
  for normalization described below.

Additionally, each clinical vocabulary term will be assigned a unique identifier when
it is registered in the LSP Sample Tracker, and will contain references to external
identifiers.


### Note on key generation using normalization

The `key` for the clinical property and vocabulary entries are formed by 
lowercasing the title, removing non-alphanumeric characters, and replacing spaces with underscores.
- Example: "Tumor Grade" -> "tumor_grade"
- Example: "Histologic Grade (WHO/ISUP)" -> "histologic_grade_who_isup"
- Example: "% Dedifferentiated" -> "percent_dedifferentiated"
- allowed characters: [a-z0-9_]
- may not start or end with an underscore
- may not contain two consecutive underscores
- certain terms may be manually normalized to avoid conflicts
- other conventions may be used, such are replacing symbols such as
  `%` by `percent` or `#` by `number`

## Workflow and tools 

Three workflows are envisioned for updating the specification:
1. **Direct update**: edit and validate `clinical_properties.csv` and `clinical_vocabularies.csv`
2. **Merge summary file data**: new fields and vocabularies from a `summary.xlsx` file.
3. **Merge external data**: new fields, vocabularies, and updates from externally generated `clinical_properties.csv`
   and `clinical_vocabularies.csv` files, e.g. from the LSP sample tracking
   database.

The **brttools** package provides tools to enable these workflows.
### Install brttools package 

```sh
pip install .
```
### 1. Validation

Validate specification files and generate a summary file.

```sh
brttool -d path-to-files/
```

Requires the complete set of specification files in the BDRT repository:
- `clinical_properties.csv`, `clinical_vocabularies.csv` 
- `metadata_properties.csv`, `resources.csv`, and `namespaces.csv`.

Actions:
1. Validate: column structure and data type using fields defined in 
   `metadata_properties.csv`
2. Validate resources and namespaces.
3. Update the ordinal column
4. Enforce unique constraints: using key columns and alternate key columns 
   (titles instead of keys)
5. Verify vocabulary terms are matched with property terms
6. Output a `summary` (xlsx) file that lists each 
   (resource, namespace) set of properties in separate sheets. Each property is 
   listed as a column header, and each vocabulary is listed in the column values.


### 2. Merge summary data

Import, merge and validate a summary file.

```sh
brttool -d path-to-files/ -s path-to-summary-file
```

Actions:
1. Read summary file
2. Merge with existing `clinical_properties.csv` and `clinical_vocabularies.csv`
3. Validate
4. Output to `clinical_properties_from_summary.csv` and `clinical_vocabularies_from_summary.csv`
5. Other:
   - interpret a single vocab value as a "prompt"
   - interpret extra vocab separated by a blank line at the end as a "description"
   - Set the property `data_type` to `integer`, `float`, or `string` based on 
     title name patterns (iif `data_type` not set in existing specification file).

### 3. Merge specification files

Merge data from one specification file to another using natural keys.

```sh
brt_mergetool -f1 new_base_specification_file -f2 overlay_specification_file [-cp or -cv]
```

Merge file2 specification data into file1:
1. Left join file1 to file2 on natural keys [resource, namespace, key]
2. Preserve non-null file1 values, merge non-null file2 values.

(see: [pandas.DataFrame.update](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.update.html#pandas.DataFrame.update))

