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
- [metadata_properties.csv](metadata_properties.csv) - the list of metadata properties used to annotate 
  [clinical_properties.csv](metadata_properties.csv) and [clinical_vocabularies.csv](clinical_vocabularies.csv).
- [namespaces.csv](namespaces.csv) - A list of namespaces (disease types)
  used to form the natural key namespace for each term and vocabulary.

## Metadata data model

Each **clinical property** is identified by a three-part natural key formed by the
combination (`resource`, `namespace`, `key`)

** NOTE: the combination of (`namespace`,`key`) must be unique - 
a `namespace` and `key` combination may not be reused between resources. 

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


Each **controlled vocabulary term** is identified by a four part 
natural key formed by the combination (`resource`, `namespace`, `field_key`, `key`)
where:
- `field_key` is the key of the field (term).
- `key` is an identifier for the vocabulary term that is 
  unique for the clinical property (`resource`, `namespace`, `field_key`) context. 
- a `key` is typically created by "normalizing" the `title` and follows the rules
  for normalization described below.


### Key normalization

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

## Metadata Properties

The [metadata_properties.csv](metadata_properties.csv) file specifies the properties that are used to 
annotate each term and controlled vocabulary in the [clinical_properties.csv](metadata_properties.csv) 
and [clinical_vocabularies.csv](clinical_vocabularies.csv) files. 

## Validation

```sh
pip install .
```
### Run the validation tool:

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


### Run the summary import and validation tool:

```sh
brttool -d path-to-files/ -s
```

This imports and validates the "summary.xlsx" file:

- summary file is an excel format used to collect clinical property and clinical 
  vocabulary terms from domain experts
- this separates terms by namespace and resource, with clinical properties listed 
  in the first row and vocabulary terms listed in columns

Actions:
1. Read in the specified summary.xls file and create a 
  `clinical_properties_from_summary.csv` and `clinical_vocabularies_from_summary.csv`
2. Validate the created files
3. Also interprets:
   - single vocab value as a "prompt"
   - extra vocab separated by a blank line at the end as a "description"
4. Set the property `data_type` to `integer`, `float`, or `string` based on title name patterns.
