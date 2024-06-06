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
- a `key` is a normalized identifier for the term that is unique
  for the `namespace` context. 
- a `key` is typically created by "normalizing" the `title` (see below).

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
- `key` is a normalized identifier for the vocabulary term that is 
  unique for the property (`resource`, `namespace`, `field_key`) context. 
- a `key` is typically created by "normalizing" the `title` (see below).


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

TODO: include a validation script


