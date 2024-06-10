import functools
import os
from collections import defaultdict, namedtuple
from enum import Enum

import pandas as pd
import re
import numpy as np

import logging

logger = logging.getLogger(__name__)


def cleanup_dataframe(dataframe):
    """
    Cleanup a dataframe read from Excel for use in data input:
    - strip whitespace in cells
    - convert empty cells to NaN
    - enforce that a blank line signals the end of data
    - Strip "(read only)" suffix from column titles
    """
    dataframe.rename(
        columns=lambda x: re.sub(r"^(.*)\s*\(read only\)\s*$", r"\1", x).strip(),
        inplace=True,
    )
    # Strip every cell of leading and trailing whitespace
    dataframe.replace(r"^\s*(.*?)\s*$", r"\1", regex=True, inplace=True)
    # Replace empty cells with NaN
    dataframe.replace(r"^\s*$", np.nan, regex=True, inplace=True)

    # Blank line signals the end of the data
    final_row = len(dataframe)
    for i, row in dataframe.iterrows():
        if row.dropna().empty:
            final_row = i
            break
    dataframe = dataframe[0:final_row]
    # Cleanup: convert nan to None
    return pd.DataFrame(dataframe).replace({np.nan: None})


def normalize(token):
    token = str(token)
    token = token.lower()
    mapping = {
        "%": " percent ",
        "#": " number ",
        "^\s?[-]{1}": " minus ",
    }
    for k, v in mapping.items():
        token = re.sub(rf"{k}", v, token)
    token = token.strip()
    token = re.sub("[^a-zA-Z0-9_]", "_", token)
    token = re.sub("[_]+", "_", token)
    token = re.sub("^[_]+", "", token)
    token = re.sub("[_]+$", "", token)

    return token


class SCHEMA:
    """
    Models and Fields that are referenced in the code
    """

    _metadata_property_fields = [
        "name",
        "resource_name",
        "ordinal",
        "title",
        "description",
        "data_type",
    ]
    MetadataProperty = namedtuple("MetadataProperty", _metadata_property_fields)(
        *_metadata_property_fields
    )

    _metadata_types = ["clinical_property", "clinical_vocabulary"]
    MetadataType = namedtuple("MetadataType", _metadata_types)(*_metadata_types)

    _clinical_property_fields = [
        "key",
        "namespace",
        "resource",
        "title",
        "data_type",
        "ordinal",
        "prompt",
    ]
    ClinicalProperty = namedtuple("ClinicalProperty", _clinical_property_fields)(
        *_clinical_property_fields
    )

    _clinical_vocabulary_fields = [
        "key",
        "field_key",
        "namespace",
        "resource",
        "title",
        "ordinal",
        "prompt",
    ]
    ClinicalVocabulary = namedtuple("ClinicalVocabulary", _clinical_vocabulary_fields)(
        *_clinical_vocabulary_fields
    )

    _data_types = [
        "string",
        "integer",
        "float",
        "boolean",
        "date",
        "arraystring",
        "arrayint",
    ]
    DataType = namedtuple("DataType", _data_types)(*_data_types)


def _validate_table(
    df, *, tablename, field_schema_map, key_fields, required_fields, alternate_keys=None
):
    missing_columns = set(field_schema_map) - set(df.columns)
    if missing_columns:
        raise ValueError(f"{tablename} is missing columns: {missing_columns}")

    # check for required fields
    null_set = df[list(required_fields)].isnull()
    if not df[null_set.any(axis=1)].empty:
        raise Exception(
            f"required fields are null in {tablename}: {required_fields}, "
            f"\n {df[null_set.any(axis=1)][list(required_fields)]}"
        )

    # parse known field data types
    errors = defaultdict(dict)

    def get_int_val(row_keys, field_name, row):
        raw_val = row[field_name]
        if not raw_val:
            return None
        try:
            return int(float(raw_val))
        except:
            errors[tuple(row[k] for k in key_fields)][
                field_name
            ] = f"int parse error: {raw_val}"

    def get_boolean_val(row_keys, field_name, row):
        raw_val = row[field_name]
        if not raw_val:
            return False
        if raw_val.lower() in {"true", "t", "1"}:
            return True
        else:
            return False

    parsers = {
        SCHEMA.DataType.integer: get_int_val,
        SCHEMA.DataType.boolean: get_boolean_val,
    }
    for field_name, field_info in field_schema_map.items():
        parser = parsers.get(field_info["data_type"])
        if parser:
            logger.info("apply parser for field: %r", field_name)
            df[field_name] = df.apply(
                functools.partial(parser, key_fields, field_name), axis=1
            )
    if errors:
        raise Exception(f"errors in {tablename}: \n{errors}")

    # check that keys are unique
    duplicated = df.duplicated(
        subset=list(key_fields), keep=False
    )  # keep=False, show all
    if duplicated.any():
        raise Exception(
            f"Duplicates found in {tablename}: \n{df[duplicated][list(key_fields)]}"
        )
    if alternate_keys:
        for alt_keyset in alternate_keys:
            duplicated = df.duplicated(subset=list(alt_keyset), keep=False)
            if duplicated.any():
                raise Exception(
                    f"Duplicates found in {tablename}: \n{df[duplicated][list(alt_keyset)]}"
                )


def read_metadata_properties_file(filename):
    """
    Read the field metadata from the properties file
    :param filename: path to the "metadata_properties.csv" file
    :return: field schema map { field_key: {field_properties: }} for
             clinical property fields, clinical vocab fields
    """

    FIELDS = SCHEMA.MetadataProperty
    # read the definitions of the fields
    metadata_properties_df = cleanup_dataframe(
        pd.read_csv(filename, sep="\t")
    )
    # Use ones-based index to match spreadsheet row
    metadata_properties_df.index += 2

    metadata_type_field_map = defaultdict(dict)
    for field_info in metadata_properties_df.to_dict(orient="records"):
        metadata_type_field_map[field_info[FIELDS.resource_name]][
            field_info[FIELDS.name]
        ] = field_info
    missing_metadata_definitions = set(SCHEMA.MetadataType) - set(
        metadata_type_field_map.keys()
    )
    if missing_metadata_definitions:
        raise ValueError(
            f"Missing metadata definitions: {missing_metadata_definitions}"
        )
    return (
        metadata_type_field_map[SCHEMA.MetadataType.clinical_property],
        metadata_type_field_map[SCHEMA.MetadataType.clinical_vocabulary],
    )


def validate_specification(
    property_field_schema, vocab_field_schema, property_df, vocab_df
):
    """
    Read and validate the clinical properties and vocabs; using field definitions from
    the metadata_properties_file
    """
    FIELDS = SCHEMA.MetadataProperty
    CP = SCHEMA.ClinicalProperty
    CV = SCHEMA.ClinicalVocabulary
    print("wooop woop!!!")
    # Validate
    _validate_table(
        property_df,
        tablename=SCHEMA.MetadataType.clinical_property,
        field_schema_map=property_field_schema,
        key_fields=(CP.resource, CP.namespace, CP.key),
        alternate_keys=((CP.namespace, CP.key), (CP.namespace, CP.title)),
        required_fields=(CP.resource, CP.namespace, CP.key, CP.title, CP.data_type),
    )

    # Set the ordinal:
    for rn, props in property_df.groupby([CP.resource, CP.namespace]):
        # set the ordinal if not set, to the alphabetical ordering
        if props[CP.ordinal].isnull().any():
            max_ordinal = props[CP.ordinal].max()
            max_ordinal = max_ordinal if not np.isnan(max_ordinal) else 0
            # - make an explicit copy, so pandas knows we know it is a copy
            to_update = props[props[CP.ordinal].isnull()].copy()
            to_update[CP.ordinal] = (
                to_update["title"].apply(normalize).rank(method="first", ascending=True)
                + max_ordinal
            )
            property_df.update(to_update)
    property_df[CP.ordinal] = property_df[CP.ordinal].astype(int)
    property_df = property_df.sort_values(by=[CP.resource, CP.namespace, CP.ordinal])

    # Validate the clincal vocabulary definitions
    _validate_table(
        vocab_df,
        tablename=SCHEMA.MetadataType.clinical_vocabulary,
        field_schema_map=vocab_field_schema,
        key_fields=(CV.resource, CV.namespace, CV.field_key, CV.key),
        alternate_keys=((CV.resource, CV.namespace, CV.field_key, CV.title),),
        required_fields=(CV.resource, CV.namespace, CV.field_key, CV.key, CV.title),
    )

    # Set the ordinal:
    for rnp, vocabs in vocab_df.groupby([CV.resource, CV.namespace, CV.field_key]):
        # set the ordinal if not set, to the alphabetical ordering
        if vocabs[CP.ordinal].isnull().any():
            max_ordinal = vocabs[CP.ordinal].max()
            max_ordinal = max_ordinal if not np.isnan(max_ordinal) else 0
            to_update = vocabs[vocabs[CP.ordinal].isnull()].copy()
            to_update[CP.ordinal] = max_ordinal + to_update["title"].apply(
                normalize
            ).rank(method="first", ascending=True)
            vocab_df.update(to_update)
    vocab_df[CV.ordinal] = vocab_df[CV.ordinal].astype(int)
    vocab_df = vocab_df.sort_values(
        by=[CV.resource, CV.namespace, CV.field_key, CV.ordinal]
    )

    # use a left merge on vocabulary dataframe to find null matches
    vocab_join_to_properties = pd.merge(
        vocab_df[[CV.resource, CV.namespace, CV.field_key, CV.key]],
        property_df[[CP.resource, CP.namespace, CP.key]],
        left_on=(CV.resource, CV.namespace, CV.field_key),
        right_on=(CP.resource, CP.namespace, CP.key),
        how="left",
        suffixes=("_v", "_p"),
    )
    if vocab_join_to_properties[f"{CP.key}_p"].isnull().any():
        unmatched = vocab_df.reset_index()[
            vocab_join_to_properties[f"{CP.key}_p"].isnull()
        ]
        unmatched.index += 2
        raise Exception(f"Missing vocabulary matches in {vocab_file}: \n{unmatched}")

    # use an inner join, then groupby to aggregate properties for vocabs
    property_to_vocabs = pd.merge(
        property_df[
            [CP.resource, CP.namespace, CP.key, CP.title, CP.ordinal, CP.prompt]
        ],
        vocab_df[
            [CV.resource, CV.namespace, CV.field_key, CV.key, CV.title, CV.ordinal]
        ],
        left_on=(CP.resource, CP.namespace, CP.key),
        right_on=(CV.resource, CV.namespace, CV.field_key),
        how="outer",
        suffixes=("_p", "_v"),
    )

    # Create output summary dataframes, one per namespace
    # - summary dataframes list properties in columns to vocabularies in rows
    output_resource_namespace_dataframes = {}
    for resource_namespace, rn_property_vocabs in property_to_vocabs.sort_values(
        by=[
            CP.resource,
            CP.namespace,
            f"{CP.ordinal}_p",
            f"{CV.ordinal}_v",
        ]
    ).groupby([CP.resource, CP.namespace], sort=False):
        property_vocab_map = {}
        for property_name, p_vocabs in rn_property_vocabs.groupby(
            [f"{CP.title}_p"], sort=False
        ):
            pn = "".join(property_name)
            if len(p_vocabs) == 1:
                # List the prompt as the vocab instead if there are no vocabs
                property_vocab_map[pn] = p_vocabs[CP.prompt].tolist()
            else:
                property_vocab_map[pn] = p_vocabs[f"{CV.title}_v"].tolist()

        df = pd.DataFrame.from_dict(
            property_vocab_map,
            orient="index",
        ).transpose()
        output_resource_namespace_dataframes[resource_namespace] = df


    return property_df, vocab_df, output_resource_namespace_dataframes
