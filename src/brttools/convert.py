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

    # # Blank line signals the end of the data
    # final_row = len(dataframe)
    # for i, row in dataframe.iterrows():
    #     if row.dropna().empty:
    #         final_row = i
    #         break
    # dataframe = dataframe[0:final_row]
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
        "description",
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
        "description",
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
            errors[tuple(row[k] for k in row_keys)][
                field_name
            ] = f"int parse error: {raw_val}"

    def get_boolean_val(_row_keys, field_name, row):
        raw_val = row[field_name]
        if not raw_val:
            return False
        if raw_val is True:
            return raw_val
        elif raw_val is False:
            return raw_val
        raw_val = str(raw_val)
        if raw_val.lower() in {"true", "t", "1"}:
            return True
        else:
            return False

    # Perform a very minimal parsing (int and bool only; all that is needed)
    parsers = {
        SCHEMA.DataType.integer: get_int_val,
        SCHEMA.DataType.boolean: get_boolean_val,
    }
    for field_name, field_info in field_schema_map.items():
        parser = parsers.get(field_info["data_type"])
        if parser:
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
    metadata_properties_df = cleanup_dataframe(pd.read_csv(filename, sep="\t"))
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
    PROP = SCHEMA.MetadataProperty
    CP = SCHEMA.ClinicalProperty
    CV = SCHEMA.ClinicalVocabulary
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
        raise Exception(f"Missing vocabulary matches: \n{unmatched}")

    property_df = property_df[
        sorted(
            property_field_schema.keys(),
            key=lambda k: property_field_schema[k][PROP.ordinal],
        )
    ]
    vocab_df = vocab_df[
        sorted(
            vocab_field_schema.keys(), key=lambda k: vocab_field_schema[k][PROP.ordinal]
        )
    ]
    return property_df, vocab_df


def create_summary(property_df, vocab_df):
    CP = SCHEMA.ClinicalProperty
    CV = SCHEMA.ClinicalVocabulary

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

    return output_resource_namespace_dataframes


def read_summary(summary_dataframes, property_field_schema, vocab_field_schema):
    """
    Read in the summary format and create a property_df and vocab_df that
    conform to the specification
    :param summary_df:
    :return: properties_df, vocabs_df
    """
    PROP = SCHEMA.MetadataProperty
    CP = SCHEMA.ClinicalProperty
    CV = SCHEMA.ClinicalVocabulary
    DT = SCHEMA.DataType

    def determine_type_from_title(title):
        title_pattern_to_data_types = {
            "^number of": DT.integer,
            "(cm)": DT.float,
            "(mm)": DT.float,
            "year": DT.integer,
            "percent": DT.float,
        }
        for pattern, data_type in title_pattern_to_data_types.items():
            if re.search(pattern, title, flags=re.IGNORECASE):
                return data_type
            elif re.search(pattern, normalize(title), flags=re.IGNORECASE):
                return data_type

        return DT.string

    logger.info("resoure-namespaces: %r", summary_dataframes.keys())

    property_df = pd.DataFrame(
        columns=sorted(
            property_field_schema.keys(),
            key=lambda k: property_field_schema[k][PROP.ordinal],
        )
    )
    vocab_df = pd.DataFrame(
        columns=sorted(
            vocab_field_schema.keys(), key=lambda k: vocab_field_schema[k][PROP.ordinal]
        )
    )

    def trim_nulls_from_end_of_series(pd_series):
        if not pd_series.size:
            return pd_series
        last_index = pd_series.last_valid_index()
        if last_index is not None:
            return pd_series.iloc[: last_index + 1]
        # otherwise, all values are None, so trim them
        return pd_series.dropna()

    for resource_name, df in summary_dataframes.items():
        rn = [x.strip() for x in resource_name.split("-")]
        assert (
            len(rn) == 2
        ), f"resource name must be in the format 'resource-namespace': {resource_name}"

        resource = rn[0]
        namespace = rn[1]

        rn_property_df = pd.DataFrame(columns=property_field_schema.keys())
        rn_property_df[CP.title] = df.columns
        rn_property_df[CP.resource] = resource
        rn_property_df[CP.namespace] = namespace
        rn_property_df[CP.ordinal] = rn_property_df.reset_index().index + 1
        rn_property_df[CP.data_type] = rn_property_df[CP.title].map(
            determine_type_from_title
        )

        for property_name in df.columns:
            # df columns are null padded to longest vocab list
            col_series = trim_nulls_from_end_of_series(df[property_name])
            if not col_series.size:
                continue
            if col_series.size > 2 and col_series.iloc[-2] is None:
                logger.info(
                    "Detected a description for property: %s: %r",
                    property_name,
                    col_series.iloc[-1],
                )
                rn_property_df.loc[
                    rn_property_df[CP.title] == property_name, CP.description
                ] = col_series.iloc[-1]
                col_series = col_series.iloc[:-2]
                logger.debug("new col series: \n%s", col_series)
            if col_series.size == 1:
                logger.info(
                    "Detected a prompt for property: %s: %r",
                    property_name,
                    col_series.iloc[0],
                )
                rn_property_df.loc[
                    rn_property_df[CP.title] == property_name, CP.prompt
                ] = col_series.iloc[0]
                continue
            property_vocab_df = pd.DataFrame(columns=vocab_field_schema.keys())
            property_vocab_df[CV.title] = col_series
            property_vocab_df[CV.key] = col_series.apply(normalize)
            property_vocab_df[CV.field_key] = normalize(property_name)
            property_vocab_df[CV.resource] = resource
            property_vocab_df[CV.namespace] = namespace
            property_vocab_df[CV.ordinal] = col_series.index + 1

            vocab_df = pd.concat([vocab_df, property_vocab_df])

        property_df = pd.concat([property_df, rn_property_df])

    vocab_df[CV.key] = vocab_df[CV.title].apply(normalize)
    property_df[CP.key] = property_df[CP.title].apply(normalize)

    property_df = property_df.replace({np.nan: None})
    vocab_df = vocab_df.replace({np.nan: None})

    return property_df, vocab_df
