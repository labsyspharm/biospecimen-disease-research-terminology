# Utilities for working wiht the biospecimen-disease-research-terminology
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
    Cleanup a dataframe for use in data input:

    - strip whitespace in cells
    - convert empty cells to None
    - convert np.nan to None
    - Strip "(read only)" suffix from column titles

    :param dataframe: a dataframe with dtype=str, keep_default_na=False
    :return: a new dataframe
    """
    dataframe.rename(
        columns=lambda x: re.sub(r"^(.*)\s*\(read only\)\s*$", r"\1", x).strip(),
        inplace=True,
    )
    # Strip every cell of leading and trailing whitespace
    dataframe.replace(r"^\s*(.*?)\s*$", r"\1", regex=True, inplace=True)
    # Replace empty cells with NaN
    dataframe.replace(r"^\s*$", None, regex=True, inplace=True)

    # # Blank line signals the end of the data
    # NOTE: because summary dataframes will be ragged, this will not be used
    # see: def trim_nulls_from_end_of_series(colname, pd_series) instead
    # final_row = len(dataframe)
    # for i, row in dataframe.iterrows():
    #     if row.dropna().empty:
    #         final_row = i
    #         break
    # dataframe = dataframe[0:final_row]
    # Cleanup: convert nan to None
    new_df = pd.DataFrame(dataframe).replace({np.nan: None})
    return new_df


def trim_nulls_from_end_of_series(colname, pd_series):
    """
    Shorten series to remove None values from the end

    :return: a new series
    """
    if not pd_series.size:
        return pd_series
    last_index = 0
    for ix, val in pd_series.iloc[::-1].items():
        # logger.debug("ix, colname, val: %r", (ix, colname, val))
        if val is not None:
            last_index = ix
            break
    logger.debug("last_index: %r", (colname, last_index))
    if last_index is not None:
        return pd_series.iloc[: last_index + 1]
    # otherwise, all values are None, so trim them
    return pd_series[~pd_series.isnull()]


def normalize(token):
    """
    Normalize a token for use as a machine-readable identifier

    :param token:
    :return: lowercased, snake-cased token
    """
    token = str(token)
    token = token.lower()
    mapping = {
        "%": " percent ",
        "#": " number ",
        r"^\s?[-]{1}": " minus ",  # minus at beginning of line
    }
    for k, v in mapping.items():
        token = re.sub(rf"{k}", v, token)
    token = token.strip()
    token = re.sub("[^a-z0-9_]", "_", token)
    token = re.sub("_+", "_", token)
    token = re.sub("^_+", "", token)
    token = re.sub("_+$", "", token)

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
    MetadataPropertyFields = namedtuple(
        "MetadataPropertyFields", _metadata_property_fields
    )(*_metadata_property_fields)

    _metadata_vocabulary_fields = ["type", "value", "ordinal", "details"]
    MetadataVocabularyFields = namedtuple(
        "MetadataVocabularyFields", _metadata_vocabulary_fields
    )(*_metadata_vocabulary_fields)

    _metadata_types = [
        "clinical_property",
        "clinical_vocabulary",
        "metadata_resource",
        "metadata_namespace",
    ]
    MetadataType = namedtuple("MetadataType", _metadata_types)(*_metadata_types)

    _clinical_property_fields = [
        "entity_id",
        "key",
        "namespace",
        "resource",
        "title",
        "description",
        "data_type",
        "ordinal",
        "prompt",
        "is_deprecated",
        "is_provisional",
    ]
    ClinicalProperty = namedtuple("ClinicalProperty", _clinical_property_fields)(
        *_clinical_property_fields
    )

    _clinical_vocabulary_fields = [
        "entity_id",
        "key",
        "field_key",
        "namespace",
        "resource",
        "title",
        "description",
        "ordinal",
        "prompt",
        "is_deprecated",
        "is_provisional",
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


PF = SCHEMA.MetadataPropertyFields
VF = SCHEMA.MetadataVocabularyFields
CP = SCHEMA.ClinicalProperty
CV = SCHEMA.ClinicalVocabulary
DT = SCHEMA.DataType

# Bootstrap the metadata property schema
METADATA_PROPERTY_SCHEMA = {
    PF.name: {PF.name: PF.name, PF.data_type: DT.string},
    PF.resource_name: {PF.name: PF.resource_name, PF.data_type: DT.string},
    PF.ordinal: {PF.name: PF.ordinal, PF.data_type: DT.integer},
    PF.title: {PF.name: PF.title, PF.data_type: DT.string},
    PF.description: {PF.name: PF.description, PF.data_type: DT.string},
    PF.data_type: {PF.name: PF.data_type, PF.data_type: DT.string},
}


def determine_type_from_title(title):
    """
    Match `title` or `normalize(title)` to SCHEMA.DataType for certain patterns.

    :return: match SCHEMA.DataType or SCHEMA.DataType.string
    """
    title_pattern_to_data_types = {
        r"^number of": DT.integer,
        r"\(cm\)": DT.float,
        r"\(mm\)": DT.float,
        r"\byear(s?)\b": DT.integer,
        r"\bday(s?)\b": DT.integer,
        r"percent": DT.float,
        r"%": DT.float,
        r"\bnumber\b": DT.integer,
        r"#": DT.integer,
        r"select all": DT.arraystring,
    }
    found_type = DT.string
    for pattern, data_type in title_pattern_to_data_types.items():
        if re.search(pattern, title, flags=re.IGNORECASE):
            found_type = data_type
        elif re.search(pattern, normalize(title), flags=re.IGNORECASE):
            found_type = data_type

    logger.debug("determine_type_from_title: %r -> %r", title, found_type)
    return found_type


def _parse_dataframe(tablename, key_fields, field_schema_map, df):
    """
    Parse dataframe cells as defined by the field_schema_map `data_type` (modify in place)
    :param string tablename:
    :param list key_fields:
    :param field_schema_map: field definitions
    :param df:
    """
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

    # Perform a very minimal parsing (int and bool only)
    parsers = {
        SCHEMA.DataType.integer: get_int_val,
        SCHEMA.DataType.boolean: get_boolean_val,
    }
    for field_name, field_info in field_schema_map.items():
        parser = parsers.get(field_info[PF.data_type])
        if parser:
            df.loc[:, field_name] = df.apply(
                functools.partial(parser, key_fields, field_name), axis=1
            )
    if errors:
        raise Exception(f"Parse errors in {tablename}: \n{errors}")


def _validate_table(
    df,
    *,
    tablename,
    field_schema_map,
    key_fields,
    required_fields,
    alternate_keys=None,
    namespaces,
    resources,
):
    """
    Validate and parse the dataframe (in place)

    - check required_fields
    - parse columns specified in field_schema_map
    - check uniqueness of key_fields and alternate_keys
    - check valid namespace and resource

    """
    missing_columns = set(field_schema_map.keys()) - set(df.columns)
    if missing_columns:
        logger.info(f"columns in table: {df.columns!r}")
        logger.info(f"field_schema_map: {field_schema_map.keys()}")
        raise ValueError(f"{tablename} is missing columns: {missing_columns}")

    # check for required fields
    null_set = df[list(required_fields)].isnull()
    if not df[null_set.any(axis=1)].empty:
        raise Exception(
            f"Some of the required fields {required_fields} are null in {tablename}: "
            f"\n {df[null_set.any(axis=1)][list(required_fields)]}"
        )

    _parse_dataframe(tablename, key_fields, field_schema_map, df)

    # check that namespace and resource are valid
    namespace_values = {ns[VF.value] for ns in namespaces}
    unknown_namespaces = df[CP.namespace].isin(namespace_values) == False
    if unknown_namespaces.any():
        raise Exception(
            f"Unknown namespaces found in {tablename}: \n{df[unknown_namespaces][list(key_fields)]}"
        )

    resource_values = {r[VF.value] for r in resources}
    unknown_resources = df[CP.resource].isin(resource_values) == False
    if unknown_resources.any():
        raise Exception(
            f"Unknown resources found in {tablename}: \n{df[unknown_resources][list(key_fields)]}"
        )

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
            # NOTE: for alternate keys, don't consider null values (should be checked by required_fields)
            # - allows for case of null entity_id
            df_filtered = df[df.index.isin(df[list(alt_keyset)].dropna().index)]
            duplicated = df_filtered[list(alt_keyset)].duplicated(
                subset=list(alt_keyset), keep=False
            )
            if duplicated.any():
                raise Exception(
                    f"{(duplicated[duplicated]).size} Duplicates found in {tablename},"
                    f" for alternate keyset {alt_keyset}: "
                    f"\n{df_filtered[duplicated==True][list(set(key_fields + alt_keyset))]}"
                )


def read_metadata_properties_file(filename):
    """
    Read the field metadata from a file
    :param filename: path to the "metadata_properties.csv" file
    :return: dict{ field_key: {field_properties: }} a field schema map for
             clinical property fields, clinical vocab fields
    """
    # read the definitions of the fields
    # TODO: validate and parse using a schema
    # Using dtype=[unspecified] allows pandas to determine data types
    # in the future, use dtype=str and parse according
    metadata_properties_df = cleanup_dataframe(
        pd.read_csv(filename, sep="\t", dtype=str, keep_default_na=False)
    )
    _parse_dataframe(
        filename,
        [PF.name, PF.resource_name, PF.ordinal],
        METADATA_PROPERTY_SCHEMA,
        metadata_properties_df,
    )
    # Use ones-based index to match spreadsheet row
    metadata_properties_df.index += 2

    metadata_type_field_map = defaultdict(dict)
    for field_info in metadata_properties_df.to_dict(orient="records"):
        metadata_type_field_map[field_info[PF.resource_name]][
            field_info[PF.name]
        ] = field_info
    missing_metadata_definitions = set(SCHEMA.MetadataType) - set(
        metadata_type_field_map.keys()
    )
    if missing_metadata_definitions:
        raise ValueError(
            f"Missing metadata definitions: {missing_metadata_definitions}"
        )
    return metadata_type_field_map


def read_namespaces_file(filename, field_schema_map):
    """
    Read namespace definitions from a file

    :return list of the namespace definition dicts
    """
    namespaces_df = cleanup_dataframe(
        pd.read_csv(filename, sep="\t", dtype=str, keep_default_na=False)
    )
    _parse_dataframe(filename, [VF.value], field_schema_map, namespaces_df)
    logger.debug(
        "namespaces_df: \n %r",
        namespaces_df,
    )
    return namespaces_df.sort_values(by=[VF.ordinal, VF.value]).to_dict(
        orient="records"
    )


def read_resources_file(filename, field_schema_map):
    """
    Read resource definitions from a file

    :return list of the resource definition dicts
    """
    resources_df = cleanup_dataframe(
        pd.read_csv(filename, sep="\t", dtype=str, keep_default_na=False)
    )
    _parse_dataframe(filename, [VF.value], field_schema_map, resources_df)
    return resources_df.sort_values(by=[VF.ordinal, VF.value]).to_dict(orient="records")


def validate_specification(
    resources,
    namespaces,
    property_field_schema,
    vocab_field_schema,
    property_df,
    vocab_df,
):
    """
    Read, parse and validate the clinical properties and vocabulary dataframes


    :param list resources: Resource definitions to check
    :param list namespaces: Namespace definitions to check
    :param list property_field_schema: Property field definitions
    :param list vocab_field_schema: Vocabulary field definitions
    :param dataframe property_df: Clinical properties dataframe to validate
    :param dataframe vocab_df: Clinical vocabularies dataframe to validate
    :return: new property and vocabulary dataframes
    """
    resource_ordering = [r[VF.value] for r in resources]
    namespace_ordering = [n[VF.value] for n in namespaces]

    # Validate
    _validate_table(
        property_df,
        tablename=SCHEMA.MetadataType.clinical_property,
        field_schema_map=property_field_schema,
        key_fields=(CP.resource, CP.namespace, CP.key),
        alternate_keys=(
            (CP.entity_id,),
            (CP.namespace, CP.key),
            (CP.namespace, CP.title),
        ),
        required_fields=(CP.resource, CP.namespace, CP.key, CP.title, CP.data_type),
        namespaces=namespaces,
        resources=resources,
    )
    logger.debug("validated property_df: \n %r", property_df)
    # Set the ordinal:
    for rn, props in property_df.groupby([CP.resource, CP.namespace]):
        # set the ordinal if not set, to the alphabetical ordering
        if props[CP.ordinal].isnull().any():
            max_ordinal = props[CP.ordinal].max()
            max_ordinal = max_ordinal if not np.isnan(max_ordinal) else 0
            # - make an explicit copy, so pandas knows we know it is a copy
            to_update = props[props[CP.ordinal].isnull()].copy()
            to_update[:, CP.ordinal] = (
                to_update["title"].apply(normalize).rank(method="first", ascending=True)
                + max_ordinal
            )
            property_df.update(to_update)
    property_df[CP.ordinal] = property_df[CP.ordinal].astype(int)
    property_df[f"{CP.resource}_ord"] = property_df[CP.resource].map(
        lambda x: resource_ordering.index(x)
    )
    property_df[f"{CP.namespace}_ord"] = property_df[CP.namespace].map(
        lambda x: namespace_ordering.index(x)
    )

    # Validate the clincal vocabulary definitions
    _validate_table(
        vocab_df,
        tablename=SCHEMA.MetadataType.clinical_vocabulary,
        field_schema_map=vocab_field_schema,
        key_fields=(CV.resource, CV.namespace, CV.field_key, CV.key),
        alternate_keys=(
            (CV.entity_id,),
            (CV.resource, CV.namespace, CV.field_key, CV.title),
        ),
        required_fields=(CV.resource, CV.namespace, CV.field_key, CV.key, CV.title),
        resources=resources,
        namespaces=namespaces,
    )
    logger.debug("validated vocab_df: \n %r", vocab_df)

    # Set the ordinal:
    for rnp, vocabs in vocab_df.groupby([CV.resource, CV.namespace, CV.field_key]):
        # set the ordinal if not set, to the alphabetical ordering
        if vocabs[CP.ordinal].isnull().any():
            max_ordinal = vocabs[CP.ordinal].max()
            max_ordinal = max_ordinal if not np.isnan(max_ordinal) else 0
            to_update = vocabs[vocabs[CP.ordinal].isnull()].copy()
            to_update[:, CP.ordinal] = max_ordinal + to_update["title"].apply(
                normalize
            ).rank(method="first", ascending=True)
            vocab_df.update(to_update)
    vocab_df[CV.ordinal] = vocab_df[CV.ordinal].astype(int)
    vocab_df[f"{CV.resource}_ord"] = vocab_df[CV.resource].map(
        lambda x: resource_ordering.index(x)
    )
    vocab_df[f"{CV.namespace}_ord"] = vocab_df[CV.namespace].map(
        lambda x: namespace_ordering.index(x)
    )
    # use a left merge on vocabulary dataframe to find null matches
    vocab_join_to_properties = pd.merge(
        vocab_df[[CV.resource, CV.namespace, CV.field_key, CV.key]],
        property_df[[CP.resource, CP.namespace, CP.key, CP.ordinal]],
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
        raise Exception(
            f"Missing vocabulary matches: \n {unmatched[list((CV.resource, CV.namespace, CV.field_key, CV.key))]}"
        )
    vocab_df = vocab_df.reset_index()
    vocab_df[f"property_{CP.ordinal}"] = vocab_join_to_properties[CP.ordinal]
    property_df = property_df.sort_values(
        by=[f"{CP.resource}_ord", f"{CP.namespace}_ord", CP.ordinal]
    )
    property_df = property_df[
        sorted(
            property_field_schema.keys(),
            key=lambda k: property_field_schema[k][PF.ordinal],
        )
    ]
    vocab_df = vocab_df.sort_values(
        by=[
            f"{CV.resource}_ord",
            f"{CV.namespace}_ord",
            f"property_{CP.ordinal}",
            CV.field_key,
            CV.ordinal,
            CV.key,
        ]
    )
    vocab_df = vocab_df[
        sorted(
            vocab_field_schema.keys(), key=lambda k: vocab_field_schema[k][PF.ordinal]
        )
    ]
    return property_df, vocab_df


def create_summary(
    resources,
    namespaces,
    property_df,
    vocab_df,
    include_deprecated=False,
    include_provisional=False,
):
    """
    Create the specification summary workbook

    - one sheet per resource - namespace
    - clinical property terms (titles) listed in row one
    - clinical vocabulary terms (titles) listed under each property

    :param list resources: Resource definitions
    :param list namespaces: Namespace definitions
    :param dataframe property_df: Validated clinical properties dataframe
    :param dataframe vocab_df: Validated clinical vocabularies dataframe
    :return: a dict of dataframes representing a summary workbook
    """
    resource_ordering = [r[VF.value] for r in resources]
    namespace_ordering = [n[VF.value] for n in namespaces]

    if not include_deprecated:
        property_df = property_df[property_df[CP.is_deprecated] == False]
        vocab_df = vocab_df[vocab_df[CV.is_deprecated] == False]
    if not include_provisional:
        property_df = property_df[property_df[CP.is_provisional] == False]
        vocab_df = vocab_df[vocab_df[CV.is_provisional] == False]

    # use an inner join, then groupby to aggregate properties for vocabs
    property_to_vocabs = pd.merge(
        property_df[
            [
                CP.resource,
                CP.namespace,
                CP.key,
                CP.title,
                CP.ordinal,
                CP.prompt,
                CP.description,
            ]
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
    property_to_vocabs[f"{CP.resource}_ord"] = property_to_vocabs[CP.resource].map(
        lambda x: resource_ordering.index(x)
    )
    property_to_vocabs[f"{CP.namespace}_ord"] = property_to_vocabs[CP.namespace].map(
        lambda x: namespace_ordering.index(x)
    )
    for resource_namespace, rn_property_vocabs in property_to_vocabs.sort_values(
        by=[
            f"{CP.resource}_ord",
            f"{CP.namespace}_ord",
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

            if p_vocabs.reset_index().at[0, CP.description]:
                if len(property_vocab_map[pn]) > 1:
                    property_vocab_map[pn].append("")
                property_vocab_map[pn].append(
                    p_vocabs.reset_index().at[0, CP.description]
                )
        df = pd.DataFrame.from_dict(
            property_vocab_map,
            orient="index",
        ).transpose()
        output_resource_namespace_dataframes[resource_namespace] = df

    return output_resource_namespace_dataframes


def read_summary(
    namespaces, summary_dataframes, property_field_schema, vocab_field_schema
):
    """
    Parse a specification summary workbook into valid clinical property and
    vocabulary dataframes

    - each "resource - namespace" is a sheet title
    - each clinical property title is listed as a column header
    - each clinical vocabulary title for the property is listed in the column
    - special cases:
      - if only one vocab is listed, it is a "prompt" value
      - if the vocab column has a blank line and then a value at the end of the
        column, this is a "description"

    :param list resources: Resource definitions
    :param dict summary_dataframes: a dict of dataframes representing the summary workbook
    :param list property_field_schema: property field definitions
    :param list vocab_field_schema: vocabulary field definitions
    :return: properties_df, vocabs_df
    """

    logger.info("resoure-namespaces: %r", summary_dataframes.keys())

    property_df = pd.DataFrame(
        columns=sorted(
            property_field_schema.keys(),
            key=lambda k: property_field_schema[k][PF.ordinal],
        )
    )
    vocab_df = pd.DataFrame(
        columns=sorted(
            vocab_field_schema.keys(), key=lambda k: vocab_field_schema[k][PF.ordinal]
        )
    )

    def find_namespace_from_value_or_detail(raw_val):
        normalized_val = normalize(raw_val)
        for namespace_info in namespaces:
            namespace = namespace_info[VF.value]
            normalized_detail = normalize(namespace_info[VF.details])
            if normalized_val == namespace or normalized_val == normalized_detail:
                return namespace
            if normalized_val in namespace or normalized_val in normalized_detail:
                # allow inexact match; because openpyxl is only reading the first 31 chars
                logger.info(
                    f"matching truncated sheet name '{raw_val}' to '{namespace}'"
                )
                return namespace

    for resource_name, df in summary_dataframes.items():
        df = cleanup_dataframe(df)
        logger.info("processing sheet %r", resource_name)
        rn = [x.strip() for x in resource_name.split("-")]
        assert (
            len(rn) >= 2
        ), f"resource name must be in the format 'resource-namespace': value: '{resource_name}'"

        resource = rn[0]
        namespace = find_namespace_from_value_or_detail("_".join(rn[1:]))
        if not namespace:
            raise Exception(
                "Namespace for sheet name could not be found: %r", resource_name
            )
        rn_property_df = pd.DataFrame(columns=property_field_schema.keys())
        rn_property_df[CP.title] = df.columns
        rn_property_df[CP.resource] = resource
        rn_property_df[CP.namespace] = namespace
        rn_property_df[CP.ordinal] = rn_property_df.reset_index().index + 1

        for property_name in df.columns:
            # df columns are null padded to longest vocab list

            col_series = trim_nulls_from_end_of_series(property_name, df[property_name])
            if not col_series.size:
                continue
            if col_series.size > 1 and pd.isnull(col_series.iloc[-2]):
                logger.info(
                    "Detected a description for property: %r: %r",
                    (resource, namespace, property_name),
                    col_series.iloc[-1],
                )
                val = col_series.iloc[-1]
                if val:
                    rn_property_df.loc[
                        rn_property_df[CP.title] == property_name, CP.description
                    ] = val
                col_series = col_series.iloc[:-2]

            col_series = col_series[~col_series.astype(str).isin(["NaN", "nan"])]
            col_series = col_series[~col_series.isnull()]
            if col_series.size == 1:
                logger.info(
                    "Detected a prompt for property: %r: %r",
                    (resource, namespace, property_name),
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
        # logger.info("resource, namespace, rn_property_df: %r", (resource, namespace, rn_property_df))
        property_df = pd.concat([property_df, rn_property_df])

    vocab_df[CV.key] = vocab_df[CV.title].apply(normalize)
    property_df[CP.key] = property_df[CP.title].apply(normalize)

    property_df = property_df.replace({np.nan: None})
    vocab_df = vocab_df.replace({np.nan: None})

    return property_df, vocab_df


def merge_original_specification_into_new(key_fields, original_df, new_df):
    """
    Merge the values from the original specification dataframe to the new dataframe:

    - perform a left join on new_dataframe (ignore unmatched rows in original dataframe),
    - overwrite null values in the new dataframe (ignore overlapping values in original dataframe)

    """
    columns = new_df.columns
    original_df = original_df.reset_index().set_index(key_fields, verify_integrity=True)
    new_df = new_df.reset_index().set_index(key_fields, verify_integrity=True)
    new_df.update(original_df, overwrite=False)
    final = new_df.reset_index()[columns]

    return final


def merge_property_dataframes(original_df, new_df):
    merged = merge_original_specification_into_new(
        [CP.resource, CP.namespace, CP.key], original_df, new_df
    )
    logger.debug(
        "null data types after merge:\n%r",
        merged.loc[
            merged[CP.data_type].isnull(),
            [CP.resource, CP.namespace, CP.key, CP.data_type],
        ],
    )
    merged.loc[merged[CP.data_type].isnull(), CP.data_type] = merged.loc[
        merged[CP.data_type].isnull(), CP.title
    ].map(determine_type_from_title)
    logger.debug("merged property_df:\n%r", merged)
    return merged


def merge_vocab_dataframes(original_df, new_df):
    merged = merge_original_specification_into_new(
        [CV.resource, CV.namespace, CV.field_key, CV.key], original_df, new_df
    )
    logger.debug("merged vocab_df:\n%r", merged)
    return merged
