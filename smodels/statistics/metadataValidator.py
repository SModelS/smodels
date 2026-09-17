"""
.. module:: metadataValidator
   :synopsis: A simple facility to check if the meta data of the onnx
   file is looking correctly

.. moduleauthor:: OLLL Collaboration

"""

__all__ = [ "validateMetaData" ]

import json
from jsonschema import validate, ValidationError

def metadata_to_dict(metadata_props : "metadata" ) -> dict:
    """ convert onnx meta data to a dictionary """
    return {entry.key: entry.value for entry in metadata_props}

# Define a JSON Schema for the metadata
METADATA_SCHEMA = {
    "type": "object",
    "required": ["preprocessing", "channels", "obs_yields", "bkg_yields",
                 "run_config" ],
    "properties": {
        "preprocessing": {
            "type": "object", # after JSON parsing
            "required": ["features_pipeline"],
            "properties": {
                "features_pipeline": {"type": "array", "items": {"type": "string"}}
            }
        },
        "channels": {
            "type": "array",
            "items": {"type": "object"}
        },
        "obs_yields": {
            "type": "array",
            "items": {"type": "array"}
        },
        "bkg_yields": {
            "type": "array",
            "items": {"type": "array"}
        },
        "run_config": {
            "type": "string",
        },
    }
}

def validateMetaData(metadata_props : "onnx metadata" ) -> dict:
    """ validate the onnx metadata, check that all ingredients are there

    :raises ValueError: If something goes wrong
    :returns: python dictionary with content
    """
    # Convert to dict
    raw = metadata_to_dict(metadata_props)

    # Parse JSON string values into Python objects
    parsed = {}
    for key, value in raw.items():
        try:
            parsed[key] = json.loads(value)
        except (json.JSONDecodeError, TypeError):
            parsed[key] = value  # keep as string if not JSON

    # Validate against schema
    try:
        validate(instance=parsed, schema=METADATA_SCHEMA)
    except ValidationError as e:
        raise ValueError(f"Metadata validation failed: {e.message}") from e

    return parsed
