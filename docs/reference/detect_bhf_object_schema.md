# Detect the Schema of a BHF Object

Classifies current, legacy 0.3.0, unsupported, and non-BHF objects
without mutating or adapting them. Legacy classification is inferred
from the frozen 0.3.0 structural signature because those objects predate
explicit schema markers.

## Usage

``` r
detect_bhf_object_schema(x)
```

## Arguments

- x:

  An object to inspect.

## Value

A `bhf_schema_detection` list with object type, schema version, status,
contract identifier, inference flag, and reason.
