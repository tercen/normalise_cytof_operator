# Normalise CyTOF

##### Description

The `Normalise CyTOF operator` is an operator to perform the Normalisation of mass
cytometry data.

##### Usage

Input projection|.
---|---
`row`        | channels 
`column`           | observations (rowId) + filename
`y-axis`        | measurement values


Input parameters|.
---|---
`Normalization`| Select normalisation beads ("dvs" or "beta")

Output relations|.
---|---
`value`        | Normalised measurements
`Diagnostic plot`        | Computed tables include a graph of the estimated compensation matrix.

##### Details

This operator uses the normalisation approach described in the [CATALYST R package](https://www.bioconductor.org/packages/devel/bioc/vignettes/CATALYST/inst/doc/preprocessing.html#normalization).

##### See Also

[debarcode_cytof_operator](https://github.com/tercen/debarcode_cytof_operator)
, [compensate_cytof_operator](https://github.com/tercen/compensate_cytof_operator)
