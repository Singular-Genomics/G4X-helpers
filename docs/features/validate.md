<br>

# `validate`
#### Validate G4X-data file and folder structure

Checks whether a G4X-data directory contains the expected files, folders, and table schemas for the detected assay type. By default, `validate` reports on both required raw inputs and secondary outputs such as single-cell and viewer artifacts.

---

## Usage
![`g4x-helpers validate --help`](../img/validate-help.svg)

--8<-- "_partials/args_optns.md"

---

--8<-- "_partials/arg_g4x_data.md"

---

### `--raw-only` / `-r`
_type_ : <span class="acc-2-code">`flag`</span>  
_default_ : `not set`

> Validate only the minimum raw data required for G4X-helpers operation.
> Use this when checking whether a sample can be processed before downstream outputs have been generated.

<br>
--8<-- "_core/_partials/end_cap.md"
