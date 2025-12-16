### `--in-place`
_type_ : <span class="acc-2-code">`flag`</span>  
_default_  : `not set`

> By default, commands write their outputs to `<G4X-DATA>/g4x_helpers/<command>`, leaving the original data untouched. Using `--in-place` instead writes outputs directly into the specified `G4X-DATA` directory, which may be required when chaining certain commands.  

> :fontawesome-solid-warning: **Note:** this will _override_ any existing artifacts modified by the command.  
> Please refer to each feature’s documentation to see which files may be updated.
