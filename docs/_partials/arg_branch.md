### `--branch` / `-b`
_type_ : <span class="acc-2-code">`string`</span>  
_default_ : command name

> Selects the processed-data branch used for command output.
> By default, some commands  i.e `redemux` write to a command-specific branch:
> `<G4X-DATA>/g4x-helpers/<command>`.

> Use a custom branch name to write to `<G4X-DATA>/g4x-helpers/<branch>`.
> This is useful when you want to keep multiple alternative outputs, or when you want to chain commands through the same branch.

```bash
g4x-helpers redemux G4X-DATA --branch custom_process ...
g4x-helpers resegment G4X-DATA --branch custom_process ...

```

> Set `--branch _src_` to write directly into the original `G4X-DATA` directory.

!!! warning
    `--branch _src_` edits the source data in place and may overwrite files produced by the command.
