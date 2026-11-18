
# `make-base-config.R`: generate `base_config.yml` from the settings registry

## What it does

`base-config/settings-registry.yml` is the single source of the GLOWanalyses settings: one entry
per setting, with its section, type, default and description. This script writes
`base-config/base_config.yml` from it: every setting with its default, under its description as
comment lines. The file is what a user reads to learn the settings, and what
`GLOWpipeline::read_glow_config()` reads for the defaults. Because it is generated, the two
cannot differ.

The output depends on the registry alone. It holds no date and no path of the machine, so two
runs give the same bytes.

## How to run it

Run it from the GLOWanalyses directory, in an R environment that has the GLOWpipeline package
with `write_base_config()` (GLOWpipeline 0.3.0 or later).

```bash
Rscript base-config/make-base-config.R            # write base-config/base_config.yml
Rscript base-config/make-base-config.R --check    # compare the committed file with the registry
```

`--check` writes nothing. It prints `base_config.yml is current` and exits 0 when the committed
file equals what the registry generates. Otherwise it prints the first differing line and exits
1.

## When to run it

After every edit of `settings-registry.yml`: regenerate the base, then stage both files. Never
edit `base_config.yml` by hand, because the next run overwrites it.

## What can go wrong

- `Run from the GLOWanalyses root`: the script was started in another directory.
- `Registry ... has N invalid default(s)`: a default in the registry does not fit the type of its
  setting. The message names each one.
- `Registry ...: ... unknown key(s)`: an entry of the registry has a misspelled field.
