# SCAN paper analysis scripts

If you are using this repo please cite:

```
Gordon, E.M., Chauvin, R.J., Van, A.N., Rajesh, A., Nielsen, A., Newbold, D.J., Lynch, C.J., Seider, N.A., 
Krimmel, S.R., Scheidter, K.M., et al. (2023). A somato-cognitive action network alternates with effector regions 
in motor cortex. Nature. 10.1038/s41586-023-05964-2.
```

## Usage

This repo is provided *as-is* and may not work out-of-the-box for your data. Some scripts have 
hard-coded paths and may need to be adjusted to work on your system. For assistance, please contact
the authors.

## Dependencies

This repo uses [wbsurfer](https://gitlab.com/vanandrew/wbsurfer), which is provided as a submodule under the `extern` directory for convenience. To download it simply call:

```
git submodule update --init
```

This will clone the repo under `extern/wbsurfer`.
