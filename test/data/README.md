# gospn fixtures

`*.npz` are gospn `mark` output, used by `test_gospn_crosscheck.jl` to compare gospn's
generator against this package's. The `.spn` next to each one is the definition it was
produced from, kept so the file can be regenerated:

```sh
gospn mark -i spnp_example1.spn -o spnp_example1.npz   # gospn 0.20.0 or later
```

The `mrspn_*.spn` fixtures are written for this crosscheck rather than taken from
gospn's `example/`: the two MRSPN nets bundled there (`raid6`, `raid10`) have
marking-dependent rates (`#Pn * lambda`), and `PetriStructure.jl` takes a constant rate,
so neither can be transcribed. Each `mrspn_*.spn` is the same net as the builder of the
matching name in `test_mrspn.jl`.

0.20.0 is the first release that writes `place` and `mark<G>`, which is what lets a row
of a matrix be keyed on its marking -- the two implementations enumerate the state space
in different orders, so nothing can be compared position by position.
