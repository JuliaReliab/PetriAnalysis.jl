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

`spnp_example2_tangible.npz` is the same net through `gospn mark -t`, which runs a
different search: it vanishes the immediate markings during the reachability walk rather
than leaving them for the elimination. Fewer states reach the file (10 rather than 11,
and no `I0I0I` block), and the generator over the tangible markings must be the same.

The `mrspn_*.npz` were produced with **gospn 0.22.0 or later**, which is the first that
writes `gentrans` and `groupgen` — which general transition each `P<k>` block is, and
what governs each group. Without them a wrong distribution is invisible: a general block
is a 0/1 jump matrix.

0.20.0 is the first release that writes `place` and `mark<G>`, which is what lets a row
of a matrix be keyed on its marking -- the two implementations enumerate the state space
in different orders, so nothing can be compared position by position.
