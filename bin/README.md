# Windows binaries (64-bit, static, built with mingw-w64 from the sources in `scripts/`)

- `wall_recursion128.exe`: the wall by the label recursion with 128-bit labels (`scripts/wall_recursion128.c`).
- `edge_machine.exe`: the edge automaton in C (`scripts/edge_machine.c`).
- `wall_recursion.exe`: the 64-bit-label recursion (`scripts/wall_recursion.c`).

Run them from PowerShell in the repository folder. To continue the seed's wall from the checkpoint in
`results/wall_recursion_logs/frontier_checkpoint.txt` toward 10¹² (about eight hours; a checkpoint line every
few seconds), open a PowerShell window you can leave open and run

    cmd /c 'bin\wall_recursion128.exe bb28d1da3719d3cdbb28d1da3719d3cd 4efb996adbd644c74efb996adbd644c7 150800000000 1000000000000 "" 128 - > wall128_to_1e12.log'

(`cmd /c` so the log is plain text). Check progress from another window with `Get-Content wall128_to_1e12.log -Tail 3`.
Self-test of the edge machine: `bin\edge_machine.exe --selftest` (must print 0 mismatches).
Quick check of the recursion (should print the same "final labels" line within a second):

    bin\wall_recursion128.exe 88d1fb8488d1fb84 7cb53b6b7cb53b6b 90000 800000 "" 0
    -> final labels: a=c6306a59c6306a59c6306a59c6306a59 b=b44b9dd3b44b9dd3b44b9dd3b44b9dd3
