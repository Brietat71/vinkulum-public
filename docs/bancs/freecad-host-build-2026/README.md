# Automated Linux FreeCAD host build

The [record.zip](record.zip) archive retains the complete build log, first and
final installer sources, build records, generated launcher and five native GUI
qualification runs. Its SHA-256 is in [sha256.json](sha256.json); each archived
file also has an internal SHA-256 entry. Cache, installed extension copies and
isolated user preferences are excluded; the packaged extension is retained.

FreeCAD 1.1.3 commit `145529fe741292ff0b3977a01195bf0247425794` was built from a
fresh clone with both Vinkulum persistence patches and the locked Linux SDK.
The final installer resumed that build after adding exclusive ownership and
surviving-session checks. CMake regenerated resources and rebuilt GUI libraries,
so all five qualification scenarios were repeated after that resume.

| Native scenario | Passed checks |
| --- | ---: |
| Exact placements across three save/reopen cycles | 36 |
| Geometry identity, edits and legacy migration | 39 |
| Rotated cylinder: real solve, persistence and invalidation | 3 |
| Independent-process reopening of that result | 3 |
| Static task interaction and worker lifecycle | 18 |

The final run hashes 1,605 executable/library/Python files before and after the
99 checks and confirms they remained unchanged. This manifest covers files in
the host build's `bin`, `lib` and `Mod` directories, not every external SDK file.
The retained SDK lock and build records identify those dependencies separately.

Ownership evidence includes the original resume counterexample, twelve isolated
process checks and the Pixi/CMake/Ninja counterexample: Ninja commands can change
process group while remaining in the recorded session. The final guard checks
the whole Linux session. Commands deliberately detaching into another session
are outside that contract. `ci/test_freecad_host.py` adds SDK-free regressions
to the mandatory local CI; it does not replace the actual SDK experiment.

The runner and captured paths describe the experiment machine and need adapting
for replay. This archive is evidence, not an installable or relocatable FreeCAD
distribution. The [host guide](../../FREECAD_HOST.md) supplies the build and
qualification commands. These bounded results do not certify all FreeCAD
workbenches, arbitrary mechanics or the complete Vinkulum kernel.
