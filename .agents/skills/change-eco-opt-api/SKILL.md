---
name: change-eco-opt-api
description: Implement or review changes to eco_opt package behavior, public classes, result objects, simulation, phases, constraints, or transcription. Use for maintaining the library, not for writing ordinary user models.
---

# Change the eco_opt API

Treat the public exports listed in `AGENTS.md` as the compatibility boundary.
Preserve that API unless the user explicitly requests a change.

Before introducing a new abstraction, inspect the existing `System`, `Ocp`,
`Phase`, results, simulation, constraint, and collocation implementations
relevant to the request. Prefer extending an existing abstraction over adding
a parallel one.

For behavioral changes:

- add focused regression tests;
- preserve existing public behavior unless the request changes it;
- update public documentation and examples when usage changes;
- avoid encoding behavior for one scientific example into the library; and
- keep user-facing code dependent only on public `eco_opt` imports.

Keep internal implementation details out of user-facing examples and
documentation. Follow the lineage restrictions in `AGENTS.md`.

Run `python -m pytest -q`. Run `python -m mkdocs build` when public
documentation changes.
