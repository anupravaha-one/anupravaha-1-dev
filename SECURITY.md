# Security Policy

## Supported versions

The `main` branch is actively maintained. Tagged releases (`v2018` and any future stable tags) are immutable snapshots and do **not** receive security patches retroactively — if you depend on a release tag, plan to upgrade to a newer tag when a security fix lands on `main`.

## Reporting a vulnerability

**Please do not open a public issue for security-sensitive bugs.** Use one of the channels below instead:

1. **Private security advisory (preferred):** open a new draft advisory at
   <https://github.com/anupravaha-one/anupravaha-1-dev/security/advisories/new>.
   Only project maintainers see the report; GitHub coordinates the disclosure workflow, including CVE assignment if applicable, and lets us cut a fix in a private fork before publishing.

2. **Email fallback:** if for any reason the advisory form is unavailable to you, write to **me12m1030@gmail.com** with `[anupravaha security]` in the subject line.

We aim to:

- **Acknowledge** a security report within **7 days** of receipt.
- **Investigate and reproduce** within **14 days**.
- **Publish a fix** within **30 days** of confirmation, coordinated with the reporter where appropriate.

## Scope

Anupravaha-I is an offline CFD solver invoked from the command line. Realistic concerns within scope:

- Crafted `.cgns` mesh files that crash the solver (denial of service) or trigger memory corruption
- Buffer overruns in input parsing, the interactive prompt handler, or user-defined-function (UDF) interfaces
- Privilege issues when Anupravaha is integrated into automated pipelines (e.g. unsafe paths written by `Settings/` or `SOLUTIONS/`)
- Vulnerabilities in the vendored CGNS sources under `Cgns/`

**Out of scope:**

- Vulnerabilities in upstream system libraries (HDF5, the system `libcgns`) — please report those upstream.
- Issues that require an attacker to already have shell access on the machine running Anupravaha.

## Disclosure preferences

We'll coordinate with you on a disclosure timeline before publishing details. The default is private fix → public release notes → 30-day embargo lifted.
