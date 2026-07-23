# Model Selection

Select models through path sampling/stepping stone analysis

## Tutorials

[Path sampling](http://www.beast2.org/path-sampling/)

[Path sampling with a GUI](http://www.beast2.org/2014/07/14/path-sampling-with-a-gui.html)

## For developer

Build the package, skipping tests:

```bash
mvn clean package -DskipTests
```

Run an example XML through BEAST:

```bash
mvn exec:exec -Dbeast.args="src/test/resources/modelselection/examples/normalTest-1.xml"
```

## Releasing

### Jar folder structure

The release module jar, bundled into the BEAST package ZIP's `lib/`, contains:

```
META-INF/
├── MANIFEST.MF
└── maven/io.github.beast2-dev/model-selection/
    ├── pom.xml
    └── pom.properties
modelselection/
├── app/tools/          (PathSampler, PathSampleAnalyser, PairedPathSampler*, GeneralisedSteppingStone, ps.png)
├── core/                (CPOLogger)
├── cpo/                 (CPOAnalyser, BEASTRunAnalyser + inner classes)
├── fxtemplates/
│   └── ModelSelection.xml
├── gss/                 (GeneralisedSteppingStone*, MCMC2GSS, MCMC2IS, GSSFromFile, TraceLog, ...)
│   ├── coalescent/
│   └── distribution/
└── inference/           (PathSampler, PathSamplingStep, PairedPathSampler*, AICMAnalyser, DiffLogger, ...)
version.xml
module-info.class
```

`modelselection/fxtemplates/ModelSelection.xml` is loaded via module path resource
lookup at runtime — `src/assembly/beast-package.xml` does not copy it anywhere else in
the ZIP.

### 1. Maven Central release (JARs)

Push a `v*` tag to trigger `.github/workflows/ci-publish.yml`, which sets the Maven
version from the tag, builds, tests, GPG-signs, and publishes to Maven Central:

```bash
git tag v1.7.0
git push origin v1.7.0
```

Monitor the run at:
https://github.com/BEAST2-Dev/model-selection/actions/workflows/ci-publish.yml

### 2. GitHub release (BEAST package ZIP)

First remove `-SNAPSHOT` from `<version>` in `pom.xml` so it matches the release
(e.g. `1.7.0-SNAPSHOT` -> `1.7.0`), then build the installable BEAST package ZIP locally:

```bash
mvn clean package -DskipTests
```

**Note:** if you skip the manual edit above, the build still succeeds, but the module
jar bundled inside the ZIP (`lib/model-selection-<version>.jar`) will carry the
`-SNAPSHOT` suffix — that's meant for dev/testing builds, not an official release.
Alternatively, instead of hand-editing `pom.xml`, run:

```bash
mvn versions:set -DnewVersion=<version> -DgenerateBackupPoms=false
```

The ZIP is written to `target/MODEL_SELECTION.v<version>.zip`. Then manually:

1. Go to https://github.com/BEAST2-Dev/model-selection/releases/new
2. Choose the matching tag (e.g. `v1.7.0`), fill in the release title/notes
3. Upload `target/MODEL_SELECTION.v<version>.zip` as a release asset
4. Publish
