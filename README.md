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

Build the installable BEAST package ZIP locally:

```bash
mvn clean package -DskipTests
```

The ZIP is written to `target/MODEL_SELECTION.v<version>.zip`. Then manually:

1. Go to https://github.com/BEAST2-Dev/model-selection/releases/new
2. Choose the matching tag (e.g. `v1.7.0`), fill in the release title/notes
3. Upload `target/MODEL_SELECTION.v<version>.zip` as a release asset
4. Publish
