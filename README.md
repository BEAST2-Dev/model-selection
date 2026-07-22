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
