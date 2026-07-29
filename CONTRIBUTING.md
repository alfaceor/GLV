# TODO

* [x] How to set up the dev environment
* [x] How to run tests
* [ ] Coding standards / linting
* [ ] Branching strategy
* [ ] How to submit a PR
* [ ] Commit message conventions (if any)

# Environment setup

## VScode devcontainer

## Install developer version

Install the dependencies to run pytest and benchmarks

```bash
pip install -e .[dev]
```

Unit test are written using pytest

```bash
pytest test
```

To run performance test

```bash
pytest benchmarks
```


