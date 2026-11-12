# FairRoot

[![license](https://alfa-ci.gsi.de/shields/badge/license-LGPL--3.0-orange.svg)](COPYRIGHT) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3896282.svg)](https://doi.org/10.5281/zenodo.3896282) [![OpenSSF Best Practices](https://bestpractices.coreinfrastructure.org/projects/7131/badge)](https://bestpractices.coreinfrastructure.org/projects/7131) [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8F-yellow)](https://fair-software.eu)

A C++ framework for the simulation, reconstruction and analysis of particle physics
experiments, based on [ROOT](https://root.cern/).

**Documentation: <https://dennisklein.github.io/FairRoot/>**

FairRoot gives the base classes for the software of an experiment:

* **Simulation.** Describe the detectors with *modules* and *detectors*. The Virtual
  Monte Carlo (VMC) transports the particles with Geant3 or Geant4, and your detector
  code does not depend on the engine.
* **Reconstruction and analysis.** *Tasks* read branches of the event tree and add new
  branches. A *run* calls the tasks for each event.
* **Parameters.** The *runtime database* reads and writes the parameter containers for
  each run.
* Time-based processing, an event display, track propagation, and online processing
  with FairMQ.

## Install

FairRoot needs ROOT, VMC, FairLogger and some more packages. The distribution
[FairSoft](https://github.com/FairRootGroup/FairSoft) contains all of them. Build FairRoot on top of
it with CMake:

```bash
export SIMPATH=/path/to/fairsoft
git clone -b dev https://github.com/FairRootGroup/FairRoot.git
cd FairRoot
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$HOME/fairroot_install"
cmake --build build --parallel
cmake --install build
```

The [install page](https://dennisklein.github.io/FairRoot/start/install/) lists the dependencies,
the build options and how to run the tests.

## Get started

* [Your first simulation](https://dennisklein.github.io/FairRoot/start/first-simulation/): run the simulation of Tutorial 1 and look at the output.
* [Your first task](https://dennisklein.github.io/FairRoot/start/first-task/): write a task, build it with a ROOT dictionary and run it.
* [Start a project from a template](https://dennisklein.github.io/FairRoot/start/new-project/): copy a project template in [`templates/`](templates) and rename it.
* [Examples](https://dennisklein.github.io/FairRoot/examples/overview/): the programs and macros in [`examples/`](examples).
* [API reference](https://dennisklein.github.io/FairRoot/api/) and [source browser](https://dennisklein.github.io/FairRoot/source/).

## Releases and support

Find all releases [on GitHub](https://github.com/FairRootGroup/FairRoot/releases) and the
changes in the [release notes](https://dennisklein.github.io/FairRoot/project/release-notes/).
Report bugs and ask questions in the [issues](https://github.com/FairRootGroup/FairRoot/issues).
Read [Contributing](CONTRIBUTING.md) before you open a pull request.

## License

FairRoot is free software under the [LGPL-3.0](LICENSE). See [COPYRIGHT](COPYRIGHT) and
[how to cite FairRoot](https://dennisklein.github.io/FairRoot/project/license/).
