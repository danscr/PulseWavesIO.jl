# PulseWavesIO.jl
Julia package for reading the PulseWaves lidar data format.

## Overview

This is a native Julia implementation to read the full waveform lidar data format [PulseWaves](https://github.com/PulseWaves/Specification/blob/master/specification.rst).

## Example usage

The PulseWaves format is a pair of files, one with ending `.pls`, and one with `.wvs`.
The `load` function takes the path to the `.pls` file, and throws an error if there is no `.wvs` file in the same directory with the same name (except the different file extension).

`load` returns the Pulse header, pulses, the (vestigial) waves header, the waveforms, and the sampling type.
The sampling type is not strictly part of the data, but it is returned as a convenience to distinguish system pulse samplings from return samplings.

Example code:
```
Using PulseWavesIO, FileIO

header, pulses, waves_header, waves, samplingType = load("path/to/pls/file")
```

## Roadmap

A goal of this project is to [register the implementation for PulseWaves with FileIO.jl](https://juliaio.github.io/FileIO.jl/stable/registering/) to make it available to users without any extra effort.
Before that happens, the code needs to be thoroughly tested and refined to ensure it works smoothly.

Therefore, __please report any issues you find or submit a pull request__. Thanks!
