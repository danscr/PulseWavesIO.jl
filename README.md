# PulseWavesIO.jl
Julia package for reading the PulseWaves lidar data format.

<img src="logo.svg" width="200">

## Overview

This is a native Julia implementation to read the full waveform lidar data format [PulseWaves](https://github.com/PulseWaves/Specification/blob/master/specification.rst).

## Example usage

The PulseWaves format is a pair of files, one with ending `.pls`, and one with `.wvs`.
The `load` function takes the path to the `.pls` file, and throws an error if there is no `.wvs` file in the same directory with the same name (except the different file extension).

`load` returns the Pulse header, pulses, the (vestigial) waves header, the waveforms, the sampling type, and a single mandatory appended variable length record.
The sampling type is not strictly part of the data, but it is returned as a convenience to distinguish system pulse samplings from return samplings.
Unless you want to write PulseWaves files, you can ignore the waves header and the appended variable length record (in the example below they are ignored by using `_` instead of naming the variable).

Example code:
```
Using PulseWavesIO, FileIO

header, pulses, _, waves, samplingType, _ = load("path/to/pls/file")
```

The pulse geometry (anchor and target XYZ) are read raw, and the function `PulsesToRays( pulses, header)` transforms coordinates into target values. The function returns an array of size `Nx2x3`, where N is the number of pulses, then the next dimension is anchor (1) or target (2) coordinates, and the last dimension is the Cartesian coordinates.
The target coordinates are a direction vector with the length of one sampling unit.

If the timestamps are needed, use the function `RescalePulseTimestamps( pulses, header)` to obtain GPS weektime data from the raw data. Note that there is a difference between GPS time and UTC of about 18-19 seconds, depending on what time your data was captured (search for "leap seconds" to learn more).

## Roadmap

A goal of this project is to [register the implementation for PulseWaves with FileIO.jl](https://juliaio.github.io/FileIO.jl/stable/registering/) to make it available to users without any extra effort.
Before that happens, the code needs to be thoroughly tested and refined to ensure it works smoothly.

Therefore, __please report any issues you find or submit a pull request__. Thanks!
