[![Build Status](https://travis-ci.org/ThomasdenH/casimir-fdfd.svg?branch=master)](https://travis-ci.org/ThomasdenH/casimir-fdfd)
[![](http://meritbadge.herokuapp.com/casimir-fdfd)](https://crates.io/crates/casimir-fdfd)
[![Coverage Status](https://coveralls.io/repos/github/ThomasdenH/casimir-fdfd/badge.svg?branch=master)](https://coveralls.io/github/ThomasdenH/casimir-fdfd?branch=master)

# casimir-fdfd
An implementation of a stress-tensor based FDFD method for computing Casimir forces. This program was created for my bachelor thesis. The corresponding report, which includes relevant theory and results, can be read [here](https://gitlab.com/denhollander-thomas/casimir-report/-/jobs/artifacts/master/raw/Intersecting%20Discs.pdf?job=compile_pdf). The program is largely based on [Virtual photons in imaginary time: Computing exact Casimir forces via standard numerical-electromagnetism techniques](https://arxiv.org/abs/0705.3661).

## Usage
To use, install the program via `cargo install casimir-fdfd`. Alternatively, you can build it yourself via the source
code. To run, pass it the path to a configuration file. Examples of configurations can be found in
[/worlds/](https://github.com/ThomasdenH/casimir-fdfd/tree/master/worlds). To display a progress bar while running, use
the flag `--progressbar`/`-p`. You can always run `casimir-fdfd --help` for usage information.

Example usage:

`casimir-fdfd "worlds/plates_medium.json" -p`

## License
MIT License

Copyright (c) 2018 Thomas den Hollander

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
