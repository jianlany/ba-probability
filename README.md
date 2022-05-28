This program calculates the bond-angle probability density distributions with kernel density estimates. Accelerated with GPU.

## Compilation
Replaced compiler path in `Makefile` with the installed Nvidia c++ compiler path and run `make`.
Note that it compiles the code with GPU acceleration by default, make sure the cuda version is compatible with the installed GPU.

## Usage
### command:
`ba-probability [path to ba file] (options)`.
### options:

`--num_theta`: number of bins on angle direction. Default: 100

`--num_length`: number of bins on bond direction. Default: 100

`--d_width_factor`: number of bins that the gaussion distribution will spread over on bond direction. Default: 1.5

`--q_width_factor`: same but on angle direction. Default: 1.5

`--d_range`: the upper and lower limit of bond. Default: 2.0  2.9

`--theta_range`: the upper and lower limit of angle. Default: 0  180

`--output-file`: the path of the output badf file. Default: p.txt

`--entropy`: enable entropy scaling. Default: false

`--renormalize`: renormalize the p so that it integrates to 1. Default: false

`--help`: print this message.

The input BA file should contain the bond angle samples, each line represents one sample and must be formatted as `l1 l2 angle`
