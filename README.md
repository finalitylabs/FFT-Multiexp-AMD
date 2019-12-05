# FFT-Multiexp-AMD
A simple C repo for compilation and execution of GPU Multiexp and FFT algorithms for EC crypto.

This needs gcc v6 to compile, if using ubuntu 16.04 this guide can get it installed.
https://gist.github.com/zuyu/7d5682a5c75282c596449758d21db5ed

#### Build
``` bash
./build.sh
```

#### Generate Inputs
``` bash
./generate_inputs
```

### Run
``` bash
./main compute inputs outputs
```