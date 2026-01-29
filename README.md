## Faster Boomerangs in Feistel Ciphers: WARP, LILLIPUT, LBlock-s

Source code of the paper *Towards Faster Boomerangs in Feistel Ciphers*.

## Model
MILP model for truncated boomerang characteristic.

## Experiments
Experimental evaluation of the distinguisher probability.
### Usage (running WARP_Boom.cpp as an example)

```
g++ -O3 -fopenmp WARP_Boom.cpp -o warp_boom
export OMP_NUM_THREADS=256
./warp_boom
```

## Acknowledgements
This project is based on the tools developed in [[HNE22](https://tosc.iacr.org/index.php/ToSC/article/view/9858)] and [[ZTX25](https://tosc.iacr.org/index.php/ToSC/article/view/12084)], with extensions.
