This repository contains code used for generation of initial structures for simulation.

With the following command:

```
BYP_num=80
TMA_num=80

bash assemble.sh $BYP_num $TMA_num
```

A structure file `full_structure160.pdb` should be produced. 

Isosahedral core was constructed with *isosahedral builder* from  [OpenMD](http://openmd.org/openmd-version-2-5/).

