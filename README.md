# MC-GCMC-Package-for-Simple-LJ-Particle-Systems
## Overview
The package named **GNM** enables MC simulations and GCMC simulations for multiple types of simple LJ particles, where **GNM** means **Grand-canonical N-particle Monte Carlo simulation**. It is appropriate as a basis for further MC simulation program writing or as a teaching tool for beginners and undergraduates. It is freely choose GCMC or MC simulation, the way to use it is similar with LAMMPS which needs to write an input file to set up all the simulation details, and only this input file is needed. There are only short-range LJ 12-6 potentials between particles, and the upper limit of particle types is very large at 500. The package is only available for Linux environments. All simulation processes have periodic boundaries. And the header files of random number are from [Numerical Recipes. 3rd Edition](https://books.google.com/books?hl=zh-CN&lr=&id=1aAOdzK3FegC&oi=fnd&pg=PA33&dq=Numerical+Recipes.+3rd+Edition&ots=3mPnDcFqpk&sig=4DGBPCGKsBsv8SgGqtlbiMQL3qQ).

---
## Usage
After the compilation, use the following command：

```./GNM {-input} {-output}```\
```--{-input} -in/-i input_file_name```\
```--{-output} -out/-o output_file_name```

For example: ```./GNM -i ./mc/in.dat -o ./mc/out```

NOTE: only the input part is required, and others will default as: 1) MC simulation; 2) output in the current path; 3) the file name starts with 'output'; 4) the configuration won't output.

---
## Input file writing guide
**Keywords:**
- **Ntype** {type} {atoms} * {chemical_potential}
- **Pair** {pair1} {pair2} {epsilon} {sigma} {*R<sub>cut</sub>*}
- **Initial_density** {Initial_density}
- **Boxlength** {Boxlength}
- **Temp** {Temp}
- **MCsteps** {MCsteps}
- **Outputsteps_atoms** {Outputsteps_atoms}`  `*#default -1 means won't output the configuration*
- **Outputfreq_atoms** {Outputfreq_atoms}`  `*#default 1 means output the configuration every steps*
- **Outputsteps_E** {Outputsteps_E}`  `*#default 0 means won't output the energy*
- **Outputfreq_E** {Outputfreq_E}`  `*#default 1 means output the energy every steps*
- **Maxdisplace** {Maxdisplace}`  `*#default is 0.5×boxlength*
- **Cellsize** {Cellsize}`  `*#default is 1.0, means the 1.0×R<sub>cut</sub>*
- **Randseed** {Randseed}`  `*#default is 32*
- **Fraction_GCMC** {Fraction_GCMC}`  `*#default is -1.0*

NOTE: You should write down these keywords line by line, and then write down your settings after each keyword, The keywords and setting values can be separated by these symbols:' :/%#'. The keywords marked with * are mandatory. If you don't set the rest keywords, the default value will be used. It should be noted that the 'Initial_density' and 'Boxlength' can't exist at the same time. And it is better not write down the output keywords if you also set it in command, because if these two conflict, they won't output anything. But if you decide to output, the step to start output must smaller than the simulation steps, and the output interval of configuration must be an integer multiple of the output interval of energy. Finally, there is no order requirement for the placement of allkeywords, and there is a complete example file for reference(**`example.in`**).

---

