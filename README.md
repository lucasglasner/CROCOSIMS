## CEAZA CROCO SIMULATIONS CONFIGURATION FILES, ROUTINES AND POSTPROCESS

---

The purpose of this repository is to maintain a backup of my progress with the ocean simulations developed with the CROCO model. 
In the root directory i will put some essential functions and routines that i have developed for the use and postprocess of typical croco simulations. 
Given that routines, the subdirectories are related to specific projects and sometimes they will have their own specific routines, but mainly they will use common python packages like xarray, scipy, scikit-learn, etc. 
Symlinks to data are removed now

The current repository status (tree -L 2) and architecture is as follows: 
```
. 
├── SIMSEQUIA
│   ├── data
│   ├── plots
│   └── *.ipynb
├── TESTSIM
│   ├── data
│   ├── plots
│   └── *.ipynb
├── MOREPROJECTS...
├── simulation_status.sh
├── utils.py
├── load.py
└── numerics.py
```
---
