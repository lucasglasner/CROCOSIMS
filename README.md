## CEAZA CROCO SIMULATIONS CONFIGURATION FILES, ROUTINES AND POSTPROCESS

---

The purpose of this repository is to maintain a backup of my progress with the ocean simulations developed with the CROCO model. 
In the root directory i will put some essential functions and routines that i have developed for the use and postprocess of typical croco simulations. 
Given that routines, the subdirectories are related to specific projects and sometimes they will have their own specific routines, but mainly they will use common python packages like xarray, scipy, scikit-learn, etc. 


The current repository status (tree -L 2) is as follows: 
```
. 
├── load.py 
├── numerics.py 
├── README.md 
├── SIMSEQUIA 
│   ├── data 
│   ├── exploration_rund0control.ipynb 
│   ├── exploration_rund0rivers.ipynb 
│   ├── exploration_rund1control.ipynb 
│   ├── grids.ipynb 
│   ├── grids_r2r.ipynb 
│   └── plots 
├── simulation_status.sh 
├── TESTSIM 
│   ├── create_zarrstore.ipynb 
│   ├── croco_testsim_atmosphericforcing.ipynb 
│   ├── croco_testsim_exp1.ipynb 
│   ├── croco_testsim_exp2.ipynb 
│   ├── croco_testsim_exp3.ipynb 
│   ├── data 
│   └── plots 
└── utils.py 
```
---