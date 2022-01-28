![Logo](plots/k4Clue_logo.png)

# k4Clue: The CLUE Algorithm in key4hep

## Standalone CLUE Algorithm on GPU and CPU

Algorithm gitlab repository [here](https://gitlab.cern.ch/kalos/clue).

Original authors:
Z.Chen, A. Di Pilato, F. Pantaleo, M. Rovere, C. Seez

If you encounter any error when compiling or running the k4Clue project, please contact:
* E. Brondolin, erica.brondolin@cern.ch
* M. Rovere, marco.rovere@cern.ch
* F. Pantaleo, felice.pantaleo@cern.ch

## Setup

### On a lxplus machine:

Source the key4hep environment:
```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh

# then setup this project
git clone --recurse-submodules https://gitlab.cern.ch/ebrondol/clue.git
cd clue
cmake -S . -B build
cmake --build build

# if installation is needed
cmake --install build/
```
If GPU/nvcc are available in the machine, the GPU version of CLUE will also be installed.
The path to the nvcc compiler will be automatically taken from the machine. nvcc <= 11.2 is needed. You can source it with:

```sh
# Get nvcc 11.2
source /cvmfs/sft.cern.ch/lcg/releases/cuda/11.2-5cee1/x86_64-centos7-gcc8-opt/setup.sh
```

## Run CLUE standalone
CLUE needs three parameters: `dc`, `rhoc` and `outlierDeltaFactor` (in the past four parameters were needed: `dc`, `deltao`, `deltac` and `rhoc`)

_dc_ is the critical distance used to compute the local density.
_rhoc_ is the minimum local density for a point to be promoted as a Seed.
_outlierDeltaFactor_ is  a multiplicative constant to be applied to `dc`.

( _deltao_ is the maximum distance for a point to be linked to a nearest higher
point.
 _deltac_ is the minimum distance for a local high density point to be promoted
as a Seed. )

### Standalone CLUE

If the projects compiles without errors, you can go run the CLUE algorithm by
```bash
cd build/src/clue
# ./main [fileName] [dc] [rhoc] [outlierDeltaFactor] [useParallel] [verbose]
./main ../../../data/input/aniso_1000.csv 20 25 2 0 1 1
```

The input files are `data/input/*.csv` with columns 
* x, y, layer, weight

The output files are `data/output/*.csv` with columns
* x, y, layer, weight, rho, delta, nh, isSeed, clusterId

### CLUE as Gaudi algorithm

If the projects compiles without errors, you can go run the CLUE algorithm by
```bash
cd build/
./run gaudirun.py ../gaudi_opts/clue_gaudi_wrapper.py
```

CLUE parameters and input/output file name are contained in `clue_gaudi_wrapper.py`.
The input files must contain data in the EDM4HEP format 
* `ECALBarrel` and `ECALEndcap` CalorimeterHit collections are required

The output file contains `ClueClusters` (currently also transformed as CaloHits).

## Run CLUE during the CLIC reconstruction

Here a simple recipe (from beginning to end):
```
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
git clone --recurse-submodules https://github.com/key4hep/k4Clue.git

git clone git@github.com:key4hep/k4MarlinWrapper.git
git clone https://github.com/iLCSoft/CLICPerformance

cd CLICPerformance/clicConfig
ddsim --steeringFile clic_steer.py --compactFile $LCGEO/CLIC/compact/CLIC_o3_v14/CLIC_o3_v14.xml --enableGun --gun.distribution uniform --gun.particle gamma --gun.energy 10*GeV --outputFile gamma_10GeV_edm4hep.root --numberOfEvents 10

cp ../../k4MarlinWrapper/test/gaudi_opts/clicRec_e4h_input.py .
k4run clicRec_e4h_input.py --EventDataSvc.input gamma_10GeV_edm4hep.root

#You can still visualise the output in slcio with:
ced2go -d ../Visualisation/CLIC_o3_v06_CED/CLIC_o3_v06_CED.xml -s 1 Output_REC.slcio

#Run CLUE in CLIC reconstruction
cp ../../k4Clue/gaudi_opts/clicRec_e4h_input_clue.py .
k4run clicRec_e4h_input_clue.py

#Run CLUE standalone
cp ../../k4Clue/gaudi_opts/clue_gaudi_wrapper.py .
k4run clue_gaudi_wrapper.py --EventDataSvc.input my_output.root --out.filename output_clue_standalone.root
```

In case you have changed something from the original repo and you have rebuild the package, you should use `source build/clueenv.sh` to make `k4run` aware of your new changes.
