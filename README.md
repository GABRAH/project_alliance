# project_alliance
Collaboration between Dr. Elisa Fadda and Dr. Jon Agirre groups

# [CAUTION: Cleaner repository accompanying the Nature Communications Demo, can be accessed in a seperate branch - https://github.com/GABRAH/project_alliance/tree/NatureCommunicationsDemo](https://github.com/GABRAH/project_alliance/tree/NatureCommunicationsDemo)

## Example commands for new version of [glycam2pdb.py](utility_scripts/glycam2pdb.py)

```sh
(privateerpython) harold@victoria:~/Dev/privateer_python/project_alliance/utility_scripts$ python glycam2pdb.py -h
```

```sh
(privateerpython) harold@victoria:~/Dev/privateer_python/project_alliance/utility_scripts$ python glycam2pdb.py -input ../glycampdbfiles/Volume -validate
```
Can run on an entire folder
```sh
(privateerpython) harold@victoria:~/Dev/privateer_python/project_alliance/utility_scripts$ python glycam2pdb.py -input ../glycampdbfiles/omannose/man9/cluster1.pdb -validate
```
Or can run on a single structure
```sh
(privateerpython) harold@victoria:~/Dev/privateer_python/project_alliance/utility_scripts$ python glycam2pdb.py -input ../glycampdbfiles/Volume 
```
If -validate flag is removed, the script will then not run Privateer after conversion to validate structures.

## Script [privateer_quick_validate.py](utility_scripts/privateer_quick_validate.py) that only does validation

```sh
(privateerpython) harold@victoria:~/Dev/privateer_python/project_alliance/utility_scripts$ python privateer_quick_validate.py -h
```

```sh
(privateerpython) harold@victoria:~/Dev/privateer_python/project_alliance/utility_scripts$ python privateer_quick_validate.py -input ../glycampdbfiles/omannoseConvertedPDB/
```
Can run on an entire folder

```sh
(privateerpython) harold@victoria:~/Dev/privateer_python/project_alliance/utility_scripts$ python privateer_quick_validate.py -input ../glycampdbfiles/omannoseConvertedPDB/man9/cluster1.pdb
```
Or on a single structure

