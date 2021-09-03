<!-- [![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url] -->



<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/GABRAH/project_alliance/tree/NatureCommunicationsDemo">
    <img src="images/figure1_panels_hr.png" alt="Logo" width="600" height="300">
  </a>

  <h3 align="center">The case for post-predictional modifications in the AlphaFold Protein Structure Database</h3>

  <p align="center">
    Supplementary Data
    <br />
    <a href="https://github.com/glycojones/privateer/issues">Report Bug</a>
    ·
    <a href="https://github.com/glycojones/privateer/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-grafting-script">About The Grafting Script</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>



## About The Grafting Script

The grafting capabilities showcased here are currently a prototype and are expected to be significantly improved upon in the future as part of Privateer software. Therefore, feedback is most welcome.

The demo [grafter.py](privateer_grafting_demo/grafter.py) script is a user-friendly rendition to get Privateer to graft N-glycans onto AlphaFoldDB models with a few command line arguments.

In order to run this script, a **development version** of [Privateer](https://github.com/glycojones/privateer/tree/privateerpython) software must be obtained that is not yet officially released through [CCP4](https://www.ccp4.ac.uk/) and [CCP-EM](https://www.ccpem.ac.uk/).

**Supported operating systems by the [development version of Privateer](https://github.com/glycojones/privateer/tree/privateerpython):**
* **Linux**(tested on Ubuntu Linux 20.04)
* **MacOS**(tested on MacOS "Catalina" 10.15.5)

**Currently unsupported operating systems by the [development version of Privateer](https://github.com/glycojones/privateer/tree/privateerpython):**
* **Windows**

## Getting Started

To get a local copy up and running follow the steps below:

### Prerequisites

Assuming a fresh installation of **Ubuntu Linux 20.04**, the following system packages need to be installed via:
```sh
sudo apt-get install cmake
sudo apt-get install bzr
sudo apt-get install virtualenv
sudo apt-get install gfortran
sudo apt-get install m4
```

On MacOS, this can be achieved via [Homebrew](https://brew.sh/) package manager and running the following commands afterwards:
```sh
brew install wget
brew install cmake
brew install bzr
brew install virtualenv
brew install gfortran
brew install m4
```

After the required system packages are installed, **[development version of Privateer](https://github.com/glycojones/privateer/tree/privateerpython)** needs to be compiled and built from scratch. For the grafting demo of AlphaFoldDB models, the installation procedure is composed of 2 steps.

### **1.) Installation of Privateer**

1. Clone Privateer repo into *privateer_python* directory:
   ```sh
   git clone https://github.com/glycojones/privateer.git privateer_python
   ```
2. Change directory to *privateer_python*:
   ```sh
   cd privateer_python
   ```
3. Switch branch to [privateerpython](https://github.com/glycojones/privateer/tree/privateerpython) to obtain the latest developments of Privateer software:
   ```sh
   git checkout privateerpython
   ```
4. Get Privateer dependencies that are hosted on GitHub([gemmi](https://github.com/project-gemmi/gemmi), [nlohmann::json](https://github.com/nlohmann/json), [pybind11](https://github.com/pybind/pybind11), [pybind11_json](https://github.com/pybind/pybind11_json)):
   ```sh
   git submodule update --init --recursive
   ```
5. Create a virtualenv for Python3 interpreter that will contain Privateer's [pybind11](https://github.com/pybind/pybind11) bindings for **C++11 backend**:
   ```sh
   virtualenv privateerpython
   ```
6. Source [CCP4](https://www.ccp4.ac.uk/) environment variables that are used by Privateer's [CCP4](https://www.ccp4.ac.uk/) dependencies:
   ```sh
   source ccp4.envsetup-sh
   ```
    **CAUTION: THIS STEP HAS TO BE REPEATED EVERY TIME A NEW TERMINAL TAB/WINDOW IS OPENED**<br/>
    Secondary caution: The output of this step may appear error-like. This is normal behaviour, as shown in the following picture:
    <p align="left">
    <img src="images/ccp4sourcescriptoutput.png" alt="CCP4Source" width="923" height="183"></p>
7. Source Python3 interpreter from virtualenv that contains [pybind11](https://github.com/pybind/pybind11) bindings to Privateer's **C++11 backend**:
   ```sh
   source privateerpython/bin/activate
   ```
    **CAUTION: THIS STEP HAS TO BE REPEATED EVERY TIME A NEW TERMINAL TAB/WINDOW IS OPENED**
8. Install necessary Python modules through **pip3** to local **privateerpython** interpreter
   ```sh
   pip install -r requirements.txt
   ```
9. Compile/build development version of Privateer along with its dependencies:
   ```sh
   python setup.py install
   ```
    Successfull installation of Privateer and its dependencies in the Terminal should appear in the following way:
    <p align="left">
    <img src="images/successfullinstall.png" alt="SuccessfulInstall" width="723" height="537"></p>

### **2.) Setup of the AlphaFoldDB glycan grafting demo script**

1. In *privateer_python* directory from **substep 2** of previous installation step, clone [project_alliane](https://github.com/GABRAH/project_alliance.git) GitHub repository:
   ```sh
   git clone https://github.com/GABRAH/project_alliance.git grafting_demo
   ```
2. Change directory to *grafting_demo*:
   ```sh
   cd grafting_demo
   ```
3. Change git branch to **[NatureCommunicationsDemo](https://github.com/GABRAH/project_alliance/tree/NatureCommunicationsDemo)**:
   ```sh
   git checkout NatureCommunicationsDemo
   ```

**After following these steps, the installation should be complete!**


<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_



<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/github_username/repo_name/issues) for a list of proposed features (and known issues).



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Your Name - [@twitter_handle](https://twitter.com/twitter_handle) - email

Project Link: [https://github.com/github_username/repo_name](https://github.com/github_username/repo_name)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* []()
* []()
* []()





<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/github_username